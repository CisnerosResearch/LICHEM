/*
###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                 By: Eric G. Kratz, Hatice Gokcan, Alice Walker,             #
#                     Erik Montelongo Vazquez, G. Andres Cisneros             #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################
*/


/*
  ###############################################################################
  #                                                                             #
  # Hatice GOKCAN                                                               #
  #                                                                             #
  # Functions for QSM Optimizations.                                            #
  # Includes:                                                                   #
  #                                                                             #
  #        Initialize hessian                : void init_Hess                   #
  #                                                                             #
  #        DBFGS algorithm for Hessian update: void DBFGS                       #
  #                                                                             #
  #        Check if spaceout is necessary    : void eqconst                     #
  #                                                                             #
  #        Spaceout tools                    :                                  #
  #                                            void spaceoutcubic               #
  #                                            void fullNCS                     #
  #                                            void getcoeff                    #
  #                                                                             #
  #                                                                             #
  #        ODE solver                        : void ODESolve                    #
  #                                                                             #
  #        Cash-Karp 4/5 Order               : void RKStep                      #
  #                                                                             #
  #        Quadratic approximation to                                           #
  #        energy and gradient               : void quad_app                    #
  #                                                                             #
  #        Update forces                     : void funupwind                   #
  #                                                                             #
  #        Calculate tangents                : void calc_tangents               #
  #                                                                             #
  #        Calculate trust radius            : void updateTR                    #
  #                                                                             #
  #        Update coordinates                : void updatepath                  #
  #                                                                             #
  #        Aligned?                          : void alignment                   #
  #        Construct path btw react and prod : void construct_path              #
  #        Linear interpolation              : void linearly_interpolate        #
  #                                                                             #
  ###############################################################################
*/


//START 
//INITIALIZE THE MATRIX THAT CONTAINS ALL HESSIANS IF IT IS THE FIRST TIME
void init_Hess(MatrixXd& Hessmat, int beadsize, int Nimages){

   MatrixXd subHess(beadsize,beadsize);
   subHess = subHess.setIdentity();

   //#pragma omp parallel for schedule(dynamic)
   for (int k = 0; k < Nimages+2; k++){
        int srow=k*beadsize;//starting row of current image
        int rsize=beadsize;//number of rows
        int csize=beadsize;//number of columns 
        Hessmat.block(srow,0,rsize,csize)=subHess; 
   }
   //#pragma omp barrier
   

}
//END

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//START: DAMPED BROYDEN-FLETCHER-GOLDFARB-SHANNO (DBFGS)
//
//Solves Eqn. 14 in 
//Burger and Yang,J. Chem. Phys. 124 054109 (2006)
//
//NOTE:  Eqn. 14 is normal BFGS.
//       In DBFGS "r" vector (Eqn. 19a and 19b) 
//       is used instead of "gamma" vector.				    			
void DBFGS(MatrixXd& Hess, VectorXd& pathdiff, VectorXd& graddiff, int ndim){

        double theta;
        double shs;
        VectorXd r((ndim));//*3);
        MatrixXd Ann(ndim,ndim);//((ndim)*3,(ndim)*3);
        MatrixXd Ann2(ndim,ndim);//((ndim)*3,(ndim)*3);
        MatrixXd Ht(ndim,ndim);//((ndim)*3,(ndim)*3);

        if(pathdiff.norm() > 1e-6)
	{
           if(pathdiff.dot(graddiff)>=0.2*(pathdiff.dot(Hess*pathdiff)))
	   {
           	theta = 1.0;
           }
           else
	   {
                theta = pathdiff.dot(Hess*pathdiff) - pathdiff.dot(graddiff);
                theta = 0.8*pathdiff.dot(Hess*pathdiff)/theta;
           }

           r = theta*graddiff + ((1-theta)*(Hess*pathdiff));
           Ann = r*(r.transpose());
           Ann2 = Ann/(pathdiff.dot(r));
           Ann = pathdiff*(pathdiff.transpose());
           Ht = Hess*Ann;
           shs = pathdiff.dot(Hess*pathdiff);
           Ann = (Ht*Hess)/shs;
           Hess = Hess + Ann2 - Ann;
	}

}
//END: DAMPED BROYDEN-FLETCHER-GOLDFARB-SHANNO (DBFGS)

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//START:: EQCONST

void eqconst(VectorXd& path, int Nimages, int beadsize, double &eqcons){
  
     VectorXd pathdiff1(beadsize);
     VectorXd pathdiff2(beadsize);
     //#pragma omp parallel for schedule(dynamic)
     for(int k=1; k<Nimages+1; k++){
         pathdiff1 = path.segment(beadsize*k,beadsize)    
                    -path.segment(beadsize*(k-1),beadsize);
         pathdiff2 = path.segment(beadsize*k,beadsize)     
                    -path.segment(beadsize*(k+1),beadsize);

         //eqcons = eqcons+abs( pathdiff1.norm() - pathdiff2.norm());
         eqcons = eqcons + 
                  abs(sqrt((pathdiff1.array().square()).sum()) - 
                      sqrt((pathdiff2.array().square()).sum()));
     }
     //#pragma omp barrier
}
//END:: EQCONST

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//START:: SPACEOUTCUBIC
void spaceoutcubic(VectorXd& wholepath,bool nebatoms[], \
                     int N, int Ndim, VectorXd& weight){
 
   //N: Nimages+2
   //Ndim: beadsize=QMdim*3
   //Ndim: beadsize=Natoms*3

   VectorXd x(N*Ndim);
   int natm=Ndim/3;
   VectorXd dists(N);//dists(N+2);
   
   VectorXd a(N); //matrix(:,i) --> i. column  
   VectorXd b(Ndim);//matrix(i,:)  --> i. row
   MatrixXd coeffmat(N,4);
   MatrixXd tmpcoeff(N,4);
   MatrixXd coeffmat2(N*Ndim,4);
   double space;
   int crd;
   int index;

   dists(0)=0;
   x=wholepath.cwiseProduct(weight);   
   
   VectorXd xcurr(Ndim);
   VectorXd xprev(Ndim);
   //#pragma omp parallel for schedule(dynamic)
   for(int k=1;k<N;k++){//image
      xcurr = x.segment(Ndim*k,Ndim); 
      xprev = x.segment(Ndim*(k-1),Ndim);
      dists(k) = dists(k-1) + (xcurr-xprev).norm(); //pathdiff.norm();
 
   }
   //#pragma omp barrier


   space=dists(N-1)/(N-1);
   x=wholepath;

   index=0;
   //#pragma omp parallel for schedule(dynamic)
   for(int i=0;i<Ndim;i++){//Ndim=beadsize
      for(int k=0;k<N;k++){//N=NBeads
         a(k) = x(i+k*Ndim);  
      } 
      getcoeff(dists,a,N,Ndim,coeffmat);//for ith atom

      //if the atom is not nebatom 
      //set the last three columns to zero
      //otherwise it will stay same
      int curratom = floor(i/3);//current atom ID
 
      if(!nebatoms[curratom]){
            tmpcoeff.col(0) = coeffmat.col(0);
            tmpcoeff.col(1) = coeffmat.col(1)*0.0;
            tmpcoeff.col(2) = coeffmat.col(2)*0.0;
            tmpcoeff.col(3) = coeffmat.col(3)*0.0;
      }
      else{
            tmpcoeff.col(0) = coeffmat.col(0);
            tmpcoeff.col(1) = coeffmat.col(1);
            tmpcoeff.col(2) = coeffmat.col(2);
            tmpcoeff.col(3) = coeffmat.col(3);
      }


      int srow=i*N;//starting row of current image
      int rsize=N;//number of rows
      int csize=4;//number of columns
      
      coeffmat2.block(srow,0,rsize,csize)=tmpcoeff;

   }
   //#pragma omp barrier

   //#pragma omp parallel for schedule(dynamic)
   //in f90
   //do i=1,n-2
   for(int k=1;k<N-1;k++){//loop images
      //in f90
      //call fullNCS(coeffs,dists,space*i,ndim,n,x(i+1,:)) 
      //x(1) reactant, x(n) reactant
      //i=1 x(2) and space*1
      //i=n-2 x(n-1) and space*(n-2)
      //
      //start from first image not react
      //x(0) is reactant, x(n-1) is product
      //start from x(1)
      b = x.segment(k*Ndim,Ndim);
      //k=1 x(1) and space*1
      //k=n-2 x(n-2) space*(n-2)
      fullNCS(coeffmat2,dists,(space*k),Ndim,N,b);
      x.segment(k*Ndim,Ndim) = b;
   }
   //#pragma omp barrier
 
   /* did not work eigen gave error
   VectorXf xx;
   VectorXd xxx;
   for(int i=0;i<x.size();i++){
       xx(i) = (float) x(i);
       xxx(i) = (double) xx(i);
   }
   wholepath = xxx;   */
   wholepath=x;  
  

}
//END:: SPACEOUTCUBIC

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void fullNCS(MatrixXd& coeffs, VectorXd& dists,
               double t, int ndim, int n, VectorXd& x)
{

//x: out
     int i;
     double crd;

     i=0;//i=1 
     //#pragma omp parallel for schedule(dynamic)
     for(int j=0; j<ndim; j++){
        //do{
        //  i=i+1;
        //}while(t>x(i) and (i<(n-1)));
        //}while(t>dists(i) and (i<(n-1)));
        while(t>dists(i) and (i<(n-1))){
            i=i+1;
        }
        i=i-1;
        //crd = x(i);
        crd = dists(i);
        x(j)= coeffs(j*n+i,0) + (t-crd)*coeffs(j*n+i,1) + \
              coeffs(j*n+i,2)*pow(t-crd,2) + \
              coeffs(j*n+i,3)*pow(t-crd,3);
     }
     //#pragma omp barrier
   
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void getcoeff(VectorXd& x, VectorXd a, int n, int Ndim, MatrixXd& coeffmat){
//n=Nimages+2
//a: coords of atoms in different images

  VectorXd c = VectorXd::Zero(n);
  VectorXd u = VectorXd::Zero(n);
  VectorXd z = VectorXd::Zero(n);
  VectorXd l = VectorXd::Ones(n);
  VectorXd d = VectorXd::Zero(n);
  VectorXd b = VectorXd::Zero(n);
  VectorXd h(n-1);
  VectorXd alpha(n-1);

  //#pragma omp parallel for schedule(dynamic)
  for(int k=0;k<(n-1);k++){//image i=1,n-1
     h(k)=x(k+1)-x(k);
     if(k>0){
       alpha(k-1) = 3/h(k)*(a(k+1)-a(k)) - 3/h(k-1)*(a(k)-a(k-1));
     }    
  }
  //#pragma omp barrier

  //#pragma omp parallel for schedule(dynamic)
  for(int k=1;k<(n-1);k++){
     l(k)=2*(x(k+1)-x(k-1))-h(k-1)*u(k-1);
     u(k)=h(k)/l(k);
     z(k)=(alpha(k-1)-h(k-1)*z(k-1))/l(k);
  } 
  //#pragma omp barrier

  //#pragma omp parallel for schedule(dynamic)
  for(int k=n-2;k>-1;k--){
     c(k)=z(k)-u(k)*c(k+1);
     b(k)=(a(k+1)-a(k))/h(k)-h(k)*(c(k+1)+2*c(k))/3;
     d(k)=(c(k+1)-c(k))/(3*h(k));
  }
  //#pragma omp barrier

  coeffmat.col(0) = a;
  coeffmat.col(1) = b; 
  coeffmat.col(2) = c;
  coeffmat.col(3) = d;


}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//START:: ODE SOLVER
void ODESolve(VectorXd& path, MatrixXd& H, VectorXd& g, VectorXd& energy, \
       int Nimages, int beadsize, VectorXd& trs, VectorXd& cons, int wholesize,\
       double ftol, double &dftot,VectorXd& weight,fstream& logFile){

   double h=0.001;  // initial h
   double hnew = h;

   double error0 = 1e-7;
   double error1;
   int maxiter = 2000;       //  Maximum number of integration steps allowed
   int iter;

   VectorXd initialpath(wholesize);
   VectorXd newpath(wholesize);
   VectorXd pathdiff(beadsize);
   VectorXd htot = VectorXd::Zero(Nimages);// total time of integration

   bool currswitch[Nimages];
   bool prevswitch[Nimages];
   for(int i=0; i<Nimages;i++){
       currswitch[i]=true;
   }

   VectorXd moved(Nimages);
   VectorXd df(Nimages);

   VectorXd equad(Nimages);
   VectorXd gquad(Nimages*beadsize);  
   VectorXd force(Nimages*beadsize);
   VectorXd tmpforce(beadsize);
 
   bool done, allswitch;
   double dftol;
   dftol = ftol/10;      // assumed that L<100 where |x - y| < L|y(x) - y(y)|
   initialpath = path;

 
   for(iter=0; iter < maxiter; iter++) {
      
      error1 = 1;

      for(int ii=0;ii<Nimages;ii++){
         prevswitch[ii]=currswitch[ii];
      }
      do{

        h = hnew;

        for(int ii=0;ii<Nimages;ii++){
           currswitch[ii]=prevswitch[ii];
        }

        RKStep(path,initialpath,H,g,energy,Nimages,beadsize,
               h,error0,error1,newpath,hnew,currswitch,cons,weight);
        
        if(h<1.e-11){
           logFile << "                 ";
           logFile << "Step size has become too small h = " << h << endl;
           break; //exit(1);
        }

      } while(error1 > error0);


      //Check status of all points

      path=newpath; //newpath comes from RKstep
      
      quad_app(path,initialpath,H,g,energy,Nimages,beadsize,equad,gquad); 
  
      funupwind(path,gquad,equad,Nimages,beadsize,force,weight); 

      VectorXd tmpdiff(wholesize);
      tmpdiff = initialpath - path;
      
      //#pragma omp parallel for schedule(dynamic)
      for(int k=1; k<Nimages+1; k++){
         pathdiff = initialpath.segment(beadsize*k,beadsize)
                  - path.segment(beadsize*k,beadsize);

         tmpforce = force.segment(beadsize*(k-1),beadsize);

         double pdsss;
         //pdsss = pathdiff.norm();
         pdsss = sqrt((pathdiff.array().square()).sum());
         moved(k-1) = pdsss; 
         double fsss;
         //fsss = tmpforce.norm();
         fsss = sqrt((tmpforce.array().square()).sum());
         df(k-1) = fsss;
         //if(currswitch(k-1)){
         if(currswitch[k-1]){
           htot(k-1) = htot(k-1) + h;
         }
      }
      //#pragma omp barrier

      done = false;
      allswitch = true;

      //#pragma omp parallel for schedule(dynamic)
      for(int i=0; i<Nimages; i++){
         if(moved(i)>trs(i)){
            done = true;
         }
         //allswitch = allswitch and (!currswitch(i));
         allswitch = (allswitch and (not currswitch[i]));
         
      }
      //#pragma omp barrier

      dftot = df.maxCoeff();
      if(done or allswitch or (dftot<dftol)){
        //logFile << "                 ";
        //logFile << "Finished integration at step " << iter << "." << endl;

        return;
      }
     
   }

   if(iter == maxiter){
          logFile << "                 ";
          logFile << "ODESolve - Exceed number of integration steps"<< endl;
   }
   /*
   logFile << "                   RKStep          : " << "FINISHED"<< endl;
   logFile << "                     ";
   logFile << "ITER          : ";
   logFile << iter << endl;
   logFile << "trs           : \n";
   logFile << trs << endl;
   logFile << "total movement: \n";
   logFile  << moved << endl;
   logFile << "df            : \n"; 
   logFile << df << endl;
   logFile << "dftot         : \n"; 
   logFile << dftot << endl; 
   logFile << "---------------------------------\n"<<endl; 
   */

   //modify const for the next run
   //#pragma omp parallel for schedule(dynamic)  
   for(int i=0;i<Nimages;i++){
      if(!currswitch[i]){
        cons(i) = cons(i)*max(0.5,htot(i)/(htot.maxCoeff()));  
      }
      else{
        if( (moved(i)/trs(i))>0.5 ){
          cons(i) = cons(i)*1.0/(moved(i)/trs(i)); 
        }  
        else{
          cons(i) = cons(i) * 2; 
        }
      }

   }
   //#pragma omp barrier
   
} 
//END:: ODE SOLVER
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//START:: CASH-KARP 4/5 ORDER
void RKStep(VectorXd& path,VectorXd& epath,MatrixXd& H,VectorXd& g,VectorXd& energy,\
              int Nimages,int beadsize,double h,double error0,double &error1,VectorXd& x5,\
              double &hnew,bool currswitch[],VectorXd& cons,VectorXd& weight){


   int wholesize = (Nimages+2)*beadsize;
   MatrixXd ks = MatrixXd::Zero(6,(Nimages+2)*beadsize);  

   VectorXd sumks(wholesize); 
   VectorXd path1(wholesize);
   VectorXd x4(wholesize);
   VectorXd equad(Nimages);
   VectorXd gquad(Nimages*beadsize);
   VectorXd force(Nimages*beadsize);
   VectorXd force1(Nimages*beadsize);
   VectorXd pathdiff(wholesize);

   if(h==0){
     x5=path;
     error1=0;
     return; 
   }

   VectorXd Gforce((Nimages+2)*beadsize);
   VectorXd hGforce((Nimages+2)*beadsize);
   VectorXd tmpf(Nimages*beadsize);
   VectorXd Gcons(Nimages*beadsize);
   //create global cons array Gcons
   for(int k=0;k<Nimages; k++){
       for(int j=0;j<beadsize;j++){
         Gcons(k*beadsize+j) = cons(k);
       }
   }
   
   // CASH-KARP 4/5 ORDER
   //STEP::1   
   quad_app(path,epath,H,g,energy,Nimages,beadsize,equad,gquad);
   funupwind(path,gquad,equad,Nimages,beadsize,force,weight);


   Gforce.setZero();
   hGforce.setZero();
   tmpf.setZero();
   tmpf=force.cwiseProduct(Gcons);
   Gforce.segment(beadsize,Nimages*beadsize) = tmpf;
   
   hGforce = -h*Gforce;
   //hGforce = -h*(Gforce.cwiseProduct(weight));
   ks.row(0) = hGforce;

   sumks = 0.2*ks.row(0);
   path1 = path + sumks;
   //----------------------------------------------------------
   //STEP::2
   quad_app(path1,epath,H,g,energy,Nimages,beadsize,equad,gquad);
   funupwind(path1,gquad,equad,Nimages,beadsize,force,weight); 
   
   Gforce.setZero();
   hGforce.setZero();
   tmpf.setZero();
   tmpf=force.cwiseProduct(Gcons);
   Gforce.segment(beadsize,Nimages*beadsize) = tmpf;
   hGforce = -h*Gforce;
   //hGforce = -h*(Gforce.cwiseProduct(weight));
   ks.row(1) = hGforce;

   sumks = 0.075*ks.row(0) + 0.225*ks.row(1);
   path1 = path + sumks;
   //----------------------------------------------------------
   //STEP::3
   quad_app(path1,epath,H,g,energy,Nimages,beadsize,equad,gquad);
   funupwind(path1,gquad,equad,Nimages,beadsize,force,weight);

   Gforce.setZero();
   hGforce.setZero();
   tmpf.setZero();
   tmpf=force.cwiseProduct(Gcons);
   Gforce.segment(beadsize,Nimages*beadsize) = tmpf;
   hGforce = -h*Gforce;
   //hGforce = -h*(Gforce.cwiseProduct(weight));
   ks.row(2) = hGforce;

   sumks = 3.0/10.0*ks.row(0) - 9.0/10.0*ks.row(1) + 6.0/5.0*ks.row(2);
   path1 = path + sumks;
   //----------------------------------------------------------
   //STEP::4
   quad_app(path1,epath,H,g,energy,Nimages,beadsize,equad,gquad);
   funupwind(path1,gquad,equad,Nimages,beadsize,force,weight);

   Gforce.setZero();
   hGforce.setZero();
   tmpf.setZero();
   tmpf=force.cwiseProduct(Gcons);
   Gforce.segment(beadsize,Nimages*beadsize) = tmpf;
   hGforce = -h*Gforce;
   //hGforce = -h*(Gforce.cwiseProduct(weight));
   ks.row(3) = hGforce;

   sumks = - 11.0/54.0*ks.row(0) + 5.0/2.0*ks.row(1) 
           - 70.0/27.0*ks.row(2) + 35.0/27.0*ks.row(3);
   path1 = path + sumks;

   //----------------------------------------------------------
   //STEP::5
   quad_app(path1,epath,H,g,energy,Nimages,beadsize,equad,gquad);
   funupwind(path1,gquad,equad,Nimages,beadsize,force,weight);

   Gforce.setZero();
   hGforce.setZero();
   tmpf.setZero();
   tmpf=force.cwiseProduct(Gcons);
   Gforce.segment(beadsize,Nimages*beadsize) = tmpf;
   hGforce = -h*Gforce;
   //hGforce = -h*(Gforce.cwiseProduct(weight));
   ks.row(4) = hGforce;

   sumks = 1631.0/55296.0*ks.row(0) + 175.0/512.0*ks.row(1) + 575.0/13824.0*ks.row(2) + 44275.0/110592.0*ks.row(3) + 253/4096.0*ks.row(4);
   path1 = path + sumks;

   //----------------------------------------------------------
   //STEP::6
   quad_app(path1,epath,H,g,energy,Nimages,beadsize,equad,gquad);
   funupwind(path1,gquad,equad,Nimages,beadsize,force,weight);

   Gforce.setZero();
   hGforce.setZero();
   tmpf.setZero();
   tmpf=force.cwiseProduct(Gcons);
   Gforce.segment(beadsize,Nimages*beadsize) = tmpf;
   hGforce = -h*Gforce;
   //hGforce = -h*(Gforce.cwiseProduct(weight));
   ks.row(5) = hGforce;
   //----------------------------------------------------------
   sumks = ks.row(0)*37.0/378.0 + ks.row(2)*250.0/621.0 
           + ks.row(3)*125.0/594.0 + ks.row(5)*512.0/1771.0;

   x5 = path + sumks;


   sumks = ks.row(0)*2825.0/27648.0 + ks.row(2)*18575.0/48384.0 + 
           ks.row(3)*13525.0/55296.0 + ks.row(4)*277.0/14336.0 + 
           ks.row(5)*1.0/4.0;
   x4 = path + sumks;
   //----------------------------------------------------------

   //COMPUTE ERROR & ADJUST STEPSIZE
   error1=0;
   pathdiff = x5 - x4;

   //#pragma omp parallel for schedule(dynamic)
   for(int k=1; k<Nimages+1; k++){
      for (int i = 0; i < beadsize; i++){
          if( (abs(pathdiff(i+k*beadsize))) > error1 ){
              error1 = abs(pathdiff(i+k*beadsize));
          } 
      }
   }
   //#pragma omp barrier

   if(error1==0){
     hnew=2*h;
   }
   else{
     hnew = 0.9*h*pow(abs(error0/error1),0.2);
   }
    
   quad_app(x5,epath,H,g,energy,Nimages,beadsize,equad,gquad);

   funupwind(x5,gquad,equad,Nimages,beadsize,force,weight);

   funupwind(epath,g,energy,Nimages,beadsize,force1,weight);

   //CHECK TO SEE IF MOVED OVER VALLEY
   VectorXd localf(beadsize);
   VectorXd localf1(beadsize);

   //#pragma omp parallel for schedule(dynamic)
   for(int i=0; i<Nimages; i++){
          localf = force.segment(beadsize*i,beadsize);
          localf1 = force1.segment(beadsize*i,beadsize);

	  currswitch[i] = currswitch[i] and (localf.dot(localf1)>=0);
   }   
   //#pragma omp barrier

   VectorXd localks(beadsize);
   VectorXd ksrowvec((Nimages+2)*beadsize);
   //#pragma omp parallel for schedule(dynamic)
   for(int j=1; j<Nimages+1; j++){
          for(int i=0; i<6; i++){
             localf1 = force1.segment(beadsize*(j-1),beadsize);
             ksrowvec = ks.row(i);//all images in row i
             localks = ksrowvec.segment(beadsize*j,beadsize);
             double dotksf;
             dotksf = localks.dot(localf1);
             //currswitch[i] = currswitch[i] and (dotksf<=0);
             currswitch[j] = currswitch[j] and (dotksf<=0);
	  }      
   }   
   //#pragma omp barrier

}
//END:: CASH-KARP 4/5 ORDER
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//START
//QUADRATIC APPROXIMATION TO ENERGY AND GRADIENT
void quad_app(VectorXd& wholepath, VectorXd& oldpath, MatrixXd& Hessmat, VectorXd& force,\
                VectorXd& energy, int Nimages, int beadsize, VectorXd& equad, VectorXd& gquad){

   VectorXd dx((Nimages+2)*beadsize);
   MatrixXd Hc(beadsize,beadsize);//Hessian of current image
   VectorXd dxc(beadsize);//difference between current and oldpath for the current image
   VectorXd Hdx(beadsize);//Hessian*dxci
   VectorXd Fc(beadsize);//force of current image
   VectorXd gquadc(beadsize);//force of current image

   double Fcdx;//dot product of hessian and path difference of the current image (Fc.dxc) 
   int srow=0;//starting row of current image
   int rsize=beadsize;//number of rows
   int csize=beadsize;//number of columns

   dx = wholepath - oldpath;
   //#pragma omp parallel for schedule(dynamic)
   for(int k=1; k<Nimages+1; k++){
      //coord of current image      
      dxc = dx.segment(beadsize*k,beadsize);

      //force of current image
      Fc = force.segment(beadsize*(k-1),beadsize);
      //get hessian of the current image
      srow=k*beadsize;//starting row of current image
      Hc=Hessmat.block(srow,0,rsize,csize);

      Hdx = Hc*dxc;//matmul(Hnew,hdx)
      Fcdx = Fc.dot(dxc);//dot_product(glast(i,:),dx(i,:))
      gquadc = Fc + Hdx;//gquad of the current image
      gquad.segment(beadsize*(k-1),beadsize) = gquadc;

      equad(k-1) = energy(k) + Fcdx + 0.5*(dxc.dot(Hdx));

   }
   //#pragma omp barrier
   
}
//END

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//START:: UPDATE FORCES
void funupwind(VectorXd& path, VectorXd& g, VectorXd& energy, int Nimages, \
                 int beadsize, VectorXd& gtan, VectorXd& weight){

   VectorXd tangents = VectorXd::Zero((Nimages)*beadsize);
 
   VectorXd g_curr(beadsize);//current image gradiend
   VectorXd t_curr(beadsize);//current image tangent
   VectorXd gt_curr(beadsize);// g_curr*t_curr
   VectorXd gtan_curr(beadsize);//current image gtan array
   double gt_sum; //sum(gt_curr)

   calc_tangents(path,Nimages,beadsize,tangents,energy,weight,0);

   //#pragma omp parallel for schedule(dynamic)
   for(int k=0; k<Nimages; k++){
         //get values from global arrays for the current image
         //start from beadsize*(k-1)
         //and get beadsize elements
         g_curr = g.segment(beadsize*k,beadsize);
         t_curr = tangents.segment(beadsize*k,beadsize);

         gt_curr=g_curr.cwiseProduct(t_curr);
         gt_sum=gt_curr.sum();
         gtan_curr = g_curr - gt_sum*t_curr;

         //update global array 
         gtan.segment(beadsize*k,beadsize) = gtan_curr;
   
    } 
    //#pragma omp barrier


}
//END:: UPDATE FORCES

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//START:: CALCULATE TANGENTS
void calc_tangents(VectorXd& path, int Nimages, int beadsize, VectorXd& tangents,\
                     VectorXd& energy, VectorXd& weight,int choice){

   VectorXd path1((Nimages+2)*beadsize);
   VectorXd w((Nimages+2)*beadsize);
   VectorXd currtan(beadsize);
   VectorXd pathdiff(beadsize);
   VectorXd tp(beadsize);
   VectorXd tm(beadsize);

   VectorXd V = VectorXd::Ones(Nimages+2);
   double Vmax;
   double Vmin;

   path1= path.cwiseProduct(weight);//wholepath*weighs;
   //path1 = path;

   if(choice==1){
      //#pragma omp parallel for schedule(dynamic)
      for(int k=1; k<Nimages+1; k++){
         double tsss = 0.0;
         //get tangents for current image   
         pathdiff = path1.segment(beadsize*(k+1),beadsize) 
                  - path1.segment(beadsize*(k-1),beadsize);
         currtan = pathdiff;

         //tsss = currtan.norm();//sqrt((currtan.array().square()).sum());
         tsss = sqrt((currtan.array().square()).sum());
         //compute the tangent values in global array
         tangents.segment(beadsize*(k-1),beadsize) = currtan/tsss;
      }
      //#pragma omp barrier
   }
   else{
       V(0)= -1*hugeNum;  //-1*V;
       V(Nimages+1) = -1*hugeNum;
       //#pragma omp parallel for schedule(dynamic)
       for(int k=1; k<Nimages+1; k++){
          V(k)=energy(k-1);
       }
       //#pragma omp barrier
 
       //#pragma omp parallel for schedule(dynamic)
       for(int k=1; k<Nimages+1; k++){
         
          Vmax=max( abs(V(k+1)-V(k)), abs(V(k-1)-V(k)) );
          Vmin=min( abs(V(k+1)-V(k)), abs(V(k-1)-V(k)) );

          tp = path1.segment(beadsize*(k+1),beadsize)
             - path1.segment(beadsize*k,beadsize);

          tm = path1.segment(beadsize*k,beadsize)
             - path1.segment(beadsize*(k-1),beadsize);

          //compute for the current image
          if( (V(k+1)>=V(k)) and (V(k)>=V(k-1)) ){
              currtan=tp;
          }
          else if( (V(k+1)<=V(k)) and (V(k)<=V(k-1)) ){
              currtan=tm;
          }
          else if( V(k+1)>V(k-1) ){
              currtan = tp*Vmax + tm*Vmin;
          }
          else{
              currtan = tp*Vmin + tm*Vmax;
          }
          
          //currtan = currtan/currtan.norm();//(sqrt((currtan.array().square()).sum()));
          currtan = currtan/(sqrt((currtan.array().square()).sum()));
          //update whole tangents with currtan (added new)
          tangents.segment(beadsize*(k-1),beadsize) = currtan;

       }
       //#pragma omp barrier

   }


}
//END:: CALCULATE TANGENTS


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//START
//CALCULATE TRUST RADIUS FOR EACH IMAGE
void updateTR (MatrixXd& Hessmat, VectorXd& glast, VectorXd& oldpath, VectorXd& wholepath,\
                 VectorXd& Eqmmm_images, VectorXd& lastenergy, VectorXd& trs, VectorXd& maxtr,\
                 int beadsize,int Nimages){

     MatrixXd Hc(beadsize,beadsize); //hessian of the current image
     VectorXd dx((Nimages+2)*beadsize);
     VectorXd dxc(beadsize);//difference between current and oldpath for the current image
     VectorXd gc(beadsize);//gradient of the current image
     double Ec, Ep; //current and previous energy of the corresponding image
     double Ediff; //change in energy
     double newTR;
     double rho;
 
     dx = wholepath - oldpath;
   
     //#pragma omp parallel for schedule(dynamic)
     for(int k=1; k<Nimages+1; k++){
        //collect energy values of previous and current runs
        Ec = Eqmmm_images(k);
        Ep = lastenergy(k);
        Ediff = Ec-Ep;
        //path difference from previous run for current image
        dxc = dx.segment(beadsize*k,beadsize);
        //current image gradient
        gc = glast.segment(beadsize*(k-1),beadsize);
        //get hessian of the current image
             int srow=k*beadsize;//starting row of current image
             int rsize=beadsize;//number of rows
             int csize=beadsize;//number of columns
             Hc=Hessmat.block(srow,0,rsize,csize);
             
        if( (abs(Ediff)>1e-10) and (dxc.norm() > 1e-7))
        {
           rho=Ediff/( 0.5*dxc.dot((Hc*dxc)) + dxc.dot(gc) );
           newTR=trs(k-1);
           if((rho<0.3) or (rho>3)){
	      	newTR=0.25*dxc.norm();
           }
           //else if((rho>0.9) and (rho<1.5) and (trs(k-1)<= ((dxc.norm())*1.4))){
           else if((rho>0.9) and (rho<1.5) and (trs(k-1)<= (sqrt((dxc.array().square()).sum())*1.4))){
                newTR = min(2*trs(k-1),maxtr(k-1));
           }
           else if ((rho < 0.75) or (rho> 1.5)){
	 	  newTR = min(dxc.norm(),trs(k-1));
           }             
           //else if ((rho < 0.85) and ( (dxc.norm()/trs(k-1))<0.8)){  
           else if ((rho < 0.85) and ((sqrt((dxc.array().square()).sum())/trs(k-1))<0.8)){ 
                newTR = trs(k-1)/2;
           }
           //else if ((rho < 0.95) and ( (dxc.norm()/trs(k-1)) < 0.5)){ 
           else if ((rho < 0.95) and ((sqrt((dxc.array().square()).sum())/trs(k-1))<0.5)){
                newTR = trs(k-1)/2;
           }
        }
        else{
            rho=1.0;
	    newTR = trs(k-1);	
        }
        trs(k-1) = newTR;    
     }
     //#pragma omp barrier


}
//END

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void updatepath(VectorXd& wholepath, vector<QMMMAtom>& QMMMData,
                  QMMMSettings& QMMMOpts,int beadsize, int Natoms,
                  bool struct_to_path){

    int counter;
   
    if(struct_to_path){//path=struct

      //#pragma omp parallel for schedule(dynamic)
      for (int p=0;p<QMMMOpts.NBeads;p++){ 
          counter = 0;
          for (int i=0;i<Natoms;i++){

            if (QMMMData[i].QMRegion or QMMMData[i].PBRegion){
               wholepath[(p*beadsize)+(counter*3)] = QMMMData[i].P[p].x;
               wholepath[(p*beadsize)+(counter*3)+1] = QMMMData[i].P[p].y;
               wholepath[(p*beadsize)+(counter*3)+2] = QMMMData[i].P[p].z;
               counter = counter+1;
            }
          }
      }
      //#pragma omp barrier
      
    }
    else{//struct=path
    
      //#pragma omp parallel for schedule(dynamic)
      for (int p=1;p<(QMMMOpts.NBeads-1);p++){ 
         counter = 0;
         for (int i=0;i<Natoms;i++){     
           if (QMMMData[i].QMRegion or QMMMData[i].PBRegion){
             QMMMData[i].P[p].x = wholepath[(p*beadsize)+(counter*3)];
             QMMMData[i].P[p].y = wholepath[(p*beadsize)+(counter*3)+1];
             QMMMData[i].P[p].z = wholepath[(p*beadsize)+(counter*3)+2];
             counter = counter+1;
           }
         }
         
      }
      //#pragma omp barrier

    }

}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//START:: ALIGNED?
void alignment(VectorXd& path,int Nimages,int beadsize,VectorXd& aligned){
     VectorXd a(beadsize);
     VectorXd b(beadsize);

     //#pragma omp parallel for schedule(dynamic)
     for(int k=1; k<Nimages+1; k++){
        for(int i=0; i<beadsize; i++){
           a(i) = path(i+k*beadsize)-path(i+(k-1)*beadsize);
           b(i) = path(i+(k+1)*beadsize)-path(i+k*beadsize);
        }
        a.normalize();
        b.normalize();
        aligned(k-1)=a.dot(b);
     }
     //#pragma omp barrier

}
//END: ALIGNED?


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



//START 
//CONSTRUCT PATH BETWEEN REACTANT AND PRODUCT USING LINEAR INTERPOLATION
void construct_path(VectorXd& wholepath,MatrixXd& Geom1, MatrixXd& Geom2,\
                      int beadsize, int Nimages, int Natoms,\
                      vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts)
{

   VectorXd reactvec(beadsize);
   VectorXd prodvec(beadsize);
      
   int reactant=0;
   int product=Nimages+1;

   //#pragma omp parallel for schedule(dynamic)
   for(int k=0; k<Nimages+2; k++){
      for(int i=0; i<Natoms; i++) {
         for (int j = 0; j < 3; j++){
           if(k==reactant){
              wholepath(k*beadsize+3*i+j) = Geom1(i,j);
              reactvec(3*i+j) = Geom1(i,j);
           }
           else if(k==product){
             wholepath(k*beadsize+3*i+j) = Geom2(i,j);
             //prodvec(3*i+j) = Geom2(i,j);
 
             // to make sure that interpolation
             // is performed only for QM and PB atoms
             // and not MM atoms
             if(QMMMData[i].QMRegion or QMMMData[i].PBRegion){          
                prodvec(3*i+j) = Geom2(i,j);
             }
             else{
                //if MM atom set coords to react for the time being
                //so that the difference (react-prod) will be 0
                //this will ensure that during iterpolation
                //MM and PB atoms of all images btw react and prod
                //will have same coords with MM atoms of react 
                prodvec(3*i+j) = Geom1(i,j);
             }

           }
         }
      }
   }
   //#pragma omp barrier

   linearly_interpolate(wholepath,reactvec,prodvec,beadsize,Nimages);

}
//END 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//START 
//LINEAR INTERPOLATION
void linearly_interpolate(VectorXd& wholepath, VectorXd& reactvec, VectorXd& prodvec, \
                            int beadsize, int Nimages){

  VectorXd RPdiff(beadsize);//reactant product difference
  RPdiff= (prodvec-reactvec)/(Nimages+1.0);

  //#pragma omp parallel for schedule(dynamic)
  for(int k=1; k<(Nimages+1); k++){
     for(int i=0; i<beadsize; i++){
        wholepath(beadsize*k+i) = reactvec(i)+k*RPdiff(i);
     }
  }
  //#pragma omp barrier
  
}
//END


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//END: Hatice GOKCAN
