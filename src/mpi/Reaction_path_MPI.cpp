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
  # Function for reaction path optimizations in parallel                        #
  # Includes:                                                                   #
  #                                                                             #
  #       QSM optimization        : void LICHEMQSM                              #
  #                                                                             #
  ###############################################################################
*/
//START: QSM optimization
void LICHEMQSMMPI(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts, \
               VectorXd& wholepath, int Nimages, int QMdim, bool &QMDone, 
               VectorXd& force, double spaceout_dist, VectorXd& Eqmmm_images,
               int macroiter,fstream& logFile)
{

  fstream ofile,qmfile,initfile,pathfile,hessfile,gfile,p0file; //Generic file stream
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(logFile); //Save settings
  string dummy; //Generic string

  /* for MPI */
  int Worldrank, Worldsize,root;
  bool master;
  MPI_Comm_rank(MPI_COMM_WORLD, &Worldrank);
  MPI_Comm_size(MPI_COMM_WORLD, &Worldsize);
  root=0;
  master=false;
  if(Worldrank==0){
    master=true;
  }
  int stepct = 0; //Counter for optimization steps
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Initialize trajectory file
  call.str("");
  call << "LICHMNEBOpt.xyz";
  qmfile.open(call.str().c_str(),ios_base::out);
  //Set end points for the optimization
  int PathStart = 0;
  int PathEnd = QMMMOpts.NBeads;
  if (QMMMOpts.frznEnds)
  {
    //Change the start and end points
    PathStart = 1;
    PathEnd = QMMMOpts.NBeads-1;
  }
  //Initialize optimization variables
  double Eold = 0; //Previous total energy
  double SumE = 0; //Current total energy
  double dftol;
  double ftol; 
  ftol = 0.0001;//QMMMOpts.QMOptTol; //??? 20* or 10*  
  double dftot;
  //to check if converged
  //double spaceout_dist;     
  double RMSdiff = 0;
  double RMSforce = 0;
  double MAXforce = 0;
  double eqcons;
  double max_dfval;

  int counter=0;
  int beadsize = 3*(Nqm+Npseudo);
  //number of elements with coords of all atoms along whole path
  int wholesize = QMMMOpts.NBeads*beadsize; 
  int index=0;
  int srow;//starting row of current image
  int rsize;//number of rows
  int csize;//number of columns
  int qsmiter = 0;//for QM part
  //int maxiter = 50;//should be 50
  int maxiter = QMMMOpts.MaxQMSteps;

  bool first_time = true;
  bool PathDoneQM = 0;

  //alignment of path
  VectorXd aligned(Nimages);

  bool before_qsm = false;
  bool struct_to_path;

  bool nebatoms[QMdim];
  for(int i=0;i<QMdim;i++){
     nebatoms[i]=false;
  }
  VectorXd weight(wholesize);
  VectorXd weighted_path(wholesize);
  // path difference between two image
  VectorXd imagediff(beadsize); 
  // force difference between two image
  VectorXd forcediff(beadsize); 
  // matrix that contains all hessians
  MatrixXd Hessmat(beadsize*(Nimages+2),beadsize);  
  //hessian of the current image
  MatrixXd tmpH(beadsize,beadsize); 
  
  // ENERGIES
  VectorXd E_images(Nimages+2);
  VectorXd Emm_images(Nimages+2);
  VectorXd Eqm_images(Nimages+2);
  //VectorXd Eqmmm_images(Nimages+2);
  E_images.setZero(); 
  Emm_images.setZero();
  Eqm_images.setZero();
  //Eqmmm_images.setZero();    

  //quadratic approximation to energy and gradient
  VectorXd equad(Nimages);
  VectorXd gquad(Nimages*beadsize);

  //previous path, gradient, and energy
  VectorXd forcefrz(Nimages*beadsize);
  VectorXd energy(Nimages);
  VectorXd prevE(Nimages);
  VectorXd glast(Nimages*beadsize);
    
  //Forces
  //getting it in LICHEMQSM from LICHEM.cpp
  VectorXd Forces(Ndof); //Local forces
  //Create array to store stats and check convergence
  MatrixXd ForceStats(QMMMOpts.NBeads,2);
  ForceStats.setZero();

  VectorXd oldpath(wholesize);
  VectorXd lastenergy(Nimages+2);
  VectorXd gtan(Nimages*beadsize);
  VectorXd gtan_curr(beadsize);

  //trust radius
  VectorXd trs(Nimages);
  trs.setConstant(0.03);

  //max trust radius
  VectorXd maxtr(Nimages);
  maxtr.setConstant(0.15);
  
  VectorXd cons(Nimages);
  cons = VectorXd::Ones(Nimages);
  VectorXd dfvals(Nimages);
 

  vector<QMMMAtom> OldQMMMData;

  weight.setZero(); 
  if(master){ 
    for(int i=0;i<Natoms;i++){
       if (QMMMData[i].QMRegion or QMMMData[i].PBRegion){ //???
          if(QMMMData[i].NEBActive){
              nebatoms[index]=true;
              for(int k=0; k<QMMMOpts.NBeads; k++){ 
                   weight(k*beadsize + index*3) = 1.0; //x
                   weight(k*beadsize + index*3 + 1) = 1.0;//y
                   weight(k*beadsize + index*3 + 2) = 1.0;//z
              }
          }
          index=index+1;  
       }
    }
  }
  VectorXd reactCoord(QMMMOpts.NBeads); //Reaction coordinate
  reactCoord.setZero();

  struct_to_path = true;
  updatepath(wholepath,QMMMData,QMMMOpts,
             beadsize,Natoms,struct_to_path);

  if(QMMMOpts.debug){
    if(Worldrank==0){

       call.str("");
       call << "rm -rf debug_" << macroiter;
       globalSys = system(call.str().c_str());

       call.str("");
       call << "mkdir debug_" << macroiter;
       globalSys = system(call.str().c_str());

       call.str("");
       call << "p0file_" << macroiter << ".txt";
       p0file.open(call.str().c_str(),ios_base::out);
       p0file << setprecision(17) << wholepath << endl;
       p0file.flush();
       p0file.close();

    }
  }

  while ( (!PathDoneQM) and (qsmiter < maxiter)){

      if(master){
        logFile << '\n';    
        logFile << "            ";
        logFile << "| QM step : ";
        logFile << qsmiter+1;
        logFile << '\n';
        logFile.flush(); //Print progress
        
        /*calculate reaction coordinate*/
        calc_react_coord(QMMMOpts, QMMMData,reactCoord);
        if(QMMMOpts.debug){
           call.str("");
           call << "Initialpath_" << macroiter << "_" << qsmiter << ".txt";
           initfile.open(call.str().c_str(),ios_base::out);
           
           call.str("");
           call << "Wholepath_" << macroiter << "_" << qsmiter << ".txt";
           pathfile.open(call.str().c_str(),ios_base::out);
           
           call.str("");
           call << "Hessian_" << macroiter << "_" << qsmiter << ".txt";
           hessfile.open(call.str().c_str(),ios_base::out);
           
           call.str("");
           call << "grad_" << macroiter << "_" << qsmiter << ".txt";
           gfile.open(call.str().c_str(),ios_base::out);
        }

      }//end if master

      //qsm iter start from 1
      if(macroiter==1 and qsmiter==0)
      //do not compute forces since 
      //it is already computed before
      //starting optimization
      {

        if(master){
           logFile << '\n';
           logFile << "               ";
           logFile << "Using forces from initial step."<< endl;
           //print bead energies
           print_progress(QMMMOpts, 0,Eqmmm_images,
                          RMSdiff, MAXforce, RMSforce,reactCoord,logFile);

           //print TS React and Prod and barriers
           getTSbead(QMMMOpts,Eqmmm_images);
           print_progress(QMMMOpts, 1,Eqmmm_images,
                          RMSdiff, MAXforce, RMSforce,reactCoord,logFile);

        }//end if master

      }  
      else
      {
        if(master){
           logFile << '\n';
           logFile << "               ";
           logFile << "Computing forces. " << endl;

           //Run optimization
           logFile << '\n';
           logFile.flush(); //Print progress
         }//end if master

         CalcForcesMPI(QMMMData,QMMMOpts,Eqm_images, Emm_images,
                       Eqmmm_images,force,beadsize,QMdim,false,logFile);

         if(master){
           //print bead energies
           print_progress(QMMMOpts, 0,Eqmmm_images,
                          RMSdiff, MAXforce, RMSforce,reactCoord,logFile);
         
           //print TS React and Prod and barriers
           getTSbead(QMMMOpts,Eqmmm_images);
           print_progress(QMMMOpts, 1,Eqmmm_images,
                          RMSdiff, MAXforce, RMSforce,reactCoord,logFile);

           //force *= -1;
           //convert to gradient. 
           //only beads btw react and prod
           //since react and prod was send when we entered
           for (int k = 1; k < Nimages+1; k++)
           {
             force.segment(beadsize*k,beadsize)=-1*force.segment(beadsize*k,beadsize);
           }
           

           /*if(QMMMOpts.KeepFiles and ((qsmiter%QMMMOpts.perQM)=0)){
             save_files(0,qsmiter,logFile);
           }*/
         }/*end if master*/

      }/* end if(macroiter==1 and qsmiter==0)*/

      if(master){ 
        //#pragma omp parallel for
        for (int k = 1; k < Nimages+1; k++){
           for (int i = 0; i < beadsize; i++){
            //forces involves all beads
            //get only beads btw react and prod
             forcefrz(i+(k-1)*beadsize) = force(i+k*beadsize);
           }
           //Eqmmm_images involves all beads
           // get only beads btw react and prod
           energy(k-1) = Eqmmm_images(k);
        }
        //#pragma omp barrier
      }/* end if master */

      //start:: if first time
      if(first_time)
      {
         if(master){ 
           //initialize hessians as identity matrices
           init_Hess(Hessmat,beadsize,Nimages);
           //Hessmat = forcefrz.squaredNorm()*Hessmat; 
           Hessmat = forcefrz.norm()*Hessmat;
           
           for (int k = 1; k < Nimages+1; k++){
           
               //update hessian using current and previous images
               imagediff = wholepath.segment(beadsize*k,beadsize)
                         - wholepath.segment(beadsize*(k-1),beadsize);
               forcediff = force.segment(beadsize*k,beadsize)
                         - force.segment(beadsize*(k-1),beadsize);
               //get hessian of the current image
               srow=k*beadsize;//starting row of current image
               rsize=beadsize;//number of rows
               csize=beadsize;//number of columns
               tmpH=Hessmat.block(srow,0,rsize,csize);    
           
               //Use DBFGS algorithm to update current hessian   
               DBFGS(tmpH,imagediff,forcediff,QMdim*3);
           
               //update hessian using current and next images
               imagediff = wholepath.segment(beadsize*k,beadsize)
                         - wholepath.segment(beadsize*(k+1),beadsize);
               forcediff = force.segment(beadsize*k,beadsize)
                         - force.segment(beadsize*(k+1),beadsize);
           
               //Use DBFGS algorithm to update current hessian
               DBFGS(tmpH,imagediff,forcediff,QMdim*3);
           
               //update the matrix (Hessmat) that contains all hessians 
               Hessmat.block(srow,0,rsize,csize)=tmpH;   
           
           }
         }/* end if master*/
         first_time = false;

      }//end:: first time

      else{

        if(master){

           //computes equad&gquad of images btw react and prod
           quad_app(wholepath,oldpath,Hessmat,forcefrz,Eqmmm_images,
                    Nimages,beadsize,equad,gquad);
           
           updateTR(Hessmat, glast, oldpath, wholepath, Eqmmm_images,
                    lastenergy, trs, maxtr, beadsize, Nimages);
           
           for (int k = 1; k < Nimages+1; k++){
           
               //update hessian using current and previous runs
               imagediff = wholepath.segment(beadsize*k,beadsize)
                         - oldpath.segment(beadsize*k,beadsize);
               forcediff = forcefrz.segment(beadsize*(k-1),beadsize)
                         - glast.segment(beadsize*(k-1),beadsize);
           
               //get hessian of the current image
               srow=k*beadsize;//starting row of current image
               rsize=beadsize;//number of rows
               csize=beadsize;//number of columns
               tmpH=Hessmat.block(srow,0,rsize,csize); 
           
               // Use DBFGS algorithm to update current hessian   
               DBFGS(tmpH,imagediff,forcediff,QMdim*3);
           
               //update the matrix (Hessmat) that contains all hessians 
               Hessmat.block(srow,0,rsize,csize)=tmpH;   
           
           }
           
           funupwind(wholepath,forcefrz,energy,Nimages,beadsize,gtan,weight);
           
           //#pragma omp parallel for
           //Frozen ends
           for(int k=0; k<Nimages; k++){
              gtan_curr=gtan.segment(beadsize*k,beadsize);
              //frozen ends: fill force stats for 
              //images between react and product
              MAXforce = abs(gtan_curr.maxCoeff());
              if (abs(gtan_curr.minCoeff()) > MAXforce)
              {
                 //Update max
                 MAXforce = abs(gtan_curr.minCoeff());
              } 
              //MAXforce = abs(Forces.maxCoeff());
              /*
              ForceStats(k+1,0) = MAXforce;
              ForceStats(k+1,1) = (gtan_curr).squaredNorm(); //RMS force           
              */
              ForceStats(k+1,0) = MAXforce*bohrRad;
              ForceStats(k+1,1) = (gtan_curr*bohrRad).squaredNorm(); //RMS force
 
              dfvals(k) =  gtan_curr.norm();
                          //sqrt((gtan_curr.array().square()).sum());
           }
           //#pragma omp barrier
           
           max_dfval=dfvals.maxCoeff();
           
           //Check convergence
           logFile << "\n";
           logFile << "               ";
           logFile << "Checking convergence:" << endl;
           logFile << "\n";

           PathDoneQM = QMConverged(QMMMData,OldQMMMData,ForceStats,
                                    qsmiter,QMMMOpts,
                                    RMSdiff, RMSforce, MAXforce);

           if(QMMMOpts.KeepFiles  and
              (((qsmiter%QMMMOpts.perQM)==0) or
              PathDoneQM or
              qsmiter==QMMMOpts.MaxQMSteps))
           {
             save_files(0,qsmiter,logFile);
           }

           //print conv. criterias and stats
           print_progress(QMMMOpts, 2, Eqmmm_images,
                          RMSdiff, MAXforce, RMSforce,reactCoord,logFile);

         }/* end if master */

         /* bcast PathDone to everyone*/
         MPI_Bcast(&PathDoneQM,1,MPI::BOOL,root,MPI_COMM_WORLD);//group_comm);

         MPI_Bcast(&max_dfval,1,MPI_DOUBLE,root,MPI_COMM_WORLD);//group_comm);

         /* wait, since there is an update and return below*/
         MPI_Barrier(MPI_COMM_WORLD);
         

         //if(PathDoneQM or (max_dfval<ftol)){
         if(PathDoneQM or (max_dfval*bohrRad<ftol)){  
            if(master){ 
              logFile << "\n";
              logFile << "\n";
              logFile << "               ";
              logFile << "Convergence achieved at step ";
              logFile << qsmiter+1;
              logFile << "\n";
              logFile << "               ";
              logFile << "Finishing QM Iterations..." << endl;
            
              struct_to_path = false;
              updatepath(wholepath,QMMMData,QMMMOpts,
                         beadsize,Natoms,struct_to_path);
              
              //Clean up files
              call.str("");
              call << "rm -f LICHMNEBOpt.xyz";
              globalSys = system(call.str().c_str());
              if(QMMMOpts.debug){              
                 for (int k = 1; k < Nimages+1; k++){
                      int srow=k*beadsize;/*starting row of current image*/
                      int rsize=beadsize;/*number of rows*/
                      int csize=beadsize;/*number of columns */
                      hessfile << setprecision(17) << Hessmat.block(srow,0,rsize,csize) << endl;
                 }
                 hessfile.flush();
                 
                 gfile << setprecision(17) << force << endl;
                 gfile.flush();
                 
                 initfile << setprecision(17) << wholepath << endl;
                 initfile.flush();

                 hessfile.close();
                 initfile.close();
                 pathfile.close();
                 gfile.close();
                 
                 call.str("");
                 call << "mv *.txt debug_" << macroiter << "/ ";
                 globalSys = system(call.str().c_str());

              }
 
            }//end if master
            QMDone = PathDoneQM;


            Send_qmmmdata(QMMMData,QMMMOpts.NBeads,0,master,Natoms);
            MPI_Barrier(MPI_COMM_WORLD);
            //Finish and return
            return;//break;
         }

      }//end::not first time

      //alignment(wholepath,Nimages,beadsize,aligned);
      //Integrate to TRs or until finished ###
      //Copy old structure and forces
      if(master){
         OldQMMMData = QMMMData; //Save structure
         oldpath=wholepath;
         glast=forcefrz;//force;
         lastenergy=Eqmmm_images;
         prevE=energy;
         if(QMMMOpts.debug){
            for (int k = 1; k < Nimages+1; k++){
                 int srow=k*beadsize;//starting row of current image
                 int rsize=beadsize;//number of rows
                 int csize=beadsize;//number of columns 
                 hessfile << setprecision(17) << Hessmat.block(srow,0,rsize,csize) << endl;
            }
            hessfile.flush();
            
            gfile << setprecision(17) << force << endl;
            gfile.flush();
         }
      }//end if master

 
      for(int iter=0; iter<4; iter++){
        
        if(master){ 


          if(iter>0){
             wholepath=oldpath;
          }
          if(QMMMOpts.debug){
             initfile << "ODE iter=" << iter  << endl;
             initfile << setprecision(17) << wholepath << endl;
             initfile.flush();
          }

          ODESolve(wholepath,Hessmat,forcefrz,Eqmmm_images,Nimages,
                   beadsize,trs,cons,wholesize,ftol,dftol,weight,logFile);

          if(QMMMOpts.debug){
             pathfile << "ODE iter=" << iter << endl;
             pathfile << setprecision(17) << wholepath << endl;
             pathfile.flush();
          }

        }//end if master

        MPI_Bcast(&dftol,1,MPI_DOUBLE,root,MPI_COMM_WORLD); 
        MPI_Barrier(MPI_COMM_WORLD);

        if(dftol < ftol/10){
           if(master){
             logFile << "                   ";
             logFile << "dftol < ftol/10." << endl;
           // wholepath is changed. Update the Struct
             struct_to_path = false;
             updatepath(wholepath,QMMMData,QMMMOpts,
                        beadsize,Natoms,struct_to_path);
           }//end if master

           //update QMMMData
           Send_qmmmdata(QMMMData,QMMMOpts.NBeads,0,master,Natoms);
           MPI_Barrier(MPI_COMM_WORLD);
           break;//return;
        }
        
      } 

      if(master){

        //Space out if necessary
        eqcons=0.0;
        
        eqconst(wholepath,Nimages,beadsize,eqcons);
        eqcons=eqcons/Nimages;
       
        /* Start: Aug 17 2018*/
        /* do not spaceout if it is the last qsmiter */
        if(eqcons > spaceout_dist){
        //if((eqcons > spaceout_dist) and ((qsmiter+1) < maxiter)){
        /* End: Aug 17 2018*/
 
          //if(QMMMOpts.debug){
              logFile << '\n' << '\n';
              logFile << "                ";         
              logFile << "spaceout distance = " << spaceout_dist << "\n";
              logFile << "                ";
              logFile << "eqconst = " << eqcons << "\n";
              logFile << "                ";
              logFile << "Spacing out... " << "\n" ;
              logFile << "\n";
          //}
          //computes equad&gquad of images btw react and prod
          quad_app(wholepath,oldpath,Hessmat,forcefrz,Eqmmm_images,
                   Nimages,beadsize,equad,gquad);
        
          funupwind(wholepath,gquad,equad,Nimages,beadsize,gtan,weight);
        
          dftot=0.0;
          
          //#pragma omp parallel for
          for (int k=0; k<Nimages;k++){
             for(int i=0;i<beadsize; i++){
                gtan_curr(i) = gtan(i+k*beadsize); 
             }
             dftot = max( dftot, gtan_curr.norm());
          }
          //#pragma omp barrier
        
          //spaceoutcubic     
          for (int k=0; k<3;k++){
              spaceoutcubic(wholepath,nebatoms,QMMMOpts.NBeads,beadsize,weight);
          }
          quad_app(wholepath,oldpath,Hessmat,forcefrz,Eqmmm_images,
                   Nimages,beadsize,equad,gquad);
        
          funupwind(wholepath,gquad,equad,Nimages,beadsize,gtan,weight);
          dftot=0.0;
        
          //#pragma omp parallel for
          for (int k=0; k<Nimages;k++){
             for(int i=0;i<beadsize; i++){
                gtan_curr(i) = gtan(i+k*beadsize);          
             }
             dftot = max( dftot, gtan_curr.norm()); 
          }
          //#pragma omp barrier
        
        }//end:: Space out if necessary
        logFile << '\n';
     
      }/*end if master*/       

      if(master){

        struct_to_path = false;
        updatepath(wholepath,QMMMData,QMMMOpts,
                   beadsize,Natoms,struct_to_path);

      }

      qsmiter += 1;
      if(QMMMOpts.debug){
         if(master){
          hessfile.close();
          initfile.close();
          pathfile.close();
          gfile.close();
         }
      }
      Send_qmmmdata(QMMMData,QMMMOpts.NBeads,0,master,Natoms);
      MPI_Barrier(MPI_COMM_WORLD);

      
  }//end while loop. qm part finished

  if(master){

    if(QMMMOpts.debug){
/*      call.str("");
      call << "mkdir debug_" << macroiter;
      globalSys = system(call.str().c_str());

      call.str("");
      call << "mv " << hessfile << " ";
      call << "debug_" << macroiter << "/ ";
      globalSys = system(call.str().c_str());

      call.str("");
      call << "mv " << pathfile << " ";
      call << "debug_" << macroiter << "/ ";
      globalSys = system(call.str().c_str());

      call.str("");
      call << "mv " << gfile << " ";
      call << "debug_" << macroiter << "/ ";
      globalSys = system(call.str().c_str());

      call.str("");
      call << "mv " << p0file << " ";
      call << "debug_" << macroiter << "/ ";
      globalSys = system(call.str().c_str());
*/
      call.str("");
      call << "mv *.txt debug_" << macroiter << "/ ";
      globalSys = system(call.str().c_str());
    }

   //Clean up files

      call.str("");
      call << "rm -f LICHMNEBOpt.xyz";
      globalSys = system(call.str().c_str());
   }
   //Finish and return
   return;
  
}
//END: QSM optimization
