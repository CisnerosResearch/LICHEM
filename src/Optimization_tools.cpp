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
  # Functions for reaction path optimizations                                   #
  # Includes:                                                                   #
  #                                                                             #
  #       Force calculation        : double CalcForces                          #
  #                                                                             #
  #       Energy calculation       : double CalcEnergy                          #
  #                                                                             #
  #       MM Optimization without                                               #
  #       restrains                : double runMMopt                            #
  #                                                                             #
  #       MM Optimization with                                                  #
  #       restrains                : double runRestrMMopt                       #
  #                                  (calls double TINKEROptRestr)              #
  #                                                                             #
  #       MM Optimization with                                                  #
  #       restrains using TINKER   : double TINKEROptRestr                      #
  #                                                                             #
  #       Frequency calculation    : void CalcFreq                              #
  #                                                                             #
  #                                                                             #
  #       Convergence Test for QSM : bool QSMConverged                          #
  #                                                                             #
  #       Convergence Test for QM  : bool QMConverged                           #
  ###############################################################################
*/

double CalcForces(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,\
VectorXd& Eqm_images, VectorXd& Emm_images,VectorXd& Eqmmm_images,\
VectorXd& force, int beadsize, int QMdim, bool first_time,
fstream& logFile)
{

      stringstream call;
      double SumE = 0; //Current total energy
      double relAU = 0.0;
      double relkcal = 0.0;
     
      //Forces
      VectorXd Forces(QMdim*3); //Local forces

      //Set end points for the optimization
      int PathStart = 0;
      int PathEnd = QMMMOpts.NBeads;
      if (QMMMOpts.frznEnds)
      {
        //Change the start and end points
        PathStart = 1;
        PathEnd = QMMMOpts.NBeads-1;
      }

      if(first_time){
         PathStart=0;
         PathEnd=QMMMOpts.NBeads;
      }

      int NewTS = 0; //Storage for new TS ID
      double NewTSEnergy = -1*hugeNum;

      for (int p=PathStart;p<PathEnd;p++)
      {
        double E = 0;
        double Emm = 0;
        double Eqm = 0;
        double EeV = 0;
        double tempFx=0.0;
        double tempFy=0.0;
        double tempFz=0.0;

        //Create blank force array
        Forces.setZero();
        //Calculate forces (QM part)
        if (Gaussian)
        {
          int tstart = (unsigned)time(0);
          Eqm += GaussianForces(QMMMData,Forces,QMMMOpts,p);
          Eqm = Eqm/har2eV;//to a.u.
          QMTime += (unsigned)time(0)-tstart;

        }
        if (PSI4)
        {
          int tstart = (unsigned)time(0);
          Eqm += PSI4Forces(QMMMData,Forces,QMMMOpts,p);
          Eqm = Eqm/har2eV;//to a.u.
          QMTime += (unsigned)time(0)-tstart;
          //Delete annoying useless files
          globalSys = system("rm -f psi.* timer.*");
        }
        if (NWChem)
        {
          int tstart = (unsigned)time(0);
          Eqm += NWChemForces(QMMMData,Forces,QMMMOpts,p);
          Eqm = Eqm/har2eV;//to a.u.
          QMTime += (unsigned)time(0)-tstart;
        }
        E += Eqm;

        Eqm_images[p]=Eqm;//E;

        //Calculate forces (MM part)
        if (TINKER)
        {
          int tstart = (unsigned)time(0);
          EeV += TINKERForces(QMMMData,Forces,QMMMOpts,p);
          if (AMOEBA)
          {
            EeV += TINKERPolForces(QMMMData,Forces,QMMMOpts,p,logFile);
          }
          //Energy is computed and returned in eV
          //convert it back to kcal
          //then convert it to a.u. since
          //qm part returns a.u. in qsm implementation
          E = E + (EeV/har2eV);//in a.u
          Emm += TINKEREnergy(QMMMData,QMMMOpts,p,logFile);
          Emm = Emm/har2eV;//in a.u
          MMTime += (unsigned)time(0)-tstart;

        }
        //if (AMBER)
        //{
        //  int tstart = (unsigned)time(0);
        //  EeV += AMBERForces(QMMMData,Forces,QMMMOpts,p);
        //  E = E + (EeV/har2eV);//in a.u
        //  Emm += AMBEREnergy(QMMMData,QMMMOpts,p);
        //  Emm = Emm/har2eV;//in a.u
        //  MMTime += (unsigned)time(0)-tstart;
        //}
        if (LAMMPS)
        {
          int tstart = (unsigned)time(0);
          EeV += LAMMPSForces(QMMMData,Forces,QMMMOpts,p);
          E = E + (EeV/har2eV); //in a.u
          Emm += LAMMPSEnergy(QMMMData,QMMMOpts,p);
          Emm = Emm/har2eV; //in a.u
          MMTime += (unsigned)time(0)-tstart;
        }
        //Update energies
        SumE += E;

        Emm_images[p]=Emm;
        
        Eqmmm_images[p]=Eqm+Emm;

        if (p == 0)
        {
          QMMMOpts.EReact = Eqm+Emm;
        }
        if (p == QMMMOpts.TSBead)
        {
          QMMMOpts.ETrans = Eqm+Emm;
        }
        if (p == (QMMMOpts.NBeads-1))
        {
          QMMMOpts.EProd = Eqm+Emm;
        }
        //get the forces of all images from local Forces
        //#pragma omp parallel for schedule(dynamic)
        Forces = Forces/har2eV;
        for(int i=0;i<QMdim;i++){   
           force[(p*beadsize)+(i*3)] = Forces(3*i);
           force[(p*beadsize)+(i*3)+1] = Forces(3*i+1); 
           force[(p*beadsize)+(i*3)+2] = Forces(3*i+2);
        }
        //#pragma omp barrier


      }//runOPT of images is finished
} 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double CalcEnergy(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,
                  VectorXd& Eqm_images, VectorXd& Emm_images,
                  VectorXd& Eqmmm_images,fstream& logFile)
{

    double SumE,SumEeV,SumEau, Eqm, Emm;
    
    QMMMOpts.ETrans = -1*hugeNum; //Locate the initial transition state
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      SumE = 0; //Clear old energies
      SumEeV = 0; //Clear old energies
      Eqm = 0;
      Emm = 0;
      //Calculate QM energy
      if (Gaussian)
      {
        int tstart = (unsigned)time(0);
        Eqm += GaussianEnergy(QMMMData,QMMMOpts,p);
        //turn it to a.u.
        Eqm = Eqm/har2eV;
        QMTime += (unsigned)time(0)-tstart;
      }
      if (PSI4)
      {
        int tstart = (unsigned)time(0);
        Eqm += PSI4Energy(QMMMData,QMMMOpts,p);
        //turn it to a.u.
        Eqm = Eqm/har2eV;
        QMTime += (unsigned)time(0)-tstart;
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tstart = (unsigned)time(0);
        Eqm += NWChemEnergy(QMMMData,QMMMOpts,p);
        //turn it to a.u.
        Eqm = Eqm/har2eV;
        QMTime += (unsigned)time(0)-tstart;
      }
    
      Eqm_images[p]=Eqm;
 
      //Calculate MM energy
      if (TINKER)
      {
        int tstart = (unsigned)time(0);
        SumEeV += TINKEREnergy(QMMMData,QMMMOpts,p,logFile);
        Emm = SumEeV/har2eV; //turn it to a.u.
        MMTime += (unsigned)time(0)-tstart;
      }
      //if (AMBER)
      //{
      //  int tstart = (unsigned)time(0);
      //  SumEeV += AMBEREnergy(QMMMData,QMMMOpts,p);
      //  SumEau = SumEeV/har2eV; //turn it to a.u.
      //  SumE = SumE + SumEau;
      //  MMTime += (unsigned)time(0)-tstart;
      //}
      if (LAMMPS)
      {
        int tstart = (unsigned)time(0);
        SumEeV += LAMMPSEnergy(QMMMData,QMMMOpts,p);
        Emm = SumEeV/har2eV; //turn it to a.u.
        MMTime += (unsigned)time(0)-tstart;
      }
      Emm_images[p]= Emm;
      SumE = Eqm + Emm;
      Eqmmm_images[p]= SumE;

      if (p == 0)
      {
        //Save reactant energy
        QMMMOpts.EReact = SumE;
      }
      else if (p == (QMMMOpts.NBeads-1))
      {
        //Save product energy
        QMMMOpts.EProd = SumE;
      }


      //Update transition state
      if (SumE > QMMMOpts.ETrans)
      {
        //Save new properties
        QMMMOpts.TSBead = p;
        QMMMOpts.ETrans = SumE;
      }
      //Copy checkpoint data to speed up first step
      stringstream call; //Stream for system calls and reading/writing files
      if (p != (QMMMOpts.NBeads-1))
      {
        if (Gaussian)
        {
          call.str("");
          call << "cp LICHM_" << (p);
          call << ".chk";
          call << " LICHM_" << (p+1);
          call << ".chk";
          globalSys = system(call.str().c_str());
        }
        if (PSI4)
        {
          call.str("");
          //call << "cp LICHM_" << (p);//hatice
          //call << ".32";//hatice
          //call << " LICHM_" << (p+1);//hatice
          //call << ".32; ";//hatice
          call << "cp LICHM_" << (p);
          call << ".180";
          call << " LICHM_" << (p+1);
          call << ".180";
          globalSys = system(call.str().c_str());
        }
      }
    }
    //==============================================================

}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double runMMopt(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,fstream& logFile)
{
       double SumE,SumEeV;
       //Set end points for the optimization
       int PathStart = 0;
       int PathEnd = QMMMOpts.NBeads;
       if (QMMMOpts.frznEnds)
       {
         //Change the start and end points
         PathStart = 1;
         PathEnd = QMMMOpts.NBeads-1;
       }
       logFile << "             ";
       logFile << "Starting MM optimization ";
       logFile << "without restraints."<< endl;
       //Run MM optimization
       for (int p=PathStart;p<PathEnd;p++)
       {
          SumEeV = 0; //Clear old energies
          SumEeV = 0;
          if (TINKER)
          {
            int tstart = (unsigned)time(0);
            int mystat=0;
            SumEeV = TINKEROpt(QMMMData,QMMMOpts,p,logFile,mystat);
            if(mystat!=0){
               exit(0);
            }
            SumE = SumEeV/har2eV; //turn it to a.u.
            MMTime += (unsigned)time(0)-tstart;
          }
          //if (AMBER)
          //{
          //  int tstart = (unsigned)time(0);
          //  SumEeV = AMBEROpt(QMMMData,QMMMOpts,p);
          //  SumE = SumEeV/har2eV; //turn it to a.u.
          //  MMTime += (unsigned)time(0)-tstart;
          //}
          if (LAMMPS)
          {
            int tstart = (unsigned)time(0);
            SumEeV = LAMMPSOpt(QMMMData,QMMMOpts,p);
            SumE = SumEeV/har2eV; //turn it to a.u.
            MMTime += (unsigned)time(0)-tstart;
          }
        }
        logFile << '\n';
        logFile << "             ";
        logFile << "MM optimization complete.";
        logFile << '\n';
        logFile.flush();
        logFile << '\n';


}
//---------------------------------------------------------------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double runRestrMMopt(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,double restr,fstream& logFile)
{
 
    logFile << '\n';
    logFile << "             ";
    logFile << "Starting MM optimization ";
    logFile << "with restraint positions.\n";
    logFile << "                  ";
    logFile << "Force constant = ";
    logFile << restr << endl;

    int PathStart = 0;
    int PathEnd = QMMMOpts.NBeads;
    if (QMMMOpts.frznEnds)
    {
       //Change the start and end points
       PathStart = 1;
       PathEnd = QMMMOpts.NBeads-1;
    }
    for(int p=PathStart;p<PathEnd;p++)
    {
       double SumE=0.0;
       int tstart = (unsigned)time(0);
       int mystat=0;
       SumE = TINKEROptRestr(QMMMData,QMMMOpts,p,restr,logFile,mystat);
       if(mystat!=0){
          exit(0);
       }
       MMTime += (unsigned)time(0)-tstart;
    }
    logFile << '\n';
    logFile << "             "; 
    logFile << "MM optimizations with ";
    logFile << "restraints are complete. \n";
    logFile << "             ";
    logFile << "Starting new QSM iteration. \n" << endl;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double TINKEROptRestr(vector<QMMMAtom>& QMMMData,
                      QMMMSettings& QMMMOpts, 
                      int Bead, double restr,fstream& logFile,int& mystat)
{
  //Runs TINKER MM optimization
  fstream ofile,ifile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  double E = 0;
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  call.str("");
  call << "cp tinker.key LICHM_";
  call << Bead << ".key";

  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << Bead << ".key";
  ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
  ofile << '\n';
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';

  if (QMMMOpts.useMMCut)
  {
    /*Apply cutoff*/
    if (QMMMOpts.useEwald and PBCon)
    {
      /*Use Ewald and truncate vdW forces*/
      ofile << "cutoff " << LICHEMFormFloat(QMMMOpts.MMOptCut,12);
      ofile << '\n';
      ofile << "ewald" << '\n';
    }
    else
    {
      /*Use smoothing functions*/
      ofile << "cutoff " << LICHEMFormFloat(QMMMOpts.MMOptCut,12);
      ofile << '\n';
      ofile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.MMOptCut,12);
      ofile << '\n';
    }
  }

  else if (QMMMOpts.useLREC)
  {
    /*Apply cutoff*/
    if (QMMMOpts.useEwald and PBCon)
    {
      /*Use Ewald or PME*/
      ofile << "ewald" << '\n';
    }
    else if (!QMMMOpts.useImpSolv)
    {
      /*Use smoothing functions*/
      ofile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      ofile << '\n';
      ofile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      ofile << '\n';
    }
  }

  if (QMMMOpts.useImpSolv)
  {
    /*Add the implicit solvation model*/
    ofile << "solvate " << QMMMOpts.solvModel;
    ofile << '\n';
  }
  ofile << "openmp-threads " << Ncpus << '\n';
  ofile << "digits 12" << '\n';

  if (PBCon)
  {
    /*PBC defined twice for safety*/
    ofile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    ofile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    ofile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    ofile << "alpha 90.0" << '\n';
    ofile << "beta 90.0" << '\n';
    ofile << "gamma 90.0" << '\n';
  }

  ct = 0; //Generic counter
  if (QMMM or (Nfreeze > 0))
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if (QMMMData[i].MMRegion or QMMMData[i].BARegion)
      {
       
       /*Start:Hatice GOKCAN*/
       /*if( (Nmm + Nbound) != Nfreeze ){
           logFile << "             ";
           logFile << "Error:\n";
           logFile << "             ";
           logFile << "Total number of MM atoms and boundary atoms\n";
           logFile << "             ";
           logFile << "are not equal to the number of frozen atoms.\n";
           mystat=1;
           return 0;
       }*/
       /*End: Hatice GOKCAN*/

        if (!QMMMData[i].frozen)
        {
          if (ct == 0)
          {
            //Start a new active line
            ofile << "active ";
          }
          else
          {
            //Place a space to separate values
            ofile << " ";
          }
          ofile << (QMMMData[i].id+1);
          ct += 1;
          if (ct == 10)
          {
            //terminate an active line
            ct = 0;
            ofile << '\n';
          }
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      ofile << '\n';
    }
  }
  if (CHRG)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMRegion)
      {
        //New charges are only needed for QM atoms
        ofile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        ofile << QMMMData[i].MP[Bead].q;
        ofile << '\n';
      }
      if (QMMMData[i].PBRegion)
      {
        //Modify the charge to force charge balance with the boundaries
        vector<int> Boundaries;
        //Boundaries = TraceBoundary(QMMMData,i);
        int mystat=0;
        Boundaries = TraceBoundary(QMMMData,i,mystat,logFile);
        if(mystat!=0){
          exit(0);
        }

        double qnew = QMMMData[i].MP[Bead].q;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qnew -= QMMMData[Boundaries[j]].MP[Bead].q;
        }
        ofile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        ofile << qnew;
        ofile << '\n';
      }
    }
  }
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMRegion)
      {
        //Write new multipole definition for the atom ID
        WriteTINKMPole(QMMMData,ofile,i,Bead);
        ofile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
      if (QMMMData[i].PBRegion)
      {
        //Modify the charge to force charge balance with the boundaries
        double qi = QMMMData[i].MP[Bead].q; //Save a copy
        vector<int> Boundaries;
        //Boundaries = TraceBoundary(QMMMData,i);
        int mystat=0;
        Boundaries = TraceBoundary(QMMMData,i,mystat,logFile);
        if(mystat!=0){
          exit(0);
        }

        double qnew = qi;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qnew -= QMMMData[Boundaries[j]].MP[Bead].q;
        }
        QMMMData[i].MP[Bead].q = qnew; //Save modified charge
        WriteTINKMPole(QMMMData,ofile,i,Bead);
        QMMMData[i].MP[Bead].q = qi; //Return to unmodified charge
        ofile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
      if (QMMMData[i].BARegion)
      {
        ofile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
    }
  }

//Start: Hatice
//input restraints

  if (QMMM or (Nfreeze > 0))
  {
    ofile << '\n';
    for (int i=0;i<Natoms;i++)
    {
      /*Add active atoms*/
      if (QMMMData[i].MMRegion or QMMMData[i].BARegion)
      {
        if (!QMMMData[i].frozen)
        {

            /*Start a new active line*/
            /*ofile << "active ";*/
            ofile << "restrain-position ";
            ofile << setw(6) << (QMMMData[i].id+1);
            ofile << " ";
            ofile << LICHEMFormFloat(QMMMData[i].P[Bead].x,16);
            ofile << " ";
            ofile << LICHEMFormFloat(QMMMData[i].P[Bead].y,16);
            ofile << " ";
            ofile << LICHEMFormFloat(QMMMData[i].P[Bead].z,16);
            ofile << " ";
            ofile << setw(6) << restr;
            ofile << '\n';
        }
      }
    }

  }  

//End: Hatice
  ofile.flush();
  ofile.close();

  /*Create TINKER xyz file from the structure*/
  call.str("");
  call << "LICHM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  /*Write atoms to the xyz file*/
  ofile << Natoms << '\n';
  if (PBCon)
  {
    /*Write box size*/
    ofile << LICHEMFormFloat(Lx,12) << " ";
    ofile << LICHEMFormFloat(Ly,12) << " ";
    ofile << LICHEMFormFloat(Lz,12) << " ";
    ofile << "90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; /*Counter for QM atoms*/
  for (int i=0;i<Natoms;i++)
  {
    ofile << setw(6) << (QMMMData[i].id+1);
    ofile << " ";
    ofile << setw(3) << QMMMData[i].MMTyp;
    ofile << " ";
    ofile << LICHEMFormFloat(QMMMData[i].P[Bead].x,16);
    ofile << " ";
    ofile << LICHEMFormFloat(QMMMData[i].P[Bead].y,16);
    ofile << " ";
    ofile << LICHEMFormFloat(QMMMData[i].P[Bead].z,16);
    ofile << " ";
    ofile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      ofile << " "; /*Avoids trailing spaces*/
      ofile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    ofile << '\n';
  }


  ofile.flush();
  ofile.close();
  //Run optimization
  call.str("");
  call << "minimize LICHM_";
  call << Bead << ".xyz ";
  call << QMMMOpts.MMOptTol << " > LICHM_";
  call << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Read new structure
  call.str("");
  call << "LICHM_" << Bead << ".xyz_2";
  ifile.open(call.str().c_str(),ios_base::in);
  getline(ifile,dummy); //Discard number of atoms
  if (PBCon)
  {
    //Discard PBC information
    getline(ifile,dummy);
  }
  for (int i=0;i<Natoms;i++)
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    //Read new positions
    line >> dummy >> dummy; //Discard atom ID and type
    line >> QMMMData[i].P[Bead].x;
    line >> QMMMData[i].P[Bead].y;
    line >> QMMMData[i].P[Bead].z;
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " LICHM_" << Bead << ".xyz";
  call << " LICHM_" << Bead << ".log";
  call << " LICHM_" << Bead << ".xyz_*";
  call << " LICHM_" << Bead << ".key";
  call << "cp LICHM_"<< Bead << ".xyz ";
  //call << "LICHM_TINKERRest_"<< Bead << "_" << restr <<".xyz";
  //call << "; cp LICHM_"<< Bead << ".log ";
  //call << "LICHM_TINKERRest_"<< Bead << "_" << restr <<".log";
  globalSys = system(call.str().c_str());
  //Change units
  E = (E/627.51);//turn it to to a.u.;
  return E;
};
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* Start: Aug 28 2018 */

void CalcFreq(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,fstream& logFile)
{

    /*if (QMMMOpts.NEBFreq)
    {*/
      /*Calculate TS frequencies*/
      int remCt = 0; /*Number of deleted translation and rotation modes*/
      int Ndof = 3*(Nqm+Npseudo); /*Number of degrees of freedom*/
      MatrixXd QMMMHess(Ndof,Ndof);
      VectorXd QMMMFreqs(Ndof);
      logFile << '\n'; /*Print blank line*/
      logFile << "TS frequencies:";
      logFile << '\n';
      logFile.flush(); /*Print progress*/
      /*Calculate QMMM frequencies*/
      QMMMHess.setZero(); /*Reset Hessian*/
      QMMMFreqs.setZero(); /*Reset frequencies*/
      /*Calculate QM Hessian*/
      if (Gaussian)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += GaussianHessian(QMMMData,QMMMOpts,QMMMOpts.TSBead);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (PSI4)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += PSI4Hessian(QMMMData,QMMMOpts,QMMMOpts.TSBead);
        QMTime += (unsigned)time(0)-tStart;
        /*Delete annoying useless files*/
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += NWChemHessian(QMMMData,QMMMOpts,QMMMOpts.TSBead);
        QMTime += (unsigned)time(0)-tStart;
      }
      /*Calculate MM Hessian*/
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += TINKERHessian(QMMMData,QMMMOpts,QMMMOpts.TSBead);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += LAMMPSHessian(QMMMData,QMMMOpts,QMMMOpts.TSBead);
        MMTime += (unsigned)time(0)-tStart;
      }
      /*Calculate frequencies*/
      QMMMFreqs = LICHEMFreq(QMMMData,QMMMHess,QMMMOpts,QMMMOpts.TSBead,remCt);
      /*Print the frequencies*/
      if (remCt > 0)
      {
        logFile << "  | Identified " << remCt;
        logFile << " translation/rotation modes";
        logFile << '\n';
      }
      logFile << "  | Frequencies:" << '\n' << '\n';
      logFile << "   ";
      remCt = 0; /*Reuse as a counter*/
      for (int i=0;i<Ndof;i++)
      {
        if (abs(QMMMFreqs(i)) > 0)
        {
          logFile << " ";
          logFile << LICHEMFormFloat(QMMMFreqs(i),10);
          remCt += 1;
          if (remCt == 3)
          {
            /*Start a new line*/
            logFile << '\n';
            logFile << "   "; /*Add extra space*/
            remCt = 0;
          }
        }
      }
      if (remCt != 0)
      {
        /*Terminate trailing line*/
        logFile << '\n';
      }
      logFile << '\n';
    /*}*/
    logFile.flush();

}
/* End: Aug 28 2018 */
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool QSMConverged(vector<QMMMAtom>& QMMMData,
                  vector<QMMMAtom>& OldQMMMData,
                  int stepct, QMMMSettings& QMMMOpts, 
                  VectorXd& Eqmmm_images,fstream& logFile)
{
  
  //Check convergence of QMMM optimizations
  stringstream call; //Stream for system calls and reading/writing files
  string dummy; //Generic string
  call.str("");
  //Initialize stats variables
  bool PathDone = 0;
  double RMSdiff = 0;
  double RMSforce = 0;
  double MAXforce = 0;
  double SumE = 0;
  double Eqm = 0;
  double Emm = 0;
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
 
    //Check if the MM region changed and gather statistics
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      SumE = 0;
      Eqm=0;
      Emm=0;
      //SumEeV=0;
      //Calculate QM energy
      if (Gaussian)
      {
        int tstart = (unsigned)time(0);
        Eqm += GaussianEnergy(QMMMData,QMMMOpts,p);
        Eqm = Eqm/har2eV;//to a.u.
        QMTime += (unsigned)time(0)-tstart;
      }
      if (PSI4)
      {
        int tstart = (unsigned)time(0);
        Eqm += PSI4Energy(QMMMData,QMMMOpts,p);
        Eqm = Eqm/har2eV;//to a.u.
        QMTime += (unsigned)time(0)-tstart;
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tstart = (unsigned)time(0);
        Eqm += NWChemEnergy(QMMMData,QMMMOpts,p);
        Eqm = Eqm/har2eV;//to a.u.
        QMTime += (unsigned)time(0)-tstart;
      }
      //Calculate MM energy
      if (TINKER)
      {
        int tstart = (unsigned)time(0);
        Emm += TINKEREnergy(QMMMData,QMMMOpts,p,logFile);
        Emm = Emm/har2eV; //turn it to a.u.
        MMTime += (unsigned)time(0)-tstart;
      }
      //if (AMBER)
      //{
      //  int tstart = (unsigned)time(0);
      //  Emm += AMBEREnergy(QMMMData,QMMMOpts,p);
      //  Emm = Emm/har2eV; //turn it to a.u.
      //  MMTime += (unsigned)time(0)-tstart;
      //}
      if (LAMMPS)
      {
        int tstart = (unsigned)time(0);
        Emm += LAMMPSEnergy(QMMMData,QMMMOpts,p);
        Emm = Emm/har2eV; //turn it to a.u.
        MMTime += (unsigned)time(0)-tstart;
      }
      //Calculate RMS displacement
      #pragma omp parallel for schedule(dynamic) reduction(+:RMSdiff)
      for (int i=0;i<Natoms;i++)
      {
        double RMStmp = 0; //Store a local sum
        for (int j=0;j<i;j++)
        {
          double Rnew = 0;
          double Rold = 0;
          Rnew = CoordDist2(QMMMData[i].P[p],QMMMData[j].P[p]).vecMag();
          Rold = CoordDist2(OldQMMMData[i].P[p],OldQMMMData[j].P[p]).vecMag();
          Rnew = sqrt(Rnew);
          Rold = sqrt(Rold);
          //Update local sum
          RMStmp += (Rnew-Rold)*(Rnew-Rold);
        }
        //Update sum
        RMSdiff += RMStmp;
      }
      #pragma omp barrier
      SumE=Eqm+Emm;
      Eqmmm_images[p] = SumE;
    }
    int AdjustedBeads; //Number of moving beads
    if (QMMMOpts.frznEnds)
    {
      //End points do not count
      AdjustedBeads = QMMMOpts.NBeads-2;
    }
    else
    {
      //All beads
      AdjustedBeads = QMMMOpts.NBeads;
    }
    RMSdiff /= (Natoms-Nfreeze)*(Natoms-Nfreeze-1)/2;
    RMSdiff /= AdjustedBeads; //Adjust for multiple replicas
    RMSdiff = sqrt(RMSdiff);
    //Update energies
    SumE = 0; //Reusing this variable to avoid making a new one
    //SumE = Es.maxCoeff();
    SumE = Eqmmm_images.maxCoeff();
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      if (p == 0)
      {
        //Using an if statement incase ends are frozen
        QMMMOpts.EReact = Eqmmm_images[p];
      }
      if (Eqmmm_images[p] == SumE)
      {
        //Save new energy
        QMMMOpts.TSBead = p;
        QMMMOpts.ETrans = Eqmmm_images[p];
      }
      if (p == (QMMMOpts.NBeads-1))
      {
        //Using an if statement incase ends are frozen
        QMMMOpts.EProd = Eqmmm_images[p];
      }
    }

    //Check convergence
    if (RMSdiff <= QMMMOpts.MMOptTol)
    {
      PathDone = 1;
    }

  return PathDone;
};
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Convergence test functions
bool QMConverged(vector<QMMMAtom>& QMMMData, vector<QMMMAtom>& OldQMMMData,
     MatrixXd& ForceStats, int stepct, QMMMSettings& QMMMOpts, 
     double &rmsdiff, double &rmsforce, double &maxforce)
{
  //Check convergence of QMMM optimizations
  stringstream call; //Stream for system calls and reading/writing files
  string dummy; //Generic string
  call.str("");
  //Initialize stats variables
  double RMSdiff = 0;
  double RMSforce = 0;
  double MAXforce = 0;
  bool PathDone = 0;

  double SumE = 0;

  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Check progress of all beads

    //Check if a QM calculation is converged
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      //Find max forces and RMSforce
      if (MAXforce < ForceStats(p,0))
      {
        MAXforce = ForceStats(p,0);
      }
      RMSforce += ForceStats(p,1);
      
      //Find RMS deviation for the whole path
      #pragma omp parallel for schedule(dynamic) reduction(+:RMSdiff)
      for (int i=0;i<Natoms;i++)
      {
        //Calculate RMS displacement
        double RMStmp = 0; //Store a local sum
        if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
        {
          for (int j=0;j<i;j++)
          {
            if (QMMMData[j].QMRegion or QMMMData[j].PBRegion)
            {
              double Rnew = 0;
              double Rold = 0;
              Rnew = CoordDist2(QMMMData[i].P[p],QMMMData[j].P[p]).vecMag();
              Rold = CoordDist2(OldQMMMData[i].P[p],OldQMMMData[j].P[p]).vecMag();
              Rnew = sqrt(Rnew);
              Rold = sqrt(Rold);
              //Update local sum
              RMStmp += (Rnew-Rold)*(Rnew-Rold);
            }
          }
        }
        //Update sum
        RMSdiff += RMStmp;
      }
      #pragma omp barrier
   
    }
    int AdjustedBeads; //Number of moving beads
    if (QMMMOpts.frznEnds)
    {
      //End points do not count
      AdjustedBeads = QMMMOpts.NBeads-2;
    }
    else
    {
      //All beads
      AdjustedBeads = QMMMOpts.NBeads;
    }

    RMSdiff /= (Nqm+Npseudo)*(Nqm+Npseudo-1)/2;
    RMSdiff /= AdjustedBeads; //Adjust for multiple replicas
    RMSdiff = sqrt(RMSdiff);

    RMSforce /= Ndof;
    RMSforce /= AdjustedBeads; //Adjust for multiple replicas
    RMSforce = sqrt(RMSforce);
    /* Forcestats is in bohrRad
    RMSforce *=bohrRad;//convert to Hartree/bohr
    MAXforce *=bohrRad;//convert to Hartree/bohr
    */
    //Check convergence criteria
    if ((RMSdiff <= QMMMOpts.QMOptTol) and
        (RMSforce <= QMMMOpts.QMRMSForceTol) and
        (MAXforce <= QMMMOpts.QMMaxForceTol))
    {
        //Finish the optimization
        PathDone = 1;
     
    }
    rmsdiff = RMSdiff;
    rmsforce = RMSforce;
    maxforce = MAXforce;

  return PathDone;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

