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


 LICHEM is licensed under GPLv3, for more information see GPL_LICENSE

 References for the LICHEM package:
 Kratz et al., J. Comput. Chem., 37, 11, 1019, (2016)

*/

//Primary LICHEM header
#include "LICHEM_headers.h"

int main(int argc, char* argv[])
{

  //Misc. initialization
  startTime = (unsigned)time(0); //Time the program starts
  srand((unsigned)time(0)); //Serial only random numbers
  //End of section

  //Initialize local variables
  string dummy; //Generic string
  double sumE,sumE2,denAvg,LxAvg,LyAvg,LzAvg,Ek; //Energies and properties
  fstream xyzFile,connectFile,regionFile,outFile; //Input and output files
  fstream logFile, errFile;
  vector<QMMMAtom> QMMMData; //Atom list
  vector<QMMMAtom> OldQMMMData; //A copy of the atoms list
  QMMMSettings QMMMOpts; //QM and MM wrapper settings
  int randNum; //Random integer
  //End of section
  int stat=0;
  //Read arguments and look for errors
  ReadArgs(argc,argv,xyzFile,connectFile,regionFile,outFile,logFile,errFile,stat);
  if(stat!=0){
    logFile.close();
    errFile.close();
    exit(0);
  }
  //End of section

  //Output stream settings
  //NB: The streams should always be returned to these settings
  logFile.precision(16);
  errFile.precision(16);
  //End of section

  //Print title and compile date
  PrintFancyTitle(logFile);
  logFile << '\n';
  logFile << "Last modification: ";
  logFile << __TIME__ << " on ";
  logFile << __DATE__ << '\n';
  logFile << '\n';
  logFile.flush();
  //End of section

  //Print early messages
  logFile << "Reading input..." << '\n';
  logFile << '\n';
  logFile.flush();
  //End of section

  //Read input and check for errors
  ReadLICHEMInput(xyzFile,connectFile,regionFile,QMMMData,QMMMOpts,logFile,stat);
  if(stat!=0){
    logFile.close();
    errFile.close();
    exit(0);
  }

  LICHEMErrorChecker(QMMMOpts,logFile,stat);
  if(stat!=0){
    logFile.close();
    errFile.close();
    exit(0);
  }

  LICHEMPrintSettings(QMMMData,QMMMOpts,logFile);
  //End of section
  //Fix PBC
  if (PBCon)
  {
    //Relatively safe PBC correction
    if (!TINKER)
    {
      PBCCenter(QMMMData,QMMMOpts); //Center the atoms in the box
    }
  }
  //End of section

  //Create backup directories
  if (CheckFile("BACKUPQM"))
  {
    stringstream call;
    call.str("");
    //Delete old files
    call << "rm -rf " << QMMMOpts.backDir << "; ";
    //Create new directory
    call << "mkdir " << QMMMOpts.backDir;
    globalSys = system(call.str().c_str());
  }
  //End of section

  /*Set keep file per step */
  QMMMOpts.perOpt=QMMMOpts.maxOptSteps;
  QMMMOpts.perQM=QMMMOpts.MaxQMSteps;
  /*
    NB: All optional simulation types should be wrapped in comments and
    else-if statements. The first comment should define what calculation is
    going to be performed, then the simulation should be enclosed in an
    else-if statement. After the else-if, an "//End of section" comment
    should be added to mark where the next simulation type begins.
  */


  //Calculate single-point energy
  if (SinglePoint)
  {
    double Eqm; //QM energy
    double Emm; //MM energy
    if (QMMMOpts.NBeads == 1)
    {
      logFile << "Single-point energy:";
    }
    if (QMMMOpts.NBeads > 1)
    {
      logFile << "Multi-point energies:";
    }
    logFile << '\n' << '\n';
    logFile.flush(); //Print progress
    //Loop over all beads
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      //Calculate QMMM energy
      Eqm = 0; //Reset QM energy
      Emm = 0; //Reset MM energy
      //Calculate QM energy
      if (QMMMOpts.NBeads > 1)
      {
        logFile << " Energy for bead: " << p << '\n';
        logFile.flush();
      }
      if (Gaussian)
      {
        int tStart = (unsigned)time(0);
        Eqm += GaussianEnergy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (PSI4)
      {
        int tStart = (unsigned)time(0);
        Eqm += PSI4Energy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tStart = (unsigned)time(0);
        Eqm += NWChemEnergy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (QMMM or QMonly)
      {
        //Print QM partial energy
        logFile << "  QM energy: " << LICHEMFormFloat(Eqm,16) << " eV";
        logFile << '\n';
        //Print progress
        logFile.flush();
      }
      //Calculate MM energy
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        Emm += TINKEREnergy(QMMMData,QMMMOpts,p,logFile);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        Emm += LAMMPSEnergy(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      //Print the rest of the energies
      if (QMMM or MMonly)
      {
        //Print MM partial energy
        logFile << "  MM energy: " << LICHEMFormFloat(Emm,16) << " eV";
        logFile << '\n';
      }
      sumE = Eqm+Emm; //Total energy
      if (QMMM)
      {
        //Print total energy
        logFile << "  QMMM energy: ";
        logFile << LICHEMFormFloat(sumE,16) << " eV";
        logFile << " ";
        logFile << LICHEMFormFloat(sumE/har2eV,16) << " a.u.";
        logFile << '\n';
      }
      logFile << '\n';
      logFile.flush(); //Print output
    }
  }
  //End of section

  //LICHEM frequency calculation
  else if (FreqCalc)
  {
    int remCt = 0; //Number of deleted translation and rotation modes
    int Ndof = 3*(Nqm+Npseudo); //Number of degrees of freedom
    MatrixXd QMMMHess(Ndof,Ndof);
    VectorXd QMMMFreqs(Ndof);
    if (QMMMOpts.NBeads == 1)
    {
      logFile << "Single-point frequencies:";
    }
    if (QMMMOpts.NBeads > 1)
    {
      logFile << "Multi-point frequencies:";
    }
    logFile << '\n';
    logFile.flush(); //Print progress
    //Loop over all beads
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      //Calculate QMMM frequencies
      QMMMHess.setZero(); //Reset frequencies
      QMMMFreqs.setZero();
      //Calculate QM energy
      if (QMMMOpts.NBeads > 1)
      {
        logFile << '\n';
        logFile << " Frequencies for bead: " << p << '\n';
        logFile.flush();
      }
      if (Gaussian)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += GaussianHessian(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (PSI4)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += PSI4Hessian(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += NWChemHessian(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      //Calculate MM energy
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += TINKERHessian(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += LAMMPSHessian(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      //Calculate frequencies
      QMMMFreqs = LICHEMFreq(QMMMData,QMMMHess,QMMMOpts,p,remCt);
      //Print the frequencies
      if (remCt > 0)
      {
        logFile << "  | Identified " << remCt;
        logFile << " translation/rotation modes";
        logFile << '\n';
      }
      logFile << "  | Frequencies:" << '\n' << '\n';
      logFile << "   ";
      remCt = 0; //Reuse as a counter
      for (int i=0;i<Ndof;i++)
      {
        if (abs(QMMMFreqs(i)) > 0)
        {
          logFile << " ";
          logFile << LICHEMFormFloat(QMMMFreqs(i),10);
          remCt += 1;
          if (remCt == 3)
          {
            //Start a new line
            logFile << '\n';
            logFile << "   "; //Add extra space
            remCt = 0;
          }
        }
      }
      if (remCt != 0)
      {
        //Terminate trailing line
        logFile << '\n';
      }
      logFile << '\n';
      logFile.flush(); //Print output
    }
  }
  //End of section

  //Optimize structure (native QM and MM package optimizers)
  else if (OptSim)
  {
    //NB: Currently only Gaussian works with this option
    VectorXd forces; //Dummy array needed for convergence tests
    int optCt = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(QMMMData,outFile,QMMMOpts);
    logFile << "Optimization:" << '\n';
    logFile.flush(); //Print progress
    //Calculate initial energy
    sumE = 0; //Clear old energies
    //Calculate QM energy
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      sumE += GaussianEnergy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      sumE += PSI4Energy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      sumE += NWChemEnergy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      sumE += TINKEREnergy(QMMMData,QMMMOpts,0,logFile);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      sumE += LAMMPSEnergy(QMMMData,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tStart;
    }
    logFile << " | Opt. step: ";
    logFile << optCt << " | Energy: ";
    logFile << LICHEMFormFloat(sumE,16) << " eV";
    logFile << '\n';
    logFile.flush(); //Print progress
    //Run optimization
    bool optDone = 0;
    if (QMMMOpts.maxOptSteps == 0)
    {
      optDone = 1;
    }
    while (!optDone)
    {
      //Copy structure
      OldQMMMData = QMMMData;
      //Run MM optimization
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        int mystat=0;
        sumE = TINKEROpt(QMMMData,QMMMOpts,0,logFile,mystat);
        if(mystat!=0){
           exit(0);
        }
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        sumE = LAMMPSOpt(QMMMData,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (QMMM)
      {
        logFile << "    MM optimization complete.";
        logFile << '\n';
        logFile.flush();
      }
      logFile << '\n';
      //Run QM optimization
      if (Gaussian)
      {
        int tStart = (unsigned)time(0);
        sumE = GaussianOpt(QMMMData,QMMMOpts,0);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (PSI4)
      {
        int tStart = (unsigned)time(0);
        sumE = PSI4Opt(QMMMData,QMMMOpts,0);
        QMTime += (unsigned)time(0)-tStart;
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tStart = (unsigned)time(0);
        sumE = NWChemOpt(QMMMData,QMMMOpts,0);
        QMTime += (unsigned)time(0)-tStart;
      }
      //Print Optimized geometry
      Print_traj(QMMMData,outFile,QMMMOpts);
      //Check convergence
      optCt += 1;
      optDone = OptConverged(QMMMData,OldQMMMData,forces,optCt,QMMMOpts,0,0,logFile);
    }
    logFile << '\n';
    logFile << "Optimization complete.";
    logFile << '\n' << '\n';
    logFile.flush();
  }
  //End of section

  //Steepest descent optimization
  else if (SteepSim)
  {
    VectorXd forces; //Dummy array needed for convergence tests
    int optCt = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(QMMMData,outFile,QMMMOpts);
    logFile << "Steepest descent optimization:" << '\n';
    logFile.flush(); //Print progress
    //Calculate initial energy
    sumE = 0; //Clear old energies
    //Calculate QM energy
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      sumE += GaussianEnergy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      sumE += PSI4Energy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.* ");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      sumE += NWChemEnergy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      sumE += TINKEREnergy(QMMMData,QMMMOpts,0,logFile);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      sumE += LAMMPSEnergy(QMMMData,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tStart;
    }
    logFile << " | Opt. step: ";
    logFile << optCt << " | Energy: ";
    logFile << LICHEMFormFloat(sumE,16) << " eV";
    logFile << '\n';
    logFile.flush(); //Print progress
    //Run optimization
    bool optDone = 0;
    while (!optDone)
    {
      //Copy structure
      OldQMMMData = QMMMData;
      //Run MM optimization
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        int mystat=0;
        sumE = TINKEROpt(QMMMData,QMMMOpts,0,logFile,mystat);
        if(mystat!=0){
           exit(0);
        }
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        sumE = LAMMPSOpt(QMMMData,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (QMMM)
      {
        logFile << "    MM optimization complete.";
        logFile << '\n';
        logFile.flush();
      }
      logFile << '\n';
      //Run QM optimization
      LICHEMSteepest(QMMMData,QMMMOpts,0,logFile);
      //Print Optimized geometry
      Print_traj(QMMMData,outFile,QMMMOpts);
      //Check convergence
      optCt += 1;
      optDone = OptConverged(QMMMData,OldQMMMData,forces,optCt,QMMMOpts,0,0,logFile);
    }
    logFile << '\n';
    logFile << "Optimization complete.";
    logFile << '\n' << '\n';
    logFile.flush();
  }
  //End of section

  //DFP minimization
  else if (DFPSim)
  {
    VectorXd forces; //Dummy array needed for convergence tests
    int optCt = 0; //Counter for optimization steps
    //Change optimization tolerance for the first step
    double savedQMOptTol = QMMMOpts.QMOptTol; //Save value from input
    double savedMMOptTol = QMMMOpts.MMOptTol; //Save value from input
    if (QMMMOpts.QMOptTol < 0.005)
    {
      QMMMOpts.QMOptTol = 0.005; //Speedy convergance on the first step
    }
    if (QMMMOpts.MMOptTol < 0.25)
    {
      QMMMOpts.MMOptTol = 0.25; //Speedy convergance on the first step
    }
    //Print initial structure
    Print_traj(QMMMData,outFile,QMMMOpts);
    logFile << "DFP optimization:" << '\n';
    logFile.flush(); //Print progress
    //Calculate initial energy
    sumE = 0; //Clear old energies
    //Calculate QM energy
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      sumE += GaussianEnergy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      sumE += PSI4Energy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      sumE += NWChemEnergy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      sumE += TINKEREnergy(QMMMData,QMMMOpts,0,logFile);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      sumE += LAMMPSEnergy(QMMMData,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tStart;
    }
    logFile << " | Opt. step: ";
    logFile << optCt << " | Energy: ";
    logFile << LICHEMFormFloat(sumE,16) << " eV";
    logFile << '\n';
    logFile.flush(); //Print progress
    //Run optimization
    bool optDone = 0;
    while (!optDone)
    {
      //Copy structure
      OldQMMMData = QMMMData;
      //Run MM optimization
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        int mystat=0;
        sumE = TINKEROpt(QMMMData,QMMMOpts,0,logFile,mystat);
        if(mystat!=0){
           exit(0);
        }
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        sumE = LAMMPSOpt(QMMMData,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (QMMM)
      {
        logFile << "    MM optimization complete.";
        logFile << '\n';
        logFile.flush();
      }
      logFile << '\n';
      //Run QM optimization
      LICHEMDFP(QMMMData,QMMMOpts,0,logFile);
      //Reset tolerance before optimization check
      QMMMOpts.QMOptTol = savedQMOptTol;
      QMMMOpts.MMOptTol = savedMMOptTol;
      //Print Optimized geometry
      Print_traj(QMMMData,outFile,QMMMOpts);
      //Check convergence
      optCt += 1;
      optDone = OptConverged(QMMMData,OldQMMMData,forces,optCt,QMMMOpts,0,0,logFile);
      if (optCt == 1)
      {
        //Avoid terminating restarts on the loose tolerance step
        optDone = 0; //Not converged
      }
    }
    logFile << '\n';
    logFile << "Optimization complete.";
    logFile << '\n' << '\n';
    logFile.flush();
  }
  //End of section

  //Run Monte Carlo
  else if (PIMCSim)
  {
    //Change units
    QMMMOpts.press *= atm2eV; //Pressure in eV/Ang^3
    //Adjust probabilities
    if (Natoms == 1)
    {
      //Remove atom centroid moves
      centProb = 0.0;
      beadProb = 1.0;
    }
    if (QMMMOpts.ensemble == "NVT")
    {
      //Remove volume changes
      volProb = 0.0;
    }
    //Initialize local variables
    sumE = 0; //Average energy
    sumE2 = 0; //Average squared energy
    denAvg = 0; //Average density
    LxAvg = 0; //Average box length
    LyAvg = 0; //Average box length
    LzAvg = 0; //Average box length
    Ek = 0; //PIMC kinietic energy
    if (QMMMOpts.NBeads > 1)
    {
      //Set kinetic energy
      Ek = 3*Natoms*QMMMOpts.NBeads/(2*QMMMOpts.beta);
    }
    int Nct = 0; //Step counter
    int ct = 0; //Secondary counter
    double Nacc = 0; //Number of accepted moves
    double Nrej = 0; //Number of rejected moves
    double Emc = 0; //Monte Carlo energy
    double Et = 0; //Total energy for printing
    bool acc; //Flag for accepting a step
    //Find the number of characters to print for the step counter
    int simCharLen;
    simCharLen = QMMMOpts.NEq+QMMMOpts.NSteps;
    simCharLen = LICHEMCount(simCharLen);
    //Start equilibration run and calculate initial energy
    logFile << "Monte Carlo equilibration:" << '\n';
    logFile.flush();
    QMMMOpts.EOld = 0;
    QMMMOpts.EOld += Get_PI_Epot(QMMMData,QMMMOpts,logFile);
    QMMMOpts.EOld += Get_PI_Espring(QMMMData,QMMMOpts);
    if (volProb > 0)
    {
      //Add PV term
      QMMMOpts.EOld += QMMMOpts.press*Lx*Ly*Lz;
    }
    Emc = QMMMOpts.EOld; //Needed if equilibration is skipped
    Nct = 0; //Reset counter to zero
    while (Nct < QMMMOpts.NEq)
    {
      Emc = 0;
      //Check step size
      if(ct == acc_Check)
      {
        if ((Nacc/(Nrej+Nacc)) > QMMMOpts.accRatio)
        {
          //Increase step size
          double randVal; //Use random values to keep from cycling up and down
          randVal = (((double)rand())/((double)RAND_MAX));
          randVal /= 10.0;
          mcStep *= 1.001+randVal;
        }
        if ((Nacc/(Nrej+Nacc)) < QMMMOpts.accRatio)
        {
          //Decrease step size
          double randVal; //Use random values to keep from cycling up and down
          randVal = (((double)rand())/((double)RAND_MAX));
          randVal /= 10.0;
          mcStep *= 0.999-randVal;
        }
        if (mcStep < stepMin)
        {
          //Set to minimum
          mcStep = stepMin;
        }
        if (mcStep > stepMax)
        {
          //Set to maximum
          mcStep = stepMax;
        }
        //Statistics
        logFile << " | Step: " << setw(simCharLen) << Nct;
        logFile << " | Step size: ";
        logFile << LICHEMFormFloat(mcStep,6);
        logFile << " | Accept ratio: ";
        logFile << LICHEMFormFloat((Nacc/(Nrej+Nacc)),6);
        logFile << '\n';
        logFile.flush(); //Print stats
        //Reset counters
        ct = 0;
        Nacc = 0;
        Nrej = 0;
      }
      //Continue simulation
      ct += 1;
      acc = MCMove(QMMMData,QMMMOpts,Emc,logFile);
      if (acc)
      {
        Nct += 1;
        Nacc += 1;
      }
      else
      {
        Nrej += 1;
      }
    }
    logFile << " Equilibration complete." << '\n';
    //Start production run
    Nct = 0; //Reset counter to zero
    Nacc = 0; //Reset counter to zero
    Nrej = 0; //Reset counter to zero
    logFile << '\n';
    logFile << "Monte Carlo production:" << '\n';
    logFile.flush();
    //Print starting conditions
    Print_traj(QMMMData,outFile,QMMMOpts);
    Et = Ek+Emc; //Calculate total energy using previous saved energy
    Et -= 2*Get_PI_Espring(QMMMData,QMMMOpts);
    logFile << " | Step: " << setw(simCharLen) << 0;
    logFile << " | Energy: " << LICHEMFormFloat(Et,12);
    logFile << " eV";
    if (QMMMOpts.ensemble == "NPT")
    {
      double rho;
      rho = LICHEMDensity(QMMMData,QMMMOpts);
      logFile << " | Density: ";
      logFile << LICHEMFormFloat(rho,8);
      logFile << " g/cm\u00B3";
    }
    logFile << '\n';
    logFile.flush(); //Print results
    //Continue simulation
    while (Nct < QMMMOpts.NSteps)
    {
      Emc = 0; //Set energy to zero
      acc = MCMove(QMMMData,QMMMOpts,Emc,logFile);
      //Update averages
      Et = 0;
      Et += Ek+Emc;
      Et -= 2*Get_PI_Espring(QMMMData,QMMMOpts);
      denAvg += LICHEMDensity(QMMMData,QMMMOpts);
      LxAvg += Lx;
      LyAvg += Ly;
      LzAvg += Lz;
      sumE += Et;
      sumE2 += Et*Et;
      //Update counters and print output
      if (acc)
      {
        //Increase counters
        Nct += 1;
        Nacc += 1;
        //Print trajectory and instantaneous energies
        if ((Nct%QMMMOpts.NPrint) == 0)
        {
          //Print progress
          Print_traj(QMMMData,outFile,QMMMOpts);
          logFile << " | Step: " << setw(simCharLen) << Nct;
          logFile << " | Energy: " << LICHEMFormFloat(Et,12);
          logFile << " eV";
          if (QMMMOpts.ensemble == "NPT")
          {
            double rho;
            rho = LICHEMDensity(QMMMData,QMMMOpts);
            logFile << " | Density: ";
            logFile << LICHEMFormFloat(rho,8);
            logFile << " g/cm\u00B3";
          }
          logFile << '\n';
          logFile.flush(); //Print results
        }
      }
      else
      {
        Nrej += 1;
      }
    }
    if ((Nct%QMMMOpts.NPrint) != 0)
    {
      //Print final geometry if it was not already written
      Print_traj(QMMMData,outFile,QMMMOpts);
    }
    sumE /= Nrej+Nacc; //Average energy
    sumE2 /= Nrej+Nacc; //Variance of the energy
    denAvg /= Nrej+Nacc; //Average density
    LxAvg /= Nrej+Nacc; //Average box size
    LyAvg /= Nrej+Nacc; //Average box size
    LzAvg /= Nrej+Nacc; //Average box size
    //Print simulation details and statistics
    logFile << '\n';
    if (QMMMOpts.NBeads > 1)
    {
      logFile << "PI";
    }
    logFile << "MC statistics:" << '\n';
    if (QMMMOpts.ensemble == "NPT")
    {
      //Print simulation box information
      logFile << " | Density: ";
      logFile << LICHEMFormFloat(denAvg,8);
      logFile << " g/cm\u00B3" << '\n';
      logFile << " | Average box size (\u212B): " << '\n';
      logFile << "  "; //Indent
      logFile << " Lx = " << LICHEMFormFloat(LxAvg,12);
      logFile << " Ly = " << LICHEMFormFloat(LyAvg,12);
      logFile << " Lz = " << LICHEMFormFloat(LzAvg,12);
      logFile << '\n';
    }
    logFile << " | Average energy: ";
    logFile << LICHEMFormFloat(sumE,16);
    logFile << " eV | Variance: ";
    logFile << LICHEMFormFloat((sumE2-(sumE*sumE)),12);
    logFile << " eV\u00B2";
    logFile << '\n';
    logFile << " | Acceptance ratio: ";
    logFile << LICHEMFormFloat((Nacc/(Nrej+Nacc)),6);
    logFile << " | Optimized step size: ";
    logFile << LICHEMFormFloat(mcStep,6);
    logFile << " \u212B";
    logFile << '\n';
    logFile << '\n';
    logFile.flush();
  }
  //End of section

  //Force-bias NEB Monte Carlo
  else if (FBNEBSim)
  {
    //Initialize local variables
    int Nct = 0; //Step counter
    int ct = 0; //Secondary counter
    double Nacc = 0; //Number of accepted moves
    double Nrej = 0; //Number of rejected moves
    int acc; //Number of steps accepted along the path
    vector<VectorXd> allForces; //Stores forces between MC steps
    VectorXd sumE(QMMMOpts.NBeads); //Average energy array
    VectorXd sumE2(QMMMOpts.NBeads); //Average squared energy array
    VectorXd Emc(QMMMOpts.NBeads); //Current MC energy
    sumE.setZero();
    sumE2.setZero();
    Emc.setIdentity(); //Initial energies should not be zero
    Emc *= hugeNum; //Forces the first step to be accepted
    //Initialize force arrays
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      //Zero force vectors make the first move an energy calculation
      VectorXd tmp(3*Natoms);
      tmp.setZero();
      allForces.push_back(tmp);
    }
    //Find the number of characters to print for the step counter
    int simCharLen;
    simCharLen = QMMMOpts.NEq+QMMMOpts.NSteps;
    simCharLen = LICHEMCount(simCharLen);
    //Start equilibration run
    logFile << "Monte Carlo equilibration:" << '\n';
    logFile.flush();
    int savedNPrint = QMMMOpts.NPrint;
    if (QMMMOpts.NPrint < 100)
    {
      //Prevent the print rate from breaking the tuning
      QMMMOpts.NPrint = 100; //Minimum value
    }
    Nct = 0; //Reset counter to zero
    while (Nct < QMMMOpts.NEq)
    {
      //Check step size
      if (ct == QMMMOpts.NPrint)
      {
        if ((Nacc/(Nrej+Nacc)) > QMMMOpts.accRatio)
        {
          //Increase step size
          double randVal; //Use random values to keep from cycling up and down
          randVal = (((double)rand())/((double)RAND_MAX));
          randVal /= 10.0;
          mcStep *= 1.001+randVal;
        }
        if ((Nacc/(Nrej+Nacc)) < QMMMOpts.accRatio)
        {
          //Decrease step size
          double randVal; //Use random values to keep from cycling up and down
          randVal = (((double)rand())/((double)RAND_MAX));
          randVal /= 10.0;
          mcStep *= 0.999-randVal;
        }
        if (mcStep < stepMin)
        {
          //Set to minimum
          mcStep = stepMin;
        }
        if (mcStep > stepMax)
        {
          //Set to maximum
          mcStep = stepMax;
        }
        //Statistics
        logFile << " | Accepted: " << setw(simCharLen) << Nct;
        logFile << " | Step size: ";
        logFile << LICHEMFormFloat(mcStep,6);
        logFile << " | Accept ratio: ";
        logFile << LICHEMFormFloat((Nacc/(Nrej+Nacc)),6);
        logFile << '\n';
        logFile.flush(); //Print stats
        //Reset counters for the next round of tuning
        ct = 0;
        Nacc = 0;
        Nrej = 0;
      }
      //Continue simulation
      ct += 1; //Equilibration counts cycles instead of steps
      acc = FBNEBMCMove(QMMMData,allForces,QMMMOpts,Emc,logFile);
      Nct += acc; //Equilibration counts acceptances instead of steps
      Nacc += acc;
      Nrej += QMMMOpts.NBeads-acc;
    }
    QMMMOpts.NPrint = savedNPrint; //Restore user defined sample rate
    logFile << " Equilibration complete." << '\n';
    //Start production run
    Nct = 0; //Reset counter to zero
    Nacc = 0; //Reset counter to zero
    Nrej = 0; //Reset counter to zero
    logFile << '\n';
    logFile << "Monte Carlo production:" << '\n';
    logFile.flush();
    //Print starting conditions
    Print_traj(QMMMData,outFile,QMMMOpts);
    logFile << " | Steps: " << setw(simCharLen) << 0;
    logFile << " | Accepted: " << setw(simCharLen) << 0;
    logFile << '\n';
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      logFile << "    Bead: ";
      logFile << setw(3) << p << " | Energy: ";
      logFile << LICHEMFormFloat(Emc(p),16) << " eV" << '\n';
    }
    logFile.flush(); //Print results
    //Continue simulation
    while (Nacc < QMMMOpts.NSteps)
    {
      acc = FBNEBMCMove(QMMMData,allForces,QMMMOpts,Emc,logFile);
      //Update statistics
      #pragma omp parallel for schedule(dynamic)
      for (int p=0;p<QMMMOpts.NBeads;p++)
      {
        sumE(p) += Emc(p);
        sumE2(p) += Emc(p)*Emc(p);
      }
      //Update counters
      Nct += QMMMOpts.NBeads;
      Nacc += acc;
      Nrej += QMMMOpts.NBeads;
      //Print output
      if ((((Nct/QMMMOpts.NBeads)%QMMMOpts.NPrint) == 0) or
         (Nacc == QMMMOpts.NSteps))
      {
        //Print progress
        Print_traj(QMMMData,outFile,QMMMOpts);
        logFile << " | Steps: " << setw(simCharLen) << Nct;
        logFile << " | Accepted: " << setw(simCharLen) << Nacc;
        logFile << '\n';
        for (int p=0;p<QMMMOpts.NBeads;p++)
        {
          logFile << "    Bead: ";
          logFile << setw(3) << p << " | Energy: ";
          logFile << LICHEMFormFloat(Emc(p),16) << " eV" << '\n';
        }
        logFile.flush(); //Print results
      }
    }
    if (((Nct/QMMMOpts.NBeads)%QMMMOpts.NPrint) != 0)
    {
      //Print final geometry if it was not already written
      Print_traj(QMMMData,outFile,QMMMOpts);
    }
    //Calculate statistics
    double scaleBy = 0.0; //Avoid integer division
    scaleBy += QMMMOpts.NBeads; //Adjust for the number of replicas
    scaleBy /= Nct; //Adjust for the total number of samples
    sumE *= scaleBy; //Scale the sum to calculate the average
    sumE2 *= scaleBy; //Scale the sum to calculate the average
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      //Standard deviation
      sumE2(p) = sqrt(sumE2(p)-sumE(p)*sumE(p));
    }
    //Print simulation details and statistics
    logFile << '\n';
    logFile << "Monte Carlo statistics:" << '\n';
    logFile << " | Acceptance ratio: ";
    logFile << LICHEMFormFloat((Nacc/Nct),6);
    logFile << " | Optimized step size: ";
    logFile << LICHEMFormFloat(mcStep,6);
    logFile << " \u212B";
    logFile << '\n';
    logFile << " | Average energies:" << '\n';
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      logFile << "    Bead: ";
      logFile << setw(3) << p << " | Energy: ";
      logFile << LICHEMFormFloat(sumE(p),16);
      logFile << " +/- " << LICHEMFormFloat(sumE2(p),16);
      logFile << " eV" << '\n';
    }
    logFile << '\n';
    logFile.flush();
  }
  //End of section

  //NEB optimization
  else if (NEBSim)
  {
    MatrixXd forceStats; //Dummy array needed for convergence tests
    int optCt = 0; //Counter for optimization steps
    //Check number of beads for Climbing image
    if (QMMMOpts.NBeads < 4)
    {
      //If the system is well behaved, then the TS bead is known.
      QMMMOpts.climb = 1; //Turn on climbing image NEB
    }
    //Change optimization tolerance for the first step
    double savedQMOptTol = QMMMOpts.QMOptTol; //Save value from input
    double savedMMOptTol = QMMMOpts.MMOptTol; //Save value from input
    //Start: Hatice
    double SavedOptTol2 =  QMMMOpts.QMRMSForceTol;
    double SavedOptTol3 = QMMMOpts.QMMaxForceTol;
    /*if (QMMMOpts.QMOptTol < 0.005)
    {
      QMMMOpts.QMOptTol = 0.005; //Speedy convergance on the first step
    }*/
    QMMMOpts.QMOptTol *= 10; //Speedy convergance on the first step
    QMMMOpts.QMRMSForceTol *= 10;
    QMMMOpts.QMMaxForceTol *= 10;
    //End: Hatice
    if (QMMMOpts.MMOptTol < 0.25)
    {
      QMMMOpts.MMOptTol = 0.25; //Speedy convergance on the first step
    }
    //Print initial structure
    Print_traj(QMMMData,outFile,QMMMOpts);
    logFile << "Nudged elastic band optimization:" << '\n';
    if (QMMMOpts.climb)
    {
      logFile << " | Short path detected. Starting climbing image NEB.";
      logFile << '\n' << '\n';
    }
    //Start: Hatice
    //logFile << " | Opt. step: 0 | Bead energies:";
    logFile << " | Initial Energies" << '\n'; 
    //logFile << "     | Bead energies: " << '\n';
    //End: Hatice
    logFile << '\n';
    logFile.flush(); //Print progress
    //Calculate reaction coordinate positions
    VectorXd reactCoord(QMMMOpts.NBeads); //Reaction coordinate
    reactCoord.setZero();
    //Start: Hatice
    VectorXd Eqmmm(QMMMOpts.NBeads);
    Eqmmm.setZero();
    //End: Hatice
    for (int p=0;p<(QMMMOpts.NBeads-1);p++)
    {
      MatrixXd geom1((Nqm+Npseudo),3); //Current replica
      MatrixXd geom2((Nqm+Npseudo),3); //Next replica
      VectorXd disp; //Store the displacement
      //Save geometries
      int ct = 0; //Reset counter for the number of atoms
      for (int i=0;i<Natoms;i++)
      {
        //Only include QM and PB regions
        if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
        {
          //Save current replica
          geom1(ct,0) = QMMMData[i].P[p].x;
          geom1(ct,1) = QMMMData[i].P[p].y;
          geom1(ct,2) = QMMMData[i].P[p].z;
          //Save replica p+1
          geom2(ct,0) = QMMMData[i].P[p+1].x;
          geom2(ct,1) = QMMMData[i].P[p+1].y;
          geom2(ct,2) = QMMMData[i].P[p+1].z;
          ct += 1;
        }
      }
      //Calculate displacement
      disp = KabschDisplacement(geom1,geom2,(Nqm+Npseudo));
      //Remove inactive atoms
      ct = 0; //Reset counter for the number of atoms
      for (int i=0;i<Natoms;i++)
      {
        //Only include QM and PB regions
        if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
        {
          //Only include active atoms in the tangent
          if (!QMMMData[i].NEBActive)
          {
            //Delete distance components
            disp(ct) = 0;
            disp(ct+1) = 0;
            disp(ct+2) = 0;
          }
          //Advance counter
          ct += 3;
        }
      }
      //Update reaction coordinate
      reactCoord(p+1) = reactCoord(p); //Start from previous bead
      reactCoord(p+1) += disp.norm(); //Add magnitude of the displacement
    }
    reactCoord /= reactCoord.maxCoeff(); //Must be between 0 and 1
    //Calculate initial energies
    QMMMOpts.ETrans = -1*hugeNum; //Locate the initial transition state
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      sumE = 0; //Clear old energies
      //Calculate QM energy
      if (Gaussian)
      {
        int tStart = (unsigned)time(0);
        sumE += GaussianEnergy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (PSI4)
      {
        int tStart = (unsigned)time(0);
        sumE += PSI4Energy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tStart = (unsigned)time(0);
        sumE += NWChemEnergy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      //Calculate MM energy
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        sumE += TINKEREnergy(QMMMData,QMMMOpts,p,logFile);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        sumE += LAMMPSEnergy(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (p == 0)
      {
        //Save reactant energy
        QMMMOpts.EReact = sumE;
      }
      else if (p == (QMMMOpts.NBeads-1))
      {
        //Save product energy
        QMMMOpts.EProd = sumE;
      }
      //Start: Hatice
      Eqmmm[p] = sumE;
      //logFile << "   Bead: ";
      //logFile << setw(LICHEMCount(QMMMOpts.NBeads)) << p;
      //logFile << " | React. coord: ";
      //logFile << LICHEMFormFloat(reactCoord(p),5);
      //logFile << " | Energy: ";
      //logFile << LICHEMFormFloat(sumE,16) << " eV";
      //logFile << '\n';
      //End: Hatice
      logFile.flush(); //Print progress
      //Update transition state
      if (sumE > QMMMOpts.ETrans)
      {
        //Save new properties
        QMMMOpts.TSBead = p;
        QMMMOpts.ETrans = sumE;
      }
      //Copy checkpoint data to speed up first step
      if ((p != (QMMMOpts.NBeads-1)) and QMMMOpts.startPathChk)
      {
        stringstream call;
        if (Gaussian and (QMMMOpts.func != "SemiEmp"))
        {
          call.str("");
          call << "cp LICHM_" << p << ".chk ";
          call << "LICHM_" << (p+1) << ".chk";
          call << " 2> LICHM_" << (p+1) << ".trash; ";
          call << "rm -f LICHM_" << (p+1) << ".trash";
          globalSys = system(call.str().c_str());
        }
        if (PSI4)
        {
          call.str("");
          call << "cp LICHM_" << p << ".180 ";
          call << "LICHM_" << (p+1) << ".180";
          call << " 2> LICHM_" << (p+1) << ".trash; ";
          call << "rm -f LICHM_" << (p+1) << ".trash";
          globalSys = system(call.str().c_str());
        }
      }
    }
    /*Start: Hatice*/
    Eqmmm[0] = QMMMOpts.EReact;
    Eqmmm[QMMMOpts.NBeads-1] = QMMMOpts.EProd;
    QMMMOpts.EReact /=har2eV;
    QMMMOpts.EProd  /=har2eV;
    QMMMOpts.ETrans /=har2eV;
    Eqmmm /= har2eV;
    /*print bead energies*/
    print_progress(QMMMOpts, 0,Eqmmm,QMMMOpts.QMOptTol,
                   QMMMOpts.QMMaxForceTol, QMMMOpts.QMRMSForceTol,
                   reactCoord,logFile);
    logFile << "\n";
    QMMMOpts.EReact *= har2eV;
    QMMMOpts.EProd  *= har2eV;
    QMMMOpts.ETrans *= har2eV;
    Eqmmm *= har2eV;

    /*End: Hatice*/

    //Run optimization
    bool pathDone = 0;
    int pathStart = 0; //First bead to optimize
    int pathEnd = QMMMOpts.NBeads; //Last bead to optimize
    if (QMMMOpts.frznEnds)
    {
      //Change the start and end points
      pathStart = 1;
      pathEnd = QMMMOpts.NBeads-1;
    }
    //Start: Hatice
    //while (!pathDone)
    while ( (!pathDone) and (optCt <= QMMMOpts.maxOptSteps))
    //End: Hatice
    {
      //Copy structure
      OldQMMMData = QMMMData;
      //Start: Hatice
      logFile << " | Opt. step    : " << optCt+1;
      logFile << '\n'; 
      //End: Hatice
      //Run MM optimization
      for (int p=pathStart;p<pathEnd;p++)
      {
        if (TINKER)
        {
          int tStart = (unsigned)time(0);
          int mystat=0;
          sumE = TINKEROpt(QMMMData,QMMMOpts,p,logFile,mystat);
          if(mystat!=0){
           exit(0);
          }
          MMTime += (unsigned)time(0)-tStart;
        }
        if (LAMMPS)
        {
          int tStart = (unsigned)time(0);
          sumE = LAMMPSOpt(QMMMData,QMMMOpts,p);
          MMTime += (unsigned)time(0)-tStart;
        }
      }
      if (QMMM)
      {
        logFile << '\n';
        logFile << "    MM optimization complete.";
        logFile << '\n';
        logFile.flush();
      }
      logFile << '\n';
      //Run QM optimization
      LICHEMNEB(QMMMData,QMMMOpts,optCt,logFile);
      //Reset tolerance before optimization check
      QMMMOpts.QMOptTol = savedQMOptTol;
      QMMMOpts.MMOptTol = savedMMOptTol;
      //Start: Hatice
      QMMMOpts.QMRMSForceTol = SavedOptTol2;
      QMMMOpts.QMMaxForceTol = SavedOptTol3;
      //End: Hatice
      //Print optimized geometry
      Print_traj(QMMMData,outFile,QMMMOpts);
      //Check convergence
      optCt += 1;
      pathDone = PathConverged(QMMMData,OldQMMMData,forceStats,optCt,
                               QMMMOpts,0,logFile);
      if (optCt == 1)
      {
        //Avoid terminating restarts on the loose tolerance step
        pathDone = 0; //Not converged
      }
    }
    BurstTraj(QMMMData,QMMMOpts);
    logFile << '\n';
    logFile << "Optimization complete.";
    logFile << '\n' << '\n';
    //Print the reaction barriers
    double dEfor = QMMMOpts.ETrans-QMMMOpts.EReact; //Forward barrier
    double dErev = QMMMOpts.ETrans-QMMMOpts.EProd; //Reverse barrier
    logFile << "NEB Results:" << '\n';
    logFile << " | Transition state bead: " << QMMMOpts.TSBead;
    logFile << '\n' << '\n';
    logFile << " | Forward barrier: ";
    logFile << LICHEMFormFloat(dEfor,16) << " eV ," << '\n';
    logFile << " |   ";
    logFile << LICHEMFormFloat(dEfor/har2eV,16) << " a.u.";
    logFile << " , ";
    logFile << LICHEMFormFloat(dEfor/kcal2eV,16) << " kcal/mol";
    logFile << '\n' << '\n';
    logFile << " | Reverse barrier: ";
    logFile << LICHEMFormFloat(dErev,16) << " eV ," << '\n';
    logFile << " |   ";
    logFile << LICHEMFormFloat(dErev/har2eV,16) << " a.u.";
    logFile << " , ";
    logFile << LICHEMFormFloat(dErev/kcal2eV,16) << " kcal/mol";
    logFile << '\n' << '\n';
    //Check stability
    if (QMMMOpts.NEBFreq)
    {
      //Calculate TS frequencies
      int remCt = 0; //Number of deleted translation and rotation modes
      int Ndof = 3*(Nqm+Npseudo); //Number of degrees of freedom
      MatrixXd QMMMHess(Ndof,Ndof);
      VectorXd QMMMFreqs(Ndof);
      logFile << '\n'; //Print blank line
      logFile << "TS frequencies:";
      logFile << '\n';
      logFile.flush(); //Print progress
      //Calculate QMMM frequencies
      QMMMHess.setZero(); //Reset Hessian
      QMMMFreqs.setZero(); //Reset frequencies
      //Calculate QM Hessian
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
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += NWChemHessian(QMMMData,QMMMOpts,QMMMOpts.TSBead);
        QMTime += (unsigned)time(0)-tStart;
      }
      //Calculate MM Hessian
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
      //Calculate frequencies
      QMMMFreqs = LICHEMFreq(QMMMData,QMMMHess,QMMMOpts,QMMMOpts.TSBead,remCt);
      //Print the frequencies
      if (remCt > 0)
      {
        logFile << "  | Identified " << remCt;
        logFile << " translation/rotation modes";
        logFile << '\n';
      }
      logFile << "  | Frequencies:" << '\n' << '\n';
      logFile << "   ";
      remCt = 0; //Reuse as a counter
      for (int i=0;i<Ndof;i++)
      {
        if (abs(QMMMFreqs(i)) > 0)
        {
          logFile << " ";
          logFile << LICHEMFormFloat(QMMMFreqs(i),10);
          remCt += 1;
          if (remCt == 3)
          {
            //Start a new line
            logFile << '\n';
            logFile << "   "; //Add extra space
            remCt = 0;
          }
        }
      }
      if (remCt != 0)
      {
        //Terminate trailing line
        logFile << '\n';
      }
      logFile << '\n';
    }
    logFile.flush();
  }
  //End of section


  //======================================================
  //START: Hatice GOKCAN
  //QSM optimization
  else if (QSMSim)
  {

    fstream ifile; //Generic file stream
    string dummy; //Generic string
    stringstream call; //Stream for system calls and reading/writing files
    int mmstep = 0;//will be used for TINKER without restraints
    int iter = 1; //for optimization steps
    int counter; //Reset counter for the number of atoms
    double gradqsm;
    bool dostep=true;
    bool PathDone = 0;
    bool QMDone = false;//will be used if only QM region
    bool before_qsm = true;//in order to compute react and prod  
 
    int optct = 0; //Counter for optimization steps
    //int macroiter=15;   
    int Nimages = QMMMOpts.NBeads;//Nimages: number of images
    int QMdim=Nqm+Npseudo; 
    int Ndof = QMdim*3; 
    int beadsize = Ndof;  
    
    int wholesize = QMMMOpts.NBeads*beadsize; //number of elements
    //QSM is frozen ends
    int PathStart = 1;
    int PathEnd = QMMMOpts.NBeads-1;
    double restr = QMMMOpts.restrConst; //default value=0.0
    double spaceout_dist=0.0; 
    
    //Initialize stats variables
    double RMSdiff = 0;
    double RMSforce = 0;
    double MAXforce = 0; 
     
    //for convergence check
    VectorXd RMSGrad = VectorXd::Zero(QMMMOpts.NBeads);
    VectorXd oldRMSGrad = VectorXd::Zero(QMMMOpts.NBeads);
    VectorXd RMSGradDiff = VectorXd::Zero(QMMMOpts.NBeads);
    VectorXd Gradconv = VectorXd::Zero(QMMMOpts.NBeads);

    //path between reactant and product (includes react and prod),
    VectorXd wholepath(wholesize); 
    if (QMMMOpts.frznEnds)
    {
      Nimages = QMMMOpts.NBeads-2;
    }

    // ENERGIES
    VectorXd E_images(Nimages+2);
    VectorXd Emm_images(Nimages+2);
    VectorXd Eqm_images(Nimages+2);
    VectorXd Eqmmm_images(Nimages+2);
    
    //FORCES AND GRADIENTS
    VectorXd Forces(Ndof); //Local forces
    VectorXd force((Nimages+2)*beadsize); //Forces of all images
    VectorXd gradient((Nimages+2)*beadsize); //gradient of all images
    //Create array to store stats and check convergence
    MatrixXd ForceStatsQM(QMMMOpts.NBeads,2);
    force.setZero();
    gradient.setZero();
    ForceStatsQM.setZero();
    
    double SavedQMOptTol = QMMMOpts.QMOptTol; //Save value from input
    double SavedMMOptTol = QMMMOpts.MMOptTol; //Save value from input
    double SavedOptTol2 =  QMMMOpts.QMRMSForceTol;
    double SavedOptTol3 = QMMMOpts.QMMaxForceTol;

    //double SavedForceTol = QMMMOpts.QMForceTol;
    //Change optimization tolerance for the first step
    QMMMOpts.QMOptTol *= 10; //Speedy convergance on the first step
    QMMMOpts.QMRMSForceTol *= 10;
    QMMMOpts.QMMaxForceTol *= 10;
    if (QMMMOpts.MMOptTol < 0.25)
    {
      QMMMOpts.MMOptTol = 0.25; //Speedy convergance on the first step
    }

    //Print initial structure
    Print_traj(QMMMData,outFile,QMMMOpts);
    //create wholepath from struct
    bool struct_to_path = true;
    updatepath(wholepath,QMMMData,QMMMOpts,
                 beadsize,Natoms,struct_to_path);
    //for spaceout distance
    VectorXd rpath(beadsize);
    VectorXd ppath(beadsize);
    VectorXd spaceout_path(beadsize);
    
    
    rpath.segment(0,beadsize)=wholepath.segment(0,beadsize);
    ppath.segment(0,beadsize)=wholepath.segment((Nimages+1)*beadsize,beadsize);
    spaceout_path=rpath-ppath;
    spaceout_dist = (spaceout_path.norm())/(10*Nimages);

    logFile << '\n';
    logFile << "   -----------------------------";
    logFile << "-----------------------------------"<< '\n'; 
    logFile << "                              ";
    logFile << "QSM OPTIMIZATION " << '\n';
    logFile << "   ---------------------------------";
    logFile << "-------------------------------"<< '\n'; 
    logFile.flush(); //Print progress
    
    logFile << '\n';
    logFile << "     Max. number of QSM Macro iterations  \n"; 
    logFile << "     is set to ";
    logFile << QMMMOpts.maxOptSteps << '\n' << endl;
 
    logFile << "     ";
    logFile << "  Points are spaced out: ";
    logFile << spaceout_dist << "\n" << endl;
    
    logFile << "     > Calculating initial ";
    logFile << "energies and forces. < " << '\n' << endl;
    logFile << '\n';
    logFile.flush(); //Print progress//open following

    CalcForces(QMMMData,QMMMOpts,Eqm_images, Emm_images,
               Eqmmm_images,force,beadsize,QMdim,before_qsm,logFile);

    /*calculate reaction coordinate*/
    VectorXd reactCoord(QMMMOpts.NBeads); //Reaction coordinate
    reactCoord.setZero();
    calc_react_coord(QMMMOpts, QMMMData,reactCoord);

    QMMMOpts.EReact=Eqmmm_images[0];
    QMMMOpts.EProd=Eqmmm_images[QMMMOpts.NBeads-1];
    getTSbead(QMMMOpts,Eqmmm_images);

    //print bead energies
    print_progress(QMMMOpts, 0,Eqmmm_images,
                   RMSdiff, MAXforce, RMSforce,reactCoord,logFile);
    
    //print TS React and Prod and barriers
    print_progress(QMMMOpts, 1,Eqmmm_images,
                   RMSdiff, MAXforce, RMSforce,reactCoord,logFile);
    
    gradient = -1*force;//convert it to gradient
    

    //Run optimization
    if(QMMMOpts.KeepFiles){
      //save initial step files
      save_files(0,0,logFile);
    }

    logFile << '\n' << endl;
    logFile << "     > Optimization Steps < " << endl; 
  
    
    while ( (!PathDone) and (iter <= QMMMOpts.maxOptSteps)) //macroiter))
    {
          logFile << "\n"; 
          logFile << "       ";
          logFile << "| Opt. step : ";
          logFile << iter;
          logFile << '\n';
          logFile.flush(); //Print progress
          
          //Copy structure
          OldQMMMData = QMMMData;

 
          if(iter==1){
              logFile << "\n";
              logFile << "         ";
              logFile << "| First optimization step.\n";
              logFile << "         ";
              logFile << "| Loose tolerance for the QM steps.\n";
              logFile << "         ";
              logFile << "| RMS Deviation = ";
              logFile << QMMMOpts.QMOptTol;
              logFile << " \u212B\n";
              logFile << "         ";
              logFile << "| RMS force     = ";
              logFile << LICHEMFormFloat(QMMMOpts.QMRMSForceTol,8);
              logFile << " Hartrees/bohr\n";
              logFile << "         ";
              logFile << "| Max. force    = ";
              logFile << LICHEMFormFloat(QMMMOpts.QMMaxForceTol,8);
              logFile << " Hartrees/bohr\n";
          }
          else{
             QMMMOpts.QMOptTol= SavedQMOptTol;
             QMMMOpts.MMOptTol = SavedMMOptTol;
             QMMMOpts.QMRMSForceTol = SavedOptTol2;
             QMMMOpts.QMMaxForceTol = SavedOptTol3;
          }   


          LICHEMQSM(QMMMData,QMMMOpts, wholepath, Nimages, QMdim, 
                    QMDone,gradient,spaceout_dist,Eqmmm_images,iter,logFile);
          //---------------------------------------------------------------------


          if(QMMM){
            //run MM optimization
            //START:restrain
            //      works only for TINKER at the moment
            //
            if(QMMMOpts.restrMM)
            {    
               //Start: do if restrain is > 2
               if (restr>=2.0)
               {
                   runRestrMMopt(QMMMData,QMMMOpts,restr,logFile);

                   restr=restr/2;/*update restr for the next iteration*/
                   PathDone=0; 
               }
               //End: do if restrain is > 2
               //Start: do if restrain is < 2
               if(restr<2.0)
               {
                   QMMMOpts.restrMM=false;//if restr<2
                   PathDone=0;
               }//Start: do if restrain is < 2
       
               if(QMMMOpts.KeepFiles and 
                  ((iter%QMMMOpts.perOpt)==0) or
                  iter==1 )
               {
                 //save MM files
                 save_files(1,iter,logFile);
               }
       
            }//END: restrain
            //START: if !QMMMOpts.restrMM
            //       QMMMOpts.restrMM became false 
            //       when restrain is < 2
            else{
               //counter for MM without restraints
               //mmstep starts from 0
               mmstep = mmstep+1; 

               runMMopt(QMMMData,QMMMOpts,logFile);

               before_qsm=false;

               PathDone=QSMConverged(QMMMData,OldQMMMData,
                               iter,QMMMOpts,Eqmmm_images,logFile);
               //print bead energies
               getTSbead(QMMMOpts,Eqmmm_images);
               print_progress(QMMMOpts, 0,Eqmmm_images,
                              RMSdiff, MAXforce, RMSforce,reactCoord,logFile);
       
               // to ensure there is at least
               // 2 mm runs without restraints
               if(mmstep<2){
                  PathDone = 0; //Not converged
               }
              
               if(QMMMOpts.KeepFiles  and 
                  (((iter%QMMMOpts.perOpt)==0) or 
                   PathDone or 
                   iter==QMMMOpts.maxOptSteps or 
                   iter==1))
               {
                 //save MM files
                 save_files(1,iter,logFile);
               }
            }//END: if !QMMMOpts.restrMM

 
         }//end: if QMMM 
       
         else{//if only QM
            PathDone = QMDone;
            if(iter==1){
              QMDone=0;
              PathDone=0;
            }  
         }//end: if only QM
         
         //Print optimized geometry
         Print_traj(QMMMData,outFile,QMMMOpts);

         /*print TS React and Prod and barriers*/
         getTSbead(QMMMOpts,Eqmmm_images);
         print_progress(QMMMOpts, 1,Eqmmm_images,
                        RMSdiff, MAXforce, RMSforce,reactCoord,logFile);

         if(QMMMOpts.KeepFiles  and 
            (((iter%QMMMOpts.perOpt)==0) or 
             PathDone or 
             iter==QMMMOpts.maxOptSteps or
             iter==1))
         {
             //save optimization step directories    
             save_files(2,iter,logFile);
         }
      
         iter = iter+1;

         if(PathDone and QMMM)
         {
             logFile << '\n';
             logFile << "               ";
             logFile << "QMMM relaxation satisfactory.";
             logFile << '\n';
         }
         
       }



      BurstTraj(QMMMData,QMMMOpts);
      logFile << '\n';
      logFile << "     > Optimization Complete <" << '\n' << endl;
      logFile << '\n' << '\n';
      logFile.flush();
  
      /*Start: Aug 28 2018 */
      if (QMMMOpts.NEBFreq)
      {
         CalcFreq(QMMMData,QMMMOpts,logFile);

         call.str("");
         call << "mkdir Freq";
         globalSys = system(call.str().c_str());

         call.str(""); 
         call << "mv LICHM_*.* ";
         call << "Freq/ ";
         globalSys = system(call.str().c_str());

         call.str("");
         call << "mv NormModes* ";
         call << "Freq/";
         globalSys = system(call.str().c_str());
      }
      /*End: Aug 28 2018 */
  
  //End of section

  }

  //END: Hatice GOKCAN
  //=======================================================

  //Inform the user if no simulations were performed
  else
  {
    logFile << "Nothing was done..." << '\n';
    logFile << "Check the simulation type in " << regFilename;
    logFile << '\n' << '\n';
    logFile.flush();
  }
  //End of section
//Start: HATICE 

    //Clean up

    if(NEBSim or QSMSim){
      //if not keepfiles, clean
      if(!QMMMOpts.KeepFiles){
        stringstream call;
        call.str("");
        call << "rm -f LICHM*"; 
        globalSys = system(call.str().c_str());
      }
    }


    if (Gaussian)
    {
      //Clear any remaining Gaussian files
      stringstream call; //Stream for system calls and reading/writing files
      call.str("");
      call << "rm -f Gau-*"; //Produced if there is a crash
      globalSys = system(call.str().c_str());
    }
    if (PSI4)
    {
      //Clear any remaining PSI4 files
      stringstream call; //Stream for system calls and reading/writing files
      call.str("");
      call << "rm -f psi*";
      globalSys = system(call.str().c_str());
    }
    if (SinglePoint or FreqCalc)
    {
      //Clear worthless output xyz file
      stringstream call; //Stream for system calls and reading/writing files
      call.str("");
      call << "rm -f ";
      for (int i=0;i<argc;i++)
      {
        //Find filename
        dummy = string(argv[i]);
        if (dummy == "-o")
        {
          call << argv[i+1];
        }
      }
      globalSys = system(call.str().c_str());
    }
    //End of section


    //Print usage statistics
    endTime = (unsigned)time(0); //Time the program completes
    double totalHours = (double(endTime)-double(startTime));
    double totalQM = double(QMTime);
    if ((QMMMOpts.NBeads > 1) and (PIMCSim or FBNEBSim))
    {
      //Average over the number of running simulations
      totalQM /= Nthreads;
    }
    double totalMM = double(MMTime);
    if ((QMMMOpts.NBeads > 1) and (PIMCSim or FBNEBSim))
    {
      //Average over the number of running simulations
      totalMM /= Nthreads;
    }
    double otherTime = totalHours-totalQM-totalMM;
    totalHours /= 3600.0; //Convert from seconds to hours
    totalQM /= 3600.0; //Convert from seconds to hours
    totalMM /= 3600.0; //Convert from seconds to hours
    otherTime /= 3600.0; //Convert from seconds to hours
    logFile << "################# Usage Statistics #################";
    logFile << '\n';
    logFile << "  Total wall time:                     ";
    logFile << LICHEMFormFloat(totalHours,6) << " hours";
    logFile << '\n';
    logFile << "  Wall time for QM Wrappers:           ";
    logFile << LICHEMFormFloat(totalQM,6) << " hours";
    logFile << '\n';
    logFile << "  Wall time for MM Wrappers:           ";
    logFile << LICHEMFormFloat(totalMM,6) << " hours";
    logFile << '\n';
    logFile << "  Wall time for LICHEM:                ";
    logFile << LICHEMFormFloat(otherTime,6) << " hours";
    logFile << '\n';
    logFile << "####################################################";
    logFile << '\n';
    logFile.flush();
    //End of section
  
    //Print a quote
    if (JOKES)
    {
        logFile << '\n';
        //logFile << "Random quote:";
        //logFile << '\n';
        //string quote; //Random quote
        //vector<string> Quotes; //Stores all possible quotes
        //FetchQuotes(Quotes); //Fetch list of quotes
        //randNum = rand() % 1000; //Randomly pick 1 of 1000 quotes
        //logFile << Quotes[randNum]; //Print quote
        //logFile << '\n';
    }
    //End of section

    //Finish output
    logFile << '\n';
    logFile << "Done.";
    logFile << '\n';
    logFile << '\n';
    logFile.flush();
    //End of section


  //Useless but supresses unused return errors for system calls
  int retValue; //= globalSys;
  //retValue = 0; //This can be changed to error messages later
  //End of section

  //Quit
  logFile.close();
  errFile.close();

  return 0;//retValue;
//};
}

