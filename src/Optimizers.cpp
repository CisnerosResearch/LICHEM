/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 A set of optimization routines for the QM part of QMMM calculculations.
 MM optimizations are performed with the native optimizers in the MM
 wrappers.

 Reference for optimization routines:
 Press et al., Numerical Recipes 3nd Edition, (2007)

*/

//Convergence test functions
bool OptConverged(vector<QMMMAtom>& Struct, vector<QMMMAtom>& OldStruct,
     VectorXd& Forces, int stepct, QMMMSettings& QMMMOpts, int Bead,
     bool QMregion)
{
  //Check convergence of QMMM optimizations
  stringstream call; //Stream for system calls and reading/writing files
  string dummy; //Generic string
  call.str("");
  //Initialize stats variables
  bool OptDone = 0;
  double RMSdiff = 0;
  double RMSforce = 0;
  double MAXforce = 0;
  double SumE = 0;
  //Check progress
  if (QMregion)
  {
    //Check if a QM calculation is converged
    for (int i=0;i<(Nqm+Npseudo);i++)
    {
      //Calculate RMS forces
      double Fx = Forces(3*i);
      double Fy = Forces(3*i+1);
      double Fz = Forces(3*i+2);
      if (abs(Fx) > MAXforce)
      {
        MAXforce = abs(Fx);
      }
      if (abs(Fy) > MAXforce)
      {
        MAXforce = abs(Fy);
      }
      if (abs(Fz) > MAXforce)
      {
        MAXforce = abs(Fz);
      }
      RMSforce += Fx*Fx+Fy*Fy+Fz*Fz;
    }
    #pragma omp parallel for reduction(+:RMSdiff)
    for (int i=0;i<Natoms;i++)
    {
      //Calculate RMS displacement (QM-QM distance matrix)
      double RMStmp = 0; //Store a local sum
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        for (int j=0;j<i;j++)
        {
          if (Struct[j].QMregion or Struct[j].PBregion)
          {
            double Rnew = 0;
            double Rold = 0;
            Rnew = CoordDist2(Struct[i].P[Bead],Struct[j].P[Bead]);
            Rold = CoordDist2(OldStruct[i].P[Bead],OldStruct[j].P[Bead]);
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
    RMSdiff /= (Nqm+Npseudo)*(Nqm+Npseudo-1)/2;
    RMSdiff = sqrt(RMSdiff);
    RMSforce /= 3*(Nqm+Npseudo);
    RMSforce = sqrt(RMSforce);
    //Print progress
    call.copyfmt(cout); //Save settings
    cout << setprecision(12);
    cout << "    QM step: " << stepct;
    cout << " | RMS dev: " << RMSdiff;
    cout << " \u212B" << '\n';
    cout << "    Max. force: " << MAXforce;
    cout << " eV/\u212B | RMS force: " << RMSforce;
    cout << " eV/\u212B" << '\n';
    cout.copyfmt(call); //Return to previous settings
    //Check convergence criteria
    if ((RMSdiff <= QMMMOpts.QMOptTol) and
       (RMSforce <= (10*QMMMOpts.QMOptTol)) and
       (MAXforce <= (20*QMMMOpts.QMOptTol)))
    {
      OptDone = 1;
      cout << "    QM optimization complete." << '\n';
    }
    cout << '\n';
    cout.flush();
  }
  if (!QMregion)
  {
    //Check energy and convergence of the whole system
    SumE = 0;
    //Calculate QM energy
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      SumE += GaussianEnergy(Struct,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      SumE += PSIEnergy(Struct,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemEnergy(Struct,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      SumE += TINKEREnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      SumE += AMBEREnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Calculate RMS displacement (distance matrix)
    #pragma omp parallel for reduction(+:RMSdiff)
    for (int i=0;i<Natoms;i++)
    {
      double RMStmp = 0; //Store a local sum
      for (int j=0;j<i;j++)
      {
        double Rnew = 0;
        double Rold = 0;
        Rnew = CoordDist2(Struct[i].P[Bead],Struct[j].P[Bead]);
        Rold = CoordDist2(OldStruct[i].P[Bead],OldStruct[j].P[Bead]);
        Rnew = sqrt(Rnew);
        Rold = sqrt(Rold);
        //Update local sum
        RMStmp += (Rnew-Rold)*(Rnew-Rold);
      }
      //Update sum
      RMSdiff += RMStmp;
    }
    #pragma omp barrier
    RMSdiff /= (Natoms-Nfreeze)*(Natoms-Nfreeze-1)/2;
    RMSdiff = sqrt(RMSdiff);
    //Print progress
    call.copyfmt(cout); //Save settings
    cout << setprecision(12);
    cout << " | Opt. step: ";
    cout << stepct << " | Energy: ";
    cout << SumE << " eV ";
    cout << " | RMS dev: " << RMSdiff;
    cout << " \u212B" << '\n';
    cout.copyfmt(call); //Replace settings
    //Check convergence
    if (RMSdiff <= QMMMOpts.MMOptTol)
    {
      OptDone = 1;
      if (QMMM)
      {
        cout << "    QMMM relaxation satisfactory.";
        cout << '\n';
      }
    }
    //Flush output
    cout.flush();
  }
  return OptDone;
};

//Optimizer functions
void LICHEMSteepest(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //Cartesian steepest descent optimizer
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Save settings
  string dummy; //Generic string
  int stepct = 0; //Counter for optimization steps
  fstream qmfile, ifile, ofile; //Generic file names
  double Eold = 0; //Old saved energy
  //Initialize files
  call.str("");
  call << "QMOpt_" << Bead << ".xyz";
  qmfile.open(call.str().c_str(),ios_base::out);
  //Initialize charges for Gaussian
  if (AMOEBA and (Gaussian or NWChem))
  {
    if (TINKER)
    {
      //Set up current multipoles
      RotateTINKCharges(Struct,Bead);
    }
    //Calculate inverse box lengths for PBC
    double ix,iy,iz; //Inverse x,y,z
    ix = 1;
    iy = 1;
    iz = 1;
    if (PBCon and NWChem)
    {
      ix /= Lx;
      iy /= Ly;
      iz /= Lz;
    }
    //Save file
    call.str("");
    call << "MMCharges_" << Bead << ".txt";
    ofile.open(call.str().c_str(),ios_base::out);
    ofile.copyfmt(cout); //Save settings
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x1*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y1*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z1*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x2*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y2*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z2*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x3*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y3*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z3*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x4*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y4*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z4*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x5*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y5*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z5*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x6*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y6*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z6*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        ofile << '\n';
      }
    }
    ofile.copyfmt(cout); //Save settings
    ofile.flush();
    ofile.close();
  }
  //Initialize optimization variables
  double stepsize = 1;
  double VecMax = 0;
  bool OptDone = 0;
  vector<QMMMAtom> OldStruct = Struct; //Previous structure
  //Run optimization
  double StepScale = QMMMOpts.StepScale;
  StepScale *= 0.70; //Take a smaller first step
  while ((!OptDone) and (stepct < QMMMOpts.MaxOptSteps))
  {
    double E = 0;
    //Create blank force array
    VectorXd Forces(3*(Nqm+Npseudo));
    Forces.setZero();
    //Calculate forces (QM part)
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      E += GaussianForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      E += PSIForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      E += NWChemForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate forces (MM part)
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
      if (AMOEBA)
      {
        E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
      }
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      E += LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Check step size
    if (E > Eold)
    {
      //Take smaller steps if the energy does not improve
      cout << "    Energy did not decrease. Reducing the step size...";
      cout << '\n';
      StepScale *= 0.60; //Reduce step size
    }
    //Save structure and energy
    Eold = E;
    OldStruct = Struct;
    //Check optimization step size
    VecMax = abs(Forces.maxCoeff());
    if (abs(Forces.minCoeff()) > VecMax)
    {
      VecMax = StepScale*abs(Forces.minCoeff());
    }
    else
    {
      VecMax *= StepScale;
    }
    stepsize = StepScale;
    if (VecMax > QMMMOpts.MaxStep)
    {
      //Scale step size
      cout << "    Scaling step size to match the maximum...";
      cout << '\n';
      stepsize *= (QMMMOpts.MaxStep/VecMax);
    }
    //Determine new structure
    int ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        Struct[i].P[Bead].x += stepsize*Forces(ct);
        Struct[i].P[Bead].y += stepsize*Forces(ct+1);
        Struct[i].P[Bead].z += stepsize*Forces(ct+2);
        ct += 3;
      }
    }
    //Print structure
    Print_traj(Struct,qmfile,QMMMOpts);
    //Check convergence
    OptDone = OptConverged(Struct,OldStruct,Forces,stepct,QMMMOpts,Bead,1);
    stepct += 1;
    //Increase step size
    StepScale *= 1.05;
    if (StepScale > QMMMOpts.StepScale)
    {
      //Prevent step size from getting too large
      StepScale = QMMMOpts.StepScale;
    }
  }
  //Clean up files
  call.str("");
  call << "rm -f QMOpt_" << Bead << ".xyz";
  call << " MMCharges_" << Bead << ".txt";
  GlobalSys = system(call.str().c_str());
  //Finish and return
  return;
};

void LICHEMQuickMin(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //Cartesian damped Verlet optimizer (aka QuickMin)
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Save settings
  string dummy; //Generic string
  int stepct = 0; //Counter for optimization steps
  fstream qmfile, ifile, ofile; //Generic file names
  //Initialize files
  call.str("");
  call << "QMOpt_" << Bead << ".xyz";
  qmfile.open(call.str().c_str(),ios_base::out);
  //Initialize charges for Gaussian
  if (AMOEBA and (Gaussian or NWChem))
  {
    if (TINKER)
    {
      //Set up current multipoles
      RotateTINKCharges(Struct,Bead);
    }
    //Calculate inverse box lengths for PBC
    double ix,iy,iz; //Inverse x,y,z
    ix = 1;
    iy = 1;
    iz = 1;
    if (PBCon and NWChem)
    {
      ix /= Lx;
      iy /= Ly;
      iz /= Lz;
    }
    //Save file
    call.str("");
    call << "MMCharges_" << Bead << ".txt";
    ofile.open(call.str().c_str(),ios_base::out);
    ofile.copyfmt(cout); //Save settings
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x1*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y1*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z1*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x2*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y2*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z2*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x3*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y3*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z3*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x4*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y4*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z4*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x5*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y5*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z5*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x6*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y6*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z6*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        ofile << '\n';
      }
    }
    ofile.copyfmt(cout); //Save settings
    ofile.flush();
    ofile.close();
  }
  //Initialize optimization variables
  double timestep = 1;
  double VecMax = 0;
  bool OptDone = 0;
  vector<QMMMAtom> OldStruct = Struct; //Previous structure
  VectorXd QMVel(3*(Nqm+Npseudo)); //Velocity vector
  QMVel.setZero(); //Start at zero Kelvin
  //Run optimization
  double StepScale = QMMMOpts.StepScale; //Make a local copy
  while ((!OptDone) and (stepct < QMMMOpts.MaxOptSteps))
  {
    double E = 0;
    OldStruct = Struct; //Save old structure
    //Create blank force array
    VectorXd Forces(3*(Nqm+Npseudo));
    Forces.setZero();
    //Calculate forces (QM part)
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      E += GaussianForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      E += PSIForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      E += NWChemForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate forces (MM part)
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
      if (AMOEBA)
      {
        E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
      }
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      E += LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Project velocities
    timestep = StepScale;
    double VdotF = QMVel.dot(Forces);
    if (VdotF <= 0)
    {
      //Delete velocities and take a steepest descent step
      QMVel = timestep*Forces;
    }
    else
    {
      //Damp velocities
      QMVel = VdotF*Forces;
    }
    //Check optimization step size
    VecMax = abs(QMVel.maxCoeff());
    if (abs(QMVel.minCoeff()) > VecMax)
    {
      VecMax = StepScale*abs(QMVel.minCoeff());
    }
    else
    {
      VecMax *= StepScale;
    }
    if (VecMax > QMMMOpts.MaxStep)
    {
      //Scale timestep
      cout << "    Scaling timestep to match the maximum displacement...";
      cout << '\n';
      timestep *= (QMMMOpts.MaxStep/VecMax);
    }
    //Determine new structure
    int ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        Struct[i].P[Bead].x += timestep*QMVel(ct);
        Struct[i].P[Bead].y += timestep*QMVel(ct+1);
        Struct[i].P[Bead].z += timestep*QMVel(ct+2);
        ct += 3;
      }
    }
    //Print structure
    Print_traj(Struct,qmfile,QMMMOpts);
    //Check convergence
    OptDone = OptConverged(Struct,OldStruct,Forces,stepct,QMMMOpts,Bead,1);
    stepct += 1;
    //Update velocities
    QMVel += timestep*Forces;
  }
  //Clean up files
  call.str("");
  call << "rm -f QMOpt_" << Bead << ".xyz";
  call << " MMCharges_" << Bead << ".txt";
  GlobalSys = system(call.str().c_str());
  //Finish and return
  return;
};

void LICHEMDFP(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //A simple DFP optimizer, which is similar to BFGS updating
  //NB: This optimizer does not have a true line search, instead
  //a steepest descent step is performed if the energy rises
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Save settings
  string dummy; //Generic string
  int stepct = 0; //Counter for optimization steps
  fstream qmfile,ifile,ofile; //Generic file streams
  //Initialize trajectory file
  call.str("");
  call << "QMOpt_" << Bead << ".xyz";
  qmfile.open(call.str().c_str(),ios_base::out);
  //Initialize variables
  double E = 0; //Energy
  double Eold = 0; //Energy from previous step
  double VecMax = 0; //Maxium atomic displacement
  //Initialize multipoles for Gaussian and NWChem
  if (AMOEBA and (Gaussian or NWChem))
  {
    if (TINKER)
    {
      //Set up current multipoles
      RotateTINKCharges(Struct,Bead);
    }
    //Calculate inverse box lengths for PBC
    double ix,iy,iz; //Inverse x,y,z
    ix = 1;
    iy = 1;
    iz = 1;
    if (PBCon and NWChem)
    {
      ix /= Lx;
      iy /= Ly;
      iz /= Lz;
    }
    //Save file
    call.str("");
    call << "MMCharges_" << Bead << ".txt";
    ofile.open(call.str().c_str(),ios_base::out);
    ofile.copyfmt(cout); //Save settings
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x1*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y1*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z1*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x2*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y2*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z2*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x3*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y3*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z3*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x4*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y4*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z4*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x5*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y5*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z5*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        ofile << '\n';
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x6*ix);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y6*iy);
        ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z6*iz);
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        ofile << '\n';
      }
    }
    ofile.copyfmt(cout); //Save settings
    ofile.flush();
    ofile.close();
  }
  //Create DFP arrays
  VectorXd OptVec(3*(Nqm+Npseudo)); //Gradient descent direction
  VectorXd GradDiff(3*(Nqm+Npseudo)); //Change in the gradient
  VectorXd Forces(3*(Nqm+Npseudo)); //Forces
  MatrixXd IHess(3*(Nqm+Npseudo),3*(Nqm+Npseudo)); //Inverse Hessian
  //Initialize arrays
  OptVec.setZero();
  GradDiff.setZero();
  Forces.setZero();
  //Create an identity matrix as the initial Hessian
  IHess.setIdentity(); //Already an "inverse" Hessian
  //Calculate forces (QM part)
  if (Gaussian)
  {
    int tstart = (unsigned)time(0);
    E += GaussianForces(Struct,Forces,QMMMOpts,Bead);
    QMTime += (unsigned)time(0)-tstart;
  }
  if (PSI4)
  {
    int tstart = (unsigned)time(0);
    E += PSIForces(Struct,Forces,QMMMOpts,Bead);
    QMTime += (unsigned)time(0)-tstart;
    //Delete annoying useless files
    GlobalSys = system("rm -f psi.* timer.*");
  }
  if (NWChem)
  {
    int tstart = (unsigned)time(0);
    E += NWChemForces(Struct,Forces,QMMMOpts,Bead);
    QMTime += (unsigned)time(0)-tstart;
  }
  //Calculate forces (MM part)
  if (TINKER)
  {
    int tstart = (unsigned)time(0);
    E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
    if (AMOEBA)
    {
      E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
    }
    MMTime += (unsigned)time(0)-tstart;
  }
  if (AMBER)
  {
    int tstart = (unsigned)time(0);
    E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
    MMTime += (unsigned)time(0)-tstart;
  }
  if (LAMMPS)
  {
    int tstart = (unsigned)time(0);
    E += LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
    MMTime += (unsigned)time(0)-tstart;
  }
  //Save forces
  VecMax = 0; //Using this variable to avoid creating a new one
  VecMax = Forces.squaredNorm(); //Calculate initial RMS force
  //Output initial RMS force
  VecMax = sqrt(VecMax/(3*(Nqm+Npseudo)));
  call.copyfmt(cout); //Save settings
  cout << setprecision(12);
  cout << "    QM step: 0";
  cout << " | RMS force: " << VecMax;
  cout << " eV/\u212B";
  cout << '\n' << '\n';
  cout.flush();
  cout.copyfmt(call); //Return to previous settings
  //Optimize structure
  Eold = E;
  bool OptDone = 0;
  double StepScale = QMMMOpts.StepScale;
  StepScale *= 0.02; //Take a very small first step
  while ((!OptDone) and (stepct < QMMMOpts.MaxOptSteps))
  {
    E = 0; // Reinitialize energy
    //Copy old structure and delete old forces force array
    vector<QMMMAtom> OldStruct = Struct;
    #pragma omp parallel for
    for (int i=0;i<(3*(Nqm+Npseudo));i++)
    {
      //Reinitialize the change in the gradient
      GradDiff(i) = Forces(i);
    }
    #pragma omp barrier
    //Determine new structure
    OptVec = IHess*Forces;
    OptVec *= StepScale;
    //Check step size
    VecMax = abs(OptVec.maxCoeff());
    if (abs(OptVec.minCoeff()) > VecMax)
    {
      VecMax = abs(OptVec.minCoeff());
    }
    if (VecMax > QMMMOpts.MaxStep)
    {
      //Scale step size
      OptVec *= (QMMMOpts.MaxStep/VecMax);
    }
    //Update positions
    int ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        Struct[i].P[Bead].x += OptVec(ct);
        Struct[i].P[Bead].y += OptVec(ct+1);
        Struct[i].P[Bead].z += OptVec(ct+2);
        ct += 3;
      }
    }
    //Print structure
    Print_traj(Struct,qmfile,QMMMOpts);
    //Calculate forces (QM part)
    Forces.setZero();
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      E += GaussianForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      E += PSIForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      E += NWChemForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate forces (MM part)
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
      if (AMOEBA)
      {
        E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
      }
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      E += LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Update Hessian
    GradDiff -= Forces;
    if (((stepct%30) == 0) and (stepct != 0))
    {
      //Build a new Hessian after 30 steps
      cout << "    Constructing new Hessian...";
      cout << '\n';
      //Shrink step size
      if (StepScale > (0.02*QMMMOpts.StepScale))
      {
        StepScale = 0.02*QMMMOpts.StepScale;
      }
      else
      {
        //Reduce step size further
        StepScale *= 0.75;
      }
      //Create new Hessian as an identity matrix
      IHess.setIdentity(); //Already an "inverse" Hessian
    }
    else if (E < Eold)
    {
      //Update Hessian
      cout << "    Updating inverse Hessian...";
      cout << '\n';
      //Start really long "line"
      IHess = IHess+((OptVec*OptVec.transpose())/(OptVec.transpose()
      *GradDiff))-((IHess*GradDiff*GradDiff.transpose()*IHess)
      /(GradDiff.transpose()*IHess*GradDiff));
      //End really long "line"
      //Increase stepsize for the next iteration
      StepScale *= 1.20;
      if (StepScale > QMMMOpts.StepScale)
      {
        //Prevent step size from getting too large
        StepScale = QMMMOpts.StepScale;
      }
    }
    else
    {
      //Take a small steepest descent step and rebuild Hessian
      cout << "    Energy did not decrease. Constructing new Hessian...";
      cout << '\n';
      //Reduce step size
      if (StepScale > (0.02*QMMMOpts.StepScale))
      {
        StepScale = 0.02*QMMMOpts.StepScale;
      }
      else
      {
        //Reduce step size further
        StepScale *= 0.75;
      }
      IHess.setIdentity(); //Already an "inverse" Hessian
    }
    //Save energy
    Eold = E;
    //Check convergence
    stepct += 1;
    OptDone = OptConverged(Struct,OldStruct,Forces,stepct,QMMMOpts,Bead,1);
  }
  //Clean up files
  call.str("");
  call << "rm -f QMOpt_" << Bead << ".xyz";
  call << " MMCharges_" << Bead << ".txt";
  GlobalSys = system(call.str().c_str());
  //Finish and return
  return;
};

//Ensemble optimizers
void EnsembleSD(vector<QMMMAtom>& Struct, fstream& traj,
     QMMMSettings& QMMMOpts, int Bead)
{
  //Ensemble steepest descent optimizer
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Save settings
  string dummy; //Generic string
  int stepct = 0; //Counter for optimization steps
  fstream ifile, ofile; //Generic file names
  //Initialize optimization variables
  double stepsize = 1;
  double StepScale = QMMMOpts.StepScale; //Saved copy
  //Run optimization
  while (stepct < QMMMOpts.MaxOptSteps)
  {
    //Run MD
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      TINKERDynamics(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
      if (AMOEBA)
      {
        //Set up current multipoles
        RotateTINKCharges(Struct,Bead);
      }
    }
    //Perform SD step
    double E = 0;
    double SumE = 0;
    //Create blank force array
    VectorXd Forces(3*(Nqm+Npseudo));
    Forces.setZero();
    //Calculate forces and energy (QM part)
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      SumE += GaussianForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      SumE += PSIForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate forces and energy (MM part)
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
      SumE += TINKEREnergy(Struct,QMMMOpts,Bead);
      if (AMOEBA)
      {
        //Force from MM polarization
        E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
      }
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
      SumE += AMBEREnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      E += LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Determine new structure
    int ct = 0; //Counter for QM and PB atoms
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        //Check X step size
        stepsize = StepScale*Forces(ct);
        if (abs(stepsize) > QMMMOpts.MaxStep)
        {
          //Scale step
          stepsize *= QMMMOpts.MaxStep/abs(stepsize);
        }
        Struct[i].P[Bead].x += stepsize;
        //Check Y step size
        stepsize = StepScale*Forces(ct+1);
        if (abs(stepsize) > QMMMOpts.MaxStep)
        {
          //Scale step
          stepsize *= QMMMOpts.MaxStep/abs(stepsize);
        }
        Struct[i].P[Bead].y += stepsize;
        //Check Z step size
        stepsize = StepScale*Forces(ct+2);
        if (abs(stepsize) > QMMMOpts.MaxStep)
        {
          //Scale step
          stepsize *= QMMMOpts.MaxStep/abs(stepsize);
        }
        Struct[i].P[Bead].z += stepsize;
        ct += 3;
      }
    }
    //Print structure and energy
    stepct += 1;
    Print_traj(Struct,traj,QMMMOpts);
    call.copyfmt(cout); //Save settings
    cout << setprecision(16);
    cout << " | Step: " << stepct;
    cout << " | Simulation time: ";
    cout << (stepct*QMMMOpts.dt*QMMMOpts.Nsteps/1000);
    cout << " ps | Energy: " << SumE;
    cout << " eV" << '\n';
    cout.copyfmt(call); //Replace settings
    cout.flush();
  }
  //Clean up files and return
  call.str("");
  call << "rm -f LICHM_" << Bead << ".dyn";
  GlobalSys = system(call.str().c_str());
  return;
};

