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
  #############################################################################
  #                                                                           #
  # Hatice GOKCAN                                                             #
  #                                                                           #
  # Functions for Gaussian in parallel LICHEM                                 #
  # Includes:                                                                 #
  #                                                                           #
  #            for controller proc:                                           #
  #                                                                           #
  #                            Writing input files :                          #
  #                                   void GaussianForcesMPIWrite             #
  #                                   void GaussianEnergyMPIWrite             #
  #                                                                           #
  #                            Reading output files:                          #
  #                                   double GaussianForcesMPIRead            #
  #                                   double GaussianEnergyMPIRead            #
  #                                                                           #
  #            for all procs:                                                 #
  #                            Running local beads :                          #
  #                                   void GaussianForcesMPI                  #
  #                                   void GaussianEnergyMPI                  #
  #                                                                           #
  #                                                                           #
  #############################################################################
*/


void GaussianForcesMPIWrite(vector<QMMMAtom>& QMMMData,
                            QMMMSettings& QMMMOpts, int bead)
{
  // Function for calculating the forces on a set of atoms
  stringstream call; // Stream for system calls and reading/writing files
  call.copyfmt(cout); // Copy print settings
  string dummy; // Generic string
  fstream QMLog; // Generic input files

  double Eqm = 0; // QM energy
  double Eself = 0; // External field self-energy

  // Check if there is a checkpoint file
  call.str("");
  call << "LICHM_" << bead << ".chk";
  bool useCheckPoint = CheckFile(call.str());
  if (QMMMOpts.func == "SemiEmp")
  {
    // Disable checkpoints for the SemiEmp force calculations
    useCheckPoint = 0;
    // Remove SemiEmp checkpoints to avoid errors
    call.str("");
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }

  // Construct Gaussian input
  call.str("");
  if (g09)
  {
    call << "%Nprocshared=" << Ncpus << '\n';
  }
  else
  {
    call << "%CPU=0-" << Ncpus-1 << '\n';
  }
  // End: Hatice
  //
  call << "#P ";
  if (QMMMOpts.func != "SemiEmp")
  {
    // Avoids defining both a basis set and method for semi-empirical
    call << QMMMOpts.func << "/"; // Print the method
  }
  call << QMMMOpts.basis << " Force=NoStep Symmetry=None" << '\n';
  /* call << "Int=Fine SCF=Big"; // Line ended further below */
  if (useCheckPoint)
  {
    // Restart if possible
    call << " Guess=TCheck";
    call << '\n';
  }
  /* call << '\n'; */
  if (QMMM)
  {
    if ((Npseudo > 0) and (QMMMOpts.func != "SemiEmp"))
    {
      // Read pseudo potential
      call << "Pseudo=Read ";
    }
    if (Nmm > 0)
    {
      // Read charges
      call << "Charge=angstroms ";
    }
    if (QMMMOpts.func != "SemiEmp")
    {
      // Avoids calculating ESP charges for semi-empirical
      call << "Population=(MK,ReadRadii)";
    }
    call << '\n';
  }
  WriteGauInputMPI(QMMMData,call.str(),QMMMOpts,bead,1);

};

/*-------------------------------------------------------------------------*/

void GaussianEnergyMPIWrite(vector<QMMMAtom>& QMMMData,
                            QMMMSettings& QMMMOpts, int bead)
{
  // Calculates the QM energy with Gaussian
  fstream QMLog;  // Generic file streams
  string dummy;  // Generic string
  stringstream call;  // Stream for system calls and reading/writing files
  call.copyfmt(cout);  // Copy print settings
  double E = 0.0;  // QM energy
  double Eself = 0.0;  // External field self-energy
  // Check if there is a checkpoint file
  call.str("");
  call << "LICHM_" << bead << ".chk";
  bool useCheckPoint = CheckFile(call.str());
  if (QMMMOpts.func == "SemiEmp")
  {
    // Disable checkpoints for the SemiEmp force calculations
    useCheckPoint = 0;
    // Remove SemiEmp checkpoints to avoid errors
    call.str("");
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }
  // Construct Gaussian input

  call.str("");

  if (g09)
  {
    call << "%Nprocshared=" << Ncpus << '\n';
  }
  else
  {
    call << "%CPU=0-" << Ncpus-1 << '\n';
  }

  call.str("");
  call << "#P ";
  if (QMMMOpts.func != "SemiEmp")
  {
    // Avoids defining both a basis set and method for semi-empirical
    call << QMMMOpts.func << "/"; // Print the method
  }
  call << QMMMOpts.basis << " SP Symmetry=None" << '\n';
  /* call << "Int=Fine SCF=Big"; */
  if (useCheckPoint)
  {
    call << " Guess=TCheck";
    call << '\n';
  }
  /* call << '\n'; */
  if (QMMM)
  {
    if ((Npseudo > 0) and (QMMMOpts.func != "SemiEmp"))
    {
      // Read pseudo potential
      call << "Pseudo=Read ";
    }
    if (Nmm > 0)
    {
      // Read charges
      call << "Charge=angstroms ";
    }
    if (QMMMOpts.func != "SemiEmp")
    {
      // Avoids calculating ESP charges for semi-empirical
      call << "Population=(MK,ReadRadii)";
    }
    call << '\n';
  }
  WriteGauInputMPI(QMMMData,call.str(),QMMMOpts,bead,0);

}

/*-------------------------------------------------------------------------*/

double GaussianForcesMPIRead(vector<QMMMAtom>& QMMMData, VectorXd& forces,
                            QMMMSettings& QMMMOpts, int bead)
{
  // Function for calculating the forces on a set of atoms
  stringstream call; // Stream for system calls and reading/writing files
  call.copyfmt(cout); // Copy print settings
  string dummy; // Generic string
  int ct; // Generic counter
  fstream QMLog; // Generic input files

  double Eqm = 0; // QM energy
  double Eself = 0; // External field self-energy

  // Extract forces
  call.str("");
  call << "LICHM_GauForce_" << bead << ".log";
  QMLog.open(call.str().c_str(),ios_base::in);
  bool gradDone = 0;
  while ((!QMLog.eof()) and QMLog.good() and (!gradDone))
  {
    // Parse file line by line
    getline(QMLog,dummy);
    stringstream line(dummy);
    line >> dummy;
    // This only works with #P
    if (dummy == "Center")
    {
      line >> dummy >> dummy;
      if (dummy == "Forces")
      {
        gradDone = 1; // Not grad school, that lasts forever
        getline(QMLog,dummy); // Clear junk
        getline(QMLog,dummy); // Clear more junk
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
          double fX = 0;
          double fY = 0;
          double fZ = 0;
          // Extract forces; Convoluted, but "easy"
          getline(QMLog,dummy);
          stringstream line(dummy);
          line >> dummy >> dummy; // Clear junk
          line >> fX;
          line >> fY;
          line >> fZ;
          // Save forces
          forces(3*i) += fX*har2eV/bohrRad;
          forces(3*i+1) += fY*har2eV/bohrRad;
          forces(3*i+2) += fZ*har2eV/bohrRad;
        }
      }
    }
    if (dummy == "Self")
    {
      line >> dummy;
      if (dummy == "energy")
      {
        line >> dummy; // Clear junk
        line >> dummy; // Ditto
        line >> dummy; // Ditto
        line >> dummy; // Ditto
        line >> Eself; // Actual self-energy of the charges
      }
    }
    // Check for partial QMMM energy
    if (dummy == "SCF")
    {
      line >> dummy;
      if (dummy == "Done:")
      {
        line >> dummy; // Clear junk
        line >> dummy; // Ditto
        line >> Eqm; // QM energy
      }
    }
    // Check for charges
    if (dummy == "Mulliken")
    {
      // Mulliken charges (fallback)
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(QMLog,dummy); // Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
          {
            // Count through all atoms in the QM calculations
            getline(QMLog,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
    if (dummy == "ESP")
    {
      // ESP (MK) charges
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(QMLog,dummy); // Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
          {
            // Count through all atoms in the QM calculations
            getline(QMLog,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
  }
  QMLog.close();

  // Check for errors
  if (!gradDone)
  {
    cerr << "Warning: No forces recovered!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to recover...";
    cerr << '\n';
    cerr.flush(); // Print warning immediately
    // Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }

  // Clean up files
  if (!QMMMOpts.KeepFiles)
  {
    call.str("");
    call << "rm -f ";
    call << "LICHM_GauForce_" << bead;
    call << ".com ";
    call << "LICHM_GauForce_" << bead << ".log ";
    globalSys = system(call.str().c_str());
  }
  // EML add
  else
  {
    call.str("");
    call << "LICHM_GauForce_" << bead << ".com";
    ct = 0; // Start counting at the second file
    while (CheckFile(call.str()))
    {
      ct += 1; // Increase file counter
      call.str(""); // Change file name
      call << "LICHM_GauForce_bead_" << bead << "_com_";
      call << ct << ".com";
    }
    call.str("");
    // Rename com
    call << "mv LICHM_GauForce_" << bead << ".com";
    call << " LICHM_GauForce_bead_" << bead << "_com_" << ct << ".com; ";
    // Rename log
    call << " mv LICHM_GauForce_" << bead << ".log";
    call << " LICHM_GauForce_bead_" << bead << "_com_";
    call << ct << ".log;";
    globalSys = system(call.str().c_str());
  }
  // EML done
  // Change units and return
  Eqm -= Eself;

  return Eqm;
};

/*-------------------------------------------------------------------------*/

double GaussianEnergyMPIRead(vector<QMMMAtom>& QMMMData,
                            QMMMSettings& QMMMOpts, int bead)
{
  fstream QMLog; // Generic file streams
  string dummy; // Generic string
  int ct; // Generic counter
  stringstream call; // Stream for system calls and reading/writing files
  call.copyfmt(cout); // Copy print settings
  double E = 0.0; // QM energy
  double Eself = 0.0; // External field self-energy

  // Read output
  call.str("");
  call << "LICHM_GauEner_" << bead << ".log";
  QMLog.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  while ((!QMLog.eof()) and QMLog.good())
  {
    stringstream line;
    getline(QMLog,dummy);
    line.str(dummy);
    line >> dummy;
    // Search for field self-energy
    if (dummy == "Self")
    {
      line >> dummy;
      if (dummy == "energy")
      {
        line >> dummy; // Clear junk
        line >> dummy; // Ditto
        line >> dummy; // Ditto
        line >> dummy; // Ditto
        line >> Eself; // Actual self-energy of the charges
      }
    }
    // Search for energy
    if (dummy == "SCF")
    {
      line >> dummy;
      if (dummy == "Done:")
      {
        line >> dummy; // Clear junk
        line >> dummy; // Ditto
        line >> E; // QM energy
        QMFinished = 1;
      }
    }
    // Check for charges
    if (dummy == "Mulliken")
    {
      // Mulliken charges (fallback)
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(QMLog,dummy); // Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
          {
            // Count through all atoms in the QM calculations
            getline(QMLog,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
    if (dummy == "ESP")
    {
      // ESP (MK) charges
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(QMLog,dummy); // Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
          {
            // Count through all atoms in the QM calculations
            getline(QMLog,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
  }
  // Check for errors
  if (!QMFinished)
  {
    cerr << "Warning: SCF did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    E = hugeNum; // Large number to reject step
    cerr.flush(); // Print warning immediately
    /* Delete checkpoint */
    call.str("");
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }
  QMLog.close();
  /* Clean up files and save checkpoint file */
  call.str("");
  if (CheckFile("BACKUPQM"))
  {
    /* Save old files */
    call << "cp LICHM_GauEner_";
    call << bead << ".* ";
    call << QMMMOpts.backDir;
    call << "/.";
    call << " 2> LICHM_GauEner_" << bead << ".trash; ";
    call << "rm -f LICHM_GauEner_" << bead << ".trash";
    call << " "; // Extra blank space before the next command
  }

  if (!QMMMOpts.KeepFiles)
  {
    call.str("");
    call << "rm -f ";
    call << "LICHM_GauEner_" << bead;
    call << ".com ";
    call << "LICHM_GauEner_" << bead << ".log ";
    globalSys = system(call.str().c_str());
  }
  // EML add
  else
  {
    call.str("");
    call << "LICHM_GauEner_" << bead << ".com";
    ct = 0; // Start counting at the second file
    while (CheckFile(call.str()))
    {
      ct += 1; // Increase file counter
      call.str(""); // Change file name
      call << "LICHM_GauEner_bead_" << bead << "_com_";
      call << ct << ".com";
    }
    call.str("");
    // Rename com
    call << "mv LICHM_GauEner_" << bead << ".com";
    call << " LICHM_GauEner_bead_" << bead << "_com_" << ct << ".com; ";
    // Rename log
    call << " mv LICHM_GauEner_" << bead << ".log";
    call << " LICHM_GauEner_bead_" << bead << "_com_";
    call << ct << ".log;";
    globalSys = system(call.str().c_str());
  }
  // EML done

  // Change units and return
  /* cout << "QME2:" << E << " SE2:" << Eself << "\n" << endl; */
  E -= Eself;
  /* E *= har2eV; */
  return E;
};

/*-------------------------------------------------------------------------*/

void GaussianForcesMPI(vector<int> mybead_list,
                       int mysize,int pathstart, int pathend)
{
  // Function for calculating the forces on a set of atoms
  stringstream call; // Stream for system calls and reading/writing files
  call.copyfmt(cout); // Copy print settings
  string dummy; // Generic string
  fstream QMLog; // Generic input files

  int myrank,wsize;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);

  MPI_Status stat;

  int count=0;

  for (int jj=0; jj<mysize;jj++)
  {
    int p=mybead_list[jj];

    // Ensure that it starts from pathstart
    //  in root proc
    //  so that no extra calc will be performed
    if (p<pathstart)
    {
      p=pathstart;
    }
    // Ensure that it ends at pathend
    //  in last proc
    //  so that no extra calc will be performed
    if (p>pathend)
    {
      // Opt is finished, break loop
      //  and convert force to hartree
      break;
    }

    // Input file
    call.str("");
    if (g09)
    {
      call << "g09 " << "LICHM_GauForce_" << p;
    }
    else
    {
      call << "g16 " << "LICHM_GauForce_" << p;
    }
    globalSys = system(call.str().c_str());
  }

  int lastdone=0;
  if (myrank==wsize-1)
  {
    lastdone=1;
    MPI_Send(&lastdone, 0, MPI_INT, 0, 44, MPI_COMM_WORLD);
  }
  if (myrank==0)
  {
    MPI_Recv(&lastdone, 0, MPI_INT, wsize-1, 44, MPI_COMM_WORLD, &stat);
  }

  // Ensure that no one exits the loop before finishing
  int value;
  if (myrank==0)
  {
    value=1;
    for (int i = 1; i < wsize; i++)\
    {
      MPI_Send(&value, 0, MPI_INT, i, 42, MPI_COMM_WORLD);
    }

  }
  else
  {
    value=-1;
    MPI_Recv(&value, 0, MPI_INT, 0, 42, MPI_COMM_WORLD, &stat);
  }
  MPI_Barrier(MPI_COMM_WORLD);

};

/*-------------------------------------------------------------------------*/

void GaussianEnergyMPI(vector<int> mybead_list,
                       int mysize,int pathstart, int pathend)
{
  // Function for calculating the forces on a set of atoms
  stringstream call; // Stream for system calls and reading/writing files
  call.copyfmt(cout); // Copy print settings
  string dummy; // Generic string
  fstream QMLog; // Generic input files

  int myrank,wsize;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  
  MPI_Status stat;

  // Synch
  MPI_Barrier(MPI_COMM_WORLD);

  for (int jj=0; jj<mysize;jj++)
  {
    int p=mybead_list[jj];

    // Ensure that it starts from pathstart
    //  in root proc
    //  so that no extra calc will be performed 
    if (p<pathstart)
    {
      p=pathstart;
    }
    // Ensure that it ends at pathend
    //  in last proc
    //  so that no extra calc will be performed
    if (p>pathend)
    {
      // Opt is finished, break loop
      // and convert force to hartree
      break;
    }

    // Input file 
    call.str("");
    if (g09)
    {
      call << "g09 ";
    }
    else
    {
      call << "g16 ";
    }
    call << "LICHM_GauEner_" << p;

    globalSys = system(call.str().c_str());
  }

  int lastdone=0;
  if (myrank==wsize-1)
  {
    lastdone=1;
    MPI_Send(&lastdone, 0, MPI_INT, 0, 44, MPI_COMM_WORLD);
  }
  if (myrank==0)
  {
    MPI_Recv(&lastdone, 0, MPI_INT, wsize-1, 44, MPI_COMM_WORLD, &stat);
  }

  // Ensure that no one exits the loop before finishing
  int value;
  if (myrank==0)
  {
    value=1;
    for (int i = 1; i < wsize; i++)
    {
      MPI_Send(&value, 0, MPI_INT, i, 42, MPI_COMM_WORLD);
    }
  }
  else
  {
    value=-1;
    MPI_Recv(&value, 0, MPI_INT, 0, 42, MPI_COMM_WORLD, &stat);
  }
  MPI_Barrier(MPI_COMM_WORLD);

};
