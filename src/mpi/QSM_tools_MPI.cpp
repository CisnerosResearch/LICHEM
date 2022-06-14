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
  # Functions for reaction path optimizations in parallel                     #
  # Includes:                                                                 #
  #                                                                           #
  #       Force calculation        : double CalcForcesMPI                     #
  #                                                                           #
  #       MM Optimization without                                             #
  #       restrains                : double runMMoptMPI                       #
  #                                                                           #
  #       MM Optimization with                                                #
  #       restrains                : double runRestrMMoptMPI                  #
  #                                  (calls double TINKEROptRestr)            #
  #                                                                           #
  #       Convergence Test for QSM : bool QSMConverged                        #
  #                                                                           #
  #############################################################################
*/

double CalcForcesMPI(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                     VectorXd& Eqm_images, VectorXd& Emm_images,
                     VectorXd& Eqmmm_images, VectorXd& force, int beadsize, 
                     int QMdim, bool first_time, fstream& logFile)
{

  stringstream call;
  double SumE = 0; // Current total energy
  double relAU = 0.0;
  double relkcal = 0.0;

  int wsize,wrank;

  MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  MPI_Status stat;

  int mysize;
  bool controller=false;
  if (wrank==0)
  {
    controller=true;
  }
  // Set end points for the optimization
  int PathStart = 0;
  int PathEnd = QMMMOpts.NBeads;
  if (QMMMOpts.frznEnds)
  {
    // Change the start and end points
    PathStart = 1;
    PathEnd = QMMMOpts.NBeads-1;
  }
  if (first_time)
  {
    PathStart=0;
    PathEnd=QMMMOpts.NBeads;
  }

  // Get local bead list
  vector<int> mybead_list;

  // Create local bead lists
  // Bead info
  // root will always have 0th and 1st beads
  if (wrank==0)
  {
    if (first_time)
    {
      mybead_list.push_back(0);
    }
    mybead_list.push_back(1);
  }
  
  int owner;
  // First and second beads are on root
  /* for (int j=2; j<QMMMOpts.NBeads;j++) */
  for (int j=2; j<PathEnd;j++)
  {
    owner = (j%wsize)-1;
    // Product is on last proc
    if (j==QMMMOpts.NBeads-1)
    {
      // Last proc
      owner = wsize - 1;
    }
    if (owner==-1)
    {
      // Last proc
      owner = wsize - 1;
    }
    
    if (wrank==owner)
    {
      // Push back bead j to mybead_list
      mybead_list.push_back(j);
    }
  }// End bead info
  mysize=mybead_list.size();

  // Forces
  VectorXd Forces(QMdim*3); // Local forces

  /* force.setZero(); */

  int NewTS = 0; // Storage for new TS ID
  double NewTSEnergy = -1*hugeNum;

  // START: WRITING GAUSSIAN
  if (wrank==0)
  {
    if (Gaussian)
    {
      for (int j=PathStart; j<PathEnd;j++)
      {
        GaussianForcesMPIWrite(QMMMData,QMMMOpts,j);
      }
    }
  }// End if wrank==0

  // END: WRITING GAUSSIAN
  MPI_Barrier(MPI_COMM_WORLD);


  // START: RUN GAUSSIAN
  // Calculate forces and energies (QM part)
  if (Gaussian)
  {
    int tstart = (unsigned)time(0);
    GaussianForcesMPI(mybead_list,mysize,PathStart,PathEnd);
    QMTime += (unsigned)time(0)-tstart;

  }// End if Gaussian
  // END: RUN GAUSSIAN

  // START: READ GAUSSIAN OUT
  // Need to read Gaussian output
  //  since we need charges for polarization
  if (wrank==0)
  {

    for (int j=PathStart; j<PathEnd;j++)
    {

      Forces.setZero();
      double Eqm=0.0;
      double Emm=0.0;

      if (Gaussian)
      {
        Eqm = GaussianForcesMPIRead(QMMMData,Forces,QMMMOpts,j);
      }
      Eqm_images[j]=Eqm;

      Forces = Forces/har2eV;
      
      // Get global forces
      for (int ii=0;ii<QMdim;ii++)
      {
        
        force[(j*beadsize)+(ii*3)] = Forces(3*ii);
        force[(j*beadsize)+(ii*3)+1] = Forces(3*ii+1);
        force[(j*beadsize)+(ii*3)+2] = Forces(3*ii+2);
        
        /*
          force[(j*beadsize)+(ii*3)] = ceil(Forces(3*ii)* 1.0e9) / 1.0e9;
          force[(j*beadsize)+(ii*3)+1] = ceil(Forces(3*ii+1)* 1.0e9) / 1.0e9;
          force[(j*beadsize)+(ii*3)+2] = ceil(Forces(3*ii+2)* 1.0e9) / 1.0e9;
        */
      }
    }
  }
  // END: READ GAUSSIAN OUT

  // Charges are updated after QM
  //  send QMMData to everyone
  Send_qmmmdata(QMMMData,QMMMOpts.NBeads,0,controller,Natoms);
  MPI_Barrier(MPI_COMM_WORLD);


  // START: WRITING TINKER
  int mystat=0;
  if (wrank==0)
  {
    if (TINKER)
    {
      for (int j=PathStart; j<PathEnd;j++)
      {
        // Forces
        TINKERForcesMPIWrite(QMMMData,QMMMOpts, j,mystat,logFile);
        // Energies
        TINKEREnergyMPIWrite(QMMMData, QMMMOpts,j);

        if (AMOEBA)
        {
          // Pol Forces
          TINKERPolForcesMPIWrite(QMMMData,QMMMOpts,j,mystat,logFile);
        }
        // End if AMOEBA

        if ((AMOEBA or GEM or QMMMOpts.useImpSolv) and QMMM)
        {
          // Pol Energies
          TINKERPolEnergyMPIWrite(QMMMData, QMMMOpts,j,mystat,logFile);
        }

      }
      // End for beads

    }
    // End if TINKER

    if (mystat!=0)
    {
      logFile.close();
    }

  }
  // End if wrank==0
  // END: WRITING TINKER

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Bcast(&mystat,1,MPI_INT,0,MPI_COMM_WORLD);
  if (mystat!=0)
  {
    MPI_Finalize();
    exit(0);
  }

  // START: RUN TINKER
  // Calculate forces (MM part)
  if (TINKER)
  {

    int tstart = (unsigned)time(0);
    
    // START: FORCES
    TINKERForcesMPI(mybead_list,mysize,PathStart,PathEnd);
    if (AMOEBA)
    {
      TINKERPolForcesMPI(mybead_list,mysize,PathStart,PathEnd);
    }
    // END: FORCES

    // START: ENERGY
    TINKEREnergyMPI(mybead_list,mysize,PathStart,PathEnd);
    if ((AMOEBA or GEM or QMMMOpts.useImpSolv) and QMMM)
    {
      TINKERPolEnergyMPI(mybead_list,mysize,PathStart,PathEnd);
    }
    // END: ENERGY

    MMTime += (unsigned)time(0)-tstart;

  }
  // End if TINKER
  // END: RUN TINKER


  // START: READ TINKER OUT
  if (wrank==0)
  {
    for (int j=PathStart; j<PathEnd;j++)
    {
    
      Forces.setZero();
      double Emm=0.0;
      if (TINKER)
      {

        // Read forces
          TINKERForcesMPIRead(QMMMData, Forces, QMMMOpts,j);
          if (AMOEBA)
          {
            TINKERPolForcesMPIRead(QMMMData,Forces,QMMMOpts,j);
          }

          // Read energies
          Emm=TINKEREnergyMPIRead(QMMMData, QMMMOpts, j);
          if ((AMOEBA or GEM or QMMMOpts.useImpSolv) and QMMM)
          {
            Emm+=TINKERPolEnergyMPIRead(QMMMData, QMMMOpts, j);
          }
          Emm = (Emm*kcal2eV)/har2eV; // In a.u.
      }

      Emm_images[j]=Emm;
      Eqmmm_images[j]=Eqm_images[j]+Emm;

      if (j == QMMMOpts.TSBead)
      {
        QMMMOpts.ETrans = Eqmmm_images[j];
      }

      // Get the forces of all images from local Forces
      //  since we already read qm forces now add Forces
      //  to global force array

      Forces = Forces/har2eV;
      // Forces = Forces*bohrRad/har2eV;
      for (int ii=0;ii<QMdim;ii++)
      {
      
        force[(j*beadsize)+(ii*3)] += Forces(3*ii);
        force[(j*beadsize)+(ii*3)+1] += Forces(3*ii+1);
        force[(j*beadsize)+(ii*3)+2] += Forces(3*ii+2);
      
        /*
          force[(j*beadsize)+(ii*3)] += ceil(Forces(3*ii)* 1.0e9) / 1.0e9;
          force[(j*beadsize)+(ii*3)+1] += ceil(Forces(3*ii+1)* 1.0e9) / 1.0e9;
          force[(j*beadsize)+(ii*3)+2] += ceil(Forces(3*ii+2)* 1.0e9) / 1.0e9;
        */
      }
    }
    // End for beads
  }
  // End if wrank==0
  // START: READ TINKER OUT


  // Update QMMMData on cores
  Send_qmmmdata(QMMMData,QMMMOpts.NBeads,0,controller,Natoms);
  MPI_Barrier(MPI_COMM_WORLD);

}

/*-------------------------------------------------------------------------*/

double runMMoptMPI(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                   bool before_qsm, fstream& logFile)
{

  double SumE,SumEeV;
  // Set end points for the optimization
  int PathStart = 0;
  int PathEnd = QMMMOpts.NBeads;
  if (QMMMOpts.frznEnds)
  {
    // Change the start and end points
    PathStart = 1;
    PathEnd = QMMMOpts.NBeads-1;
  }

  int wsize,wrank;

  MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  MPI_Status stat;

  int mysize;

  // Get local bead list
  vector<int> mybead_list;

  // Create local bead lists
  // Bead info
  // root will always have 0th and 1st beads
  if (wrank==0)
  {
    if (PathStart==0)
    {
      mybead_list.push_back(0);
    }
    mybead_list.push_back(1);
  }
  
  int owner;
  // First and second beads are on root
  /* for (int j=2; j<QMMMOpts.NBeads;j++) */
  for (int j=2; j<PathEnd;j++)
  {
    owner = (j%wsize)-1;
    // Product is on last proc
    if (j==QMMMOpts.NBeads-1)
    {
      // Last processor
      owner = wsize - 1;
    }
    if (owner==-1)
    {
      // Last processor
      owner = wsize - 1;
    }
  
    if (wrank==owner)
    {
      // Push back bead j to mybead_list
      mybead_list.push_back(j);
    }
  } // End bead info
  mysize=mybead_list.size();

  int root=0;
  bool controller=false;
  if (wrank==0)
  {
    controller=true;
  }
  if (wrank==0)
  {
    logFile << "\n";
    logFile << "             ";
    logFile << "Starting MM optimization ";
    logFile << "without restraints."<< endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  // Run MM optimization
  /* for (int p=0;p<QMMMOpts.NBeads;p++) */
  if (TINKER)
  {
    int mystat=0;
    if (wrank==0)
    {
      // Write
      for (int j=PathStart; j<PathEnd;j++)
      {
        TINKEROptMPIWrite(QMMMData, QMMMOpts,j,mystat,logFile);
      }
      if (mystat!=0)
      {
        logFile.close();
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&mystat,1,MPI_INT,0,MPI_COMM_WORLD);
    if (mystat!=0)
    {
      MPI_Finalize();
      exit(0);
    }

    // Run
    int tstart = (unsigned)time(0);
    TINKEROptMPI(mybead_list,mysize,PathStart,PathEnd,QMMMOpts);
    MMTime += (unsigned)time(0)-tstart;
    
    // Read
    if (wrank==0)
    {
      // Write
      for (int j=PathStart; j<PathEnd;j++)
      {
          TINKEROptMPIRead(QMMMData, QMMMOpts, j);
      }
    }
    // END: ENERGY
  }
  
  if (wrank==0)
  {
    logFile << '\n';
    logFile << "             ";
    logFile << "MM optimization complete.";
    logFile << '\n';
    logFile.flush();
    logFile << '\n';
  }

  Send_qmmmdata(QMMMData,QMMMOpts.NBeads,0,controller,Natoms);
  MPI_Barrier(MPI_COMM_WORLD);

}

/*-------------------------------------------------------------------------*/

double runRestrMMoptMPI(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                        double restr, fstream& logFile)
{

  double SumE,SumEeV;
  // Set end points for the optimization
  int PathStart = 0;
  int PathEnd = QMMMOpts.NBeads;
  if (QMMMOpts.frznEnds)
  {
    // Change the start and end points
    PathStart = 1;
    PathEnd = QMMMOpts.NBeads-1;
  }

  int wsize,wrank;

  MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  MPI_Status stat;

  int mysize;

  // Get local bead list
  vector<int> mybead_list;

  // Create local bead lists
  // Bead info
  // root will always have 0th and 1st beads
  if (wrank==0)
  {
    if (PathStart==0)
    {
      mybead_list.push_back(0);
    }
    mybead_list.push_back(1);
  }
  
  int owner;
  // First and second beads are on root
  /* for (int j=2; j<QMMMOpts.NBeads;j++) */
  for (int j=2; j<PathEnd;j++)
  {
    owner = (j%wsize)-1;
    // Product is on last processor
    if (j==QMMMOpts.NBeads-1)
    {
      // Last proc
      owner = wsize - 1;
    }
    if (owner==-1)
    {
      // Last proc
      owner = wsize - 1;
    }

    if (wrank==owner)
    {
      // Push back bead j to mybead_list
      mybead_list.push_back(j);
    }
  }// End bead info
  mysize=mybead_list.size();

  int root=0;
  bool controller=false;
  if (wrank==0)
  {
    controller=true;
  }
  if (wrank==0)
  {
    logFile << '\n';
    logFile << "             ";
    logFile << "Starting MM optimization ";
    logFile << "with restraint positions.\n";
    logFile << "                  ";
    logFile << "Force constant = ";
    logFile << restr << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  // Run MM optimization
  /* for (int p=0;p<QMMMOpts.NBeads;p++) */
  if (TINKER)
  {
    int mystat=0;
    if (wrank==0)
    {
      // Write
      for (int j=PathStart; j<PathEnd;j++)
      {
        TINKEROptRestrainMPIWrite(QMMMData,QMMMOpts,j,restr,mystat,logFile);
      }
      if (mystat!=0)
      {
        logFile.close();
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&mystat,1,MPI_INT,0,MPI_COMM_WORLD);
    if (mystat!=0)
    {
      MPI_Finalize();
      exit(0);
    }

    // Run
    int tstart = (unsigned)time(0);
    TINKEROptMPI(mybead_list,mysize,PathStart,PathEnd,QMMMOpts);
    MMTime += (unsigned)time(0)-tstart;
    
    // Read
    if (wrank==0)
    {
      // Write
      for (int j=PathStart; j<PathEnd;j++)
      {
        TINKEROptMPIRead(QMMMData, QMMMOpts, j);
      }
    }
    // END: ENERGY
  }
  
  if (wrank==0)
  {
    logFile << '\n';
    logFile << "             ";
    logFile << "MM optimizations with ";
    logFile << "restraints are complete. \n";
    logFile << "             ";
    logFile << "Starting new QSM iteration. \n" << endl;
  }


  // Update QMMMData
  Send_qmmmdata(QMMMData,QMMMOpts.NBeads,0,controller,Natoms);
  MPI_Barrier(MPI_COMM_WORLD);

}

/*-------------------------------------------------------------------------*/

void QSMConvergedMPI(vector<QMMMAtom>& QMMMData,
                    vector<QMMMAtom>& OldQMMMData,
                    int stepct, QMMMSettings& QMMMOpts, 
                    VectorXd& Eqmmm_images,
                    bool &PathDone, fstream& logFile)
{
  
  // Check convergence of QMMM optimizations
  stringstream call; // Stream for system calls and reading/writing files
  string dummy; // Generic string
  call.str("");
  // Initialize stats variables
  /* bool PathDone = 0; */
  double RMSdiff = 0;
  double RMSforce = 0;
  double MAXforce = 0;
  double SumE = 0;

  int Ndof = 3*(Nqm+Npseudo); // Number of QM and PB degrees of freedom

  // Set end points for the optimization
  int PathStart = 0;
  int PathEnd = QMMMOpts.NBeads;
  int wsize,wrank;

  MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  MPI_Status stat;

  int mysize;

  // Set end points for the optimization
  if (QMMMOpts.frznEnds)
  {
    // Change the start and end points
    PathStart = 1;
    PathEnd = QMMMOpts.NBeads-1;
  }

  // Get local bead list
  vector<int> mybead_list;

  // Create local bead lists
  // Bead info
  // root will always have 0th and 1st beads
  if (wrank==0)
  {
    if (PathStart==0)
    {
      mybead_list.push_back(0);
    }
      mybead_list.push_back(1);
  }
  
  int owner;
  // First and second beads are on root
  /* for (int j=2; j<QMMMOpts.NBeads;j++) */
  for (int j=2; j<PathEnd;j++)
  {
    owner = (j%wsize)-1;
    // Product is on last processor
    if (j==QMMMOpts.NBeads-1)
    {
      // Last proc
      owner = wsize - 1;
    }
    if (owner==-1)
    {
      // Last proc
      owner = wsize - 1;
    }
  
    if (wrank==owner)
    {
      // Push back bead j to mybead_list
      mybead_list.push_back(j);
    }
  }// End bead info
  mysize=mybead_list.size();


  int root=0;
  bool controller=false;
  if (wrank==0)
  {
    controller=true;
  }
  VectorXd Emm(QMMMOpts.NBeads);
  VectorXd Eqm(QMMMOpts.NBeads);
  Emm.setZero();
  Eqm.setZero();

  // START: WRITING GAUSSIAN
  if (wrank==0)
  {
    if (Gaussian)
    {
      for (int j=PathStart; j<PathEnd;j++)
      {
        GaussianEnergyMPIWrite(QMMMData,QMMMOpts,j);
      }
    }
  }
  // END: WRITING GAUSSIAN
  MPI_Barrier(MPI_COMM_WORLD);

  // START: RUN GAUSSIAN
  // Calculate Energy (QM part)
  if (Gaussian)
  {
    // If global controller, write input files
    int tstart = (unsigned)time(0);
    GaussianEnergyMPI(mybead_list,mysize,PathStart,PathEnd);
    QMTime += (unsigned)time(0)-tstart;
    
  }// End if Gaussian
  // START: RUN GAUSSIAN

  // START: READ GAUSSIAN OUT
  if (wrank==0)
  {
    for (int p=PathStart;p<PathEnd;p++)
    {
      // Read energy
      if (Gaussian)
      {
        Eqm[p]=GaussianEnergyMPIRead(QMMMData,QMMMOpts,p);
      } // End if Gaussian
    }
  }
  // END: READ GAUSSIAN OUT

  // Charges are updated after QM
  // Send QMMData to everyone
  Send_qmmmdata(QMMMData,QMMMOpts.NBeads,0,controller,Natoms);
  MPI_Barrier(MPI_COMM_WORLD);

  // START: WRITING TINKER
  int mystat=0;
  if (wrank==0)
  {
    if (TINKER)
    {
      for (int j=PathStart; j<PathEnd;j++)
      {
        TINKEREnergyMPIWrite(QMMMData, QMMMOpts,j);
      }
      if ((AMOEBA or GEM or QMMMOpts.useImpSolv) and QMMM)
      {
        for (int j=PathStart; j<PathEnd;j++)
        {
          TINKERPolEnergyMPIWrite(QMMMData, QMMMOpts,j,mystat,logFile);
        }
        if (mystat!=0)
        {
          logFile.close();
        }
      }
    }
  }
  // END: WRITING TINKER
  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Bcast(&mystat,1,MPI_INT,0,MPI_COMM_WORLD);
  if (mystat!=0)
  {
    MPI_Finalize();
    exit(0);
  }

  // START: RUN TINKER
  // Calculate energy (MM part)
  if (TINKER)
  {

    int tstart = (unsigned)time(0);
    // Run
    TINKEREnergyMPI(mybead_list,mysize,PathStart,PathEnd);
    if ((AMOEBA or GEM or QMMMOpts.useImpSolv) and QMMM)
    {
      // Run
      TINKERPolEnergyMPI(mybead_list,mysize,PathStart,PathEnd);
    }
    MMTime += (unsigned)time(0)-tstart;
  } // End if TINKER
  // END: RUN TINKER

  if (wrank==0)
  {

    for (int p=PathStart;p<PathEnd;p++)
    {

      if (TINKER)
      {
        Emm[p]=TINKEREnergyMPIRead(QMMMData, QMMMOpts, p);
        if ((AMOEBA or GEM or QMMMOpts.useImpSolv) and QMMM)
        {
          Emm[p]+=TINKERPolEnergyMPIRead(QMMMData, QMMMOpts,p);
        }
      } // End if TINKER
  
      /* Eqm[p] = Eqm[p]/har2eV; // To a.u. */
      Emm[p] = (Emm[p]*kcal2eV)/har2eV; // In a.u.
      Eqmmm_images[p] = Emm[p] + Eqm[p];

      // Calculate RMS displacement
      #pragma omp parallel for schedule(dynamic) reduction(+:RMSdiff)
      for (int i=0;i<Natoms;i++)
      {
        double RMStmp = 0; // Store a local sum
        for (int j=0;j<i;j++)
        {
          double Rnew = 0;
          double Rold = 0;
          Rnew = CoordDist2(QMMMData[i].P[p],QMMMData[j].P[p]).vecMag();
          Rold = CoordDist2(OldQMMMData[i].P[p],OldQMMMData[j].P[p]).vecMag();
          Rnew = sqrt(Rnew);
          Rold = sqrt(Rold);
          // Update local sum
          RMStmp += (Rnew-Rold)*(Rnew-Rold);
        }
        // Update sum
        RMSdiff += RMStmp;
      }
      #pragma omp barrier
    }
    
    int AdjustedBeads; // Number of moving beads
    if (QMMMOpts.frznEnds)
    {
      // End points do not count
      AdjustedBeads = QMMMOpts.NBeads-2;
    }
    else
    {
      // All beads
      AdjustedBeads = QMMMOpts.NBeads;
    }
    RMSdiff /= (Natoms-Nfreeze)*(Natoms-Nfreeze-1)/2;
    RMSdiff /= AdjustedBeads; // Adjust for multiple replicas
    RMSdiff = sqrt(RMSdiff);
    // Check convergence
    if (RMSdiff <= QMMMOpts.MMOptTol)
    {
      PathDone = 1;
    }

  }

  // Flush output
  /* logFile.flush(); */

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&PathDone,1,MPI::BOOL,root,MPI_COMM_WORLD);

  // Update QMMMData
  Send_qmmmdata(QMMMData,QMMMOpts.NBeads,0,controller,Natoms);
  MPI_Barrier(MPI_COMM_WORLD);

  /* return PathDone; */
};
