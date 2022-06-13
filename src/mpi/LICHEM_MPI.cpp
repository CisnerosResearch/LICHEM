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
  # Functions for parallel LICHEM                                             #
  # Almost same with LICHEM but;                                              #
  #                             1. Create MPI_COMM_WORLD                      #
  #                             2. Misc. initialization as in serial          #
  #                             3. Only master reads inputs                   #
  #                             4. Master shares info with workers            #
  #                             5. Perform QSM calculations                   #
  #                             5. Finalize                                   #
  #                                                                           #
  #  !!! Reading and writing to outputs are performed only by master proc     #
  #  !!! If any type of calculation except QSM, please use serial verison     #
  #############################################################################
*/

// Primary LICHEM header
#include "LICHEM_headers.h"

int main(int argc, char* argv[])
{

  // Start: MPI
  int Worldrank,Worldsize; 	/* MPI_COMM_WORLD */

  const int root=0;
  bool master=false;
  int tag, dest,source;
  int nbeads;

  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &Worldsize);  // mpirun -np all cores
                                              // (ncores*nthreads)
  MPI_Comm_rank(MPI_COMM_WORLD, &Worldrank);
  MPI_Status stat;
  MPI_Request request;

  if (Worldrank==0)
  {
    master=true;
  }

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // End: MPI_COMM_WORLD

  // Misc. initialization
  startTime = (unsigned)time(0); // Time the program starts
  srand((unsigned)time(0)); // Serial only random numbers
  // End of section

  // Initialize local variables
  string dummy; // Generic string
  double sumE,sumE2,denAvg,LxAvg,LyAvg,LzAvg,Ek; // Energies and properties
  fstream xyzFile,connectFile,regionFile,outFile,
    logFile,errFile; // Input and output files
  vector<QMMMAtom> QMMMData; // Atom list
  vector<QMMMAtom> OldQMMMData; // A copy of the atoms list
  QMMMSettings QMMMOpts; // QM and MM wrapper settings
  int randNum; // Random integer
  // End of section
  int mystat=0;

  if (master)
  {
    // Read arguments and look for errors
    ReadArgs(argc,argv,xyzFile,connectFile,regionFile,outFile,
             logFile,errFile,mystat);
    // End of section

    logFile.precision(16);
    errFile.precision(16);
    cerr.precision(16);
    if (mystat!=0)
    {
      logFile.close();
      errFile.close();
    }
  }

  MPI_Bcast(&mystat,1,MPI_INT,root,MPI_COMM_WORLD);
  if (mystat!=0)
  {
    MPI_Finalize();
    exit(0);
  }
  if (master)
  {
    // Print title and compile date
    PrintFancyTitle(logFile);
    logFile << '\n';
    logFile << "Last modification: ";
    logFile << __TIME__ << " on ";
    logFile << __DATE__ << '\n';
    logFile << '\n';
    logFile.flush();

    logFile << "Reading input..." << '\n';
    logFile << '\n';
    logFile.flush();

    // Read input and check for errors
    ReadLICHEMInput(xyzFile,connectFile,regionFile,QMMMData,QMMMOpts,
                    logFile,mystat);

    if (mystat!=0)
    {
      logFile.close();
      errFile.close();
    }

  }

  MPI_Bcast(&mystat,1,MPI_INT,root,MPI_COMM_WORLD);
  if (mystat!=0)
  {
    MPI_Finalize();
    exit(0);
  }

  if (master)
  {
    LICHEMErrorChecker(QMMMOpts,logFile,mystat);
    if (mystat!=0)
    {
      logFile.close();
      errFile.close();
    }
  }

  MPI_Bcast(&mystat,1,MPI_INT,root,MPI_COMM_WORLD);
  if (mystat!=0)
  {
    MPI_Finalize();
    exit(0);
  }

  if (master)
  {
    LICHEMPrintSettings(QMMMData,QMMMOpts,logFile);
    // End of section
    // Fix PBC
    if (PBCon)
    {
      // Relatively safe PBC correction
      if (!TINKER)
      {
        PBCCenter(QMMMData,QMMMOpts); // Center the atoms in the box
      }
    }
    // End of section

    // Create backup directories
    if (CheckFile("BACKUPQM"))
    {
      stringstream call;
      call.str("");
      // Delete old files
      call << "rm -rf " << QMMMOpts.backDir << "; ";
      // Create new directory
      call << "mkdir " << QMMMOpts.backDir;
      globalSys = system(call.str().c_str());
    }
    // End of section

    // NB: All optional simulation types should be wrapped in comments and
    //     else-if statements.
    //     The first comment should define what calculation is going to be
    //     performed, then the simulation should be enclosed in an
    //     else-if statement.
    //     After the else-if, an "// End of section" comment should be
    //     added to mark where the next simulation type begins.

    // Set keep file per step
    QMMMOpts.perOpt=QMMMOpts.maxOptSteps;
    QMMMOpts.perQM=QMMMOpts.MaxQMSteps;

  } // End if master

  // Start sharing info with workers
  Bcast_globals(root);
  MPI_Barrier(MPI_COMM_WORLD);
  Bcast_settings(QMMMOpts,root,master);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  Send_qmmmdata(QMMMData,QMMMOpts.NBeads,root,master,Natoms);
  MPI_Barrier(MPI_COMM_WORLD);

  // QSM optimization
  if (QSMSim)
  {
    /* fstream ifile; // Generic file stream */
    string dummy; // Generic string
    /* stringstream call; // Stream for system calls & reading/writing files */
    int mmstep = 0; // Will be used for TINKER without restraints
    int iter = 1; // For optimization steps
    int counter; // Reset counter for the number of atoms
    double gradqsm;
    bool dostep=true;
    bool PathDone = 0;
    bool QMDone = false; // Will be used if only QM region
    bool before_qsm = true; // In order to compute react and prod
    
    int optct = 0; // Counter for optimization steps
    /* int macroiter=15; */
    int Nimages = QMMMOpts.NBeads; // Nimages: number of images
    int QMdim=Nqm+Npseudo; 
    int Ndof = QMdim*3; 
    int beadsize = Ndof;  
    
    int wholesize = QMMMOpts.NBeads*beadsize; // Number of elements
    // QSM is frozen ends
    int PathStart = 1;
    int PathEnd = QMMMOpts.NBeads-1;
    double restr = QMMMOpts.restrConst; // Default value=0.0
    double spaceout_dist=0.0; 
    
    // Initialize stats variables
    double RMSdiff = 0;
    double RMSforce = 0;
    double MAXforce = 0; 
    
    // For convergence check
    VectorXd RMSGrad = VectorXd::Zero(QMMMOpts.NBeads);
    VectorXd oldRMSGrad = VectorXd::Zero(QMMMOpts.NBeads);
    VectorXd RMSGradDiff = VectorXd::Zero(QMMMOpts.NBeads);
    VectorXd Gradconv = VectorXd::Zero(QMMMOpts.NBeads);

    // Path between reactant and product (includes react and prod)
    VectorXd wholepath(wholesize); 
    if (QMMMOpts.frznEnds)
    {
      Nimages = QMMMOpts.NBeads-2;
    }

    // Energies
    VectorXd E_images(Nimages+2);
    VectorXd Emm_images(Nimages+2);
    VectorXd Eqm_images(Nimages+2);
    VectorXd Eqmmm_images(Nimages+2);
    
    // Forces and Gradients
    VectorXd Forces(Ndof); // Local forces
    VectorXd force((Nimages+2)*beadsize); // Forces of all images
    VectorXd gradient((Nimages+2)*beadsize); // Gradient of all images
    // Create array to store stats and check convergence
    MatrixXd ForceStatsQM(QMMMOpts.NBeads,2);
    force.setZero();
    gradient.setZero();
    ForceStatsQM.setZero();
    
    double SavedQMOptTol = QMMMOpts.QMOptTol; // Save value from input
    double SavedMMOptTol = QMMMOpts.MMOptTol; // Save value from input
    double SavedOptTol2 =  QMMMOpts.QMRMSForceTol;
    double SavedOptTol3 = QMMMOpts.QMMaxForceTol;

    /* double SavedForceTol = QMMMOpts.QMForceTol; */
    // Change optimization tolerance for the first step
    QMMMOpts.QMOptTol *= 10; // Speedy convergance on the first step
    QMMMOpts.QMRMSForceTol *= 10;
    QMMMOpts.QMMaxForceTol *= 10;
    if (QMMMOpts.MMOptTol < 0.25)
    {
      QMMMOpts.MMOptTol = 0.25; // Speedy convergance on the first step
    }
    
    // Print initial structure
    if (master)
    {
      Print_traj(QMMMData,outFile,QMMMOpts);
    }
    // Create wholepath from struct
    bool struct_to_path = true;
    updatepath(wholepath,QMMMData,QMMMOpts,
                beadsize,Natoms,struct_to_path);
    
    // For spaceout distance
    VectorXd rpath(beadsize);
    VectorXd ppath(beadsize);
    VectorXd spaceout_path(beadsize);
    
    rpath.segment(0,beadsize)=wholepath.segment(0,beadsize);
    ppath.segment(0,beadsize)=wholepath.segment((Nimages+1)*beadsize,
                                                beadsize);
    spaceout_path=rpath-ppath;
    spaceout_dist = (spaceout_path.norm())/(10*Nimages);
    /*
      spaceout_dist = ceil(((spaceout_path.norm())/(10*Nimages))* 1.0e9) /
      1.0e9;
    */

    if (master)
    { 
      logFile << '\n';
      logFile << "   -----------------------------";
      logFile << "-----------------------------------"<< '\n';
      logFile << "                              ";
      logFile << "QSM OPTIMIZATION " << '\n';
      logFile << "   ---------------------------------";
      logFile << "-------------------------------"<< '\n';
      logFile.flush(); // Print progress
      
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
      logFile.flush(); // Print progress // Open following
    }

    if (Worldrank!=0)
    {
      QMMMData.resize(Natoms);
    }

    OldQMMMData.resize(Natoms);

    CalcForcesMPI(QMMMData,QMMMOpts,Eqm_images, Emm_images,
                  Eqmmm_images,force,beadsize,QMdim,before_qsm,logFile);

    // Calculate reaction coordinate
    VectorXd reactCoord(QMMMOpts.NBeads); // Reaction coordinate
    reactCoord.setZero();
    // To tes mpi comment
    if (master)
    {
      calc_react_coord(QMMMOpts, QMMMData,reactCoord);

      QMMMOpts.EReact=Eqmmm_images[0];
      QMMMOpts.EProd=Eqmmm_images[QMMMOpts.NBeads-1];
      getTSbead(QMMMOpts,Eqmmm_images);

      // Print bead energies
      print_progress(QMMMOpts, 0,Eqmmm_images,
                    RMSdiff, MAXforce, RMSforce,reactCoord,logFile);

      // Print TS React and Prod and barriers
      print_progress(QMMMOpts, 1,Eqmmm_images,
                    RMSdiff, MAXforce, RMSforce,reactCoord,logFile);
    }
    // Send the force
    gradient = -1*force; // Convert it to gradient

    // Run optimization
    if (master)
    {
      if (QMMMOpts.KeepFiles)
      {
        // Save initial step files
        save_files(0,0,logFile);
      }
    }

    if (master)
    {
      logFile << '\n' << endl;
      logFile << "     > Optimization Steps < " << endl;
    }

    // To enter optimization together
    MPI_Barrier(MPI_COMM_WORLD);

    while ( (!PathDone) and (iter <= QMMMOpts.maxOptSteps)) /* macroiter)) */
    {
      if (master)
      {
        logFile << "\n";
        logFile << "       ";
        logFile << "| Opt. step : ";
        logFile << iter;
        logFile << '\n';
        logFile.flush(); // Print progress
      }

      // Copy structure
      // If it is just OldQMMMData = QMMMData; it gives error
      // Should use send
      if (Worldrank==0)
      {
        OldQMMMData = QMMMData;
      }

      if (iter==1)
      {
        if (master)
        {
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
      }
      else
      {
        QMMMOpts.QMOptTol= SavedQMOptTol;
        QMMMOpts.MMOptTol = SavedMMOptTol;
        QMMMOpts.QMRMSForceTol = SavedOptTol2;
        QMMMOpts.QMMaxForceTol = SavedOptTol3;
      }
      
      // Everyone needs to wait so that
      // OldQMMMData and tolerances are
      // same in every core
      MPI_Barrier(MPI_COMM_WORLD);
      
      LICHEMQSMMPI(QMMMData,QMMMOpts, wholepath, Nimages, QMdim, 
                   QMDone,gradient,spaceout_dist,Eqmmm_images,iter,logFile);

      // ---------------------------------------------------------------------

      if (QMMM)
      {
        // Run MM optimization
        // START:restrain
        //       works only for TINKER at the moment
        //
        if (QMMMOpts.restrMM)
        {    
          // Start: do if restrain is > 2
          if (restr>=2.0)
          {
            runRestrMMoptMPI(QMMMData,QMMMOpts,restr,logFile);

            MPI_Barrier(MPI_COMM_WORLD);

            restr=restr/2; // Update restr for the next iteration
            PathDone=0;
          }
          // End: do if restrain is > 2
          // Start: do if restrain is < 2
          /* else { */
          if (restr<2.0)
          {
            QMMMOpts.restrMM=false;
            PathDone=0;
            MPI_Barrier(MPI_COMM_WORLD);
          } // Start: do if restrain is < 2

          if (master)
          {
            if (QMMMOpts.KeepFiles  and
                (((iter%QMMMOpts.perOpt)==0) or
                PathDone or
                iter==QMMMOpts.maxOptSteps or
                iter==1))
            {
              // Save MM files
              save_files(1,iter,logFile);
            }
          }
        } // END: restrain
        // START: if !QMMMOpts.restrMM
        //        QMMMOpts.restrMM became false
        //        when restrain is < 2
        else
        {
          // Counter for MM without restraints
          // mmstep starts from 0
          mmstep = mmstep+1;

          runMMoptMPI(QMMMData,QMMMOpts,before_qsm,logFile);

          MPI_Barrier(MPI_COMM_WORLD);

          before_qsm=false;

          QSMConvergedMPI(QMMMData,OldQMMMData,
                          iter,QMMMOpts,Eqmmm_images,
                          PathDone,logFile);

          if (master)
          { 
            // Print bead energies
            getTSbead(QMMMOpts,Eqmmm_images);
            print_progress(QMMMOpts, 0,Eqmmm_images,
                           RMSdiff, MAXforce, RMSforce,reactCoord,logFile);
          }

          // To ensure there is at least
          // 2 mm runs without restraints
          if (mmstep<2)
          {
            PathDone = 0; // Not converged
          }
          
          if (master)
          {
            if (QMMMOpts.KeepFiles  and
                (((iter%QMMMOpts.perOpt)==0) or
                PathDone or
                iter==QMMMOpts.maxOptSteps or
                iter==1))
            {
              // Save MM files
              save_files(1,iter,logFile);
            }
          }
          MPI_Barrier(MPI_COMM_WORLD);
        } // END: if !QMMMOpts.restrMM
      } // End: if QMMM
      else
      {
        // If only QM
        PathDone = QMDone;
        if (iter==1)
        {
          QMDone=0;
          PathDone=0;
        }
      } // End: if only QM
      
      // Print optimized geometry
      if (master)
      {
        Print_traj(QMMMData,outFile,QMMMOpts);
        // Update wholepath from struct
        /*
          struct_to_path = true;
          updatepath(wholepath,QMMMData,QMMMOpts,
          beadsize,Natoms,struct_to_path);
        */
        // Print TS React and Prod and barriers
        getTSbead(QMMMOpts,Eqmmm_images);
        print_progress(QMMMOpts, 1,Eqmmm_images,
                       RMSdiff, MAXforce, RMSforce,reactCoord,logFile);
      }
      
      if (master)
      {
        if (QMMMOpts.KeepFiles  and
            (((iter%QMMMOpts.perOpt)==0) or
            PathDone or
            iter==QMMMOpts.maxOptSteps or
            iter==1))
        {
          // Save optimization step directories
          save_files(2,iter,logFile);
        }
      }

      iter = iter+1;

      if (PathDone and QMMM)
      {
        if (master)
        {
          logFile << '\n';
          logFile << "               ";
          logFile << "QMMM relaxation satisfactory.";
          logFile << '\n';
        }
      }

      MPI_Barrier(MPI_COMM_WORLD);
    }

    if (master)
    { 
      BurstTraj(QMMMData,QMMMOpts);
      logFile << '\n';
      logFile << "     > Optimization Complete <" << '\n' << endl;
      logFile << '\n' << '\n';
      logFile.flush();

      // Start: Aug 28 2018
      if (QMMMOpts.NEBFreq)
      {
        CalcFreq(QMMMData,QMMMOpts,logFile);

        stringstream call;
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
      // End: Aug 28 2018
    }

    MPI_Barrier(MPI_COMM_WORLD); // Aug 28 2018
  }
  // End of section
  // END: Hatice GOKCAN
  // =======================================================

  // Inform the user if no simulations were performed
  else
  {
    logFile << "Nothing was done..." << '\n';
    logFile << "Check the simulation type in " << regFilename;
    logFile << '\n' << '\n';
    logFile.flush();
  }
  // End of section

  // Start: HATICE
  if (master)
  {
    // Clean up

    if (NEBSim or QSMSim)
    {
      // If not keepfiles, clean
      if (!QMMMOpts.KeepFiles)
      {
        stringstream call;
        call.str("");
        call << "rm -f LICHM*";
        globalSys = system(call.str().c_str());
      }
    }

    if (Gaussian)
    {
      // Clear any remaining Gaussian files
      stringstream call; // Stream for system calls and reading/writing files
      call.str("");
      call << "rm -f Gau-*"; // Produced if there is a crash
      globalSys = system(call.str().c_str());
    }
    if (PSI4)
    {
      // Clear any remaining PSI4 files
      stringstream call; // Stream for system calls and reading/writing files
      call.str("");
      call << "rm -f psi*";
      globalSys = system(call.str().c_str());
    }
    if (SinglePoint or FreqCalc)
    {
      // Clear worthless output xyz file
      stringstream call; // Stream for system calls and reading/writing files
      call.str("");
      call << "rm -f ";
      for (int i=0;i<argc;i++)
      {
        // Find filename
        dummy = string(argv[i]);
        if (dummy == "-o")
        {
          call << argv[i+1];
        }
      }
      globalSys = system(call.str().c_str());
    }
    // End of section

    // Print usage statistics
    endTime = (unsigned)time(0); // Time the program completes
    double totalHours = (double(endTime)-double(startTime));
    double totalQM = double(QMTime);
    if ((QMMMOpts.NBeads > 1) and (PIMCSim or FBNEBSim))
    {
      // Average over the number of running simulations
      totalQM /= Nthreads;
    }
    double totalMM = double(MMTime);
    if ((QMMMOpts.NBeads > 1) and (PIMCSim or FBNEBSim))
    {
      // Average over the number of running simulations
      totalMM /= Nthreads;
    }
    double otherTime = totalHours-totalQM-totalMM;
    totalHours /= 3600.0; // Convert from seconds to hours
    totalQM /= 3600.0; // Convert from seconds to hours
    totalMM /= 3600.0; // Convert from seconds to hours
    otherTime /= 3600.0; // Convert from seconds to hours
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
    // End of section

    // Print a quote
    if (JOKES)
    {
      if (master)
      {
        logFile << '\n';
        /*
          logFile << "Random quote:";
          logFile << '\n';
          string quote; // Random quote
          vector<string> Quotes; // Stores all possible quotes
          FetchQuotes(Quotes); // Fetch list of quotes
          randNum = rand() % 1000; // Randomly pick 1 of 1000 quotes
          logFile << Quotes[randNum]; // Print quote
          logFile << '\n';
        */
      }
    }
    // End of section

    // Finish output
    logFile << '\n';
    logFile << "Done.";
    logFile << '\n';
    logFile << '\n';
    logFile.flush();
    // End of section

  }// END if master

  // Useless but supresses unused return errors for system calls
  int retValue; /* = globalSys; */
  /* retValue = 0; // This can be changed to error messages later */
  // End of section

  // Ensure that no one exits the loop before finishing
  if (Worldrank==0)
  {
    retValue=0;
    logFile.close();
    errFile.close();
    for (int i = 1; i < Worldsize; i++)
    {
      MPI_Send(&retValue, 1, MPI_INT, i, 42, MPI_COMM_WORLD);
    }
  }
  else
  {
    MPI_Recv(&retValue, 1, MPI_INT, 0, 42, MPI_COMM_WORLD, &stat);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Finalize();
  // Quit
  
  return 0; /* retValue; */
}
