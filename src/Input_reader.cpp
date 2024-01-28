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


 Routines for reading and checking the input for LICHEM.

*/

//Various input and error checking functions
void ReadArgs(int& argc, char**& argv, fstream& xyzFile,
              fstream& connectFile, fstream& regionFile, fstream& outFile,
              fstream& logFile,fstream& errFile,int& stat)
{
  //Function to read arguments
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  //fstream errFile;
  call.str("");
  call << "LICHEM" << ".err";
  errFile.open(call.str().c_str(),ios_base::out);

  //Read command line arguments
  if (argc == 1)
  {
    //Escape if there are no arguments
    logFile << "Error occured due to the missing/wrong arguments.\n";
    logFile << "Check LICHEM.err for more explanation.\n";

    errFile << '\n';
    errFile << "Missing arguments..." << '\n' << '\n';
    errFile << "Usage: lichem -n Ncpus -x Input.xyz -c Connectivity.inp ";
    errFile << "-r Regions.inp -o Output.xyz" << '\n';
    errFile << '\n';
    errFile << "Use -h or --help for detailed instructions.";
    errFile << '\n' << '\n';
    errFile.flush();
    stat=1;
    return;
  }
  dummy = string(argv[1]);
  if (dummy == "-GauExtern")
  {
    //Escape to GauExternal
    ExternalGaussian(argc,argv);
  }
  if (dummy == "-convert")
  {
    //Attempt to create LICHEM input from other formats
    dummy = string(argv[2]);
    if (dummy == "-t")
    {
      //Create LICHEM input from Gaussian input
      TINK2LICHEM(argc,argv);
    }
    if (dummy == "-b")
    {
      //Create BASIS files
      LICHEM2BASIS(argc,argv);
    }
    if (dummy == "-q")
    {
      //Create a QM connectivity file
      WriteQMConnect(argc,argv);
    }
    else
    {
      //Bad arguments
      logFile << "Error occured while executing LICHEM\n";
      logFile << "Check LICHEM.err for more explanation.\n";

      errFile << '\n';
      errFile << "Unrecognized file format.";
      errFile << '\n';
      errFile << '\n';
      errFile.flush();
    }
  }
  if (dummy == "-tinker")
  {
    //Attempt to create a TINKER XYZ file from LICHEM input
    LICHEM2TINK(argc,argv);
  }
  if (dummy == "-GlobalPoles")
  {
    //Print multipole moments in the global frame of reference
    ExtractGlobalPoles(argc,argv);
  }
  if (dummy == "-path")
  {
    //Create an initial reaction path called BeadStartStruct.xyz
    PathLinInterpolate(argc,argv);
  }
  if (dummy == "-splitpath")
  {
    //Separate a reaction path frame into a trajectory file
    SplitPathTraj(argc,argv);
  }
  if ((argc % 2) != 1)
  {
    //Check for help or missing arguments
    dummy = string(argv[1]);
    if ((dummy != "-h") and (dummy != "--help"))
    {
      //Escape if there are missing arguments
      logFile << "Error occured while executing LICHEM\n";
      logFile << "Check LICHEM.err for more explanation.\n";

      errFile << '\n';
      errFile << "Odd number of arguments..." << '\n' << '\n';
      errFile << "Usage: lichem -n Ncpus -x Input.xyz -c Connectivity.inp ";
      errFile << "-r Regions.inp -o Output.xyz" << '\n';
      errFile << '\n';
      errFile << "Use -h or --help for detailed instructions.";
      errFile << '\n' << '\n';
      errFile.flush();
      stat=1;
      return;
    }
  }
  for (int i=0;i<argc;i++)
  {
    //Read file names and CPUs
    dummy = string(argv[i]);
    if ((dummy == "-h") or (dummy == "--help"))
    {
      //Print helpful information and exit
      logFile << "Check LICHEM.err for explanation.\n";

      errFile << '\n';
      errFile << "Usage: lichem -n Ncpus -x Input.xyz -c Connectivity.inp ";
      errFile << "-r Regions.inp -o Output.xyz" << '\n';
      errFile << '\n';
      errFile << "Command line arguments:" << '\n' << '\n';
      errFile << "  -n    Number of CPUs used for the QM calculation." << '\n';
      errFile << '\n';
      errFile << "  -x    Input xyz file." << '\n' << '\n';
      errFile << "  -c    Connectivity and force field input file." << '\n';
      errFile << '\n';
      errFile << "  -r    Information about how the system is subdivided" << '\n';
      errFile << "        into QM, MM, and pseudo-atom regions." << '\n' << '\n';
      errFile << "  -o    Output xyz file for the optimized structures.";
      errFile << '\n' << '\n';
      errFile.flush();
      stat=1;
      return;
    }
    if (dummy == "-n")
    {
      //Read the number of CPUs
      Ncpus = atoi(argv[i+1]);
    }
    if (dummy == "-x")
    {
      //Read the XYZ filename
      xyzFilename = string(argv[i+1]);
      xyzFile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-c")
    {
      //Read the connectivity filename
      conFilename = string(argv[i+1]);
      connectFile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-r")
    {
      //Read the regions filename
      regFilename = string(argv[i+1]);
      regionFile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-o")
    {
      //Read the output XYZ filename
      outFile.open(argv[i+1],ios_base::out);
    }
    /*Start: Hatice GOKCAN */
    if (dummy == "-l")
    {
      //Read the output XYZ filename
      logFile.open(argv[i+1],ios_base::out);
    }
  }
  for (int i=0;i<argc;i++)
  {
    //Detect bad arguments
    dummy = string(argv[i]);
    if (dummy[0] == '-')
    {
      bool badArgs = 0; //Bad argument found
      if ((dummy != "-n") and (dummy != "-x") and
      (dummy != "-c") and (dummy != "-r") and
      (dummy != "-o") and (dummy != "-l"))
      {
        badArgs = 1;
      }
      if (badArgs)
      {
        logFile << "Error occured while executing LICHEM\n";
        logFile << "Check LICHEM.err for more explanation.\n";

        errFile << '\n';
        errFile << "Unrecognized flag..." << '\n' << '\n';
        errFile << "Usage: lichem -n Ncpus -x Input.xyz -c Connectivity.inp ";
        errFile << "-r Regions.inp -o Output.xyz" << '\n';
        errFile << '\n';
        errFile << "Use -h or --help for detailed instructions.";
        errFile << '\n' << '\n';
        errFile.flush();
        stat=1;
        return;
      }
    }
  }
  if (argc != 13)
  {
    //Escape if there are too few arguments
    logFile << "Error occured while executing LICHEM\n";
    logFile << "Check LICHEM.err for more explanation.\n";

    errFile << '\n';
    errFile << "Missing arguments..." << '\n' << '\n';
    errFile << "Usage: lichem -n Ncpus -x Input.xyz -c Connectivity.inp ";
    errFile << "-r Regions.inp -o Output.xyz -l Logfile.log" << '\n';
    errFile << '\n';
    errFile << "Use -h or --help for detailed instructions.";
    errFile << '\n' << '\n';
    errFile.flush();
    stat=1;
    return;
  }
  //Make sure input files can be read
  bool doQuit = 0;
  if (!xyzFile.good())
  {
    //Coordinate file does not exist
    logFile << "Error occured while executing LICHEM\n";
    logFile << "Check LICHEM.err for more explanation.\n";

    errFile << "Error: Could not open xyz file.";
    errFile << '\n';
    doQuit = 1;
  }
  if (!connectFile.good())
  {
    //Connectivity file does not exist
    logFile << "Error occured while executing LICHEM\n";
    logFile << "Check LICHEM.err for more explanation.\n";

    errFile << "Error: Could not open connectivity file.";
    errFile << '\n';
    doQuit = 1;
  }
  if (!regionFile.good())
  {
    //Regions file does not exist
    logFile << "Error occured while executing LICHEM\n";
    logFile << "Check LICHEM.err for more explanation.\n";

    errFile << "Error: Could not open region file.";
    errFile << '\n';
    doQuit = 1;
  }
  if (!outFile.good())
  {
    //No write permissions
    logFile << "Error occured while executing LICHEM\n";
    logFile << "Check LICHEM.err for more explanation.\n";

    errFile << "Error: Could not create output file.";
    errFile << '\n';
    doQuit = 1;
  }
  if (doQuit)
  {
    //Quit with an error
    logFile.flush(); //Print errors
    stat=1;
    return;
  }
  return;
};

void ReadLICHEMInput(fstream& xyzFile, fstream& connectFile,
                     fstream& regionFile, vector<QMMMAtom>& QMMMData,
                     QMMMSettings& QMMMOpts,fstream& logFile,int& stat)
{
  //Read input
  string dummy; //Generic string
  if (!GauExternal)
  {
    xyzFile >> Natoms;
    for (int i=0;i<Natoms;i++)
    {
      //Save atom information
      QMMMAtom tmp;
      //Set coordinates
      xyzFile >> tmp.QMTyp;
      Coord tmp2;
      xyzFile >> tmp2.x >> tmp2.y >> tmp2.z;
      tmp.P.push_back(tmp2); //Set up zeroth replica
      //Set ID and regions
      tmp.id = i;
      tmp.NEBActive = 1;
      tmp.QMRegion = 0;
      tmp.MMRegion = 1;
      tmp.PBRegion = 0;
      tmp.BARegion = 0;
      tmp.frozen = 0;
      //Set electrostatic field
      MPole tmp3; //Initialize charges and multipoles
      OctCharges tmp4; //Initialize charges and multipoles
      //Add to arrays
      tmp.MP.push_back(tmp3);
      tmp.PC.push_back(tmp4);
      //Save atomic properties
      QMMMData.push_back(tmp);
    }
  }
  for (int i=0;i<Natoms;i++)
  {
    //Save connectivity information
    int tmp;
    //id MMTyp numTyp q Nbonds [connectivity]
    connectFile >> tmp; //Atom ID
    if (tmp != QMMMData[i].id)
    {
      //Escape if connectivity errors are found
      logFile << "Error: Atoms in the connectivity file are out of order.";
      logFile << '\n';
      logFile.flush();
      stat=1;
      return;
    }
    connectFile >> QMMMData[i].MMTyp >> QMMMData[i].numTyp;
    connectFile >> QMMMData[i].m >> QMMMData[i].MP[0].q;
    connectFile >> tmp; //Number of bonds
    for (int j=0;j<tmp;j++)
    {
      //Save each bond to the atom's connectivity list
      int atomID;
      connectFile >> atomID;
      if (atomID >= Natoms)
      {
        //Search for more connectivity errors
        logFile << "Error: Atom index out of range in connectivity.";
        logFile << '\n';
        logFile << "Atom " << i << " bonded to non-existant atom ";
        logFile << atomID << '\n';
        logFile.flush();
        stat=1;
        return;
      }
      QMMMData[i].bonds.push_back(atomID); //Add bond
    }
  }
  //Read simulation keywords
  while (regionFile.good() and (!regionFile.eof()))
  {
    string keyword;
    regionFile >> keyword;
    LICHEMLowerText(keyword);
    //Check for comments
    if ((keyword[0] == '#') or (keyword[0] == '!'))
    {
      //Skip comment
    }
    //Check for simulation keywords (alphabetical)
    else if (keyword == "acceptance_ratio:")
    {
      //Read the Monte Carlo acceptance ratio
      regionFile >> QMMMOpts.accRatio;
    }
    else if (keyword == "beads:")
    {
      //Read the number of replica beads
      regionFile >> QMMMOpts.NBeads;
    }
    else if (keyword == "box_size:")
    {
      //Read the box size
      regionFile >> Lx >> Ly >> Lz;
    }
    else if (keyword == "calculation_type:")
    {
      //Set the type of calculation
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      //Single-point calculations
      if ((dummy == "single-point") or (dummy == "sp") or
         (dummy == "energy"))
      {
        //Read energy minimization options
        SinglePoint = 1;
      }
      if ((dummy == "freq") or (dummy == "frequency"))
      {
        //Read energy minimization options
        FreqCalc = 1;
      }
      //Optimizations
      if ((dummy == "opt") or (dummy == "optimize"))
      {
        //Optimize with native QM and MM optimizers
        OptSim = 1;
      }
      if ((dummy == "steep") or (dummy == "sd"))
      {
        //Optimize with the LICHEM steepest descent method
        SteepSim = 1;
      }
      if ((dummy == "dfp") or (dummy == "bfgs"))
      {
        //Optimize with the DFP optimizer
        DFPSim = 1;
        if (dummy == "bfgs")
        {
          //Print BFGS error
          logFile << "Warning: A BFGS optimizer is not implemented.";
          logFile << '\n';
          logFile << " The DFP algorithm will be used instead of BFGS.";
          logFile << '\n' << '\n';
          logFile.flush(); //Print error immediately
        }
      }
      //Reaction pathways
      if ((dummy == "neb") or (dummy == "ci-neb") or (dummy == "cineb"))
      {
        //Optimize a path with climbing image NEB
        NEBSim = 1;
      }
      //Ensemble sampling
      if (dummy == "pimc")
      {
        //Path-integral Monte Carlo
        PIMCSim = 1;
      }
      if (dummy == "fbneb")
      {
        //Force-bias Monte Carlo
        FBNEBSim = 1;
      }
      //START: HATICE
      if ((dummy == "qsm") or (dummy == "QSM") or (dummy == "Qsm"))
      {
        //Optimize a path with QSM
        QSMSim = 1;
      }
      //END: HATICE
    }
    //Start: Hatice
    else if (keyword == "nqsm:")
    {
      //number of bead in the restart file
      regionFile >> dummy;

      //all beads exist
      if(dummy=="0"){
         QMMMOpts.Nqsm = QMMMOpts.NBeads;
      }
      //only react and product exist
      else if(dummy=="1"){
         QMMMOpts.Nqsm = 2;
      }
      //wrong input
      else{
         logFile << " Error: Nqsm can be 0 or 1.\n";
         logFile << '\n';
         logFile.flush();
         //Quit
         stat=1;
         return;
      }

    }
    //End: Hatice

    else if (keyword == "electrostatics:")
    {
      //Check the type of force field
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if ((dummy == "charges") or (dummy == "charge") or
         (dummy == "point-charge"))
      {
        //Point-charge force fields
        CHRG = 1;
      }
      if (dummy == "amoeba")
      {
        //AMOEBA polarizable force field
        AMOEBA = 1;
        if (TINKER)
        {
          ExtractTINKpoles(QMMMData,0);
        }
      }
      if (dummy == "gem")
      {
        //Frozen density
        GEM = 1;
        if (TINKER)
        {
          //Collect TINKER multipoles or GEM-DM
          ExtractTINKpoles(QMMMData,0);
        }
      }
    }
    else if (keyword == "ensemble:")
    {
      //Set the thermodynamic ensemble
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if (dummy == "nvt")
      {
        //Set a consistent name for the ensemble
        QMMMOpts.ensemble = "NVT";
      }
      if (dummy == "npt")
      {
        //Set a consistent name for the ensemble
        QMMMOpts.ensemble = "NPT";
      }
    }
    else if (keyword == "eq_steps:")
    {
      //Read the number of equilibration steps
      regionFile >> QMMMOpts.NEq;
    }
    else if (keyword == "frozen_ends:")
    {
      //Check for inactive NEB end-points
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if ((dummy == "yes") or (dummy == "true"))
      {
        QMMMOpts.frznEnds = 1;
      }
    }
    else if (keyword == "init_path_chk:")
    {
      //Check for inactive NEB end-points
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if ((dummy == "no") or (dummy == "false"))
      {
        QMMMOpts.startPathChk = 0;
      }
    }
    else if (keyword == "lrec_cut:")
    {
      //Read the QMMM electrostatic cutoff for LREC
      regionFile >> QMMMOpts.LRECCut;
    }
    else if (keyword == "lrec_exponent:")
    {
      //Read the exponent for the LREC smoothing function
      regionFile >> QMMMOpts.LRECPow;
    }
    else if (keyword == "max_opt_steps:")
    {
      //Read maximum number of optimization steps
      regionFile >> QMMMOpts.maxOptSteps;
    }
    else if (keyword == "max_stepsize:")
    {
      //Read the maximum displacement during optimizations
      regionFile >> QMMMOpts.maxStep;
    }
    else if (keyword == "mm_opt_cut:")
    {
      //Read MM optimization cutoff
      regionFile >> QMMMOpts.MMOptCut;
    }
    else if (keyword == "mm_opt_tolerance:")
    {
      //Read MM optimization tolerance (RMSD value)
      regionFile >> QMMMOpts.MMOptTol;
    }
    else if (keyword == "mm_type:")
    {
      //Set MM wrapper
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if (dummy == "tinker")
      {
        TINKER = 1;
      }
      if (dummy == "lammps")
      {
        LAMMPS = 1;
      }
    }
    else if (keyword == "neb_atoms:")
    {

      if(QSMSim){
        //set all atoms to false
        for (int i=0;i<Natoms;i++){
            QMMMData[i].NEBActive = false; //false
        }      
        //Read the list of atoms to include in QSM tangents
        int numActive;
        regionFile >> numActive;

        for (int i=0;i<numActive;i++)
        {
          //Change flag
          int atomID;
          regionFile >> atomID;

          QMMMData[atomID].NEBActive = true;
        }
      }
      else{
        //Read the list of atoms to include in NEB tangents
        int numActive;
        regionFile >> numActive;
        //Temporarily mark active atoms as inactive
        for (int i=0;i<numActive;i++)
        {
          //Change flag
          int atomID;
          regionFile >> atomID;
          QMMMData[atomID].NEBActive = 0;
        }
        //Switch active and inactive groups
        #pragma omp parallel for schedule(dynamic)
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].NEBActive)
          {
            QMMMData[i].NEBActive = 0;
          }
          else
          {
            QMMMData[i].NEBActive = 1;
          }
        }
      }
    }
    else if (keyword == "opt_stepsize:")
    {
      //Read the optimization stepsize
      regionFile >> QMMMOpts.stepScale;
    }
    else if (keyword == "pbc:")
    {
      //Check for periodic boundaries
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if ((dummy == "yes") or (dummy == "true"))
      {
        PBCon = 1;
      }
    }
    //START: Hatice GOKCAN
    else if (keyword == "restrain_mm:")
    {
        //default is false
        regionFile >> dummy;
        LICHEMLowerText(dummy);
        if ((dummy == "yes") or (dummy == "true"))
        {
            QMMMOpts.restrMM = true;
        }
    }
    else if (keyword == "force_constant:")
    {
        //default is 100.0 
        regionFile >> QMMMOpts.restrConst;
    }
    else if (keyword == "qm_rms_force_tol:"){
        regionFile >> QMMMOpts.QMRMSForceTol;
    }
    else if (keyword == "qm_max_force_tol:"){
        regionFile >> QMMMOpts.QMMaxForceTol;
    }
    else if (keyword == "max_qm_steps:"){
        regionFile >> QMMMOpts.MaxQMSteps;
    }
    else if(keyword == "keep_files:"){
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if ((dummy == "yes") or (dummy == "true"))
      {
        QMMMOpts.KeepFiles=true;
      }

    }
    /* keep file per macro iter */
    else if(keyword == "per_opt_step:"){
      regionFile >> QMMMOpts.perOpt;
    }
    /* keep file per micro iter */
    else if(keyword == "per_qm_step:"){
      regionFile >> QMMMOpts.perQM;
    }
    else if(keyword == "debug:"){
        regionFile >> dummy;
        LICHEMLowerText(dummy);
        if ((dummy == "yes") or (dummy == "true"))
        {
          QMMMOpts.debug=true;
          QMMMOpts.KeepFiles=true;
        }
    }
    //END: Hatice GOKCAN
    
    else if (keyword == "potential_type:")
    {
      //Set QM, MM, and QMMM options
      regionFile >> dummy; //Potential type
      LICHEMLowerText(dummy);
      if (dummy == "qm")
      {
        //Pure QM simulation
        QMonly = 1;
        Nqm = Natoms; //Save number of QM atoms
      }
      if (dummy == "mm")
      {
        //Pure MM simulation
        MMonly = 1;
        Nmm = Natoms; //Save number of QM atoms
      }
      if (dummy == "qmmm")
      {
        //QMMM simulation
        QMMM = 1;
      }
    }
    else if (keyword == "pressure:")
    {
      //Read the pressure
      regionFile >> QMMMOpts.press;
    }
    else if (keyword == "print_normal_modes:")
    {
      //Check for inactive NEB end-points
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if ((dummy == "yes") or (dummy == "true"))
      {
        QMMMOpts.printNormModes = 1;
      }
    }
    else if (keyword == "print_steps:")
    {
      //Read the number of steps between MD and MC output
      regionFile >> QMMMOpts.NPrint;
    }
    else if (keyword == "prod_steps:")
    {
      //Read the number of production (MD or MC) steps
      regionFile >> QMMMOpts.NSteps;
    }
    else if (keyword == "qm_basis:")
    {
      //Set the basis set or semi-empirical Hamiltonian
      regionFile >> QMMMOpts.basis;
    }
    else if (keyword == "qm_charge:")
    {
      //Set the total charge on the QM region
      regionFile >> QMMMOpts.charge;
    }
    else if (keyword == "qm_memory:")
    {
      //Set the amount of memory for the QM calculations
      regionFile >> QMMMOpts.RAM;
      //Check units
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if (dummy == "mb")
      {
        //RAM is in MB
        QMMMOpts.memMB = 1;
      }
      else
      {
        //RAM is in GB
        QMMMOpts.memMB = 0;
      }
    }
    else if (keyword == "qm_method:")
    {
      //Set QM functional or method
      regionFile >> dummy;
      QMMMOpts.func = dummy; //Save name with correct case
      //Check for special methods
      LICHEMLowerText(dummy);
      if ((dummy == "semiempirical") or (dummy == "se-scf") or
         (dummy == "semi-empirical") or (dummy == "sescf") or
         (dummy == "semiemp"))
      {
        //Flag the method as a semi-empirical Hamiltonian
        QMMMOpts.func = "SemiEmp";
      }
    }
    else if (keyword == "qm_opt_tolerance:")
    {
      //Read QM optimization tolerance (RMSD value)
      regionFile >> QMMMOpts.QMOptTol;
    }
    else if (keyword == "qm_spin:")
    {
      //Set the multiplicity
      regionFile >> QMMMOpts.spin;
    }
    else if (keyword == "qm_type:")
    {
      //Set QM wrapper
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if (dummy == "psi4")
      {
        PSI4 = 1;
      }
      if (dummy == "nwchem")
      {
        NWChem = 1;
      }
/* Start: Hatice */
      //if ((dummy == "gaussian") or (dummy == "g09"))
      //{
      //  Gaussian = 1;
      //}
      if ((dummy == "gaussian") or (dummy == "g09")){
          Gaussian = 1;
          g09 = 1;
      }
      if (dummy == "g16"){
          Gaussian = 1;
      }

/* End: Hatice */
    }
    else if (keyword == "qm_units:")
    {
      //Read distance units for QM calculations
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if ((dummy == "bohr") or (dummy == "a.u."))
      {
        //Change distance units to a.u.
        QMMMOpts.unitsQM = "Bohr";
      }
    }
    else if (keyword == "solv_model:")
    {
      //Read MM implicit solvent model
      regionFile >> QMMMOpts.solvModel;
    }
    else if (keyword == "spring_constant:")
    {
      //Read the NEB spring constant
      regionFile >> QMMMOpts.kSpring;
    }
    else if (keyword == "temperature:")
    {
      //Read the temperature
      regionFile >> QMMMOpts.temp;
      //Save the inverse temperature
      QMMMOpts.beta = 1/(kBoltz*QMMMOpts.temp);
    }
    else if (keyword == "ts_freq:")
    {
      //Check for inactive NEB end-points
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if ((dummy == "yes") or (dummy == "true"))
      {
        QMMMOpts.NEBFreq = 1;
      }
    }
    else if (keyword == "use_ewald:")
    {
      //Check for MM Ewald summation
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if ((dummy == "yes") or (dummy == "true"))
      {
        //Turn on Ewald or PME
        QMMMOpts.useEwald = 1;
      }
    }
    else if (keyword == "use_lrec:")
    {
      //Turn on long-range corrections
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if ((dummy == "yes") or (dummy == "true"))
      {
        //Turn on long-range corrections
        QMMMOpts.useLREC = 1;
      }
    }
    else if (keyword == "use_mm_cutoff:")
    {
      //Check for the MM optimization cutoff
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if ((dummy == "yes") or (dummy == "true"))
      {
        //Turn on the optimization cutoff
        QMMMOpts.useMMCut = 1;
      }
    }
    else if (keyword == "use_solvent:")
    {
      //Check for MM implicit solvation
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if ((dummy == "yes") or (dummy == "true"))
      {
        //Turn on the implicit solvent
        QMMMOpts.useImpSolv = 1;
      }
    }
    //Check for region keywords
    else if (keyword == "qm_atoms:")
    {
      //Read the list of QM atoms
      regionFile >> Nqm;
      for (int i=0;i<Nqm;i++)
      {
        int atomID;
        regionFile >> atomID;
        QMMMData[atomID].QMRegion = 1;
        QMMMData[atomID].PBRegion = 0;
        QMMMData[atomID].BARegion = 0;
        QMMMData[atomID].MMRegion = 0;
      }
    }
    else if (keyword == "pseudobond_atoms:")
    {
      //Read the list of pseudobond atoms
      regionFile >> Npseudo;
      for (int i=0;i<Npseudo;i++)
      {
        int atomID;
        regionFile >> atomID;
        QMMMData[atomID].QMRegion = 0;
        QMMMData[atomID].PBRegion = 1;
        QMMMData[atomID].BARegion = 0;
        QMMMData[atomID].MMRegion = 0;
      }
    }
    else if (keyword == "boundary_atoms:")
    {
      //Read the list of boundary atoms
      regionFile >> Nbound;
      for (int i=0;i<Nbound;i++)
      {
        int atomID;
        regionFile >> atomID;
        QMMMData[atomID].QMRegion = 0;
        QMMMData[atomID].PBRegion = 0;
        QMMMData[atomID].BARegion = 1;
        QMMMData[atomID].MMRegion = 0;
      }
    }
    else if (keyword == "frozen_atoms:")
    {
      //Read the list of frozen atoms
      regionFile >> Nfreeze;
      for (int i=0;i<Nfreeze;i++)
      {
        int atomID;
        regionFile >> atomID;
        QMMMData[atomID].frozen = 1;
      }
    }
    //S:JORGE
    else if (keyword == "dispersion:")
    {
      //Adding vdW option
      regionFile >> dummy;
      QMMMOpts.dispersion = dummy;
      QMMMOpts.dispbool = 1;
    }
    else if (keyword == "gem_basis:")
    {
      //Set GEM basis
      regionFile >> dummy;
      QMMMOpts.gembasis = dummy;
    }
    else if (keyword == "gem_kexchange:")
    {
      regionFile >> QMMMOpts.kexchange;
    }
    else if (keyword == "gem_prefitted:")
    {
      //Check for prefitted coefficients
      regionFile >> dummy;
      LICHEMLowerText(dummy);
      if (dummy == "no") 
      {
        //Turn on prefitted coefficients
        QMMMOpts.prefitted = 0;
      }
    }
    //E:JORGE
    //Check for bad keywords
    else if (regionFile.good() and (!regionFile.eof()))
    {
      //Inform the user about the bad keyword
      logFile << "Error: Unrecognized keyword: ";
      logFile << keyword << '\n';
      logFile.flush();
      //Quit
      stat=1;
      return;
    }
  }
  //Reset regions for pure QM and MM
  if (QMonly)
  {
    //Reset the numbers if regions were specified in the input
    Nqm = Natoms;
    Npseudo = 0;
    Nbound = 0;
    //Redundant, but safe
    for (int i=0;i<Natoms;i++)
    {
      QMMMData[i].QMRegion = 1;
      QMMMData[i].MMRegion = 0;
      QMMMData[i].PBRegion = 0;
      QMMMData[i].BARegion = 0;
    }
    //Adjust optimization settings
    QMMMOpts.MMOptTol = QMMMOpts.QMOptTol; //Prevents early termination
  }
  if (MMonly)
  {
    //Reset the numbers if regions were specified in the input
    Nqm = 0;
    Npseudo = 0;
    Nbound = 0;
    //Redundant, but safe
    for (int i=0;i<Natoms;i++)
    {
      QMMMData[i].QMRegion = 0;
      QMMMData[i].MMRegion = 1;
      QMMMData[i].PBRegion = 0;
      QMMMData[i].BARegion = 0;
    }
  }
  Nmm = Natoms-Nqm-Npseudo-Nbound; //Set number of MM atoms
  //Replicate atoms
  if (QMMMOpts.NBeads > 1)
  {
    //Duplicate data
    for (int i=0;i<Natoms;i++)
    {
      //Create reaction-path beads
      for (int j=0;j<(QMMMOpts.NBeads-1);j++)
      {
        //Create replicas
        Coord temp = QMMMData[i].P[0];
        QMMMData[i].P.push_back(temp);
        MPole temp2 = QMMMData[i].MP[0];
        QMMMData[i].MP.push_back(temp2);
        OctCharges temp3 = QMMMData[i].PC[0];
        QMMMData[i].PC.push_back(temp3);
      }
    }
    //Set initial transition state for reaction pathways
    //Start: Hatice
    //if (NEBSim) 
    if (NEBSim or QSMSim)
    //End: Hatice
    {
      if ((QMMMOpts.NBeads%2) == 0)
      {
        //Even number of beads
        QMMMOpts.TSBead = (QMMMOpts.NBeads/2); //Slightly on the product side
      }
      else
      {
        //Odd number of beads
        QMMMOpts.TSBead = ((QMMMOpts.NBeads-1)/2); //Middle bead
      }
    }
    //Add random displacements for PIMC simulations
    if (PIMCSim)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Shift path-integral beads
        double massScale = sqrt(12.0/QMMMData[i].m); //Relative to carbon
        massScale *= 2*stepMin*centRatio; //Scale based on settings
        //Update all beads
        for (int j=0;j<(QMMMOpts.NBeads-1);j++)
        {
          //Pick random displacements
          double randX = (((double)rand())/((double)RAND_MAX));
          double randY = (((double)rand())/((double)RAND_MAX));
          double randZ = (((double)rand())/((double)RAND_MAX));
          //Place the first bead at the initial position
          if (j == 0)
          {
            randX = 0.5;
            randY = 0.5;
            randZ = 0.5;
          }
          //Update positions of active atoms
          if (!QMMMData[i].frozen)
          {
            QMMMData[i].P[j].x += (2*(randX-0.5)*massScale);
            QMMMData[i].P[j].y += (2*(randY-0.5)*massScale);
            QMMMData[i].P[j].z += (2*(randZ-0.5)*massScale);
          }
        }
      }
    }
  }
  //START: Hatice GOKCAN 
  //For QSM 
  //Read initial structures for all beads or create new ones if QSM simulation
  //if (CheckFile("QSMBeadStruct.xyz") and (!GauExternal))
  if (CheckFile("BeadStartStruct.xyz") and (!GauExternal))
  {
    //Print output
    logFile << "Reading restart information...";
    logFile << '\n' << '\n';;
    //Open file
    fstream beadfile;
    //beadfile.open("QSMBeadStruct.xyz",ios_base::in);
    beadfile.open("BeadStartStruct.xyz",ios_base::in);

    if(QSMSim){ 
      //Read and discard number of atoms
      int AtTest = 0;
      beadfile >> AtTest;
      if (AtTest != (Natoms*QMMMOpts.Nqsm))
      {
        //Print warning if the XYZ file has incorrect dimensions
        logFile << "Error: Restart file does not have the correct format!";
        logFile << '\n' << '\n';
        logFile.flush();
        //Quit
        stat=1;
        return;
      }
      else
      {
        logFile << "Restart file contains ";
        logFile << QMMMOpts.Nqsm;
        logFile << " beads.\n";
        logFile.flush();
      }

      //Read atom/bead positions
      //if Nqsm=2 then there is reactant and product
      if(QMMMOpts.Nqsm==2){
        //Read reactant XYZ coordinates
        int p=0; //reactant bead
        for (int i=0;i<Natoms;i++)
        {
            //Read atom type and discard
            beadfile >> dummy;
            //Read XYZ coordinates
            beadfile >> QMMMData[i].P[p].x;
            beadfile >> QMMMData[i].P[p].y;
            beadfile >> QMMMData[i].P[p].z;
        }
        //Read product XYZ coordinates
        p = QMMMOpts.NBeads -1; //product bead
        for (int i=0;i<Natoms;i++)
        {
            //Read atom type and discard
            beadfile >> dummy;
            //Read XYZ coordinates
            beadfile >> QMMMData[i].P[p].x;
            beadfile >> QMMMData[i].P[p].y;
            beadfile >> QMMMData[i].P[p].z;
        }
      }
     
      //if Nqsm=QMMMOpts.NBeads then whole path is present
      if(QMMMOpts.Nqsm==QMMMOpts.NBeads){
        //Read XYZ coordinates for beads
        for(int p=0;p<QMMMOpts.NBeads;p++){
          for (int i=0;i<Natoms;i++)
          {
            //Read atom type and discard
            beadfile >> dummy;
            //Read XYZ coordinates
            beadfile >> QMMMData[i].P[p].x;
            beadfile >> QMMMData[i].P[p].y;
            beadfile >> QMMMData[i].P[p].z;
          }
        }
      }
    }//endif QSMSim
    else{ //NEB etc.
      //Read and discard number of atoms
      int atTest = 0;
      beadfile >> atTest;
      if (atTest != (Natoms*QMMMOpts.NBeads))
      {
        //Print warning if the XYZ file has incorrect dimensions
        logFile << "Error: Restart file does not have the correct format!";
        logFile << '\n' << '\n';
        logFile.flush();
        //Quit
        stat=1;
        return;
      }
      //Read atom/bead positions
      //for (int i=0;i<Natoms;i++)
      for (int j=0;j<QMMMOpts.NBeads;j++)
      {
        //for (int j=0;j<QMMMOpts.NBeads;j++)
        for (int i=0;i<Natoms;i++)
        {
          //Read atom type and discard
          beadfile >> dummy;
          //Read XYZ coordinates
          beadfile >> QMMMData[i].P[j].x;
          beadfile >> QMMMData[i].P[j].y;
          beadfile >> QMMMData[i].P[j].z;
        }
      }
    } 
  }
  else if (NEBSim or QSMSim)
  {
    //Exit with an error
    logFile << "Error: No initial reaction path found in the restart file!!!";
    logFile << '\n' << '\n';
    logFile.flush();
    //Quit
    stat=1;
    return;

  }
  //END: Hatice GOKCAN

  //Collect additonal TINKER input
  if (TINKER and (!GauExternal))
  {
    //NB: Classes are not used in the QMMM
    FindTINKERClasses(QMMMData,logFile); //Finds errors
  }
  //Check if QM log files should be saved
  if (CheckFile("BACKUPQM"))
  {
    //Read backup directory
    fstream backFile;
    //Set to default value
    QMMMOpts.backDir = "Old_files";
    //Check directory
    backFile.open("BACKUPQM",ios_base::in);
    if (backFile.good())
    {
      string newName;
      backFile >> newName;
      if (!backFile.eof())
      {
        QMMMOpts.backDir = newName;
      }
    }
  }
  //Set threads based on QM CPUs and total CPUs
  if (!GauExternal)
  {
    //NB: Sanity checks and error checking are only enabled with OpenMP
    //Set default number of threads for serial builds
    Nthreads = 1;
    //Set a better more realistic number of threads for OpenMP
    #ifdef _OPENMP
      //OpenMP settings
      double Procs = double(FindMaxThreads());
      Nthreads = FindMaxThreads();
      omp_set_num_threads(Nthreads);
      //Sanity check
      if (Ncpus > Nthreads)
      {
        //Assuming only one node is used for QM
        Ncpus = Nthreads;
      }
      //Modify threads for certain multi-replica simulations
      if ((QMMMOpts.NBeads > 1) and (PIMCSim or FBNEBSim))
      {
        //Divide threads between the beads
        Nthreads = int(floor(Procs/Ncpus));
        //Set number of threads for wrappers
        omp_set_num_threads(Nthreads);
      }
    #endif
    //Set eigen threads
    setNbThreads(Nthreads);
  }
  return;
};

void LICHEMErrorChecker(QMMMSettings& QMMMOpts,fstream& logFile,int& stat)
{
  //Checks for basic errors and conflicts
  bool doQuit = 0; //Bool, quit with error
  //General errors
  if (QMMM)
  {
    //Check number of QM and MM atoms
    if ((Nqm+Npseudo) < 1)
    {
      //Make sure there are some atoms in the QM calculation
      logFile << " Error: No QM or PB atoms defined for the QMMM calculations.";
      logFile << '\n';
      doQuit = 1;
    }
    if ((Nmm+Nbound) < 1)
    {
      //Make sure there are some atoms in the MM calculations
      logFile << " Error: No MM or BA atoms defined for the QMMM calculations.";
      logFile << '\n';
      doQuit = 1;
    }
  }

  //Check LREC settings
  if (QMMMOpts.useLREC or PBCon)
  {
    //Check LREC cutoff
    if (PBCon)
    {
      //Find maximum box length
      double minLen = Lx;
      if (Ly < minLen)
      {
        minLen = Ly;
      }
      if (Lz < minLen)
      {
        minLen = Lz;
      }
      //Check cutoff
      if (QMMMOpts.useLREC and (QMMMOpts.LRECCut > (0.5*minLen)))
      {
        //Needed to make the minimum image convention safe
        QMMMOpts.LRECCut = 0.5*minLen;
        logFile << "Warning: Reducing LREC cutoff (";
        logFile << LICHEMFormFloat(QMMMOpts.LRECCut,6);
        logFile << ") due to the minimum image convention.";
        logFile << '\n' << '\n';
      }
    }
    if (QMMMOpts.useLREC and (QMMMOpts.LRECCut <= 0.10))
    {
      //Adjust cutoff to avoid divide by zero errors
      QMMMOpts.LRECCut = 0.10; //Minimum value, effectively zero
      logFile << "Warning: LREC cutoffs less than 0.1 are not allowed.";
      logFile << '\n' << '\n';
    }
    //Check LREC exponent
    if (QMMMOpts.LRECPow < 1)
    {
      //Needed to make the minimum image convention safe
      QMMMOpts.LRECPow = 3;
      logFile << "Warning: Invalid LREC exponent.";
      logFile << " LREC exponent set to 3.";
      logFile << '\n' << '\n';
    }
  }
  //Check Ewald and implicit solvation settings
  if (QMMMOpts.useEwald and (!PBCon))
  {
    //Check Ewald settings
    logFile << " Error: Ewald summation cannot be used without PBC.";
    logFile << '\n';
    doQuit = 1;
  }
  if (QMMMOpts.useImpSolv and PBCon)
  {
    //Check Ewald settings
    logFile << " Error: Implicit solvation models cannot be used with PBC.";
    logFile << '\n';
    doQuit = 1;
  }
  //Check threading
  if (Ncpus < 1)
  {
    //Checks the number of threads and continue
    logFile << " Warning: Calculations cannot run with ";
    logFile << Ncpus << " CPUs.";
    logFile << '\n';
    if (JOKES)
    {
      logFile << " Do you know how computers work?";
    }
    logFile << " Ncpus set to 1";
    logFile << '\n' << '\n';
    Ncpus = 1;
    logFile.flush(); //Print warning
  }
  //Wrapper errors
  if ((!TINKER) and (!LAMMPS) and (!QMonly))
  {
    //Check the MM wrappers
    logFile << " Error: No valid MM wrapper selected.";
    logFile << '\n';
    logFile << "  Select a wrapper if you want to run this type ";
    logFile << "of calculation.";
    logFile << '\n';
    doQuit = 1;
  }
  if ((!Gaussian) and (!PSI4) and (!NWChem) and (!MMonly))
  {
    //Check the QM wrappers
    logFile << " Error: No valid QM wrapper selected.";
    logFile << '\n';
    logFile << "  Select a wrapper if you want to run this type ";
    logFile << "of calculation.";
    logFile << '\n';
    doQuit = 1;
  }
  /* if (Gaussian and QMMM)
  {
    //Avoid options that conflict with NWChem capabilities
    if (OptSim)
    {
      //The NWChem optimizer cannot incorporate MM forces
      cout << " Error: QMMM Gaussian optimizations can only be performed";
      cout << '\n';
      cout << " with steepest descent or Davidon-Fletcher-Powell.";
      cout << '\n';
      doQuit = 1;
    }
  }*/
  if (PSI4 and QMMM)
  {
    //Avoid options that conflict with PSI4 capabilities
    if (OptSim)
    {
      //The PSI4 optimizer cannot incorporate MM forces
      logFile << " Error: QMMM PSI4 optimizations can only be performed";
      logFile << '\n';
      logFile << " with steepest descent or Davidon-Fletcher-Powell.";
      logFile << '\n';
      doQuit = 1;
    }
    if ((Npseudo != 0) or (Nbound != 0))
    {
      //PSI4 does not currently have pseudopotentials
      logFile << " Error: The PSI4 wrapper can only use QM and MM atoms.";
      logFile << '\n';
      logFile << " Remove the pseudo-bonds and boundary-atoms.";
      logFile << '\n';
      doQuit = 1;
    }
  }
  if (NWChem and QMMM)
  {
    //Avoid options that conflict with NWChem capabilities
    if (OptSim)
    {
      //The NWChem optimizer cannot incorporate MM forces
      logFile << " Error: QMMM NWChem optimizations can only be performed";
      logFile << '\n';
      logFile << " with steepest descent or Davidon-Fletcher-Powell.";
      logFile << '\n';
      doQuit = 1;
    }
  }
  if (LAMMPS and AMOEBA)
  {
    //Avoid options that conflict with LAMMPS capabilities
    logFile << " Error: LAMMPS calculations cannot be performed with";
    logFile << '\n';
    logFile << " polarizable force fields.";
    logFile << '\n';
    doQuit = 1;
  }
  //Simulation errors
  if ((QMMMOpts.ensemble == "NPT") and (!PBCon))
  {
    //Check the PBC options
    logFile << " Error: NPT simulation without PBC.";
    logFile << '\n';
    logFile << "  Turn PBC on if you want to run this type ";
    logFile << "of calculation.";
    logFile << '\n';
    doQuit = 1;
  }
  if (QMMMOpts.stepScale > 1)
  {
    //Checks the number of threads and continue
    logFile << " Warning: The optimization step scale cannot be greater";
    logFile << " than 1.";
    logFile << '\n';
    logFile << " Step scale set to 1.";
    logFile << '\n';
    QMMMOpts.stepScale = 1; //Reset step size
    logFile.flush(); //Print warning
  }
  if (doQuit)
  {
    //Quits
    logFile << '\n';
    logFile.flush();
    stat=1;
    return;
  }
  //Sarcastically continue
  logFile << "No fatal errors detected.";
  logFile << '\n';
  if (JOKES)
  {
    logFile << " And there was much rejoicing. Yay...";
    logFile << '\n';
    logFile << '\n';
    logFile.flush();
    if (CheckFile("EASTEREGG"))
    {
      PrintLapin(logFile);
    }
  }
  return;
};

void LICHEMPrintSettings(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                         fstream& logFile)
{
  //Prints out the simulation details
  logFile << "Setting up simulation..." << '\n';
  logFile << '\n';
  logFile << "Input files:" << '\n';
  logFile << " Coordinate file: " << xyzFilename << '\n';
  logFile << " Connectivity file: " << conFilename << '\n';
  logFile << " Region file: " << regFilename << '\n';
  if (CheckFile("BeadStartStruct.xyz"))
  {
    logFile << " Restart file: BeadStartStruct.xyz" << '\n';
  }
  //Start: Hatice
  if (CheckFile("QSMBeadStruct.xyz"))
  {
    logFile << " Restart file: QSMBeadStruct.xyz" << '\n';
  }
  //End: Hatice
  logFile << '\n';
  logFile << "Atoms: " << Natoms << '\n';
  if (QMonly or QMMM)
  {
    //QM regions
    logFile << " QM atoms: " << Nqm << '\n';
    logFile << "  Charge: " << QMMMOpts.charge << '\n';
    logFile << "  Spin: " << QMMMOpts.spin << '\n';
  }
  if (MMonly or QMMM)
  {
    //MM regions
    logFile << " MM atoms: " << Nmm << '\n';
    if (QMMM)
    {
      logFile << " Pseudo-atoms: " << Npseudo << '\n';
      logFile << " Boundary-atoms: " << Nbound << '\n';
    }
    if (Nfreeze > 0)
    {
      logFile << " Frozen atoms: " << Nfreeze << '\n';
    }
  }
  if (NEBSim)
  {
    //Print reaction path input for error checking
    logFile << " RP beads: " << QMMMOpts.NBeads << '\n';
    logFile << '\n';
    logFile << "Simulation mode: ";
    if (QMMM)
    {
      logFile << "QMMM";
    }
    if (QMonly)
    {
      logFile << "Pure QM";
    }
    if (MMonly)
    {
      logFile << "Pure MM";
    }
    logFile << " NEB" << '\n';
  }
  //Start: Hatice
  if (QSMSim)
  {
    //Print reaction path input for error checking
    logFile << " RP beads: " << QMMMOpts.NBeads << '\n';
    logFile << '\n';
    logFile << "Simulation mode: ";
    if (QMMM)
    {
      logFile << "QMMM";
    }
    if (QMonly)
    {
      logFile << "Pure QM";
    }
    if (MMonly)
    {
      logFile << "Pure MM";
    }
    logFile << " QSM" << '\n';
  }
  //End: Hatice
  if (PIMCSim)
  {
    //Print PIMC input for error checking
    if (QMMMOpts.NBeads > 1)
    {
      logFile << " PI beads: " << QMMMOpts.NBeads << '\n';
    }
    logFile << '\n';
    logFile << "Simulation mode: ";
    if (QMMM)
    {
      logFile << "QMMM";
    }
    if (QMonly)
    {
      logFile << "Pure QM";
    }
    if (MMonly)
    {
      logFile << "Pure MM";
    }
    logFile << " " << QMMMOpts.ensemble;
    if (QMMMOpts.NBeads > 1)
    {
      logFile << " path-integral";
    }
    logFile << " Monte Carlo" << '\n';
    logFile << " Equilibration MC steps: " << QMMMOpts.NEq << '\n';
    logFile << " Production MC steps: " << QMMMOpts.NSteps << '\n';
  }
  if (FBNEBSim)
  {
    //Print FBNEB input for error checking
    if (QMMMOpts.NBeads > 1)
    {
      logFile << " RP beads: " << QMMMOpts.NBeads << '\n';
    }
    logFile << '\n';
    logFile << "Simulation mode: ";
    if (QMMM)
    {
      logFile << "QMMM";
    }
    if (QMonly)
    {
      logFile << "Pure QM";
    }
    if (MMonly)
    {
      logFile << "Pure MM";
    }
    logFile << " NVT";
    if (QMMMOpts.NBeads > 1)
    {
      logFile << " force-bias";
    }
    logFile << " Monte Carlo" << '\n';
    logFile << " Equilibration MC steps: " << QMMMOpts.NEq << '\n';
    logFile << " Production MC steps: " << QMMMOpts.NSteps << '\n';
  }
  if (OptSim or SteepSim or DFPSim)
  {
    //Print optimization input for error checking
    logFile << '\n';
    logFile << "Simulation mode: ";
    if (QMMM)
    {
      logFile << "QMMM";
    }
    if (QMonly)
    {
      logFile << "Pure QM";
    }
    if (MMonly)
    {
      logFile << "Pure MM";
    }
    logFile << " energy minimization" << '\n';
    if (QMMM or QMonly)
    {
      logFile << " QM";
      if (QMMM)
      {
        logFile << "MM";
      }
      logFile << " minimizer: ";
      if (OptSim)
      {
        logFile << "Native QM optimizer" << '\n';
      }
      if (SteepSim)
      {
        logFile << "LICHEM steepest descent" << '\n';
      }
      if (DFPSim)
      {
        logFile << "LICHEM DFP" << '\n';
      }
    }
  }
  if (SinglePoint)
  {
    //Print single-point energy settings for error checking
    logFile << '\n';
    logFile << "Simulation mode: ";
    if (QMMM)
    {
      logFile << "QMMM";
    }
    if (QMonly)
    {
      logFile << "Pure QM";
    }
    if (MMonly)
    {
      logFile << "Pure MM";
    }
    if (QMMMOpts.NBeads == 1)
    {
      logFile << " single-point energy" << '\n';
    }
    else
    {
      logFile << " multi-point energy" << '\n';
    }
  }
  if (FreqCalc)
  {
    //Print frequency settings for error checking
    logFile << '\n';
    logFile << "Simulation mode: ";
    if (QMMM)
    {
      logFile << "QMMM";
    }
    if (QMonly)
    {
      logFile << "Pure QM";
    }
    if (MMonly)
    {
      logFile << "Pure MM";
    }
    if (QMMMOpts.NBeads == 1)
    {
      logFile << " single-point frequencies" << '\n';
    }
    else
    {
      logFile << " multi-point frequencies" << '\n';
    }
  }
  if (QMonly or QMMM)
  {
    //Print QM wrapper input for error checking
    logFile << " QM wrapper: ";
    if (PSI4)
    {
      logFile << "PSI4" << '\n';
    }
    if (Gaussian)
    {
      logFile << "Gaussian" << '\n';
    }
    if (NWChem)
    {
      logFile << "NWChem" << '\n';
    }
    logFile << " QM method: ";
    if (QMMMOpts.func != "SemiEmp")
    {
      //Avoid printing method and basis for semi-empirical
      logFile << QMMMOpts.func << "/";
    }
    logFile << QMMMOpts.basis << '\n';
  }
  if (MMonly or QMMM)
  {
    //Print MM wrapper input for error checking
    logFile << " MM wrapper: ";
    if (TINKER)
    {
      logFile << "TINKER" << '\n';
    }
    if (LAMMPS)
    {
      logFile << "LAMMPS" << '\n';
    }
    if (QMMM)
    {
      //Print QMMM wrapper input for error checking
      logFile << " MM potential: ";
      if (CHRG)
      {
        logFile << "Point-charge force field" << '\n';
      }
      if (AMOEBA)
      {
        logFile << "Polarizable force field" << '\n';
      }
      if (GEM)
      {
        logFile << "Diffuse-charge force field" << '\n';
      }
    }
    //Start: Hatice GOKCAN
    if(QMMMOpts.restrMM)
    {
      logFile << " MM restrain: YES" << '\n';
      logFile << " Force constant for restrain: " << QMMMOpts.restrConst << '\n';
    }
    //End: Hatice GOKCAN
    //Print PBC information
    if (PBCon or QMMMOpts.useLREC or QMMMOpts.useImpSolv)
    {
      logFile << '\n';
      logFile << "Simulation box settings:" << '\n';
      if (PBCon)
      {
        //Print box size and density
        double initDen = 0; //Initial density
        logFile << " Boundaries: Periodic" << '\n';
        logFile << " Box size (\u212B): ";
        logFile << LICHEMFormFloat(Lx,10) << " ";
        logFile << LICHEMFormFloat(Ly,10) << " ";
        logFile << LICHEMFormFloat(Lz,10) << '\n';
        logFile << " Density: ";
        initDen = LICHEMDensity(QMMMData,QMMMOpts);
        logFile << LICHEMFormFloat(initDen,10);
        logFile << " g/cm\u00B3" << '\n';
      }
      if (QMMMOpts.useLREC)
      {
        //Print LREC cutoff options
        logFile << " QM LREC: Yes" << '\n';
        logFile << " LREC cutoff: ";
        logFile << LICHEMFormFloat(QMMMOpts.LRECCut,8);
        logFile << " \u212B" << '\n';
        logFile << " LREC exponent: " << QMMMOpts.LRECPow << '\n';
      }
      if (QMMMOpts.useEwald)
      {
        //Print Ewald summation options
        logFile << " MM Ewald: Yes" << '\n';
      }
      if (QMMMOpts.useImpSolv)
      {
        //Print continuum solvation options
        logFile << " Implicit solvent: " << QMMMOpts.solvModel;
        logFile << '\n';
      }
    }
  }
  logFile << '\n';
  //Print parallelization settings
  logFile << "Parallelization and memory settings:" << '\n';
  logFile << " OpenMP threads: " << Nthreads << '\n';
  if (QMonly or QMMM)
  {
    logFile << " QM threads: " << Ncpus << '\n';
    logFile << " QM memory: " << QMMMOpts.RAM << " ";
    if (QMMMOpts.memMB)
    {
      logFile << "MB";
    }
    else
    {
      logFile << "GB";
    }
    logFile << '\n';
  }
  if (MMonly or QMMM)
  {
    logFile << " MM threads: " << Ncpus << '\n';
  }
  //Print Monte Carlo settings
  if (PIMCSim or FBNEBSim)
  {
    logFile << '\n';
    logFile << "Monte Carlo settings:" << '\n';
    logFile << " Temperature: " << QMMMOpts.temp;
    logFile << " K" << '\n';
    if (QMMMOpts.ensemble == "NPT")
    {
      logFile << " Pressure: " << QMMMOpts.press;
      logFile << " atm" << '\n';
    }
    if (FBNEBSim and (QMMMOpts.NBeads > 1))
    {
      logFile << " Spring constant: " << QMMMOpts.kSpring;
      logFile << " eV/\u212B\u00B2" << '\n';
    }
    logFile << " Acceptance ratio: ";
    logFile << LICHEMFormFloat(QMMMOpts.accRatio,4);
    logFile << '\n';
    logFile << " Equilibration MC steps: " << QMMMOpts.NEq;
    logFile << '\n';
    logFile << " Production MC steps: " << QMMMOpts.NSteps;
    logFile << '\n';
    logFile << " Sample every " << QMMMOpts.NPrint;
    logFile << " steps" << '\n';
  }
  //Print convergence criteria for optimizations
  //Start: Hatice
  //if (OptSim or SteepSim or DFPSim or NEBSim)
  if (OptSim or SteepSim or DFPSim or NEBSim or QSMSim)
  //End: Hatice
  {
    logFile << '\n';
    logFile << "Optimization settings:" << '\n';
    //Start: Hatice
    //if (!OptSim)
    if (!OptSim and !QSMSim)
    //End: Hatice
    {
      logFile << " Step scale factor: ";
      logFile << LICHEMFormFloat(QMMMOpts.stepScale,6);
      logFile << '\n';
    }
    //Start: Hatice
    if(!QSMSim){
      logFile << " Max. step size: ";
      logFile << LICHEMFormFloat(QMMMOpts.maxStep,6);
      logFile << " \u212B" << '\n';
      logFile << " Max. steps: " << QMMMOpts.maxOptSteps << '\n';
    }
    //End: Hatice
    if (QMMMOpts.useMMCut and (Nmm > 0))
    {
      //Print MM cutoff settings
      logFile << '\n';
      logFile << " MM cutoff: ";
      logFile << LICHEMFormFloat(QMMMOpts.MMOptCut,8);
      logFile << " \u212B" << '\n';;
    }
    if (NEBSim)
    {
      //Spring constant for the path
      logFile << '\n';
      logFile << " Spring constant: " << QMMMOpts.kSpring;
      logFile << " eV/\u212B\u00B2" << '\n';
      logFile << " End points: ";
      if (QMMMOpts.frznEnds)
      {
        logFile << "Frozen";
      }
      else
      {
        logFile << "Active";
      }
    }
    //Start: Hatice
    if (QSMSim)
    {
      //Spring constant for the path
      logFile << " End points: ";
      if (QMMMOpts.frznEnds)
      {
        logFile << "Frozen \n";
      }
      else
      {
        logFile << "Active \n";
      }
    }
    //End: hatice
    if (SteepSim or DFPSim or NEBSim or QSMSim)
    {
      logFile << " Max. opt. steps: " << QMMMOpts.maxOptSteps;
      logFile << '\n';
      logFile << " Max. QM steps per opt. step: " << QMMMOpts.MaxQMSteps;
      logFile << '\n';
      logFile << '\n';
      logFile << "QM convergence criteria:" << '\n';
      logFile << " RMS deviation: " << QMMMOpts.QMOptTol;
      logFile << " \u212B" << '\n';
      logFile << " RMS force: " << LICHEMFormFloat(QMMMOpts.QMRMSForceTol,8);
      logFile << " Hartree/bohr" << '\n';
      logFile << " Max. force: " << LICHEMFormFloat(QMMMOpts.QMMaxForceTol,8);
      logFile << " Hartree/bohr" << '\n';
    }
    
    if (Nmm > 0)
    {
      logFile << '\n';
      logFile << "MM convergence criteria:" << '\n';
      logFile << " RMS deviation: " << QMMMOpts.MMOptTol;
      logFile << " \u212B" << '\n';
      logFile << " RMS force: ";
      //Start: Hatice GOKCAN
      logFile << QMMMOpts.MMOptTol;
      logFile << " kcal/\u212B" << '\n';
      //End: Hatice GOKCAN
    }
  }
  //Print frequency analysis settings
  if (FreqCalc or QMMMOpts.NEBFreq)
  {
    logFile << '\n';
    logFile << "Frequency settings:" << '\n';
    //Always removed
    logFile << "  Remove low frequencies: Yes";
    logFile << '\n';
    //Removed for QM calculations
    logFile << "  Remove translations: ";
    if (QMMM)
    {
      logFile << "No" << '\n';
    }
    else
    {
      logFile << "Yes" << '\n';
    }
    //Removed for QM calculations
    logFile << "  Remove rotations: ";
    if (QMMM)
    {
      logFile << "No" << '\n';
    }
    else
    {
      logFile << "Yes" << '\n';
    }
  }
  logFile << '\n';
  logFile.flush(); //Flush for output being redirected to a file
  return;
};

