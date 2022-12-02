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
  # Functions to write Gaussian input in parallel LICHEM                      #
  # Includes:                                                                 #
  #                                                                           #
  #            for primary proc:                                              #
  #                                                                           #
  #                            Writing input files :                          #
  #                                   void WriteGauInputMPI                   #
  #                                                                           #
  #############################################################################
*/


void WriteGauInputMPI(vector<QMMMAtom>& QMMMData, string calcTyp,
                      QMMMSettings& QMMMOpts, int bead,
                      int gautype)
{
  // Write Gaussian input files
  stringstream call; // Stream for system calls and reading/writing files
  call.copyfmt(cout); // Copy settings from cout
  string dummy,chrgfilename; // Generic strings
  fstream inFile,outFile; // Generic file names
  // Check units
  double uConv = 1; // Units conversion constant
  if (QMMMOpts.unitsQM == "Bohr")
  {
    uConv = 1.0/bohrRad;
  }
  // Check for a charge file
  bool useChargeFile = 0;
  call.str("");
  call << "MMCharges_" << bead << ".txt";
  chrgfilename = call.str();
  useChargeFile = CheckFile(call.str());
  if (Nmm == 0)
  {
    // Skip blank charge files
    useChargeFile = 0;
  }
  // Initialize multipoles and center of mass
  bool firstCharge = 1; // Always write the first charge
  Coord QMCOM;
  if (!useChargeFile)
  {
    if (PBCon or QMMMOpts.useLREC)
    {
      QMCOM = FindQMCOM(QMMMData,QMMMOpts,bead);
    }
    if (AMOEBA)
    {
      if (TINKER)
      {
        // Set up multipoles
        RotateTINKCharges(QMMMData,bead);
      }
    }
  }
  // Construct g09 input
  call.str("");
  /*
    gautype=0 -> energy
    gautype=1 -> force 
  */
  if (gautype==0)
  {
    call << "LICHM_GauEner_" << bead << ".com";
  }
  if (gautype==1)
  {
    call << "LICHM_GauForce_" << bead << ".com";
  }
  outFile.open(call.str().c_str(),ios_base::out);
  call.str("");
  if (gautype==0)
  {
    call << "%chk=LICHM_" << bead << ".chk";
  }
  if (gautype==1)
  {
    call << "%chk=LICHM_" << bead << ".chk";
  }

  call << '\n';
  call << "%Mem=" << QMMMOpts.RAM;
  if (QMMMOpts.memMB)
  {
    call << "MB";
  }
  else
  {
    call << "GB";
  }
  call << '\n';
  // Start: Hatice
  /* call << "%NprocShared=" << Ncpus << '\n'; */
  // End: Hatice
  // Add ROUTE section
  call << calcTyp;
  // Add structure
  call << '\n'; // Blank line
  call << "QMMM" << '\n' << '\n'; // Dummy title
  call << QMMMOpts.charge << " " << QMMMOpts.spin << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (QMMMData[i].QMRegion)
    {
      call << QMMMData[i].QMTyp;
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
      call << '\n';
    }
    if (QMMMData[i].PBRegion)
    {
      call << "F";
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
      call << '\n';
    }
  }
  call << '\n'; // Blank line needed
  // Add the MM field
  if (QMMM and useChargeFile)
  {
    inFile.open(chrgfilename.c_str(),ios_base::in);
    if (inFile.good())
    {
      while (!inFile.eof())
      {
        // Copy charge file line by line
        getline(inFile,dummy);
        call << dummy << '\n';
      }
      inFile.close();
    }
  }
  else if (QMMM)
  {
    if (CHRG)
    {
      for (int i=0;i<Natoms;i++)
      {
        if (QMMMData[i].MMRegion)
        {
          // Check PBC (minimum image convention)
          Coord distCent; // Distance from QM COM
          double xShft = 0;
          double yShft = 0;
          double zShft = 0;
          if (PBCon or QMMMOpts.useLREC)
          {
            // Initialize displacements
            double dx,dy,dz; // Starting displacements
            dx = QMMMData[i].P[bead].x-QMCOM.x;
            dy = QMMMData[i].P[bead].y-QMCOM.y;
            dz = QMMMData[i].P[bead].z-QMCOM.z;
            distCent = CoordDist2(QMMMData[i].P[bead],QMCOM);
            // Calculate the shift in positions
            // NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xShft = distCent.x-dx;
              yShft = distCent.y-dy;
              zShft = distCent.z-dz;
            }
          }
          // Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.useLREC)
          {
            // Use the long-range correction
            scrq = LRECFunction(distCent,QMMMOpts);
          }
          if ((scrq > 0) or firstCharge)
          {
            // Add charges
            firstCharge = 0; // Skips writing the remaining zeros
            double tmpX,tmpY,tmpZ,tmpQ; // Temporary storage
            tmpX = (QMMMData[i].P[bead].x+xShft)*uConv;
            tmpY = (QMMMData[i].P[bead].y+yShft)*uConv;
            tmpZ = (QMMMData[i].P[bead].z+zShft)*uConv;
            tmpQ = QMMMData[i].MP[bead].q*scrq;
            call << " ";
            call << LICHEMFormFloat(tmpX,16);
            call << " ";
            call << LICHEMFormFloat(tmpY,16);
            call << " ";
            call << LICHEMFormFloat(tmpZ,16);
            call << " ";
            call << LICHEMFormFloat(tmpQ,16);
            call << '\n';
          }
        }
      }
      if (Nmm > 0)
      {
        call << '\n'; // Blank line needed
      }
    }
    if (AMOEBA)
    {
      for (int i=0;i<Natoms;i++)
      {
        if (QMMMData[i].MMRegion)
        {
          // Check PBC (minimum image convention)
          Coord distCent; // Distance from QM COM
          double xShft = 0;
          double yShft = 0;
          double zShft = 0;
          if (PBCon or QMMMOpts.useLREC)
          {
            // Initialize displacements
            double dx,dy,dz; // Starting displacements
            dx = QMMMData[i].P[bead].x-QMCOM.x;
            dy = QMMMData[i].P[bead].y-QMCOM.y;
            dz = QMMMData[i].P[bead].z-QMCOM.z;
            distCent = CoordDist2(QMMMData[i].P[bead],QMCOM);
            // Calculate the shift in positions
            // NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xShft = distCent.x-dx;
              yShft = distCent.y-dy;
              zShft = distCent.z-dz;
            }
          }
          // Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.useLREC)
          {
            // Use the long-range correction
            scrq = LRECFunction(distCent,QMMMOpts);
          }
          if ((scrq > 0) or firstCharge)
          {
            firstCharge = 0; // Skips writing the remaining zeros
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].x1+xShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].y1+yShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].z1+zShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].q1*scrq,16);
            call << '\n';
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].x2+xShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].y2+yShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].z2+zShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].q2*scrq,16);
            call << '\n';
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].x3+xShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].y3+yShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].z3+zShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].q3*scrq,16);
            call << '\n';
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].x4+xShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].y4+yShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].z4+zShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].q4*scrq,16);
            call << '\n';
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].x5+xShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].y5+yShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].z5+zShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].q5*scrq,16);
            call << '\n';
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].x6+xShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].y6+yShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].z6+zShft,16);
            call << " ";
            call << LICHEMFormFloat(QMMMData[i].PC[bead].q6*scrq,16);
            call << '\n';
          }
        }
      }
      if (Nmm > 0)
      {
        call << '\n'; // Blank line needed
      }
    }
  }
  // Add basis set information from the BASIS file
  inFile.open("BASIS",ios_base::in);
  if (inFile.good())
  {
    while (!inFile.eof())
    {
      // Copy BASIS line by line, if BASIS exists
      getline(inFile,dummy);
      call << dummy << '\n';
    }
    inFile.close();
  }
  outFile << call.str();
  outFile.flush();
  outFile.close();
  return;
};
