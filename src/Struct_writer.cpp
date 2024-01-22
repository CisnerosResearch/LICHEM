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

 Functions write input files for the QM wrappers.

*/

//QM input writers
void WriteGauInput(vector<QMMMAtom>& QMMMData, string calcTyp,
                   QMMMSettings& QMMMOpts, int bead)
{
  //Write Gaussian input files
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy,chrgfilename; //Generic strings
  fstream inFile,outFile; //Generic file names
  //Check units
  double uConv = 1; //Units conversion constant
  if (QMMMOpts.unitsQM == "Bohr")
  {
    uConv = 1.0/bohrRad;
  }
  //Check for a charge file
  bool useChargeFile = 0;
  call.str("");
  call << "MMCharges_" << bead << ".txt";
  chrgfilename = call.str();
  useChargeFile = CheckFile(call.str());
  if (Nmm == 0)
  {
    //Skip blank charge files
    useChargeFile = 0;
  }
  //Initialize multipoles and center of mass
  bool firstCharge = 1; //Always write the first charge
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
        //Set up multipoles
        RotateTINKCharges(QMMMData,bead);
      }
    }
  }
  //Construct g09 input
  call.str("");
  call << "LICHM_" << bead << ".com";
  outFile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "%chk=LICHM_" << bead << ".chk";
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
  //Start: Hatice
  if(g09){
    call << "%NprocShared=" << Ncpus << '\n';
  }
  else{
    call << "%CPU=0-" << Ncpus-1 << '\n';
  }
  //End: Hatice
  //Add ROUTE section
  call << calcTyp;
  //Add structure
  call << '\n'; //Blank line
  call << "QMMM" << '\n' << '\n'; //Dummy title
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
  call << '\n'; //Blank line needed
  //Add the MM field
  if (QMMM and useChargeFile)
  {
    inFile.open(chrgfilename.c_str(),ios_base::in);
    if (inFile.good())
    {
      while (!inFile.eof())
      {
        //Copy charge file line by line
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
          //Check PBC (minimum image convention)
          Coord distCent; //Distance from QM COM
          double xShft = 0;
          double yShft = 0;
          double zShft = 0;
          if (PBCon or QMMMOpts.useLREC)
          {
            //Initialize displacements
            double dx,dy,dz; //Starting displacements
            dx = QMMMData[i].P[bead].x-QMCOM.x;
            dy = QMMMData[i].P[bead].y-QMCOM.y;
            dz = QMMMData[i].P[bead].z-QMCOM.z;
            distCent = CoordDist2(QMMMData[i].P[bead],QMCOM);
            //Calculate the shift in positions
            //NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xShft = distCent.x-dx;
              yShft = distCent.y-dy;
              zShft = distCent.z-dz;
            }
          }
          //Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.useLREC)
          {
            //Use the long-range correction
            scrq = LRECFunction(distCent,QMMMOpts);
          }
          if ((scrq > 0) or firstCharge)
          {
            //Add charges
            firstCharge = 0; //Skips writing the remaining zeros
            double tmpX,tmpY,tmpZ,tmpQ; //Temporary storage
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
        call << '\n'; //Blank line needed
      }
    }
    if (AMOEBA)
    {
      for (int i=0;i<Natoms;i++)
      {
        if (QMMMData[i].MMRegion)
        {
          //Check PBC (minimum image convention)
          Coord distCent; //Distance from QM COM
          double xShft = 0;
          double yShft = 0;
          double zShft = 0;
          if (PBCon or QMMMOpts.useLREC)
          {
            //Initialize displacements
            double dx,dy,dz; //Starting displacements
            dx = QMMMData[i].P[bead].x-QMCOM.x;
            dy = QMMMData[i].P[bead].y-QMCOM.y;
            dz = QMMMData[i].P[bead].z-QMCOM.z;
            distCent = CoordDist2(QMMMData[i].P[bead],QMCOM);
            //Calculate the shift in positions
            //NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xShft = distCent.x-dx;
              yShft = distCent.y-dy;
              zShft = distCent.z-dz;
            }
          }
          //Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.useLREC)
          {
            //Use the long-range correction
            scrq = LRECFunction(distCent,QMMMOpts);
          }
          if ((scrq > 0) or firstCharge)
          {
            firstCharge = 0; //Skips writing the remaining zeros
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
        call << '\n'; //Blank line needed
      }
    }
  }
  //Add basis set information from the BASIS file
  inFile.open("BASIS",ios_base::in);
  if (inFile.good())
  {
    while (!inFile.eof())
    {
      //Copy BASIS line by line, if BASIS exists
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

void WriteNWChemInput(vector<QMMMAtom>& QMMMData, string calcTyp,
                      QMMMSettings& QMMMOpts, int bead)
{
  //Write NWChem input files
  fstream outFile,inFile; //Generic file streams
  string dummy,chrgfilename; //Generic strings
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  //Check units
  double uConv = 1; //Units conversion constant
  if (QMMMOpts.unitsQM == "Bohr")
  {
    uConv = 1.0/bohrRad;
  }
  //Check for a charge file
  bool useChargeFile = 0;
  call.str("");
  call << "MMCharges_" << bead << ".txt";
  chrgfilename = call.str();
  useChargeFile = CheckFile(call.str());
  if (Nmm == 0)
  {
    //Skip blank charge files
    useChargeFile = 0;
  }
  //Initialize multipoles and center of mass
  bool firstCharge = 1; //Always write the first charge
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
        //Set up multipoles
        RotateTINKCharges(QMMMData,bead);
      }
    }
  }
  //Create NWChem input
  call.str("");
  call << "LICHM_" << bead << ".nw";
  outFile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "LICHM_" << bead << ".db";
  if (CheckFile(call.str()))
  {
    outFile << "restart";
  }
  else
  {
    outFile << "start";
  }
  outFile << " LICHM_" << bead << '\n';
  outFile << "memory " << QMMMOpts.RAM;
  if (QMMMOpts.memMB)
  {
    outFile << " mb";
  }
  else
  {
    outFile << " gb";
  }
  outFile << '\n';
  outFile << "charge " << QMMMOpts.charge << '\n';
  outFile << "geometry nocenter ";
  outFile << "noautoz noautosym" << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (QMMMData[i].QMRegion)
    {
      outFile << " " << QMMMData[i].QMTyp;
      outFile << " " << (QMMMData[i].P[bead].x);
      outFile << " " << (QMMMData[i].P[bead].y);
      outFile << " " << (QMMMData[i].P[bead].z);
      outFile << '\n';
    }
    if (QMMMData[i].PBRegion)
    {
      outFile << " " << "F2pb";
      outFile << " " << (QMMMData[i].P[bead].x);
      outFile << " " << (QMMMData[i].P[bead].y);
      outFile << " " << (QMMMData[i].P[bead].z);
      outFile << '\n';
    }
  }
  outFile << "end" << '\n';
  if (CheckFile("BASIS"))
  {
    //Add basis set and ecp info
    inFile.open("BASIS",ios_base::in);
    if (inFile.good())
    {
      while (!inFile.eof())
      {
        //Copy BASIS line by line, if BASIS exists
        getline(inFile,dummy);
        stringstream line(dummy);
        if (line.str() != "")
        {
          //Avoid copying extra blank lines
          outFile << dummy << '\n';
        }
      }
    }
    inFile.close();
  }
  else
  {
    outFile << "basis" << '\n';
    outFile << " * library " << QMMMOpts.basis;
    outFile << '\n';
    outFile << "end" << '\n';
  }
  if (QMMM and useChargeFile and (Nmm > 0))
  {
    inFile.open(chrgfilename.c_str(),ios_base::in);
    if (inFile.good())
    {
      outFile << "set bq:max_nbq " << (6*(Nmm+Nbound)) << '\n';
      outFile << "bq mmchrg";
      while (!inFile.eof())
      {
        //Copy charge file line by line
        outFile << '\n'; //Avoid adding an extra blank line
        getline(inFile,dummy);
        outFile << dummy; //Print the line
      }
      inFile.close();
      outFile << "end" << '\n';
      outFile << "set bq mmchrg" << '\n';
    }
  }
  else if (QMMM and (Nmm > 0))
  {
    if (CHRG)
    {
      outFile << "set bq:max_nbq " << (Nmm+Nbound) << '\n';
      outFile << "bq mmchrg" << '\n';
      for (int i=0;i<Natoms;i++)
      {
        if (QMMMData[i].MMRegion)
        {
          //Check PBC (minimum image convention)
          Coord distCent; //Distance from QM COM
          double xShft = 0;
          double yShft = 0;
          double zShft = 0;
          if (PBCon or QMMMOpts.useLREC)
          {
            //Initialize displacements
            double dx,dy,dz; //Starting displacements
            dx = QMMMData[i].P[bead].x-QMCOM.x;
            dy = QMMMData[i].P[bead].y-QMCOM.y;
            dz = QMMMData[i].P[bead].z-QMCOM.z;
            distCent = CoordDist2(QMMMData[i].P[bead],QMCOM);
            //Calculate the shift in positions
            //NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xShft = distCent.x-dx;
              yShft = distCent.y-dy;
              zShft = distCent.z-dz;
            }
          }
          //Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.useLREC)
          {
            //Use the long-range correction
            scrq = LRECFunction(distCent,QMMMOpts);
          }
          if ((scrq > 0) or firstCharge)
          {
            //Add charges
            firstCharge = 0; //Skips writing the remaining zeros
            double tmpX,tmpY,tmpZ,tmpQ; //Temporary storage
            tmpX = (QMMMData[i].P[bead].x+xShft)*uConv;
            tmpY = (QMMMData[i].P[bead].y+yShft)*uConv;
            tmpZ = (QMMMData[i].P[bead].z+zShft)*uConv;
            tmpQ = QMMMData[i].MP[bead].q*scrq;
            outFile << " ";
            outFile << LICHEMFormFloat(tmpX,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpY,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpZ,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpQ,16);
            outFile << '\n';
          }
        }
      }
      outFile << "end" << '\n';
      outFile << "set bq mmchrg" << '\n';
    }
    if (AMOEBA)
    {
      outFile << "set bq:max_nbq " << (6*(Nmm+Nbound)) << '\n';
      outFile << "bq mmchrg" << '\n';
      for (int i=0;i<Natoms;i++)
      {
        if (QMMMData[i].MMRegion)
        {
          //Check PBC (minimum image convention)
          Coord distCent; //Distance from QM COM
          double xShft = 0;
          double yShft = 0;
          double zShft = 0;
          if (PBCon or QMMMOpts.useLREC)
          {
            //Initialize displacements
            double dx,dy,dz; //Starting displacements
            dx = QMMMData[i].P[bead].x-QMCOM.x;
            dy = QMMMData[i].P[bead].y-QMCOM.y;
            dz = QMMMData[i].P[bead].z-QMCOM.z;
            distCent = CoordDist2(QMMMData[i].P[bead],QMCOM);
            //Calculate the shift in positions
            //NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xShft = distCent.x-dx;
              yShft = distCent.y-dy;
              zShft = distCent.z-dz;
            }
          }
          //Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.useLREC)
          {
            //Use the long-range correction
            scrq = LRECFunction(distCent,QMMMOpts);
          }
          if ((scrq > 0) or firstCharge)
          {
            //Add multipoles
            firstCharge = 0; //Skips writing the remaining zeros
            double tmpX,tmpY,tmpZ,tmpQ; //Temporary storage
            //Charge 1
            tmpX = (QMMMData[i].PC[bead].x1+xShft)*uConv;
            tmpY = (QMMMData[i].PC[bead].y1+yShft)*uConv;
            tmpZ = (QMMMData[i].PC[bead].z1+zShft)*uConv;
            tmpQ = QMMMData[i].PC[bead].q1*scrq;
            outFile << " ";
            outFile << LICHEMFormFloat(tmpX,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpY,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpZ,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpQ,16);
            outFile << '\n';
            //Charge 2
            tmpX = (QMMMData[i].PC[bead].x2+xShft)*uConv;
            tmpY = (QMMMData[i].PC[bead].y2+yShft)*uConv;
            tmpZ = (QMMMData[i].PC[bead].z2+zShft)*uConv;
            tmpQ = QMMMData[i].PC[bead].q2*scrq;
            outFile << " ";
            outFile << LICHEMFormFloat(tmpX,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpY,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpZ,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpQ,16);
            outFile << '\n';
            //Charge 3
            tmpX = (QMMMData[i].PC[bead].x3+xShft)*uConv;
            tmpY = (QMMMData[i].PC[bead].y3+yShft)*uConv;
            tmpZ = (QMMMData[i].PC[bead].z3+zShft)*uConv;
            tmpQ = QMMMData[i].PC[bead].q3*scrq;
            outFile << " ";
            outFile << LICHEMFormFloat(tmpX,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpY,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpZ,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpQ,16);
            outFile << '\n';
            //Charge 4
            tmpX = (QMMMData[i].PC[bead].x4+xShft)*uConv;
            tmpY = (QMMMData[i].PC[bead].y4+yShft)*uConv;
            tmpZ = (QMMMData[i].PC[bead].z4+zShft)*uConv;
            tmpQ = QMMMData[i].PC[bead].q4*scrq;
            outFile << " ";
            outFile << LICHEMFormFloat(tmpX,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpY,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpZ,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpQ,16);
            outFile << '\n';
            //Charge 5
            tmpX = (QMMMData[i].PC[bead].x5+xShft)*uConv;
            tmpY = (QMMMData[i].PC[bead].y5+yShft)*uConv;
            tmpZ = (QMMMData[i].PC[bead].z5+zShft)*uConv;
            tmpQ = QMMMData[i].PC[bead].q5*scrq;
            outFile << " ";
            outFile << LICHEMFormFloat(tmpX,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpY,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpZ,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpQ,16);
            outFile << '\n';
            //Charge 6
            tmpX = (QMMMData[i].PC[bead].x6+xShft)*uConv;
            tmpY = (QMMMData[i].PC[bead].y6+yShft)*uConv;
            tmpZ = (QMMMData[i].PC[bead].z6+zShft)*uConv;
            tmpQ = QMMMData[i].PC[bead].q6*scrq;
            outFile << " ";
            outFile << LICHEMFormFloat(tmpX,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpY,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpZ,16);
            outFile << " ";
            outFile << LICHEMFormFloat(tmpQ,16);
            outFile << '\n';
          }
        }
      }
      outFile << "end" << '\n';
      outFile << "set bq mmchrg" << '\n';
    }
  }
  //Add DFT settings
  outFile << "dft" << '\n';
  outFile << " mult " << QMMMOpts.spin << '\n';
  outFile << " direct" << '\n';
  outFile << " grid xfine nodisk" << '\n';
  outFile << " noio" << '\n';
  outFile << " tolerances tight" << '\n';
  outFile << " xc " << QMMMOpts.func << '\n';
  //Use the checkpoint file
  call.str("");
  call << "LICHM_" << bead << ".movecs";
  if (CheckFile(call.str()))
  {
    //Tell the DFT module to read the initial vectors
    outFile << " vectors input ";
    outFile << call.str(); //Defined above
    outFile << '\n';
  }
  outFile << "end" << '\n';
  //Set calculation type
  outFile << calcTyp;
  //Print file
  outFile.flush();
  outFile.close();
  return;
};

//S:JORGE
void WritePSITHONInput(vector<QMMMAtom>& QMMMData, string calcTyp,
                    QMMMSettings& QMMMOpts, int bead)
{
  //Write PSI4 input files
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy,chrgfilename; //Generic string
  fstream inFile,outFile; //Generic file names
  //Check units
  double uConv = 1; //Units conversion constant
  uConv = 1.0/bohrRad;
  //Check for a charge file
  bool useChargeFile = 0;
  call.str("");
  call << "MMCharges_" << bead << ".txt";
  chrgfilename = call.str();
  useChargeFile = CheckFile(call.str());
  if (Nmm == 0)
  {
    //Skip blank charge files
    useChargeFile = 0;
  }
  //Initialize multipoles and center of mass
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
        //Set up multipoles
        RotateTINKCharges(QMMMData,bead);
      }
    }
  }
  //Check if there is a checkpoint file
  bool useCheckPoint;
  call.str("");
  call << "LICHM_" << bead << ".180";
  useCheckPoint = CheckFile(call.str());
  //Set up head
  call.str("");
  call << "import numpy as np" << '\n';
  call << "import psi4" << '\n';
  call << "import sys" << '\n';
  call << "import math" << '\n';
  call << "from psi4.driver import constants" << '\n';
  call << "psi4.set_memory('" << QMMMOpts.RAM; 
  if (QMMMOpts.memMB)
  {
    call << " mb";
  }
  else
  {
    call << " gb";
  }
  call << "')" << '\n';

  //Set up molecule
  call << "molA = psi4.geometry(" << '"' << '"' << '"' << '\n';
  call << " " << QMMMOpts.charge;
  call << " " << QMMMOpts.spin << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (QMMMData[i].QMRegion)
    {
      call << " " << QMMMData[i].QMTyp;
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].x*uConv,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].y*uConv,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].z*uConv,16);
      call << '\n';
    }
  }
  call << " symmetry c1" << '\n';
  call << " no_reorient" << '\n';
  call << " units bohr" << '\n';
  call << " no_com" << '\n';
  call << '"' << '"' << '"' << ")" << '\n';

  call << "molB = psi4.geometry(" << '"' << '"' << '"' << '\n';
  call << " " << QMMMOpts.charge;
  call << " " << QMMMOpts.spin << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (QMMMData[i].MMRegion)
    {
      call << " " << QMMMData[i].QMTyp;
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].x*uConv,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].y*uConv,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].z*uConv,16);
      call << '\n';
    }
  }
  call << " symmetry c1" << '\n';
  call << " no_reorient" << '\n';
  call << " units bohr" << '\n';
  call << " no_com" << '\n';
  call << '"' << '"' << '"' << ")" << '\n' << '\n';

  call << "methodname = '" << QMMMOpts.func << "'\n";
  call << "basisname = '" << QMMMOpts.basis << "'\n";
  call << "gembasisname = '" << QMMMOpts.gembasis << "'\n";
  call << "K_exchange = " << QMMMOpts.kexchange << "\n" << "\n";

  call << "epsilon = 1e-12" << '\n';
  call << "regularizer = 0.000001" << '\n';
  call << "constrain = True" << '\n' << '\n';

  call << "psi4.set_options({'puream' : False, 'print' : 1,'scf_type' : 'df'})" << '\n' << '\n';

  call << "eA, wfnA = psi4.energy(methodname+'/'+basisname, molecule=molA, return_wfn=True)" << '\n';
  if (QMMMOpts.prefitted)
    {
    call << "wfnB = psi4.core.Wavefunction.build(molB, basisname)" << '\n';
    }
  else
  {
    call << "eB, wfnB = psi4.energy(methodname+'/'+basisname, molecule=molB, return_wfn=True)" << '\n';
  }
  call << '\n';

  call << "def invert(mat):" << '\n';
  call << "    if epsilon == 0:" << '\n';
  call << "        c = np.linalg.inv(np.linalg.cholesky(mat))" << '\n';
  call << "        return np.dot(c.T,c)" << '\n';
  call << "    else:" << '\n';
  call << "        evals, evecs = np.linalg.eigh(mat)" << '\n';
  call << "        evals = np.where(evals < epsilon, 0.0, 1.0/evals)" << '\n';
  call << "        return np.einsum('ik,k,jk->ij', evecs, evals, evecs)" << '\n' << '\n';

  call << "def build_gem_field(wfn, add_coulomb, add_exchange, fit):" << '\n';
  call << "    orbital_basis = wfn.basisset()" << '\n';
  call << "    molecule = orbital_basis.molecule()" << '\n';
  call << "    aux_basis = psi4.core.BasisSet.build(molecule, 'ORBITAL', gembasisname)" << '\n';
  call << "    aux_basis.apply_hermite_normalization()" << '\n';
  call << "    zero_basis = psi4.core.BasisSet.zero_ao_basis_set()" << '\n';
  call << "    mints = psi4.core.MintsHelper(orbital_basis)" << '\n';
  call << "    field = psi4.QMMM().extern" << '\n';
  call << "    if fit:" << '\n';
  call << "        J_PQ = np.squeeze(mints.ao_eri(aux_basis, zero_basis, aux_basis, zero_basis))" << '\n';
  call << "        J_PQinv = invert(J_PQ + regularizer*np.eye(aux_basis.nbf()))" << '\n';
  call << "        J_Pmn = np.squeeze(mints.ao_eri(aux_basis, zero_basis, orbital_basis, orbital_basis))" << '\n';
  call << "        d =  2 * np.einsum('Qmn,mn->Q', J_Pmn, wfn.Da())" << '\n';
  call << "        fit_coefficients = np.einsum('PQ,Q->P', J_PQinv, d)" << '\n';
  call << "        if constrain:" << '\n';
  call << "            def dfact(n):" << '\n';
  call << "                if n <= 1:" << '\n';
  call << "                    return 1.0" << '\n';
  call << "                if n%2 == 1:" << '\n';
  call << "                    k = (n+1)//2" << '\n';
  call << "                    return math.factorial(2*k) / (2**k * math.factorial(k))" << '\n';
  call << "                else:" << '\n';
  call << "                    k = n//2" << '\n';
  call << "                    return 2**k * math.factorial(k)" << '\n';
  call << "            nshell = aux_basis.nshell()" << '\n';
  call << "            q = []" << '\n';
  call << "            for sh in range(nshell):" << '\n';
  call << "                shell = aux_basis.shell(sh)" << '\n';
  call << "                nprim = shell.nprimitive" << '\n';
  call << "                if nprim > 1:" << '\n';
  call << "                    raise('This code currently assumes uncontracted GEM basis sets')" << '\n';
  call << "                L = shell.am" << '\n';
  call << "                ex = shell.exp(0)" << '\n';
  call << "                coef = shell.coef(0)" << '\n';
  call << "                prefac = coef * np.sqrt((np.pi**3)/(2**L * ex**(L+3)))" << '\n';
  call << "                for lx in range(L,-1,-1):" << '\n';
  call << "                    for lz in range(L-lx+1):" << '\n';
  call << "                        ly = L - lx - lz" << '\n';
  call << "                        if (lx%2==0) and (ly%2==0) and (lz%2==0):" << '\n';
  call << "                            q.append(prefac*dfact(lx-1)*dfact(ly-1)*dfact(lz-1))" << '\n';
  call << "                        else:" << '\n';
  call << "                            q.append(0.0)" << '\n';
  call << "            natom = molecule.natom()" << '\n';
  call << "            Q = 0" << '\n';
  call << "            for i in range(natom):" << '\n';
  call << "                Q += molecule.Z(i)" << '\n';
  call << "            Q -= molecule.molecular_charge()" << '\n';
  call << "            old_nelec = np.dot(fit_coefficients, q)" << '\n';
  call << "            psi4.core.print_out(f'Number of electrons before constraining: {old_nelec:.6f}')" << '\n';
  call << "            qdot = np.dot(J_PQinv, q)" << '\n';
  call << "            lam = (Q - np.dot(qdot,  d)) / np.dot(qdot, q)" << '\n';
  call << "            fit_coefficients += lam * qdot" << '\n';
  call << "            new_nelec = np.dot(fit_coefficients, q)" << '\n';
  call << "            psi4.core.print_out(f'\\nNumber of electrons after constraining: {new_nelec:.6f}')" << '\n';
  call << "    else:" << '\n';
  call << "        fit_coefficients = np.loadtxt('fitcoeffs.dat')" << '\n';
  call << "    coefs_vec = psi4.core.Vector(aux_basis.nbf())" << '\n';
  call << "    for n, coef in enumerate(fit_coefficients): coefs_vec.set(n, coef)" << '\n';
  call << "    field.addBasis(aux_basis, coefs_vec)" << '\n';
  call << "    if add_exchange:" << '\n';
  call << "        K_coefs_vec = psi4.core.Vector(aux_basis.nbf())" << '\n';
  call << "        for n, coef in enumerate(fit_coefficients): K_coefs_vec.set(n, K_exchange*coef)" << '\n';
  call << "        field.addExchangeBasis(aux_basis, K_coefs_vec)" << '\n';
  call << "    for n in range(molecule.natom()):" << '\n';
  call << "        field.addCharge(molecule.Z(n), molecule.x(n), molecule.y(n), molecule.z(n))" << '\n';
  call << "    return field" << '\n' << '\n';

  if (QMMMOpts.prefitted)
  {
  call << "field = build_gem_field(wfnB, add_coulomb=True, add_exchange=False, fit=False)" << '\n';
  }
  else
  {
  call << "field = build_gem_field(wfnB, add_coulomb=True, add_exchange=False, fit=True)" << '\n';
  }
  call << "psi4.core.set_global_option_python('EXTERN', field)" << '\n';
  call << "eA_Jperturbed = psi4.energy(methodname+'/'+basisname, molecule=molA)" << '\n' << '\n';

  if (QMMMOpts.prefitted)
  {
  call << "field = build_gem_field(wfnB, add_coulomb=True, add_exchange=True, fit=False)" << '\n';
  }
  else
  {
  call << "field = build_gem_field(wfnB, add_coulomb=True, add_exchange=True, fit=True)" << '\n';
  }
  
  call << "psi4.core.set_global_option_python('EXTERN', field)" << '\n';
  call << "eA_JKperturbed = psi4.energy(methodname+'/'+basisname, molecule=molA)" << '\n' << '\n';

  call << "psi4.core.print_out('\\n Coulomb Energy: {:6} kcal/mol'.format(constants.hartree2kcalmol*(eA_Jperturbed-eA)))" << '\n';
  call << "psi4.core.print_out('\\n Exchange Energy: {:6} kcal/mol\\n'.format(constants.hartree2kcalmol*(eA_JKperturbed-eA_Jperturbed)))" << '\n';
  call << "psi4.core.print_out('\\n QM+GEM Energy: {:6} \\n'.format(eA_JKperturbed))" << '\n';

  //Create file
  dummy = call.str(); //Store file as a temporary variable
  call.str("");
  call << "newLICHM_" << bead << ".py";
  outFile.open(call.str().c_str(),ios_base::out);
  outFile << dummy << '\n';
  outFile.flush();
  outFile.close();

  return;
};
//E:JORGE

void WritePSI4Input(vector<QMMMAtom>& QMMMData, string calcTyp,
                    QMMMSettings& QMMMOpts, int bead)
{
  //Write PSI4 input files
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy,chrgfilename; //Generic string
  fstream inFile,outFile; //Generic file names
  //Check units
  double uConv = 1; //Units conversion constant
  if (QMMMOpts.unitsQM == "Bohr")
  {
    uConv = 1.0/bohrRad;
  }
  //Check for a charge file
  bool useChargeFile = 0;
  call.str("");
  call << "MMCharges_" << bead << ".txt";
  chrgfilename = call.str();
  useChargeFile = CheckFile(call.str());
  if (Nmm == 0)
  {
    //Skip blank charge files
    useChargeFile = 0;
  }
  //Initialize multipoles and center of mass
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
        //Set up multipoles
        RotateTINKCharges(QMMMData,bead);
      }
    }
  }
  //Check if there is a checkpoint file
  bool useCheckPoint;
  call.str("");
  call << "LICHM_" << bead << ".180";
  useCheckPoint = CheckFile(call.str());
  //Set up memory
  call.str("");
  call << "set_num_threads(" << Ncpus << ")" << '\n';
  call << "memory " << QMMMOpts.RAM;
  if (QMMMOpts.memMB)
  {
    call << " mb";
  }
  else
  {
    call << " gb";
  }
  call << '\n';
  //Set options
  if (QMMMOpts.spin == 1)
  {
    //Closed shell reference
    call << "set reference rhf";
  }
  else
  {
    //Open shell reference
    call << "set reference uhf";
  }
  call << '\n';
  call << "set basis ";
  call << QMMMOpts.basis << '\n';
  if (useCheckPoint)
  {
    call << "set guess read";
  }
  else
  {
    call << "set guess sad";
  }
  call << '\n';
  call << "set scf_type df" << '\n';
  if (QMMMOpts.unitsQM == "Bohr")
  {
    //Use Bohr for distances
    call << "set units bohr" << '\n';
  }
  call << '\n';
  //Keep the checkpoint files
  //NB: MOs->180
  call << "psi4_io.set_specific_path(180,'./')" << '\n';
  call << "psi4_io.set_specific_retention(180,True)" << '\n';
  call << '\n';
  //Set up molecules
  call << "molecule LICHM_";
  call << bead << " {" << '\n';
  call << " " << QMMMOpts.charge;
  call << " " << QMMMOpts.spin << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (QMMMData[i].QMRegion)
    {
      call << " " << QMMMData[i].QMTyp;
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].x*uConv,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].y*uConv,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].z*uConv,16);
      call << '\n';
    }
  }
  call << " symmetry c1" << '\n';
  call << " no_reorient" << '\n';
  call << " no_com" << '\n';
  call << "}" << '\n' << '\n';
  //Set up MM field
  if (QMMM and useChargeFile and (Nmm > 0))
  {
    inFile.open(chrgfilename.c_str(),ios_base::in);
    if (inFile.good())
    {
      call << "Chrgfield = QMMM()" << '\n';
      while (!inFile.eof())
      {
        //Copy charge file line by line
        getline(inFile,dummy);
        call << dummy << '\n';
      }
      inFile.close();
      call << "psi4.set_global_option_python('EXTERN',Chrgfield.extern)";
      call << '\n' << '\n';
    }
  }
  else if (QMMM and (Nmm > 0))
  {
    if (CHRG)
    {
      call << "Chrgfield = QMMM()" << '\n';
      for (int i=0;i<Natoms;i++)
      {
        if (QMMMData[i].MMRegion)
        {
          //Check PBC (minimum image convention)
          Coord distCent; //Distance from QM COM
          double xShft = 0;
          double yShft = 0;
          double zShft = 0;
          if (PBCon or QMMMOpts.useLREC)
          {
            //Initialize displacements
            double dx,dy,dz; //Starting displacements
            dx = QMMMData[i].P[bead].x-QMCOM.x;
            dy = QMMMData[i].P[bead].y-QMCOM.y;
            dz = QMMMData[i].P[bead].z-QMCOM.z;
            distCent = CoordDist2(QMMMData[i].P[bead],QMCOM);
            //Calculate the shift in positions
            //NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xShft = distCent.x-dx;
              yShft = distCent.y-dy;
              zShft = distCent.z-dz;
            }
          }
          //Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.useLREC)
          {
            //Use the long-range correction
            scrq = LRECFunction(distCent,QMMMOpts);
          }
          if (scrq > 0)
          {
            //Add charge
            double tmpX,tmpY,tmpZ,tmpQ; //Temporary storage
            tmpX = (QMMMData[i].P[bead].x+xShft)*uConv;
            tmpY = (QMMMData[i].P[bead].y+yShft)*uConv;
            tmpZ = (QMMMData[i].P[bead].z+zShft)*uConv;
            tmpQ = QMMMData[i].MP[bead].q*scrq;
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(tmpQ,16) << ",";
            call << LICHEMFormFloat(tmpX,16) << ",";
            call << LICHEMFormFloat(tmpY,16) << ",";
            call << LICHEMFormFloat(tmpZ,16) << ")";
            call << '\n';
          }
        }
      }
      call << "psi4.set_global_option_python('EXTERN',Chrgfield.extern)";
      call << '\n' << '\n';
    }
    if (AMOEBA)
    {
      call << "Chrgfield = QMMM()" << '\n';
      for (int i=0;i<Natoms;i++)
      {
        if (QMMMData[i].MMRegion)
        {
          //Check PBC (minimum image convention)
          Coord distCent; //Distance from QM COM
          double xShft = 0;
          double yShft = 0;
          double zShft = 0;
          if (PBCon or QMMMOpts.useLREC)
          {
            //Initialize displacements
            double dx,dy,dz; //Starting displacements
            dx = QMMMData[i].P[bead].x-QMCOM.x;
            dy = QMMMData[i].P[bead].y-QMCOM.y;
            dz = QMMMData[i].P[bead].z-QMCOM.z;
            distCent = CoordDist2(QMMMData[i].P[bead],QMCOM);
            //Calculate the shift in positions
            //NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xShft = distCent.x-dx;
              yShft = distCent.y-dy;
              zShft = distCent.z-dz;
            }
          }
          //Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.useLREC)
          {
            //Use the long-range correction
            scrq = LRECFunction(distCent,QMMMOpts);
          }
          if (scrq > 0)
          {
            //Add multipoles
            double tmpX,tmpY,tmpZ,tmpQ; //Temporary storage
            //Charge 1
            tmpX = (QMMMData[i].PC[bead].x1+xShft)*uConv;
            tmpY = (QMMMData[i].PC[bead].y1+yShft)*uConv;
            tmpZ = (QMMMData[i].PC[bead].z1+zShft)*uConv;
            tmpQ = QMMMData[i].PC[bead].q1*scrq;
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(tmpQ,16) << ",";
            call << LICHEMFormFloat(tmpX,16) << ",";
            call << LICHEMFormFloat(tmpY,16) << ",";
            call << LICHEMFormFloat(tmpZ,16) << ")";
            call << '\n';
            //Charge 2
            tmpX = (QMMMData[i].PC[bead].x2+xShft)*uConv;
            tmpY = (QMMMData[i].PC[bead].y2+yShft)*uConv;
            tmpZ = (QMMMData[i].PC[bead].z2+zShft)*uConv;
            tmpQ = QMMMData[i].PC[bead].q2*scrq;
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(tmpQ,16) << ",";
            call << LICHEMFormFloat(tmpX,16) << ",";
            call << LICHEMFormFloat(tmpY,16) << ",";
            call << LICHEMFormFloat(tmpZ,16) << ")";
            call << '\n';
            //Charge 3
            tmpX = (QMMMData[i].PC[bead].x3+xShft)*uConv;
            tmpY = (QMMMData[i].PC[bead].y3+yShft)*uConv;
            tmpZ = (QMMMData[i].PC[bead].z3+zShft)*uConv;
            tmpQ = QMMMData[i].PC[bead].q3*scrq;
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(tmpQ,16) << ",";
            call << LICHEMFormFloat(tmpX,16) << ",";
            call << LICHEMFormFloat(tmpY,16) << ",";
            call << LICHEMFormFloat(tmpZ,16) << ")";
            call << '\n';
            //Charge 4
            tmpX = (QMMMData[i].PC[bead].x4+xShft)*uConv;
            tmpY = (QMMMData[i].PC[bead].y4+yShft)*uConv;
            tmpZ = (QMMMData[i].PC[bead].z4+zShft)*uConv;
            tmpQ = QMMMData[i].PC[bead].q4*scrq;
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(tmpQ,16) << ",";
            call << LICHEMFormFloat(tmpX,16) << ",";
            call << LICHEMFormFloat(tmpY,16) << ",";
            call << LICHEMFormFloat(tmpZ,16) << ")";
            call << '\n';
            //Charge 5
            tmpX = (QMMMData[i].PC[bead].x5+xShft)*uConv;
            tmpY = (QMMMData[i].PC[bead].y5+yShft)*uConv;
            tmpZ = (QMMMData[i].PC[bead].z5+zShft)*uConv;
            tmpQ = QMMMData[i].PC[bead].q5*scrq;
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(tmpQ,16) << ",";
            call << LICHEMFormFloat(tmpX,16) << ",";
            call << LICHEMFormFloat(tmpY,16) << ",";
            call << LICHEMFormFloat(tmpZ,16) << ")";
            call << '\n';
            //Charge 6
            tmpX = (QMMMData[i].PC[bead].x6+xShft)*uConv;
            tmpY = (QMMMData[i].PC[bead].y6+yShft)*uConv;
            tmpZ = (QMMMData[i].PC[bead].z6+zShft)*uConv;
            tmpQ = QMMMData[i].PC[bead].q6*scrq;
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(tmpQ,16) << ",";
            call << LICHEMFormFloat(tmpX,16) << ",";
            call << LICHEMFormFloat(tmpY,16) << ",";
            call << LICHEMFormFloat(tmpZ,16) << ")";
            call << '\n';
          }
        }
      }
      call << "psi4.set_global_option_python('EXTERN',Chrgfield.extern)";
      call << '\n';
      call << '\n';
    }
    if (GEM)
    {
      //Add generic field field from a file (psithon)
      if (CheckFile("FIELD"))
      {
        //Read a block of psithon code
        inFile.open("FIELD",ios_base::in);
        while ((!inFile.eof()) and inFile.good())
        {
          getline(inFile,dummy);
          call << dummy << '\n';
        }
        //Save the field
        call << "psi4.set_global_option_python('EXTERN',Chrgfield.extern)";
        call << '\n';
        //Make sure the QM region is active
        call << "activate(LICHM_" << bead << ")" << '\n';
        call << '\n';
        inFile.close();
      }
    }
  }
  //Add calculation type
  call << calcTyp;
  //Create file
  dummy = call.str(); //Store file as a temporary variable
  call.str("");
  call << "LICHM_" << bead << ".dat";
  outFile.open(call.str().c_str(),ios_base::out);
  outFile << dummy << '\n';
  outFile.flush();
  outFile.close();
  return;
};

//Other QM files
void WriteQMConnect(int& argc,char**& argv)
{
  //Write the connectivity input for pure QM calculations
  fstream posfile,outFile; //File streams
  string dummy; //Generic string
  xyzFilename = "NOFILE"; //Global XYZ filename
  //Read arguments
  Natoms = 0; //For safety
  bool doQuit = 0; //Exit with an error
  cout << "Reading LICHEM input: ";
  for (int i=0;i<argc;i++)
  {
    dummy = string(argv[i]);
    //Check regions file
    if (dummy == "-q")
    {
      stringstream file;
      file << argv[i+1];
      if (!CheckFile(file.str()))
      {
        cout << "Error: Could not open XYZ file!!!";
        cout << '\n';
        doQuit = 1;
      }
      xyzFilename = file.str();
      posfile.open(argv[i+1],ios_base::in);
      cout << argv[i+1];
    }
  }
  cout << '\n' << '\n'; //Terminate output
  //Error check
  if (!CheckFile(xyzFilename))
  {
    cout << "Error: Missing XYZ file!!!";
    cout << '\n' << '\n';
    doQuit = 1;
  }
  if (!doQuit)
  {
    //Write connectivity information
    outFile.open("connect.inp",ios_base::out);
    posfile >> Natoms; //Number of atoms
    for (int i=0;i<Natoms;i++)
    {
      //Read the atom type
      string AtTyp;
      posfile >> AtTyp; //Read element
      //Clear junk position data
      posfile >> dummy >> dummy >> dummy;
      //Write connectivity line
      outFile << i << " "; //Index
      outFile << AtTyp << " "; //Element
      outFile << chemTable.revTyping(AtTyp) << " "; //Atomic number
      outFile << chemTable.getAtMass(AtTyp) << " "; //Mass
      outFile << "0.00 0" << '\n'; //Charge and bonds
    }
    cout << "Connectivity data written to connect.inp";
    cout << '\n' << '\n';
    cout.flush();
    outFile.flush();
    outFile.close();
  }
  //Quit
  posfile.close();
  exit(0);
  return;
};

