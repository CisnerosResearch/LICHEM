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
  # Functions for TINKER in parallel LICHEM                                     #
  # Includes:                                                                   #
  #                                                                             #
  #            for master proc:                                                 #
  #                                                                             #
  #                            Writing input files :                            #
  #                                   void TINKERForcesMPIWrite                 #
  #                                   void TINKERPolForcesMPIWrite              #
  #                                   double TINKEREnergyMPIWrite               #
  #                                   double TINKERPolEnergyMPIWrite            #
  #                                   void TINKEROptMPIWrite                    #
  #                                                                             #
  #                            Reading output files:                            #
  #                                   void TINKERForcesMPIRead                  #
  #                                   void TINKERPolForcesMPIRead               #
  #                                   double TINKEREnergyMPIRead                #
  #                                   double TINKERPolEnergyMPIRead             #
  #                                   void TINKEROptMPIRead                     #
  #                                                                             #
  #            for all procs:                                                   #
  #                            Running local beads :                            #
  #                                   void TINKERForcesMPI                      #
  #                                   void TINKERPolForcesMPI                   #
  #                                   double TINKEREnergyMPI                    #
  #                                   double TINKERPolEnergyMPI                 #
  #                                   void TINKEROptMPI                         #
  #                                                                             #
  ###############################################################################
*/


void TINKERForcesMPIWrite(vector<QMMMAtom>& QMMMData,
                    QMMMSettings& QMMMOpts, int bead,int& mystat,fstream& logFile)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream outFile,inFile; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double Emm = 0.0;
  int ct; //Generic counter
  //Construct MM forces input for TINKER
  call.str("");
  call << "cp tinker.key LICHM_TINKERForces_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_TINKERForces_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n';
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
  }
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
    {
      if (ct == 0)
      {
        //Start a new active line
        outFile << "active ";
      }
      else
      {
        //Place a space to separate values
        outFile << " ";
      }
      outFile << (QMMMData[i].id+1);
      ct += 1;
      if (ct == 10)
      {
        //terminate an active line
        ct = 0;
        outFile << '\n';
      }
    }
  }
  if (ct != 0)
  {
    //Terminate trailing actives line
    outFile << '\n';
  }
  outFile << "group-inter" << '\n'; //Modify interactions
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add group 1 atoms
    if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
    {
      if (ct == 0)
      {
        //Start a new group line
        outFile << "group 1 ";
      }
      else
      {
        //Place a space to separate values
        outFile << " ";
      }
      outFile << (QMMMData[i].id+1);
      ct += 1;
      if (ct == 10)
      {
        //terminate a group line
        ct = 0;
        outFile << '\n';
      }
    }
  }
  if (ct != 0)
  {
    //Terminate trailing group line
    outFile << '\n';
  }
  if (CHRG)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion or QMMMData[i].BARegion)
      {
        //New charges are needed for QM and PB atoms
        outFile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        outFile << "0.0"; //Delete charges
        outFile << '\n';
      }
    }
  }
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion or QMMMData[i].BARegion)
      {
        double qi = 0;
        //Remove charge
        qi = QMMMData[i].MP[bead].q;
        QMMMData[i].MP[bead].q = 0;
        WriteTINKMPole(QMMMData,outFile,i,bead);
        QMMMData[i].MP[bead].q += qi; //Restore charge
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
    }
  }
  outFile.flush();
  outFile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_TINKERForces_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,12);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,12);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,12);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
};

void TINKERPolForcesMPIWrite(vector<QMMMAtom>& QMMMData, 
                       QMMMSettings& QMMMOpts, int bead,int& mystat,fstream& logFile)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream outFile,inFile; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double Emm = 0.0;
  int ct; //Generic counter
  //Construct MM forces input for TINKER
  call.str("");
  call << "cp tinker.key LICHM_TINKERPolForces_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_TINKERPolForces_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n';
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else if (!QMMMOpts.useImpSolv)
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
  }
  if (AMOEBA)
  {
    //Get rid of non-polarization interactions
    outFile << "polarizeterm only" << '\n';
  }
  if (QMMMOpts.useImpSolv)
  {
    //Add the implicit solvation model
    outFile << "solvateterm" << '\n';
    outFile << "solvate " << QMMMOpts.solvModel;
    outFile << '\n';
  }
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
    {
      if (ct == 0)
      {
        //Start a new active line
        outFile << "active ";
      }
      else
      {
        //Place a space to separate values
        outFile << " ";
      }
      outFile << (QMMMData[i].id+1);
      ct += 1;
      if (ct == 10)
      {
        //terminate an active line
        ct = 0;
        outFile << '\n';
      }
    }
  }
  if (ct != 0)
  {
    //Terminate trailing actives line
    outFile << '\n';
  }
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMRegion)
      {
        WriteTINKMPole(QMMMData,outFile,i,bead);
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].PBRegion)
      {
        //Modify the charge to force charge balance with the boundaries
        double qi = QMMMData[i].MP[bead].q; //Save a copy
        vector<int> boundaries;
        boundaries = TraceBoundary(QMMMData,i,mystat,logFile);
        double qNew = qi;
        for (unsigned int j=0;j<boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qNew -= QMMMData[boundaries[j]].MP[bead].q;
        }
        QMMMData[i].MP[bead].q = qNew; //Save modified charge
        WriteTINKMPole(QMMMData,outFile,i,bead);
        QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].BARegion)
      {
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
    }
  }
  outFile.flush();
  outFile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_TINKERPolForces_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,12);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,12);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,12);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();

};/*finish write*/  

void TINKEREnergyMPIWrite(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                    int bead)
{
  //Runs TINKER MM energy calculations
  fstream outFile,inFile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  double E = 0.0;
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  call.str("");
  call << "cp tinker.key LICHM_TINKEREnergy_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_TINKEREnergy_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n';
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else if (!QMMMOpts.useImpSolv)
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
  }
  if (QMMM)
  {
    outFile << "polarizeterm none" << '\n'; //Remove polarization energy
  }
  else if (QMMMOpts.useImpSolv)
  {
    //Add the implicit solvation model for pure MM calculations
    outFile << "solvate " << QMMMOpts.solvModel;
    outFile << '\n';
  }
  ct = 0; //Generic counter
  if (QMMM)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if (QMMMData[i].MMRegion or QMMMData[i].BARegion)
      {
        if (ct == 0)
        {
          //Start a new active line
          outFile << "active ";
        }
        else
        {
          //Place a space to separate values
          outFile << " ";
        }
        outFile << (QMMMData[i].id+1);
        ct += 1;
        if (ct == 10)
        {
          //terminate an active line
          ct = 0;
          outFile << '\n';
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      outFile << '\n';
    }
  }
  if (CHRG)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion or QMMMData[i].BARegion)
      {
        //New charges are only needed for QM atoms
        outFile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        outFile << "0.0" << '\n';
      }
    }
  }
  if (AMOEBA or GEM)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add multipoles
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion or QMMMData[i].BARegion)
      {
        double qi = 0;
        //Remove charge
        qi = QMMMData[i].MP[bead].q;
        QMMMData[i].MP[bead].q = 0;
        //Write new multipole definition for the atom ID
        WriteTINKMPole(QMMMData,outFile,i,bead);
        //Restore charge
        QMMMData[i].MP[bead].q = qi;
        //Remove polarization
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
    }
  }
  outFile.flush();
  outFile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_TINKEREnergy_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,12);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,12);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,12);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
};

void TINKERPolEnergyMPIWrite(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                       int bead,int& mystat,fstream& logFile)
{
  //Function to extract the polarization energy
  fstream outFile,inFile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  double EPol = 0; //Polarization energy
  double ESolv = 0; //Solvation energy
  double E = 0; //Total energy for error checking
  int ct; //Generic counter
  //Create TINKER xyz file
  call.str("");
  call << "LICHM_TINKERPolEnergy_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,12);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,12);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,12);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
  //Create new TINKER key file
  call.str("");
  call << "cp tinker.key LICHM_TINKERPolEnergy_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_TINKERPolEnergy_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n';
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else if (!QMMMOpts.useImpSolv)
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
  }
  if (AMOEBA)
  {
    //Get rid of non-polarization interactions
    outFile << "polarizeterm only" << '\n';
  }
  if (QMMMOpts.useImpSolv)
  {
    //Add the implicit solvation model
    outFile << "solvateterm" << '\n';
    outFile << "solvate " << QMMMOpts.solvModel;
    outFile << '\n';
  }
  ct = 0; //Generic counter
  if (QMMM)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if (QMMMData[i].MMRegion or QMMMData[i].BARegion)
      {
        if (ct == 0)
        {
          //Start a new active line
          outFile << "active ";
        }
        else
        {
          //Place a space to separate values
          outFile << " ";
        }
        outFile << (QMMMData[i].id+1);
        ct += 1;
        if (ct == 10)
        {
          //terminate an active line
          ct = 0;
          outFile << '\n';
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      outFile << '\n';
    }
  }

  for (int i=0;i<Natoms;i++)
  {
    //Add nuclear charges
    if (QMMMData[i].QMRegion)
    {
      //Write new multipole definition for the atom ID
      WriteTINKMPole(QMMMData,outFile,i,bead);
      outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
      outFile << '\n';
    }
    if (QMMMData[i].PBRegion)
    {
      //Modify the charge to force charge balance with the boundaries
      double qi = QMMMData[i].MP[bead].q; //Save a copy
      vector<int> boundaries;
      boundaries = TraceBoundary(QMMMData,i,mystat,logFile);
      double qNew = qi;
      for (unsigned int j=0;j<boundaries.size();j++)
      {
        //Subtract boundary atom charge
        qNew -= QMMMData[boundaries[j]].MP[bead].q;
      }
      QMMMData[i].MP[bead].q = qNew; //Save modified charge
      WriteTINKMPole(QMMMData,outFile,i,bead);
      QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
      outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
      outFile << '\n';
    }
    if (QMMMData[i].BARegion)
    {
      outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
      outFile << '\n';
    }
  }
  outFile.flush();
  outFile.close();
  
};

void TINKEROptMPIWrite(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts, int bead,
                       int& mystat,fstream& logFile)
{

  //Runs TINKER MM optimization
  fstream outFile,inFile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  double E = 0;
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  call.str("");
  call << "cp tinker.key LICHM_TINKEROpt_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_TINKEROpt_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n';
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useMMCut)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald and truncate vdW forces
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.MMOptCut,12);
      outFile << '\n';
      outFile << "ewald" << '\n';
    }
    else
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.MMOptCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.MMOptCut,12);
      outFile << '\n';
    }
  }
  else if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else if (!QMMMOpts.useImpSolv)
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  if (QMMMOpts.useImpSolv)
  {
    //Add the implicit solvation model
    outFile << "solvate " << QMMMOpts.solvModel;
    outFile << '\n';
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
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
           return;
       }*/
       /*End: Hatice GOKCAN*/

        if (!QMMMData[i].frozen)
        {
          if (ct == 0)
          {
            //Start a new active line
            outFile << "active ";
          }
          else
          {
            //Place a space to separate values
            outFile << " ";
          }
          outFile << (QMMMData[i].id+1);
          ct += 1;
          if (ct == 10)
          {
            //terminate an active line
            ct = 0;
            outFile << '\n';
          }
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      outFile << '\n';
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
        outFile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        outFile << QMMMData[i].MP[bead].q;
        outFile << '\n';
      }
      if (QMMMData[i].PBRegion)
      {
        //Modify the charge to force charge balance with the boundaries
        vector<int> boundaries;
        boundaries = TraceBoundary(QMMMData,i,mystat,logFile);
        double qNew = QMMMData[i].MP[bead].q;
        for (unsigned int j=0;j<boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qNew -= QMMMData[boundaries[j]].MP[bead].q;
        }
        outFile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        outFile << qNew;
        outFile << '\n';
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
        WriteTINKMPole(QMMMData,outFile,i,bead);
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].PBRegion)
      {
        //Modify the charge to force charge balance with the boundaries
        double qi = QMMMData[i].MP[bead].q; //Save a copy
        vector<int> boundaries;
        boundaries = TraceBoundary(QMMMData,i,mystat,logFile);
        double qNew = qi;
        for (unsigned int j=0;j<boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qNew -= QMMMData[boundaries[j]].MP[bead].q;
        }
        QMMMData[i].MP[bead].q = qNew; //Save modified charge
        WriteTINKMPole(QMMMData,outFile,i,bead);
        QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].BARegion)
      {
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
    }
  }
  outFile.flush();
  outFile.close();
  
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_TINKEROpt_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,12);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,12);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,12);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();

}
void TINKEROptRestrainMPIWrite(vector<QMMMAtom>& QMMMData, 
                               QMMMSettings& QMMMOpts, int bead,
                               double restr,int& mystat,fstream& logFile)
{

  fstream ofile,ifile; 
  stringstream call; 
  call.copyfmt(cout);
  string dummy; 
  double E = 0;
  int ct;
  call.str("");
  /*Copy the original key file and make changes*/
  call.str("");
  call << "cp tinker.key LICHM_TINKEROpt_";
  call << bead << ".key";

  globalSys = system(call.str().c_str());
  /*Update key file*/
  call.str("");
  call << "LICHM_TINKEROpt_";
  call << bead << ".key";
  ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
  ofile << '\n';
  ofile << "#QM force field parameters"; 
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
  ct = 0; 
  if (QMMM or (Nfreeze > 0))
  {
    for (int i=0;i<Natoms;i++)
    {
      /*Add active atoms*/
      if (QMMMData[i].MMRegion or QMMMData[i].BARegion)
      {

       /*if( (Nmm + Nbound) != Nfreeze ){
           logFile << "             ";
           logFile << "Error:\n";
           logFile << "             ";
           logFile << "Total number of MM atoms and boundary atoms\n";
           logFile << "             ";
           logFile << "are not equal to the number of frozen atoms.\n";
           mystat=1;
           return;
       }*/
        if (!QMMMData[i].frozen)
        {
          if (ct == 0)
          {
            /*Start a new active line*/
            ofile << "active ";
          }
          else
          {
            /*Place a space to separate values*/
            ofile << " ";
          }
          ofile << (QMMMData[i].id+1);
          ct += 1;
          if (ct == 10)
          {
            /*terminate an active line*/
            ct = 0;
            ofile << '\n';
          }
        }
      }
    }
    if (ct != 0)
    {
      /*Terminate trailing actives line*/
      ofile << '\n';
    }
  }
  if (CHRG)
  {
    for (int i=0;i<Natoms;i++)
    {
      /*Add nuclear charges*/
      if (QMMMData[i].QMRegion)
      {
        /*New charges are only needed for QM atoms*/
        ofile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        ofile << QMMMData[i].MP[bead].q;
        ofile << '\n';
      }
      if (QMMMData[i].PBRegion)
      {
        /*Modify the charge to force charge balance with the boundaries*/
        vector<int> Boundaries;
        Boundaries = TraceBoundary(QMMMData,i,mystat,logFile);
        double qnew = QMMMData[i].MP[bead].q;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          /*Subtract boundary atom charge*/
          qnew -= QMMMData[Boundaries[j]].MP[bead].q;
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
      /*Add nuclear charges*/
      if (QMMMData[i].QMRegion)
      {
        /*Write new multipole definition for the atom ID*/
        WriteTINKMPole(QMMMData,ofile,i,bead);
        ofile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
      if (QMMMData[i].PBRegion)
      {
        /*Modify the charge to force charge balance with the boundaries*/
        double qi = QMMMData[i].MP[bead].q; 
        vector<int> Boundaries;
        Boundaries = TraceBoundary(QMMMData,i,mystat,logFile);
        double qnew = qi;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          /*Subtract boundary atom charge*/
          qnew -= QMMMData[Boundaries[j]].MP[bead].q;
        }
        QMMMData[i].MP[bead].q = qnew; /*Save modified charge*/
        WriteTINKMPole(QMMMData,ofile,i,bead);
        QMMMData[i].MP[bead].q = qi; /*Return to unmodified charge*/
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

/*Start: Hatice
 * input restraints*/
  
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
            ofile << LICHEMFormFloat(QMMMData[i].P[bead].x,12);
            ofile << " ";
            ofile << LICHEMFormFloat(QMMMData[i].P[bead].y,12);
            ofile << " ";
            ofile << LICHEMFormFloat(QMMMData[i].P[bead].z,12);
            ofile << " ";
            ofile << setw(6) << restr;
            ofile << '\n';
        }
      }
    }

  }

/*End: Hatice*/
  ofile.flush();
  ofile.close();

  /*Create TINKER xyz file from the structure*/
  call.str("");
  call << "LICHM_TINKEROpt_" << bead << ".xyz";
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
    ofile << LICHEMFormFloat(QMMMData[i].P[bead].x,12);
    ofile << " ";
    ofile << LICHEMFormFloat(QMMMData[i].P[bead].y,12);
    ofile << " ";
    ofile << LICHEMFormFloat(QMMMData[i].P[bead].z,12);
    ofile << " ";
    ofile << setw(4) << QMMMData[i].numTyp;

    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      ofile << " "; /*Avoids trailing spaces*/
      ofile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    ofile.copyfmt(cout); /*Copy settings from cout*/
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();


}
/*-----------------------*/
/*-----------------------*/

void TINKERForcesMPIRead(vector<QMMMAtom>& QMMMData, VectorXd& forces,
                    QMMMSettings& QMMMOpts, int bead)
{  
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  //Collect MM forces
  fstream MMGrad; //QMMM output
  //Open files
  call.str("");
  call << "LICHM_TINKERForces_" << bead << ".grad";
  MMGrad.open(call.str().c_str(),ios_base::in);
  //Read derivatives
  bool gradDone = 0;
  double Emm=0.0;
  
  while ((!MMGrad.eof()) and MMGrad.good() and (!gradDone))
  {
    getline(MMGrad,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Type")
    {
      line >> dummy >> dummy;
      if (dummy == "dE/dX")
      {
        gradDone = 1; //Not grad school, that lasts forever
        getline(MMGrad,dummy);
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
          double fX = 0;
          double fY = 0;
          double fZ = 0;
          //Convoluted, but "easy"
          getline(MMGrad,dummy);
          stringstream line(dummy);
          line >> dummy >> dummy; //Clear junk
          line >> fX;
          line >> fY;
          line >> fZ;
          //Change from gradient to force
          fX *= -1;
          fY *= -1;
          fZ *= -1;
          //Switch to eV/A and save forces
          forces(3*i) += fX*kcal2eV;
          forces(3*i+1) += fY*kcal2eV;
          forces(3*i+2) += fZ*kcal2eV;
        }
      }
    }
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "Energy")
      {
        //Collect partial MM energy
        line >> dummy >> Emm;
      }
    }
  }
  MMGrad.close();
  //Clean up files
  if(!QMMMOpts.KeepFiles){
    call.str("");
    call << "rm -f";
    call << " LICHM_TINKERForces_" << bead << ".xyz";
    call << " LICHM_TINKERForces_" << bead << ".key";
    call << " LICHM_TINKERForces_" << bead << ".grad";
    call << " LICHM_TINKERForces_" << bead << ".err";
    globalSys = system(call.str().c_str());
  }
  //Return
  Emm *= kcal2eV;
};

void TINKERPolForcesMPIRead(vector<QMMMAtom>& QMMMData, VectorXd& forces,
                       QMMMSettings& QMMMOpts, int bead)
{

  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  double Emm = 0.0;
 
  //Collect MM forces
  fstream MMGrad; //QMMM output
  //Open files
  call.str("");
  call << "LICHM_TINKERPolForces_" << bead << ".grad";
  MMGrad.open(call.str().c_str(),ios_base::in);
  //Read derivatives
  bool gradDone = 0;
  while ((!MMGrad.eof()) and MMGrad.good() and (!gradDone))
  {
    getline(MMGrad,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Type")
    {
      line >> dummy >> dummy;
      if (dummy == "dE/dX")
      {
        gradDone = 1; //Not grad school, that lasts forever
        getline(MMGrad,dummy);
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
          double fX = 0;
          double fY = 0;
          double fZ = 0;
          //Convoluted, but "easy"
          getline(MMGrad,dummy);
          stringstream line(dummy);
          line >> dummy >> dummy; //Clear junk
          line >> fX;
          line >> fY;
          line >> fZ;
          //Change from gradient to force
          fX *= -1;
          fY *= -1;
          fZ *= -1;
          //Switch to eV/A and change sign
          forces(3*i) += fX*kcal2eV;
          forces(3*i+1) += fY*kcal2eV;
          forces(3*i+2) += fZ*kcal2eV;
        }
      }
    }
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "Energy")
      {
        //Collect partial MM energy
        line >> dummy >> Emm;
      }
    }
  }
  MMGrad.close();
  //Clean up files
  if(!QMMMOpts.KeepFiles){
    call.str("");
    call << "rm -f";
    call << " LICHM_TINKERPolForces_" << bead << ".xyz";
    call << " LICHM_TINKERPolForces_" << bead << ".key";
    call << " LICHM_TINKERPolForces_" << bead << ".grad";
    call << " LICHM_TINKERPolForces_" << bead << ".err";
    globalSys = system(call.str().c_str());
  }

  //Return
  Emm *= kcal2eV;
};

double TINKEREnergyMPIRead(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                    int bead)  
{  
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  fstream inFile;
  double E=0.0;

  call.str("");
  call << "LICHM_TINKEREnergy_" << bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  //Read MM potential energy
  bool EFound = 0;
  while ((!inFile.eof()) and inFile.good())
  {
    inFile >> dummy;
    if (dummy == "Total")
    {
      inFile >> dummy >> dummy;
      if (dummy == "Energy")
      {
        inFile >> dummy >> E;
        EFound = 1;
      }
    }
  }
  if (!EFound)
  {
    //Warn user if no energy was found
    cerr << "Warning: No MM energy found after a calculation!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    E = hugeNum; //Large number to reject step
  }
  inFile.close();
  //Clean up files
  if(!QMMMOpts.KeepFiles){
    call.str("");
    call << "rm -f";
    call << " LICHM_TINKEREnergy_" << bead << ".xyz";
    call << " LICHM_TINKEREnergy_" << bead << ".log";
    call << " LICHM_TINKEREnergy_" << bead << ".key";
    call << " LICHM_TINKEREnergy_" << bead << ".err";
    globalSys = system(call.str().c_str());
  }

  return E;
};

double TINKERPolEnergyMPIRead(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                    int bead)  
{

  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  fstream inFile;
  double EPol = 0; //Polarization energy
  double ESolv = 0; //Solvation energy
  double E = 0; //Total energy for error checking

  //Extract polarization energy
  call.str("");
  call << "LICHM_TINKERPolEnergy_" << bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  bool EFound = 0;
  while ((!inFile.eof()) and inFile.good())
  {
    inFile >> dummy;
    if (dummy == "Total")
    {
      inFile >> dummy >> dummy;
      if (dummy == "Energy")
      {
        inFile >> dummy >> E;
        EFound = 1;
      }
    }
    if (dummy == "Polarization")
    {
      inFile >> EPol;
    }
    if (dummy == "Implicit")
    {
      inFile >> dummy;
      if (dummy == "Solvation")
      {
        inFile >> ESolv;
      }
    }
  }
  if (!EFound)
  {
    //Warn user if no energy was found
    cerr << "Warning: No MM energy found after a calculation!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    EPol = 0; //Prevents errors when polarization is off
    ESolv = 0; //Prevents errors when implicit solvation is off
  }
  inFile.close();
  //Clean up files
  
  if(!QMMMOpts.KeepFiles){
    call.str("");
    call << "rm -f";
    call << " LICHM_TINKERPolEnergy_" << bead << ".xyz";
    call << " LICHM_TINKERPolEnergy_" << bead << ".log";
    call << " LICHM_TINKERPolEnergy_" << bead << ".key";
    call << " LICHM_TINKERPolEnergy_" << bead << ".err";
    globalSys = system(call.str().c_str());
  }
  //Return polarization and solvation energy in kcal/mol
  double Emmpol=(EPol+ESolv);//*kcal2eV;
  //cout << "Emm2=" << Emmpol << "\n";
  //return EPol+ESolv;
  return Emmpol;
};

void TINKEROptMPIRead(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts, int bead)
{

  fstream outFile,inFile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  double E = 0;

  //Read new structure
  call.str("");
  call << "LICHM_TINKEROpt_" << bead << ".xyz_2";
  inFile.open(call.str().c_str(),ios_base::in);
  getline(inFile,dummy); //Discard number of atoms
  if (PBCon)
  {
    //Discard PBC information
    getline(inFile,dummy);
  }
  for (int i=0;i<Natoms;i++)
  {
    getline(inFile,dummy);
    stringstream line(dummy);
    //Read new positions
    line >> dummy >> dummy; //Discard atom ID and type
    line >> QMMMData[i].P[bead].x;
    line >> QMMMData[i].P[bead].y;
    line >> QMMMData[i].P[bead].z;
  }
  inFile.close();
  //Clean up files
  
  if(!QMMMOpts.KeepFiles){
    call.str("");
    call << "rm -f";
    call << " LICHM_TINKEROpt_" << bead << ".xyz";
    call << " LICHM_TINKEROpt_" << bead << ".log";
    call << " LICHM_TINKEROpt_" << bead << ".xyz_*";
    call << " LICHM_TINKEROpt_" << bead << ".key";
    call << " LICHM_TINKEROpt_" << bead << ".err";
    globalSys = system(call.str().c_str());
  }
  //Change units
  E *= kcal2eV;
};

/*-----------------------*/
void TINKERForcesMPI(vector<int> mybead_list,
                  int mysize,int pathstart, int pathend)
{ 

  //Function for calculating the forces on a set of atoms
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  string dummy; //Generic string
  fstream QMLog; //Generic input files

  int myrank,wsize;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  
  MPI_Status stat;


  for(int jj=0; jj<mysize;jj++)
  {

        int p=mybead_list[jj];

        /* ensure that it starts from pathstart*/
        /* in root proc */
        /* so that no extra calc will be performed */
        if(p<pathstart){
           p=pathstart;
        }
        /* ensure that it ends at pathend*/
        /* in last proc */
        /* so that no extra calc will be performed */
        if(p>pathend){
           /* opt is finished, break loop */
           /* and convert force to hartree */
           break;
        }

        /* input file */ 
        call.str("");
        call << "testgrad ";
        call << "LICHM_TINKERForces_" << p << ".xyz";
        call << " Y N N > ";
        call << "LICHM_TINKERForces_" << p << ".grad";

        globalSys = system(call.str().c_str());

  }

  int lastdone=0;
  if(myrank==wsize-1){
    lastdone=1;
    MPI_Send(&lastdone, 0, MPI_INT, 0, 44, MPI_COMM_WORLD);
  }
  if(myrank==0){
    MPI_Recv(&lastdone, 0, MPI_INT, wsize-1, 44, MPI_COMM_WORLD, &stat);
  }


  /* ensure that no one exits the loop before finishing */
  int value;
  if(myrank==0){
    value=1;
    for(int i = 1; i < wsize; i++){
       MPI_Send(&value, 0, MPI_INT, i, 42, MPI_COMM_WORLD);
    }

  }
  else{
    value=-1;
    MPI_Recv(&value, 0, MPI_INT, 0, 42, MPI_COMM_WORLD, &stat);
  }
  MPI_Barrier(MPI_COMM_WORLD);



}; 

void TINKERPolForcesMPI(vector<int> mybead_list,
                  int mysize,int pathstart, int pathend)
{ 

  //Function for calculating the forces on a set of atoms
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  string dummy; //Generic string
  fstream QMLog; //Generic input files

  int myrank,wsize;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  
  MPI_Status stat;

  for(int jj=0; jj<mysize;jj++)
  {

        int p=mybead_list[jj];

        /* ensure that it starts from pathstart*/
        /* in root proc */
        /* so that no extra calc will be performed */
        if(p<pathstart){
           p=pathstart;
        }
        /* ensure that it ends at pathend*/
        /* in last proc */
        /* so that no extra calc will be performed */
        if(p>pathend){
           /* opt is finished, break loop */
           /* and convert force to hartree */
           break;
        }

        /* input file */ 
        call.str("");
        call << "testgrad ";
        call << "LICHM_TINKERPolForces_" << p << ".xyz";
        call << " Y N N > ";
        call << "LICHM_TINKERPolForces_" << p << ".grad";

        globalSys = system(call.str().c_str());
  }

  int lastdone=0;
  if(myrank==wsize-1){
    lastdone=1;
    MPI_Send(&lastdone, 0, MPI_INT, 0, 44, MPI_COMM_WORLD);
  }
  if(myrank==0){
    MPI_Recv(&lastdone, 0, MPI_INT, wsize-1, 44, MPI_COMM_WORLD, &stat);
  }


  /* ensure that no one exits the loop before finishing */
  int value;
  if(myrank==0){
    value=1;
    for(int i = 1; i < wsize; i++){
       MPI_Send(&value, 0, MPI_INT, i, 42, MPI_COMM_WORLD);
    }

  }
  else{
    value=-1;
    MPI_Recv(&value, 0, MPI_INT, 0, 42, MPI_COMM_WORLD, &stat);
  }
  MPI_Barrier(MPI_COMM_WORLD);



};/* finish run */ 


void TINKEREnergyMPI(vector<int> mybead_list,
                     int mysize,int pathstart, int pathend)  
{  

  //Function for calculating the forces on a set of atoms
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  string dummy; //Generic string
  fstream QMLog; //Generic input files

  int myrank,wsize;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  
  MPI_Status stat;

  for(int jj=0; jj<mysize;jj++)
  {

        int p=mybead_list[jj];

        /* ensure that it starts from pathstart*/
        /* in root proc */
        /* so that no extra calc will be performed */
        if(p<pathstart){
           p=pathstart;
        }
        /* ensure that it ends at pathend*/
        /* in last proc */
        /* so that no extra calc will be performed */
        if(p>pathend){
           /* opt is finished, break loop */
           /* and convert force to hartree */
           break;
        }

        /* input file */ 
        call.str("");
        call << "analyze LICHM_TINKEREnergy_";
        call << p << ".xyz E > LICHM_TINKEREnergy_";
        call << p << ".log";


        globalSys = system(call.str().c_str());


  }

  int lastdone=0;
  if(myrank==wsize-1){
    lastdone=1;
    MPI_Send(&lastdone, 0, MPI_INT, 0, 44, MPI_COMM_WORLD);
  }
  if(myrank==0){
    MPI_Recv(&lastdone, 0, MPI_INT, wsize-1, 44, MPI_COMM_WORLD, &stat);
  }


  /* ensure that no one exits the loop before finishing */
  int value;
  if(myrank==0){
    value=1;
    for(int i = 1; i < wsize; i++){
       MPI_Send(&value, 0, MPI_INT, i, 42, MPI_COMM_WORLD);
    }

  }
  else{
    value=-1;
    MPI_Recv(&value, 0, MPI_INT, 0, 42, MPI_COMM_WORLD, &stat);
  }
  MPI_Barrier(MPI_COMM_WORLD);



};/* finish run */ 


void TINKERPolEnergyMPI(vector<int> mybead_list,
                     int mysize,int pathstart, int pathend)  
{  

  //Function for calculating the forces on a set of atoms
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  string dummy; //Generic string
  fstream QMLog; //Generic input files

  int myrank,wsize;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  
  MPI_Status stat;

  for(int jj=0; jj<mysize;jj++)
  {

        int p=mybead_list[jj];

        /* ensure that it starts from pathstart*/
        /* in root proc */
        /* so that no extra calc will be performed */
        if(p<pathstart){
           p=pathstart;
        }
        /* ensure that it ends at pathend*/
        /* in last proc */
        /* so that no extra calc will be performed */
        if(p>pathend){
           /* opt is finished, break loop */
           /* and convert force to hartree */
           break;
        }

        /* input file */ 
        call.str("");
        call << "analyze LICHM_TINKERPolEnergy_";
        call << p << ".xyz E > LICHM_TINKERPolEnergy_";
        call << p << ".log";

        globalSys = system(call.str().c_str());


  }

  int lastdone=0;
  if(myrank==wsize-1){
    lastdone=1;
    MPI_Send(&lastdone, 0, MPI_INT, 0, 44, MPI_COMM_WORLD);
  }
  if(myrank==0){
    MPI_Recv(&lastdone, 0, MPI_INT, wsize-1, 44, MPI_COMM_WORLD, &stat);
  }


  /* ensure that no one exits the loop before finishing */
  int value;
  if(myrank==0){
    value=1;
    for(int i = 1; i < wsize; i++){
       MPI_Send(&value, 0, MPI_INT, i, 42, MPI_COMM_WORLD);
    }

  }
  else{
    value=-1;
    MPI_Recv(&value, 0, MPI_INT, 0, 42, MPI_COMM_WORLD, &stat);
  }
  MPI_Barrier(MPI_COMM_WORLD);



};/* finish run */ 

void TINKEROptMPI(vector<int> mybead_list,
                  int mysize,int pathstart, int pathend,
                  QMMMSettings& QMMMOpts)
{ 

  //Function for calculating the forces on a set of atoms
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  string dummy; //Generic string
  fstream QMLog; //Generic input files

  int myrank,wsize;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &wsize);
  
  MPI_Status stat;


  for(int jj=0; jj<mysize;jj++)
  {

        int p=mybead_list[jj];

        /* ensure that it starts from pathstart*/
        /* in root proc */
        /* so that no extra calc will be performed */
        if(p<pathstart){
           p=pathstart;
        }
        /* ensure that it ends at pathend*/
        /* in last proc */
        /* so that no extra calc will be performed */
        if(p>pathend){
           /* opt is finished, break loop */
           /* and convert force to hartree */
           break;
        }

        /* input file */ 
        call.str("");
        call << "minimize LICHM_TINKEROpt_";
        call << p << ".xyz ";
        call << QMMMOpts.MMOptTol << " > LICHM_TINKEROpt_";
        call << p << ".log";

        globalSys = system(call.str().c_str());

  }


  int lastdone=0;
  if(myrank==wsize-1){
    lastdone=1;
    MPI_Send(&lastdone, 0, MPI_INT, 0, 44, MPI_COMM_WORLD);
  }
  if(myrank==0){
    MPI_Recv(&lastdone, 0, MPI_INT, wsize-1, 44, MPI_COMM_WORLD, &stat);
  }


  /* ensure that no one exits the loop before finishing */
  int value;
  if(myrank==0){
    value=1;
    for(int i = 1; i < wsize; i++){
       MPI_Send(&value, 0, MPI_INT, i, 42, MPI_COMM_WORLD);
    }

  }
  else{
    value=-1;
    MPI_Recv(&value, 0, MPI_INT, 0, 42, MPI_COMM_WORLD, &stat);
  }
  MPI_Barrier(MPI_COMM_WORLD);


};
