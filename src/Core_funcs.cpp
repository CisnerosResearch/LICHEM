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


 Primary functions for LICHEM.

*/

//Core utility functions
void PrintFancyTitle(fstream& logFile)
{
  logFile << '\n';
  logFile << "#######################################";
  logFile << "########################################";
  logFile << '\n';
  logFile << "#                                      ";
  logFile << "                                       #";
  logFile << '\n';
  logFile << "#                 ";
  logFile << "LICHEM: Layered Interacting CHEmical Models";
  logFile << "                 #";
  logFile << '\n';
  logFile << "#                                      ";
  logFile << "                                       #";
  logFile << '\n';
  logFile << "#                      ";
  logFile << "Symbiotic Computational Chemistry";
  logFile << "                      #";
  logFile << '\n';
  logFile << "#                                      ";
  logFile << "                                       #";
  logFile << '\n';
  logFile << "#######################################";
  logFile << "########################################";
  logFile << '\n' << '\n';
  logFile.flush();
  return;
};

double LICHEMFactorial(int n)
{
  //Calculate a factorial
  double val = 1;
  while (n > 1)
  {
    val *= n;
    n -= 1;
  }
  return val;
};

bool CheckFile(const string& file)
{
  //Checks if a file exists
  struct stat buffer;
  if (stat(file.c_str(),&buffer) != -1)
  {
    //Yep...
    return 1;
  }
  //Nope...
  return 0;
};

int FindMaxThreads()
{
  //Function to count the number of allowed threads
  int ct = 0; //Generic counter
  #pragma omp parallel reduction(+:ct)
  ct += 1; //Add one for each thread
  #pragma omp barrier
  //Return total count
  return ct;
};

double Bohring(double ri)
{
  //Convert ri (Bohr) to Angstroms
  double r = 0;
  r = ri/bohrRad;
  return r;
};

Coord CoordDist2(Coord& a, Coord& b)
{
  //Signed displacements
  double dx = a.x-b.x;
  double dy = a.y-b.y;
  double dz = a.z-b.z;
  //Check PBC
  if (PBCon)
  {
    bool check = 1; //Continue checking PBC
    while (check)
    {
      //Overkill, bit it checks if atoms are wrapped multiple times
      check = 0; //Stop if there are no changes
      if (abs(dx) > (0.5*Lx))
      {
        //Update X displacement
        check = 1; //The displacement was updated
        if (dx > 0)
        {
          dx -= Lx;
        }
        else
        {
          dx += Lx;
        }
      }
      if (abs(dy) > (0.5*Ly))
      {
        //Update Y displacement
        check = 1; //The displacement was updated
        if (dy > 0)
        {
          dy -= Ly;
        }
        else
        {
          dy += Ly;
        }
      }
      if (abs(dz) > (0.5*Lz))
      {
        //Update Z displacement
        check = 1; //The displacement was updated
        if (dz > 0)
        {
          dz -= Lz;
        }
        else
        {
          dz += Lz;
        }
      }
    }
  }
  //Save displacements
  Coord dispAB; //Distance between A and B
  dispAB.x = dx;
  dispAB.y = dy;
  dispAB.z = dz;
  return dispAB;
};

double LRECFunction(Coord& dist, QMMMSettings& QMMMOpts)
{
  //Calculates the LREC scale factor based on the distance between points
  double R = dist.vecMag(); //Distance between the points
  //Check if the charge is within the cutoff
  if (R <= (QMMMOpts.LRECCut*QMMMOpts.LRECCut))
  {
    //Scale the charge
    double scrq = 1;
    R = sqrt(R);
    double scrqA,scrqB; //Temporary variables
    //Calculate temp. variables
    scrqA = (QMMMOpts.LRECCut-R)/QMMMOpts.LRECCut;
    scrqB = -3*scrqA*scrqA;
    scrqA *= 2*scrqA*scrqA;
    //Combine temp. variables
    scrqA += scrqB+1;
    //Set the scale factor
    scrq -= pow(scrqA,QMMMOpts.LRECPow);
    return scrq;
  }
  //Charge is not within the cutoff
  return 0;
};

//Functions to check connectivity
//vector<int> TraceBoundary(vector<QMMMAtom>& QMMMData, int atID)
vector<int> TraceBoundary(vector<QMMMAtom>& QMMMData, int atID, int& mystat,fstream& logFile)
{
  //Function to find all boundary atoms connected to a pseudobond atom
  bool bondError = 0; //Checks if the molecular structure "breaks" the math
  vector<int> boundAtoms; //Final list of atoms
  //Add atoms bonded to atom "AtID"
  for (unsigned int i=0;i<QMMMData[atID].bonds.size();i++)
  {
    int bondID = QMMMData[atID].bonds[i];
    if (QMMMData[bondID].BARegion)
    {
      boundAtoms.push_back(bondID);
    }
    if (QMMMData[bondID].PBRegion and (bondID != atID))
    {
      //Two PBs are connected and this system will fail
      bondError = 1;
    }
  }
  //Check find other boundary atoms bonded to the initial set
  bool moreFound = 1; //More boundary atoms were found
  while (moreFound and (!bondError))
  {
    moreFound = 0; //Break the loop
    vector<int> tmp;
    for (unsigned int i=0;i<boundAtoms.size();i++)
    {
      int BAID = boundAtoms[i];
      for (unsigned int j=0;j<QMMMData[BAID].bonds.size();j++)
      {
        int bondID = QMMMData[BAID].bonds[j];
        //Check if it is on the list
        bool isThere = 0;
        for (unsigned int k=0;k<boundAtoms.size();k++)
        {
          if (bondID == boundAtoms[k])
          {
            isThere = 1;
          }
        }
        if (!isThere)
        {
          if (QMMMData[bondID].BARegion)
          {
            moreFound = 1; //Keep going
            tmp.push_back(bondID);
          }
          if (QMMMData[bondID].PBRegion and (bondID != atID))
          {
            //Two PBs are connected and this system will fail
            bondError = 1;
          }
        }
      }
    }
    //Add them to the list
    for (unsigned int i=0;i<tmp.size();i++)
    {
      bool isThere = 0;
      for (unsigned int j=0;j<boundAtoms.size();j++)
      {
        //Avoid adding the atom twice
        if (tmp[i] == boundAtoms[j])
        {
          isThere = 1;
        }
      }
      if (!isThere)
      {
        boundAtoms.push_back(tmp[i]);
      }
    }
  }
  //Check for errors
  if (bondError)
  {
    logFile << "Error: Two pseudo-bonds are connected through boudary atoms!!!";
    logFile << '\n';
    logFile << " The connections prevent LICHEM from correctly updating";
    logFile << " the charges" << '\n' << '\n';
    logFile.flush();
    //Quit to avoid unphysical results
    //exit(0);
    //return boundAtoms;
  }
  //Return list if there are no errors
  return boundAtoms;
};

bool Bonded(vector<QMMMAtom>& QMMMData, int atom1, int atom2)
{
  //Function to check if two atoms are 1-2 connected
  bool isBound = 0;
  for (unsigned int i=0;i<QMMMData[atom2].bonds.size();i++)
  {
    if (atom1 == QMMMData[atom2].bonds[i])
    {
      isBound = 1;
    }
  }
  return isBound;
};

bool Angled(vector<QMMMAtom>& QMMMData, int atom1, int atom3)
{
  //Function to check if two atoms are 1-3 connected
  bool isBound = 0;
  for (unsigned int i=0;i<QMMMData[atom1].bonds.size();i++)
  {
    int atom2 = QMMMData[atom1].bonds[i];
    for (unsigned int j=0;j<QMMMData[atom2].bonds.size();j++)
    {
      if (atom3 == QMMMData[atom2].bonds[j])
      {
        isBound = 1;
      }
    }
  }
  return isBound;
};

bool Dihedraled(vector<QMMMAtom>& QMMMData, int atom1, int atom4)
{
  //Function to check if two atoms are 1-4 connected
  bool isBound = 0;
  for (unsigned int i=0;i<QMMMData[atom1].bonds.size();i++)
  {
    int atom2 = QMMMData[atom1].bonds[i];
    for (unsigned int j=0;j<QMMMData[atom2].bonds.size();j++)
    {
      int atom3 = QMMMData[atom2].bonds[j];
      for (unsigned int k=0;k<QMMMData[atom3].bonds.size();k++)
      {
        if (atom4 == QMMMData[atom3].bonds[k])
        {
          isBound = 1;
        }
      }
    }
  }
  return isBound;
};

//Structure correction functions
void PBCCenter(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts)
{
  //Move the system to the center of the simulation box
  double avgX = 0;
  double avgY = 0;
  double avgZ = 0;
  #pragma omp parallel
  {
    #pragma omp for nowait schedule(dynamic) reduction(+:avgX)
    for (int i=0;i<Natoms;i++)
    {
      //Loop over all beads
      double centX = 0;
      for (int j=0;j<QMMMOpts.NBeads;j++)
      {
        //Update local average postion
        centX += QMMMData[i].P[j].x;
      }
      //Upate full average
      avgX += centX;
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:avgY)
    for (int i=0;i<Natoms;i++)
    {
      //Loop over all beads
      double centY = 0;
      for (int j=0;j<QMMMOpts.NBeads;j++)
      {
        //Update local average postion
        centY += QMMMData[i].P[j].y;
      }
      //Upate full average
      avgY += centY;
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:avgZ)
    for (int i=0;i<Natoms;i++)
    {
      //Loop over all beads
      double centZ = 0;
      for (int j=0;j<QMMMOpts.NBeads;j++)
      {
        //Update local average postion
        centZ += QMMMData[i].P[j].z;
      }
      //Upate full average
      avgZ += centZ;
    }
  }
  #pragma omp barrier
  //Convert sums to averages
  avgX /= Natoms*QMMMOpts.NBeads;
  avgY /= Natoms*QMMMOpts.NBeads;
  avgZ /= Natoms*QMMMOpts.NBeads;
  //Move atoms to the center of the box
  #pragma omp parallel
  {
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<Natoms;i++)
    {
      //Loop over all beads
      for (int j=0;j<QMMMOpts.NBeads;j++)
      {
        //Move bead to the center
        QMMMData[i].P[j].x -= avgX;
        QMMMData[i].P[j].x += 0.5*Lx;
      }
    }
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<Natoms;i++)
    {
      //Loop over all beads
      for (int j=0;j<QMMMOpts.NBeads;j++)
      {
        //Move bead to the center
        QMMMData[i].P[j].y -= avgY;
        QMMMData[i].P[j].y += 0.5*Ly;
      }
    }
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<Natoms;i++)
    {
      //Loop over all beads
      for (int j=0;j<QMMMOpts.NBeads;j++)
      {
        //Move bead to the center
        QMMMData[i].P[j].z -= avgZ;
        QMMMData[i].P[j].z += 0.5*Lz;
      }
    }
  }
  #pragma omp barrier
  //Return with updated structure
  return;
};

Coord FindQMCOM(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts, int bead)
{
  //Find the center of mass for the QM region
  Coord QMCOM; //Center of mass position
  double avgX = 0; //Average x position
  double avgY = 0; //Average y position
  double avgZ = 0; //Average z position
  double totM = 0; //Total mass
  #pragma omp parallel
  {
    #pragma omp for nowait schedule(dynamic) reduction(+:totM)
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
      {
        totM += QMMMData[i].m;
      }
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:avgX)
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
      {
        avgX += QMMMData[i].m*QMMMData[i].P[bead].x;
      }
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:avgY)
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
      {
        avgY += QMMMData[i].m*QMMMData[i].P[bead].y;
      }
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:avgZ)
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
      {
        avgZ += QMMMData[i].m*QMMMData[i].P[bead].z;
      }
    }
  }
  #pragma omp barrier
  //Save center of mass
  QMCOM.x = avgX/totM;
  QMCOM.y = avgY/totM;
  QMCOM.z = avgZ/totM;
  return QMCOM;
};

//Misc.
void PrintLapin(fstream& logFile)
{
  //Print a nice picture
  stringstream lapin;
  lapin.str("");
  //Add obfuscated text
  lapin << '\n';
  lapin << "                   /^\\";
  lapin << "     /^\\";
  lapin << '\n';
  lapin << "                   |  \\";
  lapin << "   /  |";
  lapin << '\n';
  lapin << "                   |  |";
  lapin << "   |  |";
  lapin << '\n';
  lapin << "                   |";
  lapin << "  |___|  |";
  lapin << '\n';
  lapin << "                  /";
  lapin << "           \\";
  lapin << '\n';
  lapin << "                 /";
  lapin << "   &";
  lapin << "     &";
  lapin << "   \\";
  lapin << '\n';
  lapin << "                 |";
  lapin << "             |";
  lapin << '\n';
  lapin << "                 |";
  lapin << "     .";
  lapin << " .";
  lapin << "     |";
  lapin << '\n';
  lapin << "                  \\";
  lapin << "   \\___/";
  lapin << "   /";
  lapin << '\n';
  lapin << "                   \\_";
  lapin << "       _/";
  lapin << '\n';
  lapin << "                     |";
  lapin << "     |";
  lapin << '\n' << '\n';
  lapin << "         ______";
  lapin << "       ______";
  lapin << "       ______";
  lapin << '\n';
  lapin << "        /";
  lapin << "      \\";
  lapin << "     /";
  lapin << "      \\";
  lapin << "     /";
  lapin << "      \\";
  lapin << '\n';
  lapin << "       /";
  lapin << "        \\";
  lapin << "   /";
  lapin << "        \\";
  lapin << "   /";
  lapin << "        \\";
  lapin << '\n';
  lapin << "      |";
  lapin << "----------";
  lapin << "|";
  lapin << " |";
  lapin << "          |";
  lapin << " |";
  lapin << "**********";
  lapin << "|";
  lapin << '\n';
  lapin << "      [";
  lapin << "          ]";
  lapin << " [";
  lapin << "%%%%%%%%%%";
  lapin << "]";
  lapin << " [";
  lapin << "          ]";
  lapin << '\n';
  lapin << "      |";
  lapin << "----------";
  lapin << "|";
  lapin << " |";
  lapin << "          |";
  lapin << " |";
  lapin << "**********";
  lapin << "|";
  lapin << '\n';
  lapin << "       \\";
  lapin << "        /";
  lapin << "   \\";
  lapin << "        /";
  lapin << "   \\";
  lapin << "        /";
  lapin << '\n';
  lapin << "        \\";
  lapin << "______";
  lapin << "/";
  lapin << "     \\";
  lapin << "______";
  lapin << "/";
  lapin << "     \\";
  lapin << "______";
  lapin << "/";
  lapin << '\n' << '\n' << '\n';
  //Print text and return
  logFile << lapin.str();
  logFile.flush();
  return;
};

void FetchQuotes(vector<string>& Quotes)
{
  //Generate random quotes
  string dummy; //Generic string
  //dummy = "\'It is difficult to prove that this quote is not random.\'";
  //dummy += '\n';
  //dummy += "                                           -Eric G. Kratz";
  for (int i=0;i<1000;i++)
  {
    //Add quotes to the list
    Quotes.push_back(dummy);
  }
  return;
};

