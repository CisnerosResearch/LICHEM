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

 Path integral functions using QM, MM, and QMMM energies. Calls to wrappers
 are parallel over the number of beads. Other functions are mostly parallel
 over the number of atoms.

*/

//Path integral Monte Carlo functions
double Get_PI_Espring(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts)
{
  //Calculate total harmonic PI ring energy
  double E = 0.0; //Final energy
  double wZero; //Mass-independent force constant
  wZero = 1/(QMMMOpts.beta*hbar);
  wZero *= wZero*toeV*QMMMOpts.NBeads;
  #pragma omp parallel for schedule(dynamic) reduction(+:E)
  for (int i=0;i<Natoms;i++)
  {
    QMMMData[i].Ep = 0.0; //Reset safed energy
    double w = wZero*QMMMData[i].m; //Mass-scaled force constant
    for (int j=0;j<QMMMOpts.NBeads;j++)
    {
      //Bead energy, one bond to avoid double counting
      int j2 = j-1;
      if (j2 == -1)
      {
        j2 = QMMMOpts.NBeads-1; //Ring PBC
      }
      //Calculate displacement with PBC
      double dr2; //Squared displacement
      dr2 = CoordDist2(QMMMData[i].P[j],QMMMData[i].P[j2]).vecMag();
      QMMMData[i].Ep += 0.5*w*dr2; //Harmonic energy
    }
    E += QMMMData[i].Ep; //Save energy
  }
  return E;
};

double Get_PI_Epot(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,fstream& logFile)
{
  //Potential for all beads
  double E = 0.0;
  //Fix parallel for classical MC
  int mcThreads = Nthreads;
  if (QMMMOpts.NBeads == 1)
  {
    mcThreads = 1;
  }
  //Calculate energy
  #pragma omp parallel for schedule(dynamic) num_threads(mcThreads) \
          reduction(+:E,QMTime,MMTime)
  for (int p=0;p<QMMMOpts.NBeads;p++)
  {
    //Run the wrappers for all beads
    double Es = 0.0;
    //Timer variables
    int t_qm_start = 0;
    int t_mm_start = 0;
    int times_qm = 0;
    int times_mm = 0;
    //Calculate QM energy
    if (Gaussian)
    {
      t_qm_start = (unsigned)time(0);
      Es += GaussianEnergy(QMMMData,QMMMOpts,p);
      times_qm += (unsigned)time(0)-t_qm_start;
    }
    if (PSI4)
    {
      t_qm_start = (unsigned)time(0);
      Es += PSI4Energy(QMMMData,QMMMOpts,p);
      times_qm += (unsigned)time(0)-t_qm_start;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      t_qm_start = (unsigned)time(0);
      Es += NWChemEnergy(QMMMData,QMMMOpts,p);
      times_qm += (unsigned)time(0)-t_qm_start;
    }
    //Calculate MM energy
    if (TINKER)
    {
      t_mm_start = (unsigned)time(0);
      Es += TINKEREnergy(QMMMData,QMMMOpts,p,logFile);
      times_mm += (unsigned)time(0)-t_mm_start;
    }
    if (LAMMPS)
    {
      t_mm_start = (unsigned)time(0);
      Es += LAMMPSEnergy(QMMMData,QMMMOpts,p);
      times_mm += (unsigned)time(0)-t_mm_start;
    }
    //Add temp variables to the totals
    E += Es;
    QMTime += times_qm;
    MMTime += times_mm;
  }
  //Calculate the average energy
  E /= QMMMOpts.NBeads;
  return E;
};

bool MCMove(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts, double& Emc,fstream& logFile)
{
  //Function to perform Monte Carlo moves and accept/reject the moves
  bool acc = 0; //Accept or reject
  //Copy QMMMData
  vector<QMMMAtom> QMMMData2;
  QMMMData2 = QMMMData;
  //Pick random move and apply PBC
  double randNum = (((double)rand())/((double)RAND_MAX));
  if (randNum > (1-centProb))
  {
    //Move a centroid
    int p;
    bool frozenAt = 1;
    while (frozenAt)
    {
      //Make sure the atom is not frozen
      p = (rand()%Natoms);
      if (QMMMData2[p].frozen == 0)
      {
        frozenAt = 0;
      }
    }
    double randX = (((double)rand())/((double)RAND_MAX));
    double randY = (((double)rand())/((double)RAND_MAX));
    double randZ = (((double)rand())/((double)RAND_MAX));
    double dx = 2*(randX-0.5)*mcStep*centRatio;
    double dy = 2*(randY-0.5)*mcStep*centRatio;
    double dz = 2*(randZ-0.5)*mcStep*centRatio;
    //Update positions
    #pragma omp parallel
    {
      #pragma omp for nowait schedule(dynamic)
      for (int i=0;i<QMMMOpts.NBeads;i++)
      {
        QMMMData2[p].P[i].x += dx;
      }
      #pragma omp for nowait schedule(dynamic)
      for (int i=0;i<QMMMOpts.NBeads;i++)
      {
        QMMMData2[p].P[i].y += dy;
      }
      #pragma omp for nowait schedule(dynamic)
      for (int i=0;i<QMMMOpts.NBeads;i++)
      {
        QMMMData2[p].P[i].z += dz;
      }
    }
    #pragma omp barrier
  }
  if (randNum < beadProb)
  {
    //Move all beads in a centroid
    int p;
    bool frozenAt = 1;
    while (frozenAt)
    {
      //Make sure the atom is not frozen
      p = (rand()%Natoms);
      if (QMMMData2[p].frozen == 0)
      {
        frozenAt = 0;
      }
    }
    for (int i=0;i<QMMMOpts.NBeads;i++)
    {
      //Randomly displace each bead
      double randX = (((double)rand())/((double)RAND_MAX));
      double randY = (((double)rand())/((double)RAND_MAX));
      double randZ = (((double)rand())/((double)RAND_MAX));
      double dx = 2*(randX-0.5)*mcStep;
      double dy = 2*(randY-0.5)*mcStep;
      double dz = 2*(randZ-0.5)*mcStep;
      QMMMData2[p].P[i].x += dx;
      QMMMData2[p].P[i].y += dy;
      QMMMData2[p].P[i].z += dz;
    }
  }
  //Initialize energies
  double EOld = QMMMOpts.EOld;
  double ENew = 0;
  //Save box lengths
  double LxSave = Lx;
  double LySave = Ly;
  double LzSave = Lz;
  //Attempt a volume move
  randNum = (((double)rand())/((double)RAND_MAX));
  if (randNum < volProb)
  {
    //Anisotropic volume change
    if (isotrop == 0)
    {
      //Assumes that MM cutoffs are safe
      randNum = (((double)rand())/((double)RAND_MAX));
      Lx += 2*(randNum-0.5)*mcStep;
      randNum = (((double)rand())/((double)RAND_MAX));
      Ly += 2*(randNum-0.5)*mcStep;
      randNum = (((double)rand())/((double)RAND_MAX));
      Lz += 2*(randNum-0.5)*mcStep;
    }
    //Isotropic volume change
    if (isotrop == 1)
    {
      //Assumes that MM cutoffs are safe
      randNum = (((double)rand())/((double)RAND_MAX));
      Lx += 2*(randNum-0.5)*mcStep;
      Ly += 2*(randNum-0.5)*mcStep;
      Lz += 2*(randNum-0.5)*mcStep;
    }
    //Decide how to scale the centroids
    bool scaleRing = 0; //Shift the ring
    randNum = (((double)rand())/((double)RAND_MAX));
    if (randNum >= 0.5)
    {
      //Evenly scale the size of the ring
      scaleRing = 1;
    }
    //Scale positions
    if (scaleRing)
    {
      #pragma omp parallel
      {
        #pragma omp for nowait schedule(dynamic)
        for (int i=0;i<Natoms;i++)
        {
          for (int j=0;j<QMMMOpts.NBeads;j++)
          {
            double shift;
            shift = QMMMData2[i].P[j].x;
            //Check PBC without wrapping the molecules
            bool check = 1; //Continue the PBC checks
            while (check)
            {
              //Check the value
              check = 0;
              if (shift > Lx)
              {
                shift -= Lx;
                check = 1;
              }
              if (shift < 0)
              {
                shift += Lx;
                check = 1;
              }
            }
            shift = ((Lx/LxSave)-1)*shift;
            QMMMData2[i].P[j].x += shift;
          }
        }
        #pragma omp for nowait schedule(dynamic)
        for (int i=0;i<Natoms;i++)
        {
          for (int j=0;j<QMMMOpts.NBeads;j++)
          {
            double shift;
            shift = QMMMData2[i].P[j].y;
            //Check PBC without wrapping the molecules
            bool check = 1; //Continue the PBC checks
            while (check)
            {
              //Check the value
              check = 0;
              if (shift > Ly)
              {
                shift -= Ly;
                check = 1;
              }
              if (shift < 0)
              {
                shift += Ly;
                check = 1;
              }
            }
            shift = ((Ly/LySave)-1)*shift;
            QMMMData2[i].P[j].y += shift;
          }
        }
        #pragma omp for nowait schedule(dynamic)
        for (int i=0;i<Natoms;i++)
        {
          for (int j=0;j<QMMMOpts.NBeads;j++)
          {
            double shift;
            shift = QMMMData2[i].P[j].z;
            //Check PBC without wrapping the molecules
            bool check = 1; //Continue the PBC checks
            while (check)
            {
              //Check the value
              check = 0;
              if (shift > Lz)
              {
                shift -= Lz;
                check = 1;
              }
              if (shift < 0)
              {
                shift += Lz;
                check = 1;
              }
            }
            shift = ((Lz/LzSave)-1)*shift;
            QMMMData2[i].P[j].z += shift;
          }
        }
      }
      #pragma omp barrier
    }
    else
    {
      #pragma omp parallel
      {
        #pragma omp for nowait schedule(dynamic)
        for (int i=0;i<Natoms;i++)
        {
          //Find centroids
          double shift = 0; //Change of position for the centroid
          for (int j=0;j<QMMMOpts.NBeads;j++)
          {
            shift += QMMMData2[i].P[j].x; //Add to the position sum
          }
          shift /= QMMMOpts.NBeads; //Average position
          //Check PBC without wrapping the molecules
          bool check = 1; //Continue the PBC checks
          while (check)
          {
            //Check the value
            check = 0;
            if (shift > Lx)
            {
              shift -= Lx;
              check = 1;
            }
            if (shift < 0)
            {
              shift += Lx;
              check = 1;
            }
          }
          //Calculate the change in position
          shift = ((Lx/LxSave)-1)*shift;
          for (int j=0;j<QMMMOpts.NBeads;j++)
          {
            //Update the position
            QMMMData2[i].P[j].x += shift;
          }
        }
        #pragma omp for nowait schedule(dynamic)
        for (int i=0;i<Natoms;i++)
        {
          //Find centroids
          double shift = 0; //Change of position for the centroid
          for (int j=0;j<QMMMOpts.NBeads;j++)
          {
            shift += QMMMData2[i].P[j].y; //Add to the position sum
          }
          shift /= QMMMOpts.NBeads; //Average position
          //Check PBC without wrapping the molecules
          bool check = 1; //Continue the PBC checks
          while (check)
          {
            //Check the value
            check = 0;
            if (shift > Ly)
            {
              shift -= Ly;
              check = 1;
            }
            if (shift < 0)
            {
              shift += Ly;
              check = 1;
            }
          }
          //Calculate the change in position
          shift = ((Ly/LySave)-1)*shift;
          for (int j=0;j<QMMMOpts.NBeads;j++)
          {
            //Update the position
            QMMMData2[i].P[j].y += shift;
          }
        }
        #pragma omp for nowait schedule(dynamic)
        for (int i=0;i<Natoms;i++)
        {
          //Find centroids
          double shift = 0; //Change of position for the centroid
          for (int j=0;j<QMMMOpts.NBeads;j++)
          {
            shift += QMMMData2[i].P[j].z; //Add to the position sum
          }
          shift /= QMMMOpts.NBeads; //Average position
          //Check PBC without wrapping the molecules
          bool check = 1; //Continue the PBC checks
          while (check)
          {
            //Check the value
            check = 0;
            if (shift > Lz)
            {
              shift -= Lz;
              check = 1;
            }
            if (shift < 0)
            {
              shift += Lz;
              check = 1;
            }
          }
          //Calculate the change in position
          shift = ((Lz/LzSave)-1)*shift;
          for (int j=0;j<QMMMOpts.NBeads;j++)
          {
            //Update the position
            QMMMData2[i].P[j].z += shift;
          }
        }
      }
      #pragma omp barrier
    }
  }
  //Update energies
  ENew += Get_PI_Epot(QMMMData2,QMMMOpts,logFile);
  ENew += Get_PI_Espring(QMMMData2,QMMMOpts);
  if (QMMMOpts.ensemble == "NPT")
  {
    //Add PV energy term
    ENew += QMMMOpts.press*Lx*Ly*Lz;
  }
  //Accept or reject
  double dE = QMMMOpts.beta*(ENew-EOld);
  if (QMMMOpts.ensemble == "NPT")
  {
    //Add Nln(V) term
    double volTerm;
    volTerm = Lx*Ly*Lz; //New volume
    volTerm /= LxSave*LySave*LzSave; //Divide by old volume
    volTerm = log(volTerm); //Take the natural logarithm
    volTerm *= Natoms*QMMMOpts.NBeads; //Scale by number of particles
    dE -= volTerm; //Subtract from the energy
  }
  double prob = exp(-1*dE);
  randNum = (((double)rand())/((double)RAND_MAX));
  if (randNum <= prob)
  {
    //Accept
    QMMMData = QMMMData2;
    Emc = ENew;
    QMMMOpts.EOld = ENew;
    acc = 1;
  }
  else
  {
    //Reject
    Emc = EOld;
    //Revert to old box sizes
    Lx = LxSave;
    Ly = LySave;
    Lz = LzSave;
  }
  //Return decision
  return acc;
};

