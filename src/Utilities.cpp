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
  # Functions for random utilities                                            #
  # Includes:                                                                 #
  #                                                                           #
  #       Reaction coordinate calculation : void calc_react_coord             #
  #                                                                           #
  #       Printing utility                : void print_progress               #
  #                                                                           #
  #       Keep files utility              : void save_files                   #
  #                                                                           #
  #       Find TS and complex beads &                                         #
  #       their energies                  : void getTSbead                    #
  #                                                                           #
  #############################################################################
*/

void calc_react_coord(QMMMSettings& QMMMOpts, vector<QMMMAtom>& QMMMData, 
                      VectorXd& reactCoord)
{
  // Calculate reaction coordinate
  for (int p=0;p<(QMMMOpts.NBeads-1);p++)
  {
    MatrixXd geom1((Nqm+Npseudo),3); // Current replica
    MatrixXd geom2((Nqm+Npseudo),3); // Next replica
    VectorXd disp; // Store the displacement
    // Save geometries
    int ct = 0; // Reset counter for the number of atoms
    for (int i=0;i<Natoms;i++)
    {
      // Only include QM and PB regions
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
      {
        // Save current replica
        geom1(ct,0) = QMMMData[i].P[p].x;
        geom1(ct,1) = QMMMData[i].P[p].y;
        geom1(ct,2) = QMMMData[i].P[p].z;
        // Save replica p+1
        geom2(ct,0) = QMMMData[i].P[p+1].x;
        geom2(ct,1) = QMMMData[i].P[p+1].y;
        geom2(ct,2) = QMMMData[i].P[p+1].z;
        ct += 1;
      }
    }
    // Calculate displacement
    disp = KabschDisplacement(geom1,geom2,(Nqm+Npseudo));
    // for NEB and QSM
    // Remove inactive atoms
    ct = 0; // Reset counter for the number of atoms
    for (int i=0;i<Natoms;i++)
    {
      // Only include QM and PB regions
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
      {
        // Only include active atoms in the tangent
        if (!QMMMData[i].NEBActive)
        {
          // Delete distance components
          disp(ct) = 0;
          disp(ct+1) = 0;
          disp(ct+2) = 0;
        }
        // Advance counter
        ct += 3;
      }
    }
    
    // Update reaction coordinate
    reactCoord(p+1) = reactCoord(p); // Start from previous bead
    reactCoord(p+1) += disp.norm(); // Add magnitude of the displacement
  }
  
  reactCoord /= reactCoord.maxCoeff(); // Must be between 0 and 1
  
}

/*-------------------------------------------------------------------------*/

void print_progress(QMMMSettings& QMMMOpts, int print_level, 
                    VectorXd& Eqmmm_images, double RMSdiff, double MAXforce, 
                    double RMSforce, VectorXd& reactCoord, fstream& logFile)
{

  double reactE = QMMMOpts.EReact;
  double prodE = QMMMOpts.EProd;
  double TSE = QMMMOpts.ETrans; 
  // Send print_level to this function
  //-------------------------------------------
  // Print bead energies
  if (print_level==0)
  {
    // Bead energies (CalcEner and CalcForce)
    stringstream call; // Stream for system calls and reading/writing files
    call.copyfmt(logFile); // Save settings
    logFile << "               ";
    logFile << "Bead energies:\n";
    logFile << "                 ";
    logFile << "Bead |  Coord |       E (a.u)        |  RE (kcal/mol)  ";
    logFile << '\n';
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      logFile << right;
      logFile << fixed;
      logFile << "                 ";
      logFile << setw(4) << p << " |  ";
      logFile << setprecision(3);
      logFile << setw(5) << reactCoord(p) << " | ";
      logFile << setprecision(8);
      logFile << setw(20) << Eqmmm_images[p];
      logFile << " | ";
      logFile << setprecision(4);
      /* reactE)*627.51;//<< " kcal/mol"; */
      logFile << setw(14) << (Eqmmm_images[p] - QMMMOpts.EReact)*627.51;
      logFile << '\n'; 
      logFile.flush(); // Print progress  
    }
    logFile.copyfmt(call); // Replace settings
  }

  // -------------------------------------------
  // End of CalcForce
  // Print TS bead and barriers
  if (print_level==1)
  {
    stringstream call; // Stream for system calls and reading/writing files
    call.copyfmt(logFile); // Save settings
    logFile << fixed;
    logFile << "\n";
    logFile << "               ";
    logFile << "| TS Bead          : ";
    logFile << QMMMOpts.TSBead << " \n";
    logFile << "               ";
    logFile << "| Reactant Energy  : ";
    logFile << setprecision(8) << reactE << " a.u \n";
    logFile << "               ";
    logFile << "| TS Energy        : ";
    logFile << setprecision(8) << TSE << " a.u \n";
    logFile << "               ";
    logFile << "| Product Energy   : ";
    logFile << setprecision(8) << prodE << " a.u \n";
    // Start: Aug 28 2018
    /*
      logFile << "               ";
      logFile << "| Forward Barrier  : ";
      logFile << setprecision(8);
      logFile << (QMMMOpts.ETrans- QMMMOpts.EReact) << " a.u ";
      logFile << "   | ";
      logFile << setprecision(4);
      logFile << (QMMMOpts.ETrans- QMMMOpts.EReact)*627.51 << " kcal/mol";
      logFile << '\n';
      logFile << "               ";
      logFile << "| Backward Barrier : ";
      logFile << setprecision(8);
      logFile << (QMMMOpts.ETrans- QMMMOpts.EProd) << " a.u ";
      logFile << "   | ";
      logFile << setprecision(4);
      logFile << (QMMMOpts.ETrans- QMMMOpts.EProd)*627.51 << " kcal/mol";
      logFile << '\n'<< endl;
    */
    if (QMMMOpts.EComplforw!=QMMMOpts.EReact)
    {
      logFile << "               ";
      logFile << "| Forward Complex  : ";
      logFile << setprecision(8) <<  QMMMOpts.EComplforw << " a.u \n"; 
    }
    if (QMMMOpts.EComplback!=QMMMOpts.EProd)
    {
      logFile << "               ";
      logFile << "| Backward Complex : ";
      logFile << setprecision(8) <<  QMMMOpts.EComplback << " a.u \n";
    }
    logFile << "               ";
    logFile << "| Forward Barrier  : ";
    logFile << setprecision(8);
    logFile << (QMMMOpts.ETrans- QMMMOpts.EComplforw) << " a.u ";
    logFile << "   | ";
    logFile << setprecision(4);
    logFile << (QMMMOpts.ETrans- QMMMOpts.EComplforw)*627.51 << " kcal/mol";
    logFile << '\n';
    logFile << "               ";
    logFile << "| Backward Barrier : ";
    logFile << setprecision(8);
    logFile << (QMMMOpts.ETrans- QMMMOpts.EComplback) << " a.u ";
    logFile << "   | ";
    logFile << setprecision(4);
    logFile << (QMMMOpts.ETrans- QMMMOpts.EComplback)*627.51 << " kcal/mol";
    logFile << '\n'<< endl;
    // End: Aug 28 2018
    logFile.copyfmt(call); // Replace settings
  }

  // -------------------------------------------
  // Print convergence criterias QM
  if (print_level==2)
  {
    // QMConverged
    // Print progress
    stringstream call; // Stream for system calls and reading/writing files
    call.copyfmt(logFile); // Save settings
    logFile << setprecision(12);
    logFile << "          ";
    logFile << "     | RMS dev. tol.  : " << QMMMOpts.QMOptTol << " \u212B\n";
    logFile << "          ";
    logFile << "     | RMS force tol. : ";
    logFile << LICHEMFormFloat(QMMMOpts.QMRMSForceTol,8); 
    /* logFile << "  Hartrees/\u212B\n"; */
    logFile << " Hartree/bohr\n";
    logFile << "          ";
    logFile << "     | Max force tol. : ";
    logFile << LICHEMFormFloat(QMMMOpts.QMMaxForceTol,8); 
    /* logFile << "  Hartrees/\u212B\n"; */
    logFile << " Hartree/bohr\n";
    logFile << "          ";
    logFile << "     | RMS dev.  : " << LICHEMFormFloat(RMSdiff,8);
    logFile << " \u212B" << '\n';
    logFile << "          ";
    logFile << "     | RMS force : " << LICHEMFormFloat(RMSforce,8);
    logFile << " Hartree/bohr" << '\n'; //" Hartrees/\u212B" << '\n';
    logFile << "          ";
    logFile << "     | Max. force: " << LICHEMFormFloat(MAXforce,8);
    logFile << " Hartree/bohr"; //" Hartrees/\u212B"; // << '\n'; 
    logFile.copyfmt(call); // Return to previous settings
  }
}

/*-------------------------------------------------------------------------*/

void save_files(int CalcType, int iter, fstream& logFile)
{
  stringstream call;

  if (CalcType==0)
  {
    // Save QM outputs 
    call.str("");
    call << "rm -rf ";
    call << "LICHM_QM_";
    call << iter;
    call << "; mkdir ";
    call << "LICHM_QM_";
    call << iter;
    globalSys = system(call.str().c_str());
    call.str("");
    call << "mv LICHM_*.* ";
    call << "LICHM_QM_";
    call << iter;
    call << "/ ";
    globalSys = system(call.str().c_str());
    call.str("");
    call << "cp LICHM_QM_";
    call << iter;
    call << "/LICHM_*.chk .";
    globalSys = system(call.str().c_str());
  }
  else if (CalcType==1)
  {
    // Save MM files
    call.str("");
    call << "rm -rf ";
    /*
      call << "LICHM_MM_";
      call << iter;
    */
    call << "LICHM_MM";
    call << "; mkdir ";
    /*
      call << "LICHM_MM_";
      call << iter;
    */
    call << "LICHM_MM";
    globalSys = system(call.str().c_str());
    call.str("");
    call << "mv LICHM_*.* ";
    /*
      call << "LICHM_MM_";
      call << iter; 
    */
    call << "LICHM_MM";
    call << "/ ";
    globalSys = system(call.str().c_str());
  }
  else if (CalcType==2)
  {
    // Save optimization step directories    
    call.str("");
    call << "rm -rf ";
    call << "LICHM_QSM_Opt_";
    call << iter;
    call << "; mkdir ";
    call << "LICHM_QSM_Opt_";
    call << iter;
    // mv qm directories
    call << "; mv LICHM_QM_* ";
    call << "LICHM_QSM_Opt_";
    call << iter;
    call << "/. ";
    // mv mm directories
    /* call << "; mv LICHM_MM_* "; */
    call << "; mv LICHM_MM ";
    call << "LICHM_QSM_Opt_";
    call << iter;
    call << "/.; ";
    globalSys = system(call.str().c_str());
  }
  else
  {
    logFile << "\n" << endl;
  }

}

/*-------------------------------------------------------------------------*/

void getTSbead(QMMMSettings& QMMMOpts, VectorXd& Eqmmm_images)
{
  double MaxE = 0;
  MaxE = Eqmmm_images.maxCoeff();
  double Ecomplex=0.0; // Update this value so that everyone has it
  

  for (int p=0;p<QMMMOpts.NBeads;p++)
  {
    if (Eqmmm_images[p] == MaxE)
    {
      QMMMOpts.TSBead = p;
      QMMMOpts.ETrans = Eqmmm_images[p];
    }
  }
  
  // Start: Aug 28 2018
  // Forward complex
  QMMMOpts.EComplforw=QMMMOpts.EReact;
  for (int p=0;p<QMMMOpts.TSBead;p++)
  { 
    if (Eqmmm_images[p] < QMMMOpts.EReact)
    {
      QMMMOpts.EComplforw=Eqmmm_images[p];
    }
  }
  // Backward complex
  QMMMOpts.EComplback=QMMMOpts.EProd;
  for (int p=QMMMOpts.TSBead+1;p<QMMMOpts.NBeads;p++)
  {
    if (Eqmmm_images[p] < QMMMOpts.EProd)
    {
      QMMMOpts.EComplback=Eqmmm_images[p];
    }
  }
  // End: Aug 28 2018

}
