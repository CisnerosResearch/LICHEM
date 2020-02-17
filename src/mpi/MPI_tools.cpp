/*
  ###############################################################################
  #                                                                             #
  # Hatice GOKCAN                                                               #
  #                                                                             #
  # Functions for data sharing in parallel LICHEM                               #
  # Includes:                                                                   #
  #                                                                             #
  #        Broadcast global variables in LICHEM_globals.h : void Bcast_globals  #
  #        Broadcast elements of class QMMMSettings       : void Bcast_settings #
  #        Send elements of class QMMMAtom,                                     #
  #                         class Coord,                                        #
  #                         class MPol,                                         #
  #                         class OctCharges              : void Send_qmmmdata  #
  #                                                                             #
  ###############################################################################
*/

void Bcast_globals(int root)
{

  /*integers*/
  MPI_Bcast(&globalSys,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&Nthreads,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&Ncpus,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&Nfreeze,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&Npseudo,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&Nbound,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&Natoms,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&Nqm,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&Nmm,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&startTime,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&endTime,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMTime,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&MMTime,1,MPI_INT,root,MPI_COMM_WORLD);

  /*doubles*/
  MPI_Bcast(&mcStep,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&Lx,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&Ly,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&Lz,1,MPI_DOUBLE,root,MPI_COMM_WORLD);

  /*booleans*/
  MPI_Bcast(&GEM,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&AMOEBA,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&CHRG,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&PSI4,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&NWChem,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&Gaussian,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&TINKER,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&LAMMPS,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&PBCon,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMM,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&MMonly,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMonly,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&OptSim,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&SteepSim,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&DFPSim,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&NEBSim,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&PIMCSim,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&FBNEBSim,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&FreqCalc,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&SinglePoint,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&GauExternal,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QSMSim,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&g09,1,MPI::BOOL,root,MPI_COMM_WORLD);
  //MPI_Bcast(&g16,1,MPI::BOOL,root,MPI_COMM_WORLD);
 
}
/*---------------------------------------------------*/
void Bcast_settings(QMMMSettings& QMMMOpts,int root,
                    bool master)
{

  /*integers*/
  MPI_Bcast(&QMMMOpts.RAM,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.charge,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.spin,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.LRECPow,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.NEq,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.NSteps,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.NBeads,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.NPrint,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.maxOptSteps,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.TSBead,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.Nqsm,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.perOpt,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.perQM,1,MPI_INT,root,MPI_COMM_WORLD);

  /*bools*/
  MPI_Bcast(&QMMMOpts.memMB,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.useLREC,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.useMMCut,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.useEwald,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.useImpSolv,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.climb,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.frznEnds,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.NEBFreq,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.printNormModes,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.startPathChk,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.restrMM,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.KeepFiles,1,MPI::BOOL,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.debug,1,MPI::BOOL,root,MPI_COMM_WORLD);

  /*doubles*/
  MPI_Bcast(&QMMMOpts.LRECCut,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.MMOptCut,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.temp,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.beta,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.press,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.accRatio,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.MMOptTol,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.QMOptTol,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.stepScale,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.maxStep,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.kSpring,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.EOld,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.EReact,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.EProd,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.ETrans,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.restrConst,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.QMRMSForceTol,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.QMMaxForceTol,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.MaxQMSteps,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  /* Start: Aug 28 2018 */
  MPI_Bcast(&QMMMOpts.EComplforw,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.EComplback,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  /* End: Aug 28 2018 */
  /*strings*/
  /* get size of strings */
  int sz1,sz2,sz3,sz4,sz5,sz6;
  if(master){
    sz1 = QMMMOpts.func.size();
    sz2 = QMMMOpts.basis.size();
    sz3 = QMMMOpts.unitsQM.size();
    sz4 = QMMMOpts.backDir.size();
    sz5 = QMMMOpts.solvModel.size();
    sz6 = QMMMOpts.ensemble.size();
  }
  /* broadcast size of strings*/
  MPI_Bcast(&sz1,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&sz2,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&sz3,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&sz4,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&sz5,1,MPI_INT,root,MPI_COMM_WORLD);
  MPI_Bcast(&sz6,1,MPI_INT,root,MPI_COMM_WORLD);

  /* resize strings in worker procs*/
  if(!master){
    QMMMOpts.func.resize(sz1);
    QMMMOpts.basis.resize(sz2);
    QMMMOpts.unitsQM.resize(sz3);
    QMMMOpts.backDir.resize(sz4);
    QMMMOpts.solvModel.resize(sz5);
    QMMMOpts.ensemble.resize(sz6);
  }
  /* broadcast strings */
  MPI_Bcast(&QMMMOpts.func[0],sz1+1,MPI_CHAR,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.basis[0],sz2+1,MPI_CHAR,root,MPI_COMM_WORLD);
  MPI_Bcast(&QMMMOpts.unitsQM[0],sz3+1,MPI_CHAR,root,MPI_COMM_WORLD) ;
  MPI_Bcast(&QMMMOpts.backDir[0],sz4+1,MPI_CHAR,root,MPI_COMM_WORLD) ;
  MPI_Bcast(&QMMMOpts.solvModel[0],sz5+1,MPI_CHAR,root,MPI_COMM_WORLD) ;
  MPI_Bcast(&QMMMOpts.ensemble[0],sz6+1,MPI_CHAR,root,MPI_COMM_WORLD);

}

/*-----------------------------------------------------*/
void Send_qmmmdata(vector<QMMMAtom>& QMMMData,
                   int Nbeads,int root,bool master,
                   int Natoms)
{

  int nprocs;
  int receiver;
  int mybead;
  int myrank;
  MPI_Status stat;
  
  if(!master){
    QMMMData.resize(Natoms);
  }
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //cout << "Gr:" << myrank << " inside send_qmmm\n" << endl;

  /* Create struct for class Mpole,coords,and octahedral charges */
  struct{
    bool chiralFlip; /*Flip y axis*/
    /*string type; */ /*Bisector, Z-then-X, Z-Only, 3-Fold, or Z-Bisect*/
    int atom1; /*Atom which defines the z axis*/
    int atom2; /*Atom which defines the x axis*/
    int atom3; /*Atom which defines the y axis (chiral only)*/
    /*Monopole moment*/
    double q;
    /*Cartesian dipole moments*/
    double Dx;
    double Dy;
    double Dz;
    /*Cartesian induced dipole moments (global frame)*/
    double IDx;
    double IDy;
    double IDz;
    /*Cartesian quadrupole moments (Q_ij = Q_ji)*/
    double Qxx;
    double Qxy;
    double Qxz;
    double Qyy;
    double Qyz;
    double Qzz;
    /* Coords */
    double cx;
    double cy;
    double cz;
    /* octahedral charges */
    double q1; /*Charge in the +x direction*/
    double q2; /*Charge in the +y direction*/
    double q3; /*Charge in the +z direction*/
    double q4; /*Charge in the -x direction*/
    double q5; /*Charge in the -y direction*/
    double q6; /*Charge in the -z direction*/
    /*Positions of the charges in the global frame*/
    double x1; /*Position of charge 1*/
    double y1; /*Position of charge 1*/
    double z1; /*Position of charge 1*/
    double x2; /*Position of charge 2*/
    double y2; /*Position of charge 2*/
    double z2; /*Position of charge 2*/
    double x3; /*Position of charge 3*/
    double y3; /*Position of charge 3*/
    double z3; /*Position of charge 3*/
    double x4; /*Position of charge 4*/
    double y4; /*Position of charge 4*/
    double z4; /*Position of charge 4*/
    double x5; /*Position of charge 5*/
    double y5; /*Position of charge 5*/
    double z5; /*Position of charge 5*/
    double x6; /*Position of charge 6*/
    double y6; /*Position of charge 6*/
    double z6; /*Position of charge 6*/
  } mystruct;
  /*last member is MPI_UB to be safe*/ 
  int blocklengths[45]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 
                        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
 
  MPI_Datatype types[45]={MPI::BOOL, 
                          MPI_INT,
                          MPI_INT, 
                          MPI_INT, 
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_DOUBLE,
                          MPI_UB};

  MPI_Aint displacements[45];
  MPI_Datatype mystruct_type;
  MPI_Aint intex,doubleex,boolex;
  
  MPI_Type_extent(MPI_INT,&intex);
  MPI_Type_extent(MPI_DOUBLE,&doubleex);
  MPI_Type_extent(MPI::BOOL,&boolex);

  displacements[0] = (MPI_Aint) 0;                     /* chiralFlip */
  displacements[1] = boolex;                           /* atom1      */
  displacements[2] = displacements[1] + intex;         /* atom2      */
  displacements[3] = displacements[2] + intex;         /* atom3      */
  displacements[4] = displacements[3] + intex;         /* q          */
  displacements[5] = displacements[4] + doubleex;      /* Dx         */
  displacements[6] = displacements[5] + doubleex;      /* Dy         */
  displacements[7] = displacements[6] + doubleex;      /* Dz         */
  displacements[8] = displacements[7] + doubleex;      /* IDx        */
  displacements[9] = displacements[8] + doubleex;      /* IDy        */
  displacements[10]= displacements[9] + doubleex;      /* IDz        */
  displacements[11]= displacements[10] + doubleex;     /* Qxx        */
  displacements[12]= displacements[11] + doubleex;     /* Qxy        */
  displacements[13]= displacements[12] + doubleex;     /* Qxz        */
  displacements[14]= displacements[13] + doubleex;     /* Qyy        */
  displacements[15]= displacements[14] + doubleex;     /* Qyz        */
  displacements[16]= displacements[15] + doubleex;     /* Qzz        */  
  displacements[17]= displacements[16] + doubleex;     /* cx         */
  displacements[18]= displacements[17] + doubleex;     /* cy         */
  displacements[19]= displacements[18] + doubleex;     /* cz         */ 
  displacements[20]= displacements[19] + doubleex;     /* q1         */
  displacements[21]= displacements[20] + doubleex;     /* q2         */
  displacements[22]= displacements[21] + doubleex;     /* q3         */
  displacements[23]= displacements[22] + doubleex;     /* q4         */
  displacements[24]= displacements[23] + doubleex;     /* q5         */
  displacements[25]= displacements[24] + doubleex;     /* q6         */
  displacements[26]= displacements[25] + doubleex;     /* x1         */
  displacements[27]= displacements[26] + doubleex;     /* y1         */
  displacements[28]= displacements[27] + doubleex;     /* z1         */
  displacements[29]= displacements[28] + doubleex;     /* x2         */
  displacements[30]= displacements[29] + doubleex;     /* y2         */
  displacements[31]= displacements[30] + doubleex;     /* z2         */
  displacements[32]= displacements[31] + doubleex;     /* x3         */
  displacements[33]= displacements[32] + doubleex;     /* y3         */
  displacements[34]= displacements[33] + doubleex;     /* z3         */
  displacements[35]= displacements[34] + doubleex;     /* x4         */
  displacements[36]= displacements[35] + doubleex;     /* y4         */
  displacements[37]= displacements[36] + doubleex;     /* z4         */
  displacements[38]= displacements[37] + doubleex;     /* x5         */
  displacements[39]= displacements[38] + doubleex;     /* y5         */
  displacements[40]= displacements[39] + doubleex;     /* z5         */
  displacements[41]= displacements[40] + doubleex;     /* x6         */
  displacements[42]= displacements[41] + doubleex;     /* y6         */
  displacements[43]= displacements[42] + doubleex;     /* z6         */
  displacements[44]= displacements[43] + doubleex;     /* MPI_UB     */

  MPI_Type_struct(44,blocklengths,displacements,types,&mystruct_type);
  MPI_Type_commit(&mystruct_type);

  /*first and second beads are on root */
  for(int j=2; j<Nbeads;j++)
  {
     /* Define receiver rank */
     /* make sure the order is correct */
     /* Frozen ends!!!
 *     
 *      i.e. : Nbeads = 9, Nimages=7
 *             for load balance 
 *             nproc should be 7
 *
 *      beads: 0 1 2 3 4 5 6 7 8
 *      procs: 0 0 1 2 3 4 5 6 6
 *
 *      i.e. : Nbeads = 8, Nimages=6
 *             for load balance 
 *             nproc should be 2, 3, 6
 *
 *      nproc: 2
 *      beads: 0 1 2 3 4 5 6 7
 *      procs: 0 0 1 0 1 0 1 1
 *
 *      nproc: 3 (not very good load balance)
 *      beads: 0 1 2 3 4 5 6 7
 *      procs: 0 0 1 2 0 1 2 2
 *
 *      nproc: 6       
 *      beads: 0 1 2 3 4 5 6 7
 *      procs: 0 0 1 2 3 4 5 5 
 * 
 *
 *   */

     receiver = (j%nprocs)-1;/*check for !FrozenEnds*/
     /*send product to last proc*/
     if(j==Nbeads-1)
     {
       receiver = nprocs - 1;/* last proc */
     }
     if(receiver==-1){
        receiver = nprocs - 1; /* last proc */
     }
     //if(master){
     //   cout << "receiver: " << receiver << "\n" << endl;
     //}
     /* Start atoms for jth bead */
     for(int i=0;i<Natoms;i++)
     {

        /* size of the strings*/
        int sz1,sz2,sz3;
      
        int bsz;//bonds size
        /* Mpole,coord and octahedral charge size */
        int mpsize, csize, ocsize;
        /* sizes */
        if(master){
           sz1 = QMMMData[i].QMTyp.size();
           sz2 = QMMMData[i].MMTyp.size();
           sz3 = QMMMData[i].MP[j].type.size();
           bsz = QMMMData[i].bonds.size();
        }

        /* Root already owns info, no need sending */
        if(receiver!=root)/*start sending*/
        {
          /* if root proc, start sending info */
          if(master){

            /*doubles*/
            MPI_Send(&QMMMData[i].Ep,1,MPI_DOUBLE,receiver,101,MPI_COMM_WORLD);
            MPI_Send(&QMMMData[i].m,1,MPI_DOUBLE,receiver,102,MPI_COMM_WORLD);
            //if(i==1){
            //  cout << "Rec:" << receiver << " b:" << j << " m:"<< QMMMData[i].m <<" ...\n" << endl;
            //}

            /*bools*/
            MPI_Send(&QMMMData[i].NEBActive,1,MPI::BOOL,receiver,103,MPI_COMM_WORLD);
            MPI_Send(&QMMMData[i].QMRegion,1,MPI::BOOL,receiver,104,MPI_COMM_WORLD);
            MPI_Send(&QMMMData[i].MMRegion,1,MPI::BOOL,receiver,105,MPI_COMM_WORLD);
            MPI_Send(&QMMMData[i].PBRegion,1,MPI::BOOL,receiver,106,MPI_COMM_WORLD);
            MPI_Send(&QMMMData[i].BARegion,1,MPI::BOOL,receiver,107,MPI_COMM_WORLD);
            MPI_Send(&QMMMData[i].frozen,1,MPI::BOOL,receiver,108,MPI_COMM_WORLD);

            /*integers*/
            MPI_Send(&QMMMData[i].numTyp,1,MPI_INT,receiver,109,MPI_COMM_WORLD);
            MPI_Send(&QMMMData[i].numClass,1,MPI_INT,receiver,110,MPI_COMM_WORLD);
            MPI_Send(&QMMMData[i].id,1,MPI_INT,receiver,111,MPI_COMM_WORLD);

            /*strings */
          
            /* send the size of strings  */
            MPI_Send(&sz1,1,MPI_INT,receiver,112,MPI_COMM_WORLD);
            MPI_Send(&sz2,1,MPI_INT,receiver,113,MPI_COMM_WORLD);

            /* now send the strings */
            MPI_Send(&QMMMData[i].QMTyp[0],sz1+1,MPI_CHAR,receiver,114,MPI_COMM_WORLD);
            MPI_Send(&QMMMData[i].MMTyp[0],sz2+1,MPI_CHAR,receiver,115,MPI_COMM_WORLD);

            /* send size of bonds */
            MPI_Send(&bsz,1,MPI_INT,receiver,116,MPI_COMM_WORLD);
            /* send bonds */

            /*bond info*/
            MPI_Send(&QMMMData[i].bonds[0],bsz,MPI_INT,receiver,121,MPI_COMM_WORLD);
    
            /* send Mpole,coord and octahedral charge size */
            mpsize = QMMMData[i].MP.size();
            csize = QMMMData[i].P.size();
            ocsize = QMMMData[i].PC.size();
            MPI_Send(&mpsize,1,MPI_INT,receiver,122,MPI_COMM_WORLD);
            MPI_Send(&csize,1,MPI_INT,receiver,123,MPI_COMM_WORLD);
            MPI_Send(&ocsize,1,MPI_INT,receiver,124,MPI_COMM_WORLD);

            /* send Mpole class as struct */
            /* set var */
            mystruct.chiralFlip= QMMMData[i].MP[j].chiralFlip; 
            mystruct.atom1 = QMMMData[i].MP[j].atom1; 
            mystruct.atom2 = QMMMData[i].MP[j].atom2; 
            mystruct.atom3 = QMMMData[i].MP[j].atom3; 
            mystruct.q = QMMMData[i].MP[j].q;
            mystruct.Dx = QMMMData[i].MP[j].Dx;
            mystruct.Dy = QMMMData[i].MP[j].Dy;
            mystruct.Dz = QMMMData[i].MP[j].Dz;
            mystruct.IDx = QMMMData[i].MP[j].IDx;
            mystruct.IDy = QMMMData[i].MP[j].IDy;
            mystruct.IDz = QMMMData[i].MP[j].IDz;
            mystruct.Qxx = QMMMData[i].MP[j].Qxx;
            mystruct.Qxy = QMMMData[i].MP[j].Qxy;
            mystruct.Qxz = QMMMData[i].MP[j].Qxz;
            mystruct.Qyy = QMMMData[i].MP[j].Qyy;
            mystruct.Qyz = QMMMData[i].MP[j].Qyz;
            mystruct.Qzz = QMMMData[i].MP[j].Qzz;
            /* Coords */
            mystruct.cx = QMMMData[i].P[j].x;
            mystruct.cy = QMMMData[i].P[j].y; 
            mystruct.cz = QMMMData[i].P[j].z;   
            /* Octahedral Charges */
            mystruct.q1 = QMMMData[i].PC[j].q1;
            mystruct.q2 = QMMMData[i].PC[j].q2;
            mystruct.q3 = QMMMData[i].PC[j].q3;
            mystruct.q4 = QMMMData[i].PC[j].q4;
            mystruct.q5 = QMMMData[i].PC[j].q5;
            mystruct.q6 = QMMMData[i].PC[j].q6;
            mystruct.x1 = QMMMData[i].PC[j].x1;
            mystruct.y1 = QMMMData[i].PC[j].y1;
            mystruct.z1 = QMMMData[i].PC[j].z1;
            mystruct.x2 = QMMMData[i].PC[j].x2;
            mystruct.y2 = QMMMData[i].PC[j].y2;
            mystruct.z2 = QMMMData[i].PC[j].z2;
            mystruct.x3 = QMMMData[i].PC[j].x3;
            mystruct.y3 = QMMMData[i].PC[j].y3;
            mystruct.z3 = QMMMData[i].PC[j].z3;
            mystruct.x4 = QMMMData[i].PC[j].x4;
            mystruct.y4 = QMMMData[i].PC[j].y4;
            mystruct.z4 = QMMMData[i].PC[j].z4;
            mystruct.x5 = QMMMData[i].PC[j].x5;
            mystruct.y5 = QMMMData[i].PC[j].y5;
            mystruct.z5 = QMMMData[i].PC[j].z5;
            mystruct.x6 = QMMMData[i].PC[j].x6;
            mystruct.y6 = QMMMData[i].PC[j].y6;
            mystruct.z6 = QMMMData[i].PC[j].z6;     
            MPI_Send(&mystruct,1,mystruct_type,receiver,1001,MPI_COMM_WORLD);
            
            MPI_Send(&sz3,1,MPI_INT,receiver,1132,MPI_COMM_WORLD);
            if(sz3!=0){/* in the test case sz3 was 0 so it was giving error*/
              MPI_Send(&QMMMData[i].MP[j].type,sz3+1,MPI_CHAR,receiver,1152,MPI_COMM_WORLD);
            }

          }/* endif master */

          /* if worker proc, start receiving info */
          if((myrank==receiver)){

            /*doubles*/
            MPI_Recv(&QMMMData[i].Ep,1,MPI_DOUBLE,root,101,MPI_COMM_WORLD,&stat);
            MPI_Recv(&QMMMData[i].m,1,MPI_DOUBLE,root,102,MPI_COMM_WORLD,&stat);

            //if(i==1){
            //  cout << "Mr:" << myrank << " b:" << j << " m:"<< QMMMData[i].m <<" ...\n" << endl;
            //}

            /*bools*/
            MPI_Recv(&QMMMData[i].NEBActive,1,MPI::BOOL,root,103,MPI_COMM_WORLD,&stat);
            MPI_Recv(&QMMMData[i].QMRegion,1,MPI::BOOL,root,104,MPI_COMM_WORLD,&stat);
            MPI_Recv(&QMMMData[i].MMRegion,1,MPI::BOOL,root,105,MPI_COMM_WORLD,&stat);
            MPI_Recv(&QMMMData[i].PBRegion,1,MPI::BOOL,root,106,MPI_COMM_WORLD,&stat);
            MPI_Recv(&QMMMData[i].BARegion,1,MPI::BOOL,root,107,MPI_COMM_WORLD,&stat);
            MPI_Recv(&QMMMData[i].frozen,1,MPI::BOOL,root,108,MPI_COMM_WORLD,&stat);

            /*integers*/
            MPI_Recv(&QMMMData[i].numTyp,1,MPI_INT,root,109,MPI_COMM_WORLD,&stat);
            MPI_Recv(&QMMMData[i].numClass,1,MPI_INT,root,110,MPI_COMM_WORLD,&stat);
            MPI_Recv(&QMMMData[i].id,1,MPI_INT,root,111,MPI_COMM_WORLD,&stat);

            /*strings*/
            /*receive size of strings */
            MPI_Recv(&sz1,1,MPI_INT,root,112,MPI_COMM_WORLD,&stat);
            MPI_Recv(&sz2,1,MPI_INT,root,113,MPI_COMM_WORLD,&stat);

            /*resize strings*/
            QMMMData[i].QMTyp.resize(sz1);
            QMMMData[i].MMTyp.resize(sz2);

            /*receive strings*/
            MPI_Recv(&QMMMData[i].QMTyp[0],sz1+1,MPI_CHAR,root,114,MPI_COMM_WORLD,&stat);
            MPI_Recv(&QMMMData[i].MMTyp[0],sz2+1,MPI_CHAR,root,115,MPI_COMM_WORLD,&stat);

            /* vectors */
            /* receive size */
            MPI_Recv(&bsz,1,MPI_INT,root,116,MPI_COMM_WORLD,&stat);
            /* resize */
            QMMMData[i].bonds.resize(bsz);            
            /* receive vector */
            MPI_Recv(&QMMMData[i].bonds[0],bsz,MPI_INT,root,121,MPI_COMM_WORLD,&stat);

            /* receive mpole,coord and octahedral charge sizes*/
            MPI_Recv(&mpsize,1,MPI_INT,root,122,MPI_COMM_WORLD,&stat);
            MPI_Recv(&csize,1,MPI_INT,root,123,MPI_COMM_WORLD,&stat);
            MPI_Recv(&ocsize,1,MPI_INT,root,124,MPI_COMM_WORLD,&stat);
            /* resize */
            QMMMData[i].MP.resize(mpsize);
            QMMMData[i].P.resize(csize);
            QMMMData[i].PC.resize(ocsize);

            //cout << "Proc-" << myrank << " start receiving mystruct\n" << endl;   
 
            /* Receive mpole class */
            MPI_Recv(&mystruct,1,mystruct_type,root,1001,MPI_COMM_WORLD,&stat);
            /* set class inst. */
            QMMMData[i].MP[j].chiralFlip = mystruct.chiralFlip;   
            QMMMData[i].MP[j].atom1 = mystruct.atom1;   
            QMMMData[i].MP[j].atom2 = mystruct.atom2;   
            QMMMData[i].MP[j].atom3 = mystruct.atom3;   
            QMMMData[i].MP[j].q = mystruct.q;
	    QMMMData[i].MP[j].Dx = mystruct.Dx;
            QMMMData[i].MP[j].Dy = mystruct.Dy;
            QMMMData[i].MP[j].Dz = mystruct.Dz;
            QMMMData[i].MP[j].IDx = mystruct.IDx;
            QMMMData[i].MP[j].IDy = mystruct.IDy;
            QMMMData[i].MP[j].IDz = mystruct.IDz;
            QMMMData[i].MP[j].Qxx = mystruct.Qxx;
            QMMMData[i].MP[j].Qxy = mystruct.Qxy;
            QMMMData[i].MP[j].Qxz = mystruct.Qxz;
            QMMMData[i].MP[j].Qyy = mystruct.Qyy;
            QMMMData[i].MP[j].Qyz = mystruct.Qyz;
            QMMMData[i].MP[j].Qzz = mystruct.Qzz;
            /* coords */
            QMMMData[i].P[j].x = mystruct.cx ;
            QMMMData[i].P[j].y = mystruct.cy ;
            QMMMData[i].P[j].z = mystruct.cz ;
            /* Octahedral charges */ 
            QMMMData[i].PC[j].q1 = mystruct.q1;
            QMMMData[i].PC[j].q2 = mystruct.q2;
            QMMMData[i].PC[j].q3 = mystruct.q3;
            QMMMData[i].PC[j].q4 = mystruct.q4;
            QMMMData[i].PC[j].q5 = mystruct.q5;
            QMMMData[i].PC[j].q6 = mystruct.q6;
            QMMMData[i].PC[j].x1 = mystruct.x1;
            QMMMData[i].PC[j].y1 = mystruct.y1;
            QMMMData[i].PC[j].z1 = mystruct.z1;
            QMMMData[i].PC[j].x2 = mystruct.x2;
            QMMMData[i].PC[j].y2 = mystruct.y2;
            QMMMData[i].PC[j].z2 = mystruct.z2;
            QMMMData[i].PC[j].x3 = mystruct.x3;
            QMMMData[i].PC[j].y3 = mystruct.y3;
            QMMMData[i].PC[j].z3 = mystruct.z3;
            QMMMData[i].PC[j].x4 = mystruct.x4;
            QMMMData[i].PC[j].y4 = mystruct.y4;
            QMMMData[i].PC[j].z4 = mystruct.z4;
            QMMMData[i].PC[j].x5 = mystruct.x5;
            QMMMData[i].PC[j].y5 = mystruct.y5;
            QMMMData[i].PC[j].z5 = mystruct.z5;
            QMMMData[i].PC[j].x6 = mystruct.x6;
            QMMMData[i].PC[j].y6 = mystruct.y6;
            QMMMData[i].PC[j].z6 = mystruct.z6; 

            //cout << "Proc-" << myrank << " finish receiving mystruct\n" << endl;

            MPI_Recv(&sz3,1,MPI_INT,root,1132,MPI_COMM_WORLD,&stat);
            QMMMData[i].MP[j].type.resize(sz3);
            if(sz3!=0){/* in the test case sz3 was 0 so it was giving error*/
              MPI_Recv(&QMMMData[i].MP[j].type,sz3+1,MPI_CHAR,root,1152,MPI_COMM_WORLD,&stat);
            }
          
          }/* end if myrank==receiver */

        }/*end if receiver!=root*/

     }/*end for Natoms*/
      
  } /*end for Nbeads*/

  //cout << "Proc-" << myrank << " finish receiving\n" << endl;

}
/*-----------------------------------------------------*/
