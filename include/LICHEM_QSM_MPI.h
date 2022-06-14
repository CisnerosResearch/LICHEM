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


 QSM header for parallel LICHEM.

*/

// Start: Hatice GOKCAN

#ifndef LICHEM_QSM
#define LICHEM_QSM

void construct_path(VectorXd& wholepath,MatrixXd& Geom1, MatrixXd& Geom2,
                    int beadsize, int Nimages, int Natoms,
                    vector<QMMMAtom>&,QMMMSettings&);

void linearly_interpolate(VectorXd& wholepath,VectorXd& reactvec,
                          VectorXd& prodvec, int beadsize, int Nimages);

void init_Hess(MatrixXd& Hessmat, int beadsize, int Nimages);

void DBFGS(MatrixXd& Hess, VectorXd& pathdiff, VectorXd& graddiff, int ndim);

void quad_app(VectorXd& wholepath, VectorXd& oldpath, MatrixXd& Hessmat,
              VectorXd& force, VectorXd& energy, int Nimages,
              int beadsize, VectorXd& equad,VectorXd& gquad);

void updateTR (MatrixXd& Hessmat, VectorXd& glast, VectorXd& oldpath,
               VectorXd& wholepath, VectorXd& Eqmmm_images,
               VectorXd& lastenergy, VectorXd& trs,
               VectorXd& maxtr, int beadsize,int Nimages);

void ODESolve(VectorXd& path, MatrixXd& H, VectorXd& g,
              VectorXd& energy,int Nimages, int beadsize,
              VectorXd& trs, VectorXd& cons, int wholesize,
              double ftol, double &dftot,VectorXd& weight,fstream&);

void RKStep(VectorXd& path1,VectorXd& path0,MatrixXd& H,VectorXd& g,
            VectorXd& energy,int Nimages,int beadsize,double h,
            double error0,double &error1,VectorXd& newpath,
            double &hnew,bool currswitch[],VectorXd& cons,VectorXd& weight);

void funupwind(VectorXd& path, VectorXd& g, VectorXd& energy,\
int Nimages, int beadsize, VectorXd& gtan, VectorXd& weight);

void calc_tangents(VectorXd& path, int Nimages, int beadsize,
                   VectorXd& tangents, VectorXd& energy,
                    VectorXd& weight,int choice);

void alignment(VectorXd& path,int Nimages,int beadsize,VectorXd& aligned);

void eqconst(VectorXd& path, int Nimages, int beadsize, double &eqcons);

void weigh_path(VectorXd& weight, bool nebatoms[],int QMdim,
                int beadsize,int Nimages);

void spaceoutcubic(VectorXd& wholepath, bool nebatoms[], int N,
                   int Ndim, VectorXd& weight);

void getcoeff(VectorXd& x, VectorXd a, int n, int Ndim,MatrixXd& coeffmat);

void fullNCS(MatrixXd& coeffs, VectorXd& dists,double t, int ndim,
             int n, VectorXd& x);

void updatepath(VectorXd& wholepath, vector<QMMMAtom>& QMMMData,
                QMMMSettings& QMMMOpts,int beadsize, int Natoms,
                bool path_to_struct);

double CalcEnergy(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,
                  VectorXd& Eqm_images, VectorXd& Emm_images,
                  VectorXd& Eqmmm_images,fstream& logFile);

double runMMopt(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts);

double TINKEROptRestr(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                      int Bead, double restr,fstream&, int&);

bool QSMConverged(vector<QMMMAtom>& QMMMData, vector<QMMMAtom>& OldQMMMData,
                  int stepct, QMMMSettings& QMMMOpts, VectorXd& Eqmmm_images,
                  fstream& logFile);

double CalcForces(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,
                  VectorXd& Eqm_images, VectorXd& Emm_images,
                  VectorXd& Eqmmm_images,VectorXd& force,
                  int beadsize, int QMdim, bool first_time,fstream& logFile);

bool QMConverged(vector<QMMMAtom>& QMMMData,vector<QMMMAtom>& OldQMMMData,
                 MatrixXd& ForceStats, int stepct, QMMMSettings& QMMMOpts,
                 double &rmsdiff, double &rmsforce, double &maxforce);

void print_progress(QMMMSettings& QMMMOpts, int print_level,
                    VectorXd& Eqmmm_images,
                    double RMSdiff, double MAXforce, double RMSforce,
                    VectorXd& reactCoord,fstream&);

double runRestrMMopt(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,
                     double restr);

void calc_react_coord(QMMMSettings& QMMMOpts, vector<QMMMAtom>& QMMMData,
                      VectorXd& reactCoord);

void save_files(int CalcType,int iter,fstream&);

void getTSbead(QMMMSettings& QMMMOpts,VectorXd& Eqmmm_images);

void CalcFreq(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,fstream&);

// MPI_Tools
void Bcast_globals(int root);
void Bcast_settings(QMMMSettings& QMMMOpts,int root,
                    bool controller);

void Send_qmmmdata(vector<QMMMAtom>& QMMMData,
                   int Nbeads,int root,bool controller,
                   int Natoms);

void LICHEMQSMMPI(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
               VectorXd& wholepath, int Nimages, int QMdim, bool &QMDone,
               VectorXd& force, double spaceout_dist, VectorXd& Eqmmm_images,
               int macroiter,fstream&);

void GaussianForcesMPI(vector<int> mybead_list,
                       int mysize,int pathstart, int pathend);

void GaussianForcesMPIWrite(vector<QMMMAtom>& QMMMData,
                      QMMMSettings& QMMMOpts, int bead);

double GaussianForcesMPIRead(vector<QMMMAtom>& QMMMData,VectorXd& forces,
                      QMMMSettings& QMMMOpts, int bead);


void GaussianEnergyMPI(vector<int> mybead_list,
                       int mysize,int pathstart, int pathend);

void GaussianEnergyMPIWrite(vector<QMMMAtom>& QMMMData,
                      QMMMSettings& QMMMOpts, int bead);

double GaussianEnergyMPIRead(vector<QMMMAtom>& QMMMData,
                      QMMMSettings& QMMMOpts, int bead);

void TINKERForcesMPIWrite(vector<QMMMAtom>& QMMMData,
                    QMMMSettings& QMMMOpts, int bead,int&,fstream&);

void TINKERPolForcesMPIWrite(vector<QMMMAtom>& QMMMData,
                       QMMMSettings& QMMMOpts, int bead,int&,fstream&);

void TINKEREnergyMPIWrite(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                    int bead);

void TINKERPolEnergyMPIWrite(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                       int bead,int&,fstream&);

void TINKEROptMPIWrite(vector<QMMMAtom>& QMMMData,
                       QMMMSettings& QMMMOpts, int bead,int&,fstream&);

void TINKEROptRestrainMPIWrite(vector<QMMMAtom>& QMMMData,
                               QMMMSettings& QMMMOpts, int bead,
                               double restr,int&,fstream& logFile);

void TINKERForcesMPIRead(vector<QMMMAtom>& QMMMData, VectorXd& forces,
                    QMMMSettings& QMMMOpts, int bead);

void TINKERPolForcesMPIRead(vector<QMMMAtom>& QMMMData, VectorXd& forces,
                            QMMMSettings& QMMMOpts, int bead);

double TINKEREnergyMPIRead(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                    int bead);

double TINKERPolEnergyMPIRead(vector<QMMMAtom>& QMMMData,
                              QMMMSettings& QMMMOpts, int bead);

void TINKEROptMPIRead(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                    int bead);

void TINKERForcesMPI(vector<int> mybead_list,
                  int mysize,int pathstart, int pathend);

void TINKERPolForcesMPI(vector<int> mybead_list,
                  int mysize,int pathstart, int pathend);

void TINKEREnergyMPI(vector<int> mybead_list,
                     int mysize,int pathstart, int pathend);

void TINKERPolEnergyMPI(vector<int> mybead_list,
                     int mysize,int pathstart, int pathend);

void TINKEROptMPI(vector<int> mybead_list,
                  int mysize,int pathstart, int pathend,
                  QMMMSettings& QMMMOpts);

void QSMConvergedMPI(vector<QMMMAtom>& QMMMData, vector<QMMMAtom>& OldQMMMData,
                  int stepct, QMMMSettings& QMMMOpts, VectorXd& Eqmmm_images,
                  bool &PathDone,fstream&);

double CalcForcesMPI(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,
                  VectorXd& Eqm_images, VectorXd& Emm_images,
                  VectorXd& Eqmmm_images,VectorXd& force,
                  int beadsize, int QMdim, bool first_time,fstream& logFile);

double runMMoptMPI(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,
                   bool before_qsm,fstream&);

void WriteGauInputMPI(vector<QMMMAtom>& QMMMData, string calcTyp,
                   QMMMSettings& QMMMOpts, int bead,
                   int gautype);

double runRestrMMoptMPI(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                        double restr,fstream&);

#include "Reaction_path_MPI.cpp"
#include "QSM_tools.cpp"
#include "Optimization_tools.cpp"
#include "Utilities.cpp"
#include "QSM_tools_MPI.cpp"
#include "MPI_tools.cpp"
#include "GaussianMPI.cpp"
#include "TINKERMPI.cpp"
#include "Struct_writerMPI.cpp"

#include <mpi.h>

#endif
// End: Hatice GOKCAN

