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

 QSM header for LICHEM.

*/

/*Start: Hatice GOKCAN*/
#ifndef LICHEM_QSM
#define LICHEM_QSM

void LICHEMQSM(vector<QMMMAtom>&,QMMMSettings&,VectorXd& wholepath,
               int,int,bool &QMDone,VectorXd& force,double,
               VectorXd& Eqmmm_images, int macroiter,fstream& logFile);


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

double runMMopt(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,fstream&);

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


void print_progress(QMMMSettings& QMMMOpts, int print_level, VectorXd& Eqmmm_images,
                    double RMSdiff, double MAXforce, double RMSforce,
                    VectorXd& reactCoord,fstream&);

double runRestrMMopt(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,double restr,fstream&);

void calc_react_coord(QMMMSettings& QMMMOpts, vector<QMMMAtom>& QMMMData,
                      VectorXd& reactCoord);

void save_files(int CalcType,int iter,fstream&);

void getTSbead(QMMMSettings& QMMMOpts,VectorXd& Eqmmm_images);

void CalcFreq(vector<QMMMAtom>& QMMMData,QMMMSettings& QMMMOpts,fstream&);

#include "QSM_tools.cpp"
#include "Optimization_tools.cpp"
#include "Utilities.cpp"

#endif
/*End: Hatice GOKCAN*/
