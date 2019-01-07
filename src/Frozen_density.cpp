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


 LICHEM wrapper functions for frozen density QMMM calculations using the
 Gaussian electrostatic model (GEM) or similar potentials

 Reference for GEM:
 

 Reference for conversion to point-charges:
 Stone, The Theory of Intermolecular Forces, (2013)
 Devereux et al., J. Chem. Theory Comp., 10, 10, 4229, (2014)
 Gao et al., Chem. Phys. Lett., 593, 165, (2014)

*/

//GEM utility functions
double GEMC6(double C6, Coord& POSi, Coord& POSj, double Rcut)
{
  //Function to calculate the LJ style dispersion
  double Eij = 0; //Dispersion energy
  double Rij6 = 0; //Distance between atoms
  //Calculate distance
  Rij6 = CoordDist2(POSi,POSj).vecMag(); //Rij^2
  //Check cutoff
  if (Rij6 <= (Rcut*Rcut))
  {
    //Update energy
    Rij6 *= Rij6*Rij6; //Rij^6
    Eij -= C6/Rij6; //Dispersion energy
  }
  //Return energy
  return Eij;
};

double GEMBuffC7(double C7, double Rmin, Coord& POSi, Coord& POSj,
                 double Rcut)
{
  //Function to calculate buffered 14-7 style dispersion
  const double A7 = -2*pow(1.07,7); //Part of C7 from the 14-7 shift
  double Eij = 0; //Dispersion energy
  double Rij = 0; //Distance between atoms
  //Calculate distance
  Rij = CoordDist2(POSi,POSj).vecMag(); //Rij^2
  //Check cutoff
  if (Rij <= (Rcut*Rcut))
  {
    //Update energy
    double Rij7;
    Rij = sqrt(Rij); //Rij
    Rij /= Rmin; //Scale distance
    Rij += 0.07; //Add shift factor
    Rij7 = Rij*Rij*Rij; //Shifted Rij^3
    Rij7 *= Rij7*Rij; //Shifted Rij^7
    Eij = (A7*C7)/Rij7; //Dispersion energy
  }
  //Return energy
  return Eij;
};

//Functions to calculate GEM energy

