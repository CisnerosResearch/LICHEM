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

 Compile time options for LICHEM.

*/

//Make including safe
#ifndef LICHEM_OPTS
#define LICHEM_OPTS

namespace LICHEMOpts
{
  //Compile options
  const bool JOKES = 1; //Print humorous comments

  //Monte Carlo options
  const bool isotrop = 1; //Force isotropic expansion in NPT Monte Carlo
  const double stepMin = 0.005; //Minimum Monte Carlo step size (Angstroms)
  const double stepMax = 1.0; //Maximum Monte Carlo step size (Angstroms)
  const double centRatio = 5.0; //Scales step size for path-integral centroids
  const int acc_Check = 2000; //Eq. Monte Carlo steps before checking accratio

  //Move Probabilities for PIMC
  /*

  NB: The Monte Carlo engine allows for multi-particle moves:

  If (BeadProb+CentProb) > 1) then there is a chance that multiple particles
  move during a step.

  If (BeadProb+CentProb) < 1) then there is a chance that no prarticles move
  during a step.

  If (VolProb > 0) and the simulation is not periodic, then VolProb does
  nothing.

  If the simulation is periodic, multi-particle moves will be common since
  changing the box size moves all particles.

  */
  double beadProb = 0.60; //Probability to move all beads for an atom
  double centProb = 0.50; //Probability to move a centroid
  double volProb = 0.35; //Volume change probability
};

#endif

