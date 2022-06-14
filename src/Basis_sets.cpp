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


 Compiled databases for looking up Cartesian or Hermite Gaussian basis sets

 Reference for basis sets:
 

*/

/*

  Includes:
    - vector<HermGau> HermBasis

*/

// SECTION: Basis set definitions

/*
  vector<HermGau> HermBasis
  -------------------------
  Generates Hermite basis functions for an atom in the system.

  Parameters
  ----------
  Typ: Element or atom type name.
  basName: Name of the basis set.

  Returns
  -------
  newBasis: Array of basis functions for the atom.
*/
vector<HermGau> HermBasis(string Typ, string basName)
{
  // Function to set specific Hermite basis sets
  bool badBasis = 0; // Flag to exit if no basis set is found
  vector<HermGau> newBasis;
  // NB: Organized by basis name, then element
  
  // Check for errors
  if (badBasis)
  {
    cerr << "Error: Basis set " << basName;
    cerr << " is not defined for atom " << Typ;
    cerr << "!!!" << '\n' << '\n';
    cerr.flush();
    exit(0);
  }
  // Return basis set if it was found in the database
  return newBasis;
};
