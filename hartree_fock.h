/*****************************************

Programmer: Alex Dyhdalo

Revision History: July 18, 2014 - created
		  July 24, 2014 - complete for s-wave only Minn

Header file that defines functions and variables used in the HF
class solver. 

The class holds general definitions for a Hartree-Fock solver
that has a 1-body term and a 2-body interaction.

*****************************************/

#ifndef HARTREE_FOCK_H // define guard
#define HARTREE_FOCK_H

#include <armadillo>
#include "phys_system.h"
#include <string>

class hartree_fock
{
  private:
   int nmax;
   int particles;
   int total_iter;
   int orbital;
   arma::mat density_mat;
   void update_density_mat (arma::mat &, phys_system &);
   double eps;
   double calc_energy (phys_system &, double);
   std::string choice;
   int sloppy (int);
  public:
   hartree_fock (int, int, std::string, int);
   void run (phys_system &, double);
};

#endif
