/*****************************************

Programmer: Alex Dyhdalo

Revision History: July 18, 2014 - created

Class to hold information about the phys_system to analyzed. Here we consider neutron drops in an external potential interacting via some 2-body force (here the minnesota potential).

The class here is derived from the base class phys_system and solves for the antisymmetric two-body matrix
elements for pure s (l = 0) states.  

*****************************************/

#ifndef SWAVE_H // define guard
#define SWAVE_H

#include <armadillo>
#include <cmath>
#include <string>

class swave : public phys_system
{
  private:
   double l_param;
   double v_r, v_t, v_s, k_r, k_t, k_s;
   arma::field<arma::mat> twobody_pot;
   double alex_radial (int, int, double, double);
   double interaction (double, double);
   void test_ortho (int, double [], double [], double, int);
   void populate_2body (double [], double [], int);
   double individ_element (int, int, int, int, double [], double [], int);
  public:
   swave (int, double, double[], double [], int, std::string);
   double pot_access_modifier (int, int, int, int, int);
   double kinetic_access (int, int, double);
   double orb_degen (int, int, std::string, int);
   int basis_extent ();
};

#endif
