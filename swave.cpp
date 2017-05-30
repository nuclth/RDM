/*********************************************************

File to define functions declared in the "swave.h" and the
"phys_system.h" files. 

*********************************************************/

#include "phys_system.h"
#include "swave.h"
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <omp.h>
#include <iomanip>
#include <string>

static const double v_r = 200.;				// parameter for the minnesota potential (the same down to k_s)
static const double v_t = -178.;
static const double v_s = -91.85;
static const double k_r = 1.487;
static const double k_t = 0.639;
static const double k_s = 0.465;




/*******************************************

Constructor for the system. Various parameters of the minnesota
potential given below. System takes nmax oscillator shells (with a
twofold degeneracy from spin --> 5 = nmax ---> 10 total states). The
final nmax multiplication by 2 in the basis_extent function reflects 
this fact.  

*******************************************/


create swave (int osc_n_input, double length_param, double x[], double w[], int mesh_size, std::string choice) : phys_system (osc_n_input, length_param, x, w, mesh_size, choice)
{
  l_param = length_param;		// length parameter for the HO
  v_r = 200.;				// parameter for the minnesota potential (the same down to k_s)
  v_t = -178.;
  v_s = -91.85;
  k_r = 1.487;
  k_t = 0.639;
  k_s = 0.465;

  nmax = osc_n_input;			// basis size input

  // test for orthonormality, if no problem occurs --> calculate two-body matrix elements
  // otherwise throw exception
  try 
  {
	test_ortho (nmax, x, w, l_param, mesh_size); 
	populate_2body (x, w, mesh_size);
  }

  catch (const char * msg)  
  {
	std::cerr << msg << std::endl;
  }

}

/***************************************************

Function to access the relevant parts of the twobody potential. 
Called in the Hartree Fock solver. 

***************************************************/

double swave::pot_access_modifier (int a, int b, int c, int d, int)
{
  double potential_factor = twobody_pot (a,b)(c,d);
  return potential_factor;
}

/***********************************************

Function to output the kinetic energy and HO 1-body potential
energy contributions. E = (2n + l + 3/2) * hbar * omega where
n is the principle quantum number and l is the orbital angular 
momentum. This term is diagonal in the principle quantum number n
e.g. <n1|H_0|n2> = E*\delta_{n1,n2}.

***********************************************/

double swave::kinetic_access (int a, int b, double hw)
{   

  if (a >= nmax && b >= nmax)			// a < nmax corresponds to spin in one direction (e.g. up)
  {						// a > nmax corresponds to spin in the other direction (e.g. down)
	a = a - nmax;
	b = b - nmax;
  }

  if (a == b)
	return ((2.*(double)a  + 1.5) * hw);

  else
	return 0.;

}

/*************************************************************

Function that computes radial wavefunction values for the harmonic
oscillator. Input is a given principle quantum number n, an orbital angular
momentum term l, a length parameter b (b = sqrt (hbar/mass*omega)), and 
the spatial coordinate point r. 

**************************************************************/

double ho_radial (int n, int l, double b, double r)
{
  double xi = r/b;
  double alpha = (double)l + 0.5;

  double result;
  double constant = sqrt((gsl_sf_fact(n)  * pow(2., (double)n + (double)l + 2.))/ (sqrt(M_PI) * gsl_sf_doublefact(2.*(double)n + 2.*(double)l + 1. )));

  result = pow(b, -1.5) * constant * exp(-xi*xi/2.) * pow(xi,(double)l) * gsl_sf_laguerre_n (n, alpha, (xi*xi));
  
  return result;
}

/**************************************************************************

Function to populate the twobody potential for the case of s-wave only. The matrix
elements are calculated by doing a double integral over the radial wavefunctions using
Gauss-Legendre quadrature.

Two-body matrix elements are stored in a armadillo field object twobody_pot (a,b) 
that holds a matrix at each separate point a,b. For example, for the values 
twobody_pot (0, 1), there exists a 2*nmax by 2*nmax matrix. This is a way for 
storing a 4-index object. The factor of 2 takes into account spin degeneracy (up,down)
for each oscillator shell. 

All terms of a given spin orientation (e.g. up) are stored in elements 0 to (nmax - 1). 
The terms of opposite spin orientation are stored in elements nmax to (2*nmax - 1).

**************************************************************************/

void populate_2body (double x[], double w[], int mesh_size)
{

  // instantiate the field object
  
  twobody_pot.set_size (2*nmax, 2*nmax);

  for (int a = 0; a < 2 * nmax; a++)
  {
	for (int b = 0; b < 2 * nmax; b++)
	{
		twobody_pot(a,b).set_size(2*nmax, 2*nmax);
	}
  }

  std::cout << "Initializing Field Elements" << std::endl;

  // set all values to zero explicitly

  for(int a = 0; a < 2*nmax; a++)
  {  
    for(int b = 0; b < 2*nmax; b++)
    {	
	for(int c = 0; c < 2*nmax; c++)
	{
	  for(int d = 0; d < 2*nmax; d++)
	  {
		   twobody_pot (a,b)(c,d) = 0.;
          }

	}
    }
  }


  std::cout << "Computing Matrix Elements" << std::endl;

  omp_set_num_threads(4);                                // NUMBER OF CORES SET HERE

  #pragma omp parallel for

  for(int a = 0; a < nmax; a++)
  {  
  for(int b = 0; b < nmax; b++)
  {	
  for(int c = 0; c < nmax; c++)
  {
  for(int d = 0; d < nmax; d++)      // a-d loops over orbital shells 
  {
    for (int s1 = 0; s1 < 2; s1++)
    {
    for (int s2 = 0; s2 < 2; s2++)   // s1, s2 loops over spin degeneracy
    {

	// matrix elements are zero for aligned spins, redundant included for code clarity
	if (s1 == s2)
	{
	  twobody_pot (a,b)(c,d) = 0.;
	  twobody_pot (a + nmax, b + nmax) (c+nmax, d+nmax) = 0.;
	}

  	if (s1 > s2)
	  twobody_pot (a + nmax,b)(c + nmax,d) = twobody_pot(a, b+nmax)(c, d + nmax);

	if (s2 > s1)
	  twobody_pot (a,b + nmax)(c,d + nmax) = individ_element(a,b,c,d, x, w, mesh_size);
	
    }
    }   // end spin loops
  }
  }
  }
  }     // end oscillator shell loops

}

/***********************************************************

Function that calculates individual matrix elements for the twobody potential. The direct and exchange
terms are calculated separately then added together. The two are different by differing in the placement of
the principal quantum number n and which coordinate we are doing integration over.

Direct: <n1n2|v|n3n4> --> integration with n1, n3 over r1 and n2, n4 over r2

Exchange: <n1n2|v|n4n3> --> integration over n1, n4 over r1 and n2, n3 over r2

***********************************************************/

double swave_ME (int a, int b, int c, int d, double x[], double w[], int mesh_size)
{

  double direct = 0.;
  double exchange = 0.;


  for (int i = 0; i < mesh_size; i++)
  {
  for (int j = 0; j < mesh_size; j++)
  {

	direct += x[i] * x[j] * w[i] * w[j] * alex_radial(a, 0, l_param, x[i]) * alex_radial(c, 0, l_param, x[i]) * alex_radial(b, 0, l_param, x[j]) * alex_radial(d, 0, l_param, x[j]) * interaction (x[i], x[j]);

  }
  }

  for (int i = 0; i < mesh_size; i++)
  {
  for (int j = 0; j < mesh_size; j++)
  {
	exchange += x[i] * x[j] * w[i] * w[j] * alex_radial(a, 0, l_param, x[i]) * alex_radial(d, 0, l_param, x[i]) * alex_radial(b, 0, l_param, x[j]) * alex_radial(c, 0, l_param, x[j]) * interaction (x[i], x[j]);
  }
  }

  return (direct + exchange);
}

/**********************************************************

Part of the double integral that contains the interaction. In the present case,
this is a Minnesota potential with gaussians and prefactors v_s & v_r. 

**********************************************************/

double potential (double i, double j) 
{
  double output = 0.;
  double numerical_prefactor = 1./8.;

  output = v_r *(exp(-k_r * (i-j)*(i-j)) - exp(-k_r *(i+j)*(i+j)) )/k_r 
  + v_s * (exp(-k_s * (i-j)*(i-j)) - exp(-k_s *(i+j)*(i+j)) )/k_s;

  output *= numerical_prefactor;

  return output;
}

/************************************************************

Function to test the orthonormality of the radial wavefunctions for the HO. 

*************************************************************/

void swave::test_ortho (int total_states, double x[], double w[], double b_param, int mesh_size)
{

  double ortho_tolerance = 1.0E-10;   // numerical cutoff for orthonormality

  double ortho_num;

  // test each radial function combination

  for (int i = 0; i < total_states; i++)
  {
	for (int j = 0; j < total_states; j++)
	{

	  ortho_num = 0.;	  

	  for (int l = 0; l < mesh_size; l++)    // sum over all mesh points
	  {
		ortho_num += x[l] * x[l] * w[l] * alex_radial (i, 0, b_param, x[l]) * alex_radial (j, 0, b_param, x[l]);
	  }

	// Error for same radial wavefunctions not being equal to 1
	if (i == j && ortho_num < 0.9)
		throw "ERROR: Radial functions not orthonormal! R_i R_i != 1";

	// Error for different radial wavefunctions not being equal to 0
	if (i != j && ortho_num > ortho_tolerance)
		throw "ERROR: Radial functions not orthonormal! R_i R_j != 0";

	}
  }


}


/*********************************************************

FUNCTION RETURNS ONLY NOT 1 FOR J-SCHEME

Function to account for the degeneracy factor in the density matrix. The
swave only hamiltonian already takes spin into account in the basis size 
(e.g. for one orbital n = 0, l = 0, the Hamiltonian will be 2x2 for spin
up and spin down).  

NOTE: THIS FUNCTION ONLY WORKS CORRECTLY FOR CLOSED-SHELL SYSTEMS

**********************************************************/

double swave::orb_degen(int, int, std::string, int)
{
  
  return 1.;

}


/**********************************************************

Function to return the size of our basis for the HF solver.
Factor of 2 accounts for spin degeneracy. 

**********************************************************/

int swave::basis_extent ()
{
  return (2*nmax);
}
