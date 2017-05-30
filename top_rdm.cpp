/*
plan:
	create mesh
	allocate space for 1 and 2-body ME
	populate 1 and 2-body matrix elements
	output ME to dat file in SDPA format
	call CSDP program	
*/
/*****************************************************************************

Programmer: Alex Dyhdalo

This is the main wrapper file that calls all the other machinery to solve the HF
equations. 

The struct defined here holds all the user specified parameters which are read in via the 
read_in_inputs function below. User inputs are read in from the file "inputs.inp". This is 
done to avoid having to recompile when the user changes input parameters. 

Then, depending on the user choice, an instance of the corresponding derived class is 
created (swave, fullm, or fullj). These derived classes are derived from the base class phys_system.
The basis size is extracted via a class function. Then, an instance of the hartree_fock class is created
 and the hartree_fock class function ".run" is called which actually outputs the ground state energy 
for the system. 

*****************************************************************************/

// necessary definitions
#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "boost/cstdlib.hpp"
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <cstdio>
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

static const double v_r = 200.;				// parameter for the minnesota potential (the same down to k_s)
static const double v_t = -178.;
static const double v_s = -91.85;
static const double k_r = 1.487;
static const double k_t = 0.639;
static const double k_s = 0.465;

struct parameters 
{
	int nodes;
	int mesh_size;
	int s_wave_basis;
	int particles;
	double hc;
	double mass;
	double hw;
};



// function prototypes
extern void gauss (int npts, int job, double a, double b, double xpts[], double weights[]);
parameters read_in_inputs ();

template <typename two_array>
void populate_1body (two_array & h1_mat, const double hw)
{

	size_t mat_length = h1_mat.size();

  size_t spatial_orbs = mat_length/2;

	for (size_t i = 0; i < mat_length; ++i)
  {
    size_t j = i;

    if (j >= spatial_orbs)
      j = i - spatial_orbs; 

		h1_mat[i][i] = (2.*(double)j + 1.5) * hw;
  }

}


template <typename Array>
void print(std::ostream& os, const Array& A)
{
  typename Array::const_iterator i;
  os << "[";
  for (i = A.begin(); i != A.end(); ++i) {
    print(os, *i);
    if (boost::next(i) != A.end())
      os << ',';
  }
  os << "]";
  os << "\n";
}


template<> void print<double>(std::ostream& os, const double& x)
{
  os << x;
}


double ho_radial (const int n, const int l, const double b, const double r)
{
  double xi = r/b;
  double alpha = (double)l + 0.5;

  double result;
  double constant = sqrt((gsl_sf_fact(n)  * pow(2., (double)n + (double)l + 2.))/ (sqrt(M_PI) * gsl_sf_doublefact(2.*(double)n + 2.*(double)l + 1. )));

  result = pow(b, -1.5) * constant * exp(-xi*xi/2.) * pow(xi,(double)l) * gsl_sf_laguerre_n (n, alpha, (xi*xi));
  
  return result;
}

double potential (double i, double j) 
{
  double output = 0.;
  double numerical_prefactor = 1./8.;

  output = v_r *(exp(-k_r * (i-j)*(i-j)) - exp(-k_r *(i+j)*(i+j)) )/k_r 
  + v_s * (exp(-k_s * (i-j)*(i-j)) - exp(-k_s *(i+j)*(i+j)) )/k_s;

  output *= numerical_prefactor;

  return output;
}

double swave_ME (const int a, const int b, const int c, const int d, const double x[], const double w[], const int mesh_size, const double l_param)
{

  double direct = 0.;
  double exchange = 0.;


  for (int i = 0; i < mesh_size; i++)
  {
  for (int j = 0; j < mesh_size; j++)
  {

  direct += x[i] * x[j] * w[i] * w[j] * ho_radial(a, 0, l_param, x[i]) * ho_radial(c, 0, l_param, x[i]) * ho_radial(b, 0, l_param, x[j]) * ho_radial(d, 0, l_param, x[j]) * potential (x[i], x[j]);

  }
  }

  for (int i = 0; i < mesh_size; i++)
  {
  for (int j = 0; j < mesh_size; j++)
  {
  exchange += x[i] * x[j] * w[i] * w[j] * ho_radial(a, 0, l_param, x[i]) * ho_radial(d, 0, l_param, x[i]) * ho_radial(b, 0, l_param, x[j]) * ho_radial(c, 0, l_param, x[j]) * potential (x[i], x[j]);
  }
  }

  return (direct + exchange);
}


template <typename four_array>
void populate_2body (four_array & h2_mat, const double x[], const double w[], const int mesh_size, const double l_param)
{

  size_t mat_length = h2_mat.size();

  size_t s_orbs = mat_length / 2;

  std::cout << mat_length << "\n";

  for(size_t a = 0; a < s_orbs; a++)
  {  
  for(size_t b = 0; b < s_orbs; b++)
  { 
  for(size_t c = 0; c < s_orbs; c++)
  {
  for(size_t d = 0; d < s_orbs; d++)      // a-d loops over orbital shells 
  {
    for (int s1 = 0; s1 < 2; s1++)
    {
    for (int s2 = 0; s2 < 2; s2++)   // s1, s2 loops over spin degeneracy
    {

  // matrix elements are zero for aligned spins, redundant included for code clarity
      if (s1 == s2)
      {
        h2_mat [a][b][c][d] = 0.;
        h2_mat [a + s_orbs][b + s_orbs][c + s_orbs][d + s_orbs] = 0.;
      }

      if (s1 > s2)
        h2_mat [a + s_orbs][b][c + s_orbs][d] = h2_mat[a][b + s_orbs][c][d + s_orbs];

      if (s2 > s1)
        h2_mat [a][b + s_orbs][c][d + s_orbs] = swave_ME(a,b,c,d, x, w, mesh_size, l_param);
  
    }
    }   // end spin loops
  }
  }
  }
  }     // end oscillator shell loops

}

template <typename two_array, typename four_array>
  void create_spda_file (const two_array h1_mat, const four_array h2_mat)
  {
    
    size_t h1_len = h1_mat.size();

    size_t s_orbs = h1_len / 2;

    for (size_t i = 0; i < h1_len; ++i)
    {
      for (size_t j = i; j < h1_len; ++j)
      {
        if (h1_mat[i][j] != 0.)
          std::cout << h1_mat[i][j] << " ";
      }
    }  

    
  }


/********************************************

BEGIN MAIN PROGRAM

********************************************/

int main ()
{
  parameters input_params;						// create a struct to hold user specified parameters
  input_params = read_in_inputs (); 					// read in the values from inputs.inp, store in input_params	

  // Now create/copy input parameters for ease of reading the code into newly defined local variables


  const int mesh_size = input_params.mesh_size;		// how many points to use for Gauss-Legendre (G-L) quadrature 
  const int bsize = 2 * input_params.s_wave_basis;	// the size of the HO basis used, ONLY ACTIVE FOR S-WAVE
  const int particles = input_params.particles;		// the number of neutrons (particles) in the trap
  const int nodes = input_params.nodes;
  const double hc = input_params.hc;					// hbar * c - given in MeV * fm
  const double mass = input_params.mass;				// Nucleon mass - given in MeV/c^2
  const double hw = input_params.hw;					// hbar * omega 
  const double l_param = hc / sqrt(mass * hw);		// the relevant length parameter for the HO


  std::cout << "Building system..." << std::endl;

  double weight[nodes], x[nodes];							// create arrays to hold gaussian weights/positions
  gauss (mesh_size, 2, 0., 4., x, weight);						// creating weights/positions for s-wave G-L quad.
/*  for (int i = 0; i < nodes; i++)
  {
	printf ("%f", x[i]);
	printf ("\n");
  }*/

  int test = 10;

  typedef boost::multi_array<double, 2> two_array;
//  typedef two_array::index index;
  two_array h1_mat(boost::extents[bsize][bsize]);


  print(std::cout, h1_mat); 

  populate_1body (h1_mat, hw);

/*  // Assign values to the elements
  int values = 0;
  for(size_t i = 0; i != bsize; ++i) 
    for(size_t j = 0; j != bsize; ++j)
       		h1_mat[i][j] = values++;
*/
  print(std::cout, h1_mat); 

//  test_pass_2 (h1_mat);


  typedef boost::multi_array<double, 4> four_array;
//  typedef four_array::index index;
  four_array h2_mat(boost::extents[bsize][bsize][bsize][bsize]);

  populate_2body (h2_mat, x, weight, mesh_size, l_param);

  print(std::cout, h2_mat);

  create_spda_file (h1_mat, h2_mat);


//  print(std::cout, h2_mat);

/*  // Verify values
  int verify = 0;
  for(size_t i = 0; i != bsize; ++i) 
    for(size_t j = 0; j != bsize; ++j)
        	assert(h1_mat[i][j] == verify++);


*/

//  arma::mat h1_mat (bsize, bsize, arma::fill::zeros);
//  arma::field <arma::mat> h2_mat (bsize);

//  h2_mat.print();

//  h1_mat (1,1) = 4.0;

//  test_pass (h1_mat);

//  h1_mat.print();

/* swave (s_wave_basis_size, b_param, x, weight, mesh_size, choice);	// create an instance of the swave class
    basis_size = neutron_drops.basis_extent();  					// output the size of the basis
    std::cout << "Constructing Hartree Fock Instance" << std::endl;
    std::cout << "Basis Size = " << basis_size << std::endl;
    hartree_fock solution (basis_size, total_iterations, choice, particles);		// create an instance of the HF class
    std::cout << "System built: Iterating Hartree-Fock" << std::endl;
    solution.run (neutron_drops, hbw);   						// self-consistently solve the HF eqns
*/
 
  return EXIT_SUCCESS;
}

/************************************************

END MAIN PROGRAM

************************************************/

/**********************************************

Function which reads in the user specified values in
inputs.inp to an instance of the struct parameters.

**********************************************/

parameters read_in_inputs ()
{
  struct parameters input;				// local struct to hold user values
 
  std::ifstream input_in ("inputs.inp");		// the relevant input file stream
  std::string dummy;
  int counter = 0;					// counter lets us know which variable we're reading in

  while (std::getline (input_in, dummy)) 		// while runs as long as there are lines to read in the input file
  {
	if (!dummy.length() || dummy[0] == '#')		// skip zero length lines and lines that start with #
	  continue;

	std::stringstream ss;

	ss << dummy;					// read in the line to the stringstream ss

	if (counter == 0)				// read in parameters by the way they're ordered in "inputs.inp"
	  ss >> input.nodes;

	if (counter == 1)
	  ss >> input.mesh_size;
	
	if (counter == 2)
	  ss >> input.s_wave_basis;

	if (counter == 3)
	  ss >> input.particles;

	if (counter == 4)
	  ss >> input.hc;
	
	if (counter == 5)
	  ss >> input.mass;

	if (counter == 6)
	  ss >> input.hw;

	counter++;					// increase the counter after reading in a parameter	
  }

  return input;						// return the struct 
}
