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
  std::string m_ref;
  std::string m_mat;
};



// function prototypes
extern void gauss (int npts, int job, double a, double b, double xpts[], double weights[]);
parameters read_in_inputs ();

template <typename two_array>
void populate_1body (const two_array & ref_m, two_array & h1_mat, const double hw)
{

	size_t mat_length = h1_mat.size();


	for (size_t i = 0; i < mat_length; ++i)
  {
    double n = ref_m [i][1];
    double l = ref_m [i][2];

		h1_mat[i][i] = (2.*n + l + 1.5) * hw;
  }

}


template <typename Array>
void print(std::ostream& os, const Array & A)
{
  typename Array::const_iterator i;
  os.precision(2);
  os << "[";
  for (i = A.begin(); i != A.end(); ++i) {
    print(os, *i);
//    if (boost::next(i) != A.end())
//      os << ',' << "\t";
  }
  os << "]";
  os << "\n";
}


template<> void print<double>(std::ostream& os, const double & x)
{
  os << x << "\t";
}

template<> void print<int>(std::ostream& os, const int & x)
{
  os << x;
}

/*double ho_radial (const int n, const int l, const double b, const double r)
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
}*/

/*double swave_ME (const int a, const int b, const int c, const int d, const double x[], const double w[], const int mesh_size, const double l_param)
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
}*/


/*template <typename four_array>
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

}*/

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




template <typename ref_array>
void read_in_reference_m_scheme (ref_array & ref_m, std::string m_ref_file)
{
  const char * m_reference_file = (m_ref_file).c_str();
  // input file stream for m_scheme
  std::ifstream m_ref_in (m_reference_file);
 
  if (m_ref_in.fail())
  throw "ERROR: Cannot open m-scheme reference file";
 
  size_t total_lines = 0;
  std::string dummy;
  double ref_num, n, l, j , m_j, tz;

  // find total number of defined reference lines
  while (std::getline (m_ref_in, dummy))
  ++total_lines;


  // clear the file stream, reset to read in the elements
  m_ref_in.clear();
  m_ref_in.seekg (0, std::ios::beg);

  size_t ref_size = ref_m.size();
  size_t ele_in = 0;

  // read in and assign the references line by line
  // to the matrix ref_m
  for (size_t i = 0; i < total_lines; i++)
  {
  std::string orbit_dummy_1;
  std::string orbit_dummy_2;
  std::getline (m_ref_in, dummy);

  if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
    continue;

  std::stringstream ss;

  ss << dummy;            // read in the line to stringstream ss

  ss >> orbit_dummy_1 >> orbit_dummy_2 >> ref_num >> n >> l >> j >> m_j >> tz;  // assign values of the line

  std::cout << ref_num << " " << n << " " << l << " " << j << " " << m_j << " " << tz << "\n";


  // only extract neutron-neutron states
  if (tz == 1.)
  {
    ref_m [ele_in][0] = ref_num;    // reference number of the line
    ref_m [ele_in][1] = n;          // principle quantum number
    ref_m [ele_in][2] = l;          // orbital angular mom.
    ref_m [ele_in][3] = j * 0.5;    // total angular mom.
    ref_m [ele_in][4] = m_j * 0.5;  // total angular mom. projection
    ref_m [ele_in][5] = tz;         // isospin projection (should all be +1.0)
    ele_in++;
  }

  if (ele_in >= ref_size)
    break;

  }
  
  
  

  print(std::cout, ref_m);            // print the resulting matrix


  
}

/**********************************************

Function to read in the actual m scheme matrix  elements 
calculated from Morten's code. 

THIS WILL BREAK IF THE INPUT FILE .dat CHANGES TITLE OR FORMAT.

**********************************************/

template <typename ref_array, typename four_array>
void read_in_matrix_m_scheme (ref_array & ref_m, four_array & h2_mat, std::string m_mat_file)
{
  // input file stream
  const char * m_matrix_file = (m_mat_file).c_str();
  std::ifstream m_matrix_in (m_matrix_file);

  size_t total_lines = 0;
  std::string dummy;
  int alpha, beta, gamma, delta;
  double value;

  if (m_matrix_in.fail())
  throw "ERROR: Cannot read m-scheme matrix file";

  // get the total number of lines in the matrix elements file
  while(std::getline(m_matrix_in, dummy))
  ++total_lines;
 
  // clear and reset the file stream
  m_matrix_in.clear();
  m_matrix_in.seekg(0, std::ios::beg);

  size_t m_size = h2_mat.size();

  // set the size of the m_scheme holder
  //m_scheme_matrix.set_size (m_size, m_size);


  int i,j,k,l;

  for (size_t n = 0; n < total_lines; n++)
  {
  std::getline (m_matrix_in, dummy);

  if (!dummy.length() || dummy[0] == '#')
    continue;

  i = -1;
  j = -1;
  k = -1;
  l = -1;

  std::stringstream ss;

  ss << dummy;

  ss >> alpha >> beta >> gamma >> delta >> value;

  // check to see if alpha, beta, gamma, delta we're reading in
  // matches the reference numbers in our reference matrix
  // if yes --> read it in, otherwise discard it

  for (size_t w = 0; w < m_size; w++)
  {
    if (alpha == ref_m[w][0]) 
      i = w;
  }
  
  for (size_t w = 0; w < m_size; w++)
  {
    if (beta == ref_m[w][0]) 
      j = w;
  }

  for (size_t w = 0; w < m_size; w++)
  {
    if (gamma == ref_m[w][0]) 
      k = w;
  }

  for (size_t w = 0; w < m_size; w++)
  {
    if (delta == ref_m[w][0]) 
      l = w;
  }

  if (i >= 0 && j >= 0 && k >= 0 && l >= 0) // only activates if all 4 "if" statements above are true
    h2_mat [i][j][k][l] = value;

  }

  // reset index term on reference matrix to span 0 to m_size instead of e.g. 3, 4, 9, 10, etc... for neutrons

  for (size_t n = 0; n < m_size; n++)
  ref_m[n][0] = n;

  print(std::cout, ref_m);


}


template <typename two_array, typename four_array>
void fullm_populate_hamiltonian (two_array & ref_m, two_array & h1_mat, four_array & h2_mat, const std::string reference_file, const std::string matrix_file, const double hw) 
{

  try
  {
    read_in_reference_m_scheme (ref_m, reference_file);
    populate_1body (ref_m, h1_mat, hw);
    read_in_matrix_m_scheme (ref_m, h2_mat, matrix_file);
  }

  catch (const char * msg)
  {
    std::cerr << msg << std::endl;  
  }
 
}

template <typename two_array, typename four_array>
void compactify_h2 (const two_array & ref_m, two_array & comp_h2, four_array & h2_mat)
{
  size_t mat_length = h2_mat.size();

  double value;

  for (size_t i = 0; i < mat_length; ++i)
  {
    for (size_t j = 0; j < mat_length; ++j)
    {
      for (size_t k = 0; k < mat_length; ++k)
      {
        for (size_t l = 0; l < mat_length; ++l)
        {
            value = h2_mat [i][j][k][l];

            if (value != 0.)
            {
            	std::cout << i << " " << j << " " << k << " " << l << " " << value << "\n";

            }


            size_t left  = j * mat_length + i;

            size_t right = k * mat_length + l; 

  /*         if (i == 0 && j == 4 && k == 0)
            {
            	std::cout << ref_m[l][1] << "\t" << ref_m[l][2] << "\n";
            }

           if (i == 0 && j == 1)
            {
            	std::cout << right << " " << l << " " << k << "\t" << value << "\n"; 
            }
*/
            comp_h2 [left][right] = value;

        }

      }

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
  const int bsize = input_params.s_wave_basis;	// the size of the HO basis used, ONLY ACTIVE FOR S-WAVE
  const int particles = input_params.particles;		// the number of neutrons (particles) in the trap
  const int nodes = input_params.nodes;
  const double hc = input_params.hc;					// hbar * c - given in MeV * fm
  const double mass = input_params.mass;				// Nucleon mass - given in MeV/c^2
  const double hw = input_params.hw;					// hbar * omega 
  const std::string m_ref = input_params.m_ref;    // single particle reference file - m scheme
  const std::string m_mat = input_params.m_mat;    // m scheme matrix elements

  const double l_param = hc / sqrt(mass * hw);		// the relevant length parameter for the HO


  std::cout << "Building system..." << std::endl;

  double weight[nodes], x[nodes];							// create arrays to hold gaussian weights/positions
  gauss (mesh_size, 2, 0., 4., x, weight);						// creating weights/positions for s-wave G-L quad.
/*  for (int i = 0; i < nodes; i++)
  {
	printf ("%f", x[i]);
	printf ("\n");
  }*/

 
  typedef boost::multi_array<double, 2> two_array;
  two_array ref_m (boost::extents[bsize][6]);
  two_array h1_mat(boost::extents[bsize][bsize]);

  typedef boost::multi_array<double, 4> four_array;
  four_array h2_mat(boost::extents[bsize][bsize][bsize][bsize]);











  fullm_populate_hamiltonian (ref_m, h1_mat, h2_mat, m_ref, m_mat, hw);

  print(std::cout, h1_mat); 

  two_array comp_h2 (boost::extents[bsize*bsize][bsize*bsize]);

  compactify_h2 (ref_m, comp_h2, h2_mat);

  print(std::cout, comp_h2);

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

  if (counter == 7)
    ss >> input.m_ref;

  if (counter == 8)
    ss >> input.m_mat;


	counter++;					// increase the counter after reading in a parameter	
  }

  return input;						// return the struct 
}


