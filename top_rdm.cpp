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

// struct to hold the input file variables

struct parameters 
{
	size_t nodes;
	size_t mesh_size;
	size_t basis;
	size_t particles;
	double hc;
	double mass;
	double hw;
  std::string m_ref;
  std::string m_mat;
};

struct con_flags
{
	bool F1_flag;
	bool F2_flag;
	bool F3_flag;
	bool F4_flag;
  bool F5_flag;
  bool F6_flag;
	bool diag_toggle;
	bool two_body_toggle;
  bool Q_flag;
  bool G_flag;
};


// function prototypes
extern void gauss (int npts, int job, double a, double b, double xpts[], double weights[]);
parameters read_in_inputs ();



/***************************************************************


 Function to populate the 1-body part of the Hamiltonian. 
 The 1-body part here is T + U for kinetic energy T and external potential U


***************************************************************/


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

/***************************************************************


 Function to print out an arbitrary array


***************************************************************/

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

/***************************************************************


 Individual print functions depending on whether template takes 
 double or int.


***************************************************************/


template<> void print<double>(std::ostream& os, const double & x)
{
  os << x << "\t";
}

template<> void print<int>(std::ostream& os, const int & x)
{
  os << x;
}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///        FUNCTIONS RELATING TO SWAVE IMPLEMENTATION
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///         END SWAVE IMPLEMENTATION
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


/***************************************************************

 
 Function to automate output of constraint matrices.


***************************************************************/

template <typename two_array>
void con_matrix_out (const two_array & m_pass, const size_t con_count, const size_t bsize, const struct con_flags flag_pass, std::ofstream & spda_out)
{

  size_t lower = 0;
  size_t upper = bsize;

  for (size_t i = lower; i < upper; i++)
  {
  for (size_t j = i;     j < upper; j++)
  {
    size_t k = i + 1;
    size_t l = j + 1;

    if (m_pass[i][j] != 0.0)
      spda_out << con_count << " " << 1 << " " << k << " " << l << " " << m_pass[i][j] << std::endl;
  }
  }

  lower += bsize;
  upper += bsize;

  for (size_t i = lower; i < upper; i++)
  {
  for (size_t j = i;     j < upper; j++)
  {
    size_t k = i + 1 - lower;
    size_t l = j + 1 - lower;

    if (m_pass[i][j] != 0.0)
      spda_out << con_count << " " << 2 << " " << k << " " << l << " " << m_pass[i][j] << std::endl;
  }
  }

  lower += bsize;
  upper += bsize*bsize;

  if (flag_pass.two_body_toggle)
  {
    for (size_t i = lower; i < upper; i++)
    {
    for (size_t j = i;     j < upper; j++)
    {
      size_t k = i + 1 - lower;
      size_t l = j + 1 - lower;

      if (m_pass[i][j] != 0.0)
        spda_out << con_count << " " << 3 << " " << k << " " << l << " " << m_pass[i][j] << std::endl;
    }
    }
  }

  // did this very quickly, check again!
  lower += bsize*bsize; 
  upper += bsize*bsize;

}



/***************************************************************

 
 Function to output the SPD file in correct format.


***************************************************************/

template <typename one_array, typename two_array, typename three_array>
void create_spda_file (const two_array c_matrix, const struct con_flags flag_pass, const two_array F1_con, const three_array F2_con, const one_array F2_val, const three_array F3_con, const one_array F3_val, const two_array F4_con, const three_array F5_con, const one_array F5_val, const int particles, const size_t bsize, std::ofstream & spda_out)
{
  size_t F2num = F2_con.size();
  size_t F3num = F3_con.size();
  size_t F5num = F5_con.size();
//  size_t cmat_extent = F1_con.size();


  // number of constraints start
  size_t num_cons = 0;  
  
  if (flag_pass.F1_flag)
  	num_cons += 1;	// N trace constraint 

  if (flag_pass.F2_flag)			
  	num_cons += F2num; //  p + q constraints  

  if (flag_pass.two_body_toggle)
  {
    if (flag_pass.F3_flag)
      num_cons += F3num;

  	if (flag_pass.F4_flag)
  		num_cons += 1; // N(N+1)/2 trace constraint

  	if (flag_pass.F5_flag)
  		num_cons += F5num;
  }

  spda_out << num_cons << std::endl;         // output number of constraints




  // number of blocks start
  size_t blocks = 2;

  if (flag_pass.two_body_toggle)
  {
    blocks ++;

    if (flag_pass.Q_flag)
      blocks++;

    if (flag_pass.G_flag)
      blocks++;
  }

  spda_out << blocks << std::endl; 						 // output number of blocks in X
 



  // block sizes start
  spda_out << bsize << " " << bsize << " ";     // output block sizes
  
  if (flag_pass.two_body_toggle)
  {
    spda_out << bsize*bsize << " ";

    if (flag_pass.Q_flag)
    	spda_out << "Q SIZE HERE" << " ";

    if (flag_pass.G_flag)
    	spda_out << "G SIZE HERE" << " ";

  }
  
  spda_out << std::endl;	             




  // constraint values start

  if (flag_pass.F1_flag)
    spda_out << particles << " "; 							 // output N trace constraint

  if (flag_pass.F2_flag)
  {
    for (size_t i = 0; i < F2num; i++)
      spda_out << F2_val[i] << " ";
  }

  if (flag_pass.F3_flag)
  {
    for (size_t i = 0; i < F3num; i++)
      spda_out << F3_val[i] << " ";
  }

  if (flag_pass.F4_flag)
    spda_out << (particles*(particles+1)/2) << " ";

  if (flag_pass.F5_flag)
  {
    for (size_t i = 0; i < F5num; i++)
      spda_out << F5_val[i] << " ";
  }

  spda_out << std::endl;






  // constraint matrix start

  size_t con_count = 0;

  con_matrix_out (c_matrix, con_count, bsize, flag_pass, spda_out);
  con_count++;

  if (flag_pass.F1_flag)
  {
  	con_matrix_out (F1_con, con_count, bsize, flag_pass, spda_out);
  	con_count++;
  }

  if (flag_pass.F2_flag)
  {
	for (size_t cnum = 0; cnum < F2num; cnum++)
	{
	   con_matrix_out (F2_con[cnum], con_count, bsize, flag_pass, spda_out);
	   con_count++;
	}
  }

  if (flag_pass.F3_flag)
  {
	for (size_t cnum = 0; cnum < F3num; cnum++)
	{
	   con_matrix_out (F3_con[cnum], con_count, bsize, flag_pass, spda_out);
	   con_count++;
	}
  }

  if (flag_pass.F4_flag)
  {
  	con_matrix_out (F4_con, con_count, bsize, flag_pass, spda_out);
  	con_count++;
  }

  if (flag_pass.F5_flag)
  {
	for (size_t cnum = 0; cnum < F5num; cnum++)
	{
	   con_matrix_out (F5_con[cnum], con_count, bsize, flag_pass, spda_out);
	   con_count++;
	}
  }

/*  for (size_t i = 0; i < bsize; i++)
  {
  for (size_t j = i; j < bsize; j++)
  {
    if (c_matrix[i][j] != 0.0)
      spda_out << con_count << " " << 1 << " " << i+1 << " " << j+1 << " " << c_matrix[i][j] << std::endl;
  }
  }

  for (size_t i = bsize; i < 2*bsize; i++)
  {
  for (size_t j = i; j < 2*bsize; j++)
  {
    if (c_matrix[i][j] != 0.0)
      spda_out << con_count << " " << 2 << " " << i+1-bsize << " " << j+1-bsize << " " << c_matrix[i][j] << std::endl;
  }
  }

  for (size_t i = 2*bsize; i < 2*bsize+bsize*bsize; i++)
  {
  for (size_t j = i; j < 2*bsize+bsize*bsize; j++)
  {
    if (c_matrix[i][j] != 0.0)
      spda_out << con_count << " " << 3 << " " << i+1-2*bsize << " " << j+1-2*bsize << " " << c_matrix[i][j] << std::endl;
  }
  }

  con_count++;*/
/*  
  for (size_t i = 0; i < bsize; i++)
  {
  for (size_t j = i; j < bsize; j++)
  {
    if (F1_con[i][j] != 0.0)
      spda_out << con_count << " " << 1 << " " << i+1 << " " << j+1 << " " << F1_con[i][j] << std::endl;
  }
  }

  for (size_t i = bsize; i < 2*bsize; i++)
  {
  for (size_t j = i; j < 2*bsize; j++)
  {
    if (F1_con[i][j] != 0.0)
      spda_out << con_count << " " << 2 << " " << i+1-bsize << " " << j+1-bsize << " " << F1_con[i][j] << std::endl;
  }
  }


  con_count++;*/

/*  for (size_t cnum = 0; cnum < F2num; cnum++)
  {
    for (size_t i = 0; i < bsize; i++)
    {
    for (size_t j = i; j < bsize; j++)
    {
      if (F2_con[cnum][i][j] != 0.0)
       spda_out << con_count << " " << 1 << " " << i+1 << " " << j+1 << " " << F2_con[cnum][i][j] << std::endl;
    }
    }

    for (size_t i = bsize; i < 2*bsize; i++)
    {
    for (size_t j = i; j < 2*bsize; j++)
    {
      if (F2_con[cnum][i][j] != 0.0)
        spda_out << con_count << " " << 2 << " " << i+1-bsize << " " << j+1-bsize << " " << F2_con[cnum][i][j] << std::endl;
    }
    }

    con_count++;

  }
*/
//  print (std::cout, F4_con);


//  std::cout << bsize*bsize << std::endl;

/*  for (size_t i = 2*bsize; i < 2*bsize+bsize*bsize; i++)
  {
  for (size_t j = i; j < 2*bsize+bsize*bsize; j++)
  {
	std::cout << F4_con[i][j] << "\t" << i << " " << j << std::endl;  	
      if (F4_con[i][j] != 0.0)
       spda_out << con_count << " " << 3 << " " << i+1-2*bsize << " " << j+1-2*bsize << " " << F4_con[i][j] << std::endl;
  }
  }

*/

}

/***************************************************************

 
 Function to read in the reference file for m scheme.


***************************************************************/

template <typename ref_array>
void read_in_reference_m_scheme (ref_array & ref_m, const std::string m_ref_file, std::ofstream & diag_out, const bool diag_toggle)
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

//  std::cout << ref_num << " " << n << " " << l << " " << j << " " << m_j << " " << tz << "\n";


  // only extract neutron-neutron states
  if (tz == 1.)
  {
    ref_m [ele_in][0] = ref_num;    // reference number of the line
    ref_m [ele_in][1] = n;          // principle quantum number
    ref_m [ele_in][2] = l;          // orbital angular mom.
    ref_m [ele_in][3] = j * 0.5;    // total angular mom.
    ref_m [ele_in][4] = m_j * 0.5;  // total angular mom. projection
    ref_m [ele_in][5] = tz;         // isospin projection (should all be +1.0)
    ref_m [ele_in][6] = ele_in;     // new index for sp orbital
    ele_in++;
  }

  if (ele_in >= ref_size)
    break;

  }
  
  
  
  if(diag_toggle)
  {
    diag_out << "Single particle orbitals pulled from REF file" << std::endl << std::endl;
    print(diag_out, ref_m);            // print the resulting matrix
    diag_out << std::endl << std::endl;
  }


  
}

/***************************************************************

Function to read in the actual m scheme matrix  elements 
calculated from Morten's code. 

THIS WILL BREAK IF THE INPUT FILE .dat CHANGES TITLE OR FORMAT.

***************************************************************/

template <typename ref_array, typename five_array>
void read_in_matrix_m_scheme (const ref_array & ref_m, five_array & h2_mat, const std::string m_mat_file)
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

  const size_t m_size = h2_mat.size();

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
    {
      h2_mat [i][j][k][l][0] = value;
      h2_mat [i][j][k][l][1] = alpha;
      h2_mat [i][j][k][l][2] = beta;
      h2_mat [i][j][k][l][3] = gamma;
      h2_mat [i][j][k][l][4] = delta;
    }
  }

  // reset index term on reference matrix to span 0 to m_size instead of e.g. 3, 4, 9, 10, etc... for neutrons

//  for (size_t n = 0; n < m_size; n++)
//  ref_m[n][0] = n;

//  print(std::cout, ref_m);


}

/***************************************************************

Wrapper function to read in single-particle m scheme reference 
file, populate 1-body Hamiltonian matrix elements, and then 
populate 2-body matrix elements.

***************************************************************/


template <typename two_array, typename five_array>
void fullm_populate_hamiltonian (two_array & ref_m, two_array & h1_mat, five_array & h2_mat, const std::string reference_file, const std::string matrix_file, const double hw, std::ofstream & diag_out, const bool diag_toggle) 
{

  try
  {
    read_in_reference_m_scheme (ref_m, reference_file, diag_out, diag_toggle);
    populate_1body (ref_m, h1_mat, hw);
    read_in_matrix_m_scheme (ref_m, h2_mat, matrix_file);
  }

  catch (const char * msg)
  {
    std::cerr << msg << std::endl;  
  }
 
}


/***************************************************************



***************************************************************/

template <typename two_array, typename five_array>
void compactify_h2 (const two_array & ref_m, two_array & comp_h2, five_array & h2_mat, std::ofstream & diag_out, const bool diag_toggle)
{
  const size_t mat_length = h2_mat.size();

  double value;

  size_t alpha, beta, gamma, delta;

  if(diag_toggle)
  {
    diag_out << "Output of nonzero 2-body ME" << std::endl << std::endl;    
  }

  for (size_t i = 0; i < mat_length; ++i)
  {
    for (size_t j = 0; j < mat_length; ++j)
    {
      for (size_t k = 0; k < mat_length; ++k)
      {
        for (size_t l = 0; l < mat_length; ++l)
        {
            value = h2_mat [i][j][k][l][0];

            alpha = h2_mat [i][j][k][l][1];
            beta  = h2_mat [i][j][k][l][2];
            gamma = h2_mat [i][j][k][l][3];
            delta = h2_mat [i][j][k][l][4];

            if (diag_toggle)
            {
              if (value != 0.)
              {
              	diag_out << alpha << " " << beta << " " << gamma << " " << delta << " " << value << " \t";

                for (size_t loop = 0; loop < mat_length; loop++)
                {
                  if (alpha == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\t";
                }


                for (size_t loop = 0; loop < mat_length; loop++)
                {
                  if (beta == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\t";
                }


                for (size_t loop = 0; loop < mat_length; loop++)
                {
                  if (gamma == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\t";
                }


                for (size_t loop = 0; loop < mat_length; loop++)
                {
                  if (delta == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\n";
                }

              }
            }


            size_t left  = i * mat_length + j;

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

  if (diag_toggle)
    diag_out << std::endl << std::endl;
}

/***************************************************************

Function to create the F0 constraint matrix for the SDP solver.
This matrix is just the Hamiltonian matrix elements for the 1RDM
and 2RDM parts (where applicable) with zeros in the other parts
of the matrix.

Note that since the SDP solver maximizes the objective function,
here we multiply all elements by -1 so that we are maximizing
the negative energy.

***************************************************************/

template <typename two_array>
void create_c_matrix (const two_array & h1_mat, const two_array & h2_mat, two_array & c_matrix, const bool two_body_toggle)
{
  size_t h1_len = h1_mat.size();
  size_t h2_len = h2_mat.size();

  for (size_t i = 0; i < h1_len; i++) 
  {
    for (size_t j = 0; j < h1_len; j++)
    {
      c_matrix[i][j] = h1_mat [i][j] * -1.0;
    }

  }

  if (two_body_toggle)
  {
    for (size_t i = 0; i < h2_len; i++) 
    {
    for (size_t j = 0; j < h2_len; j++)
    {
      size_t ic = i + 2*h1_len;
      size_t jc = j + 2*h1_len;

      c_matrix[ic][jc] = h2_mat [i][j] * -1.0;
    }
    }
  }


  print(std::cout, c_matrix);
}

/***************************************************************

Kronecker delta function. Returns 1 if i=j and 0 otherwise.

***************************************************************/

double kron_del(const size_t i, const size_t j)
{

  if (i == j)
    return 1.;

  return 0.;
}


/***************************************************************

Function to create our F1 constraint matrix to fix the total
particle number N (fixes trace of the 1RDM to be N). 

In keeping with other constraint functions, first populate the
F1_con_1 matrix (only one constraint) then populate the full F1_con
matrix. 

Note that F1_con_1 is a basis x basis size matrix while F1_con 
matches the dimensions of the F0 constraint matrix.

***************************************************************/

template <typename two_array>
void init_F1_flag (two_array & F1_con_1, two_array & F1_con, const size_t bsize)
{

//  size_t cmat_extent = F1_con.size();

  for (size_t i = 0; i < bsize; i++)
  {
    F1_con_1[i][i] = 1.;    
  }

  for (size_t i = 0; i < bsize; i++)
  {
    if (i < bsize)
      F1_con[i][i] = F1_con_1[i][i];
  }
}


/***************************************************************

Function to create our F2 constraint matrices to fix the relation
between the 1RDM and its twin q. 

First populate the F2_con_1 matrix which serves as the basis x basis
size constraint matrix for the 1RDM sector and the q sector (both 
have the same constraint matrix).

Then populate this matrix into the full F0 sized constraint term 
F2_con. 

***************************************************************/


template <typename three_array, typename one_array, typename four_array>
void init_F2_flag (four_array & F2_con_1, three_array & F2_con, one_array & F2_val, const size_t bsize)
{


//  size_t F2num = F2_con.size();           // number of F2 constraint matrices
//  size_t bsize = F2_con_1.size();         // basis size

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
    for (size_t k = 0; k < bsize; k++)    // loop over matrix row
    {
    for (size_t l = 0; l < bsize; l++)    // loop over matrix column
    {

        F2_con_1[i][j][k][l] = (1./2.)*(kron_del(i,k)*kron_del(j,l) + kron_del(i,l)*kron_del(j,k));
    }
    }

  }
  }

  size_t counter = 0;

  for (size_t i1 = 0;  i1 < bsize; i1++)
  {
  for (size_t i2 = i1; i2 < bsize; i2++)
  {
    for (size_t j = 0; j < bsize; j++)
    {
    for (size_t k = 0; k < bsize; k++)
    {
      size_t jq = j + bsize;
      size_t kq = k + bsize;

      size_t i = counter;

      F2_con[i][j][k]   = F2_con_1 [i1][i2][j][k];
      F2_con[i][jq][kq] = F2_con_1 [i1][i2][j][k];
    }
    }

    if (i1 == i2)
      F2_val[counter] = 1.;
    else
      F2_val[counter] = 0.;

    counter++;
  }
  }

//  printf("%i", F2num);
//  printf("\n");
//  printf("%i", bsize);
//  printf("\n");
}


/***************************************************************

Function to create our F3 constraint matrix to fix the relation
between the 2RDM partial trace and the 1RDM.

***************************************************************/

template <typename one_array, typename three_array, typename four_array, typename six_array>
void init_F3_flag (four_array & F3_build_1, six_array & F3_build_3, three_array & F3_con, one_array & F3_val, const size_t bsize, const size_t N)
{

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)      // loop over jth constraint matrix
  {
    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
      F3_build_1 [i][k][ip][kp] = -1.0 * (N - 1.0) / 4.0 * (kron_del(i,ip)*kron_del(k,kp) + kron_del(k,ip)*kron_del(i,kp));
    }
    }
  }
  }

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)      // loop over jth constraint matrix
  {
    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
     F3_build_3 [i][k][ip][jp][kp][lp] = 1./2. * (kron_del(i,ip)*kron_del(k,kp) + kron_del(k,ip)*kron_del(i,kp)) * kron_del(jp,lp);
    }
    }
    }
    }
  }
  }


  size_t counter = 0;

  for (size_t i1 = 0; i1 < bsize; i1++)
  {
  for (size_t i2 = i1; i2 < bsize; i2++)
  {
    for (size_t j = 0; j < bsize; j++)
    {
    for (size_t k = 0; k < bsize; k++)
    {

      size_t b = counter;

      F3_con[b][j][k] = F3_build_1 [i1][i2][j][k];


    }
    }

    counter++;
  }
  }


  counter = 0;

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t k = i; k < bsize; k++)      // loop over jth constraint matrix
  {
    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {

      size_t b = counter;

      size_t left  = ip * bsize + jp;
      size_t right = kp * bsize + lp;

      left  += 2*bsize;
      right += 2*bsize;

      F3_con[b][left][right] = F3_build_3 [i][k][ip][jp][kp][lp];

    }
    }
    }
    }

    F3_val[counter] = 0.;

    counter++;
  }
  }

}

/***************************************************************

Function to create our F4 constraint matrix to fix the 2RDM trace.

***************************************************************/

template <typename two_array>
void init_F4_flag (two_array & F4_con, const size_t bsize)
{


//  size_t cmat  = F4_con.size();           // number of F2 constraint matrices
//  size_t bsize = F4_con_3.size();         // basis size

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
    for (size_t k = 0; k < bsize; k++)    // loop over matrix row
    {
    for (size_t l = 0; l < bsize; l++)    // loop over matrix column
    {

        size_t left  = i * bsize + j;

        size_t right = k * bsize + l; 

        left  += 2*bsize;
        right += 2*bsize;

        F4_con [left][right] = kron_del(i,k) * kron_del(j,l);

    }
    }

  }
  }

}

/***************************************************************

Function to create our F5 constraint matrix to fix antisymmetry 
of the 2RDM in the upper two indices.

***************************************************************/

template <typename one_array, typename three_array, typename eight_array>
void init_F5_flag (eight_array & F5_build_3, three_array & F5_con, one_array & F5_val, const size_t bsize)
{


//  size_t cmat  = F4_con.size();           // number of F2 constraint matrices
//  size_t bsize = F4_con_3.size();         // basis size

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = 0; l < bsize; l++)    // loop over matrix column
  {
  	for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
  	{
  	for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
  	{
  	for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
  	{
  	for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
  	{

  		F5_build_3 [i][j][k][l][ip][jp][kp][lp] = (1./2.)* (
  			(kron_del(i,ip)*kron_del(j,jp) + kron_del(i,jp)*kron_del(j,ip)) * kron_del(k,kp)*kron_del(l,lp)
  			+
  			(kron_del(i,kp)*kron_del(j,lp) + kron_del(i,lp)*kron_del(j,kp)) * kron_del(k,ip)*kron_del(l,jp)
  			);
  	}
  	}
  	}
  	}
  }
  }
  }
  }

  size_t counter = 0;

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = 0; l < bsize; l++)    // loop over matrix column
  {
  	for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
  	{
  	for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
  	{
  	for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
  	{
  	for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
  	{
  		size_t b = counter;

  		size_t left  = ip * bsize + jp;
  		size_t right = kp * bsize + lp;

  		left  += 2*bsize;
  		right += 2*bsize;

  		F5_con[b][left][right] = F5_build_3 [i][j][k][l][ip][jp][kp][lp];
  	}
  	}
  	}
  	}

  	F5_val [counter] = 0.;

  	counter++;
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
  
  parameters input_params;						             // create a struct to hold user specified parameters
  input_params = read_in_inputs (); 					     // read in the values from inputs.inp, store in input_params	

  // Now create/copy input parameters for ease of reading the code into newly defined local variables


//  const int mesh_size = input_params.mesh_size;		 // how many points to use for Gauss-Legendre (G-L) quadrature 
  const size_t bsize = input_params.basis;	   // the size of the HO basis used, ONLY ACTIVE FOR S-WAVE
  const size_t particles = input_params.particles;		 // the number of neutrons (particles) in the trap
//  const int nodes = input_params.nodes;
//  const double hc = input_params.hc;					     // hbar * c - given in MeV * fm
//  const double mass = input_params.mass;			     // Nucleon mass - given in MeV/c^2
  const double hw = input_params.hw;					     // hbar * omega 
  const std::string m_ref = input_params.m_ref;    // single particle reference file - m scheme
  const std::string m_mat = input_params.m_mat;    // m scheme matrix elements

//  const double l_param = hc / sqrt(mass * hw);		 // the relevant length parameter for the HO

  std::cout << "Building system..." << std::endl;

//  const char * m_matrix_file = (m_mat_file).c_str();
  


  // define flags for all different combinations of conditions in RDM

  const bool F1_flag = true;
  const bool F2_flag = true;
  const bool F3_flag = true;
  const bool F4_flag = false;
  const bool F5_flag = true;
  const bool F6_flag = false;

  const bool two_body_toggle = true;
  const bool Q_flag = false;
  const bool G_flag = false;

  if (!two_body_toggle and (F4_flag or F5_flag or F6_flag))
  {
  	std::cerr << "ERROR: TWO-BODY AND F CONSTRAINT FLAGS INCOMPATIBLE - GIVING UP" << std::endl;
  	return EXIT_FAILURE;
  }	

  const bool diag_toggle = true;


  struct con_flags flag_pass;

  flag_pass.F1_flag = F1_flag;
  flag_pass.F2_flag = F2_flag;
  flag_pass.F3_flag = F3_flag;
  flag_pass.F4_flag = F4_flag;
  flag_pass.F5_flag = F5_flag;
  flag_pass.F6_flag = F6_flag;

  flag_pass.two_body_toggle = two_body_toggle;

  flag_pass.Q_flag = Q_flag;
  flag_pass.G_flag = G_flag;

  flag_pass.diag_toggle = true;


//  const size_t F2num = bsize*bsize;
  const size_t F2num = bsize * (bsize + 1)/2;
  const size_t F3num = bsize * (bsize + 1)/2;
//  const size_t F3num = bsize*bsize;
  const size_t F5num = bsize*bsize*bsize*bsize;

  // turn off or on two-body potential



  size_t cmat_extent = 0;

  if (two_body_toggle)
  {
    cmat_extent = 2*bsize + bsize*bsize;
  }

  else
  {
    cmat_extent = 2*bsize;
  }

  const std::string diag_file = "diagnostic_out/test_diag.dat";
  const std::string spda_file = "spd_files/test_spd.dat";

  std::ofstream diag_out (diag_file);
  std::ofstream spda_out (spda_file);

//  double weight[nodes], x[nodes];							// create arrays to hold gaussian weights/positions
//  gauss (mesh_size, 2, 0., 4., x, weight);						// creating weights/positions for s-wave G-L quad.
/*  for (int i = 0; i < nodes; i++)
  {
	printf ("%f", x[i]);
	printf ("\n");
  }*/

 
  typedef boost::multi_array<double, 1> one_array;
  typedef boost::multi_array<double, 2> two_array;
  typedef boost::multi_array<double, 3> three_array;
  typedef boost::multi_array<double, 4> four_array;
  typedef boost::multi_array<double, 5> five_array;
  typedef boost::multi_array<double, 6> six_array;
  typedef boost::multi_array<double, 7> seven_array;
  typedef boost::multi_array<double, 8> eight_array;


  two_array ref_m (boost::extents[bsize][7]);
  two_array h1_mat(boost::extents[bsize][bsize]);

  five_array h2_mat(boost::extents[bsize][bsize][bsize][bsize][5]);


  two_array  F1_con_1 (boost::extents[bsize][bsize]);
  four_array F2_con_1 (boost::extents[bsize][bsize][bsize][bsize]);
 
  four_array F3_build_1 (boost::extents[bsize][bsize][bsize][bsize]);
  six_array  F3_build_3 (boost::extents[bsize][bsize][bsize][bsize][bsize][bsize]);

  four_array F4_con_3 (boost::extents[bsize][bsize][bsize][bsize]);

  eight_array F5_build_3 (boost::extents[bsize][bsize][bsize][bsize][bsize][bsize][bsize][bsize]);

  two_array   F1_con (boost::extents[cmat_extent][cmat_extent]);
  three_array F2_con (boost::extents[F2num][cmat_extent][cmat_extent]);
  three_array F3_con (boost::extents[F3num][cmat_extent][cmat_extent]);
  two_array	  F4_con (boost::extents[cmat_extent][cmat_extent]);
  three_array F5_con (boost::extents[F5num][cmat_extent][cmat_extent]);

  one_array   F2_val (boost::extents[F2num]);
  one_array   F3_val (boost::extents[F3num]);
  one_array	  F5_val (boost::extents[F5num]);

  if (F1_flag)
  {
    init_F1_flag (F1_con_1, F1_con, bsize);
  }

  if (F2_flag)
  {
    init_F2_flag (F2_con_1, F2_con, F2_val, bsize);
  }

  if (F3_flag)
  {
    init_F3_flag (F3_build_1, F3_build_3, F3_con, F3_val, bsize, particles);
  }

  if (F4_flag)
  {
    init_F4_flag (F4_con, bsize);
  }

  if (F5_flag)
  {
    init_F5_flag (F5_build_3, F5_con, F5_val, bsize);
  }

  fullm_populate_hamiltonian (ref_m, h1_mat, h2_mat, m_ref, m_mat, hw, diag_out, diag_toggle);



  two_array comp_h2 (boost::extents[bsize*bsize][bsize*bsize]);

  compactify_h2 (ref_m, comp_h2, h2_mat, diag_out, diag_toggle);

  print(std::cout, comp_h2);

  two_array c_matrix (boost::extents[cmat_extent][cmat_extent]);

  create_c_matrix (h1_mat, comp_h2, c_matrix, two_body_toggle);

  create_spda_file (c_matrix, flag_pass, F1_con, F2_con, F2_val, F3_con, F3_val, F4_con, F5_con, F5_val, particles, bsize, spda_out);

  four_array test_h2 (boost::extents[bsize][bsize][bsize][bsize]);

  double q = 0.;

  for (size_t i = 0; i < bsize; i++)
  {
  for (size_t j = 0; j < bsize; j++)
  {
  for (size_t k = 0; k < bsize; k++)
  {
  for (size_t l = 0; l < bsize; l++)
  {
    test_h2 [i][j][k][l] = h2_mat [i][j][k][l][0];//q;//h2_mat [i][j][k][l][0];
    q++;    
  }
  }
  }
  }


  if (diag_toggle)
  {

    diag_out << "1 body Hamiltonian matrix" << std::endl << std::endl;

    print(diag_out, h1_mat); 

    diag_out << std::endl << std::endl;

    if (two_body_toggle)
    {
      diag_out << "Original 2 body Hamiltonian matrix" << std::endl << std::endl;
      print (diag_out, test_h2);
      diag_out << std::endl << std::endl;

      diag_out << "Compacted 2 body Hamiltonian matrix" << std::endl << std::endl;
      print (diag_out, comp_h2);
      diag_out << std::endl << std::endl;
    }

    else
      diag_out << "Two-body toggle flag set to false - No 2-body output" << std::endl << std::endl;



  	diag_out << "F0 Constraint Matrix" << std::endl << std::endl;

  	print(diag_out, c_matrix);

  	diag_out << std::endl << std::endl;



  	diag_out << "1 term - " << "F1 constraint matrix output" << std::endl << std::endl;

  	print(diag_out, F1_con);

  	diag_out << std::endl << std::endl;



  	diag_out << F2num << " terms - " << "F2 constraint matrix output" << std::endl << std::endl;

  	print(diag_out, F2_con);
  	
  	diag_out << std::endl << std::endl;



  	diag_out << "F2 constraint values" << std::endl << std::endl;

  	print(diag_out, F2_val);

  	diag_out << std::endl << std::endl;



    diag_out << F3num << " terms - " << "F3 constraint matrix output" << std::endl << std::endl;

    print(diag_out, F3_con);
    
    diag_out << std::endl << std::endl;



    diag_out << "F3 constraint values" << std::endl << std::endl;

    print(diag_out, F3_val);

    diag_out << std::endl << std::endl;




  	diag_out << "1 term - " << "F4 constraint matrix output" << std::endl << std::endl;

  	print(diag_out, F4_con);

  	diag_out << std::endl << std::endl;




    diag_out << F5num << " terms - " << "F5 constraint matrix output" << std::endl << std::endl;

    print(diag_out, F5_con);
    
    diag_out << std::endl << std::endl;



    diag_out << "F5 constraint values" << std::endl << std::endl;

    print(diag_out, F5_val);

    diag_out << std::endl << std::endl;

  }

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
	  ss >> input.basis;

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


