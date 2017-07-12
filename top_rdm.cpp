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
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>


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
  bool F7_flag;
  bool F8_flag;
  bool F9_flag;
  bool F10_flag;
	bool diag_toggle;
	bool two_body_toggle;
  bool redundant_check;
  bool Q_flag;
  bool G_flag;
};


// function prototypes
extern void gauss (int npts, int job, double a, double b, double xpts[], double weights[]);
parameters read_in_inputs ();


  typedef boost::multi_array<double, 1> one_array;
  typedef boost::multi_array<double, 2> two_array;
  typedef boost::multi_array<double, 3> three_array;
  typedef boost::multi_array<double, 4> four_array;
  typedef boost::multi_array<double, 5> five_array;
  typedef boost::multi_array<double, 6> six_array;
//  typedef boost::multi_array<double, 7> seven_array;
  typedef boost::multi_array<double, 8> eight_array;





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




/***************************************************************

 
 Function to automate output of constraint matrices.


***************************************************************/

template <typename two_array>
void con_matrix_out (const two_array & m_pass, const size_t con_count, const size_t bsize, const struct con_flags flag_pass, std::ofstream & spda_out)
{


  const size_t PQsize = bsize*(bsize-1)/2;

  size_t lower = 0;
  size_t upper = bsize;

  for (size_t i = lower; i < upper; i++)
  {
  for (size_t j = lower; j < upper; j++)
  {
      spda_out << m_pass[i][j] << " ";
  }
  }

  lower += bsize;
  upper += bsize;

  for (size_t i = lower; i < upper; i++)
  {
  for (size_t j = lower; j < upper; j++)
  {
      spda_out << m_pass[i][j] << " ";
  }
  }

  lower += bsize;
  upper += PQsize;

  if (flag_pass.two_body_toggle)
  {
    for (size_t i = lower; i < upper; i++)
    {
    for (size_t j = lower; j < upper; j++)
    {
      spda_out << m_pass[i][j] << " ";
    }
    }


    if (flag_pass.Q_flag)
    {

      lower += PQsize; 
      upper += PQsize;

      for (size_t i = lower; i < upper; i++)
      {
      for (size_t j = lower;     j < upper; j++)
      {
      	spda_out << m_pass[i][j] << " ";
      }
      }

    }

    if (flag_pass.G_flag)
    {

      lower += PQsize; 
      upper += bsize*bsize;

      for (size_t i = lower; i < upper; i++)
      {
      for (size_t j = lower;     j < upper; j++)
      {
      	spda_out << m_pass[i][j] << " ";
      }
      }

    }

  }

  spda_out << std::endl;

/*
  const size_t PQsize = bsize*(bsize-1)/2;

  size_t lower = 0;
  size_t upper = bsize;

  for (size_t i = lower; i < upper; i++)
  {
  for (size_t j = i;     j < upper; j++)
  {
//    size_t k = i + 1 - lower;
//    size_t l = j + 1 - lower;

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
  upper += PQsize;

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


    if (flag_pass.Q_flag)
    {

      lower += PQsize; 
      upper += PQsize;

      for (size_t i = lower; i < upper; i++)
      {
      for (size_t j = i;     j < upper; j++)
      {
        size_t k = i + 1 - lower;
        size_t l = j + 1 - lower;

        if (m_pass[i][j] != 0.0)
          spda_out << con_count << " " << 4 << " " << k << " " << l << " " << m_pass[i][j] << std::endl;
      }
      }

    }

    if (flag_pass.G_flag)
    {

      lower += PQsize; 
      upper += bsize*bsize;

      for (size_t i = lower; i < upper; i++)
      {
      for (size_t j = i;     j < upper; j++)
      {
        size_t k = i + 1 - lower;
        size_t l = j + 1 - lower;

        if (m_pass[i][j] != 0.0)
          spda_out << con_count << " " << 5 << " " << k << " " << l << " " << m_pass[i][j] << std::endl;
      }
      }

    }

  }


*/
}

/***************************************************************

 
 Function to look for redundancies in constraint terms.


***************************************************************/

template <typename three_array>
void check_constraints (const three_array & Con, const std::string & label, const bool output, const size_t num, const size_t cmat_extent)
{

  std::cout << label << " REDUNDANCY CHECK" << std::endl;

  size_t total = 0;


  for (size_t b = 0; b < num; b++)
  {
  for (size_t q = b; q < num; q++)
  {
    if (q == b)
      continue;

    bool mat_used = false;

    for (size_t l1 = 0;  l1 < cmat_extent; l1++)
    {
    for (size_t l2 = l1; l2 < cmat_extent; l2++)
    {
      if (Con[q][l1][l2] != 0.)
        mat_used = true;
    }
    }

    if (!mat_used)
    	continue;

    if (Con[b] == Con[q])
    {

      total++;
      if (output)
      {
        std::cout << "b " << b << "\t" << "q " << q << std::endl;
        if (false)
        {
          print (std::cout, Con[b]);
          std::cout << "\n";
          print (std::cout, Con[q]);
        }
      }
    }
  }
  }

  std::cout << "Total redundancies: " << total << std::endl;
}



/***************************************************************

 
 Function to output the SPD file in correct format.


***************************************************************/

template <typename one_array, typename two_array, typename three_array>
void create_spda_file (const two_array & c_matrix, const struct con_flags flag_pass, const three_array & F1_con, const one_array & F1_val, const three_array & F2_con, const one_array & F2_val, const three_array & F3_con, const one_array & F3_val, const three_array & F4_con, const one_array & F4_val, const size_t F7num, const three_array & F7_con, const one_array & F7_val, const three_array & F10_con, const one_array & F10_val, const size_t bsize, const size_t cmat_extent, std::ofstream & spda_out)
{
  const size_t F1num = F1_con.size();
  const size_t F2num = F2_con.size();
  const size_t F3num = F3_con.size();
  const size_t F4num = F4_con.size();



  const size_t F10num = F10_con.size();
  //  size_t cmat_extent = F1_con.size();

  const size_t PQsize = bsize*(bsize-1)/2;

  if(flag_pass.redundant_check)
  { 
    check_constraints (F1_con, "F1", true, F1num, cmat_extent);
    check_constraints (F2_con, "F2", true, F2num, cmat_extent);
    check_constraints (F3_con, "F3", true, F3num, cmat_extent);
    check_constraints (F4_con, "F4", true, F4num, cmat_extent);

    check_constraints (F7_con, "F7", false, F7num, cmat_extent);

    check_constraints (F10_con, "F10", true, F10num, cmat_extent);                                
  }


  // number of constraints start
  size_t num_cons = 0;  
  
  if (flag_pass.F1_flag)
  	num_cons += F1num;	// N trace constraint 

  if (flag_pass.F2_flag)			
  	num_cons += F2num; //  p + q constraints  

  if (flag_pass.two_body_toggle)
  {
    if (flag_pass.F3_flag)
      num_cons += F3num;

  	if (flag_pass.F4_flag)
  		num_cons += F4num; // N(N+1)/2 trace constraint


    if (flag_pass.F7_flag)
      num_cons += F7num;


    if (flag_pass.F10_flag)
      num_cons += F10num;
  }

  spda_out << num_cons << std::endl;         // output number of constraints




  // number of blocks start
  size_t blocks = 2;

  if (flag_pass.two_body_toggle)
  {
    blocks++;

    if (flag_pass.Q_flag and PQsize > 1)
      blocks++;

    if (flag_pass.G_flag)
      blocks++;
  }

  spda_out << blocks << std::endl; 						 // output number of blocks in X
 



  // block sizes start
  if (flag_pass.F1_flag)
    spda_out << bsize << " ";

  if (flag_pass.F2_flag)
    spda_out << bsize << " ";     // output block sizes
  
  if (flag_pass.two_body_toggle)
  {
    if (PQsize > 1)
      spda_out << PQsize << " ";

    if (!flag_pass.Q_flag and PQsize == 1)
      spda_out << -1 << " ";

    if (flag_pass.Q_flag)
    {
      if (PQsize > 1)
    	 spda_out << PQsize << " ";

      else
        spda_out << -2 << " ";
    }


    if (flag_pass.G_flag)
    	spda_out << bsize*bsize << " ";

  }
  
  spda_out << std::endl;	             




  // constraint values start

  if (flag_pass.F1_flag)
  {
    for (size_t i = 0; i < F1num; i++)
      spda_out << F1_val[i] << " ";
  }
						 // output N trace constraint

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
  {
    for (size_t i = 0; i < F4num; i++)
      spda_out << F4_val[i] << " ";
  }


  if (flag_pass.F7_flag)
  {
    for (size_t i = 0; i < F7num; i++)
      spda_out << F7_val[i] << " ";
  }


  if (flag_pass.F10_flag)
  {
    for (size_t i = 0; i < F10num; i++)
      spda_out << F10_val[i] << " ";
  }

  spda_out << std::endl;






  // constraint matrix start

  size_t con_count = 0;



  con_matrix_out (c_matrix, con_count, bsize, flag_pass, spda_out);
  con_count++;


  if (flag_pass.F1_flag)
  {
  for (size_t cnum = 0; cnum < F1num; cnum++)
  {
     con_matrix_out (F1_con[cnum], con_count, bsize, flag_pass, spda_out);
     con_count++;
  }
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
  for (size_t cnum = 0; cnum < F4num; cnum++)
  {
     con_matrix_out (F4_con[cnum], con_count, bsize, flag_pass, spda_out);
     con_count++;
  }
  }

 

  if (flag_pass.F7_flag)
  {
  for (size_t cnum = 0; cnum < F7num; cnum++)
  {
     con_matrix_out (F7_con[cnum], con_count, bsize, flag_pass, spda_out);
     con_count++;
  }
  }

 


  if (flag_pass.F10_flag)
  {
  for (size_t cnum = 0; cnum < F10num; cnum++)
  {
     con_matrix_out (F10_con[cnum], con_count, bsize, flag_pass, spda_out);
     con_count++;
  }
  }


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

    if (beta == ref_m[w][0]) 
      j = w;

    if (gamma == ref_m[w][0]) 
      k = w;

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
    std::cout << "REFERENCE READ" << std::endl;
    populate_1body (ref_m, h1_mat, hw);
    std::cout << "1 BODY POPULATED" << std::endl;
    read_in_matrix_m_scheme (ref_m, h2_mat, matrix_file);
    std::cout << "2 BODY POPULATED" << std::endl;
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
  const size_t bsize = h2_mat.size();

  double value;

  size_t alpha, beta, gamma, delta;

  if(diag_toggle)
  {
    diag_out << "Output of nonzero 2-body ME" << std::endl << std::endl;    
  }

  for (size_t i = 0; i < bsize; ++i)
  {
  for (size_t j = 0; j < bsize; ++j)
  {
    if (i >= j)
      continue;

  for (size_t k = 0; k < bsize; ++k)
  {
  for (size_t l = 0; l < bsize; ++l)
  {
    if (k >= l)
      continue;

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

                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (alpha == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\t";
                }


                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (beta == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\t";
                }


                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (gamma == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\t";
                }


                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (delta == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\n";
                }

              }
            }

            size_t ips = i + 1;
            size_t jps = j + 1;
            size_t kps = k + 1;
            size_t lps = l + 1;

            size_t left  = jps - ips + (2*bsize - ips) * (ips - 1)/2 - 1;
            size_t right = lps - kps + (2*bsize - kps) * (kps - 1)/2 - 1;


  /*         if (i == 0 && j == 4 && k == 0)
            {
            	std::cout << ref_m[l][1] << "\t" << ref_m[l][2] << "\n";
            }

           if (i == 0 && j == 1)
            {
            	std::cout << right << " " << l << " " << k << "\t" << value << "\n"; 
            }
*/

//            std::cout << value << std::endl;

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


void init_C_matrix (const con_flags flag_pass, std::ofstream & spda_out, const two_array & h1_mat, const two_array & h2_mat, size_t & con_count)
{
  size_t h1_len = h1_mat.size();
  size_t h2_len = h2_mat.size();


    for (size_t ip = 0;  ip < h1_len; ip++)
    {
    for (size_t jp = ip; jp < h1_len; jp++)
    {

	    double val1 = h1_mat [ip][jp] * -1.0;

      size_t n = ip + 1;
      size_t m = jp + 1; 

      if (val1 != 0.)
        spda_out << con_count << " " << 1 << " " << n << " " << m << " " << val1 << "\n";
    }
    }





	if (flag_pass.F3_flag)
	{
	    for (size_t ip = 0; ip < h2_len; ip++)      // loop over ith constraint matrix
	    {
	    for (size_t jp = ip; jp < h2_len; jp++)      // loop over jth constraint matrix
	    {

	      double val3 = h2_mat [ip][jp] * -1./2.;

        size_t n = ip + 1;
        size_t m = jp + 1; 

      if (val3 != 0.)
        spda_out << con_count << " " << 3 << " " << n << " " << m << " " << val3 << "\n";

	    }
	    }
	}

  con_count++;

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


void init_con_values (const con_flags flag_pass, std::ofstream & spda_out, const size_t bsize, const size_t particles)
{

//  size_t cmat_extent = F1_con.size();


  // F1 Flag - N Trace condition
	if (flag_pass.F1_flag)
	  spda_out << particles << " ";

  // F2 Flag - Linear p q relations
	if (flag_pass.F2_flag)
	{
	  for (size_t i1 = 0;  i1 < bsize; i1++)
	  {
	  for (size_t i2 = i1; i2 < bsize; i2++)
	  {

	    if (i1 == i2)
	      spda_out << 1. << " ";
	    else
	      spda_out << 0. << " ";


	  }
	  }
	 }


  // F3 Flag - P and p trace relation
	 if (flag_pass.F3_flag)
	 {
	  	for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
	  	{
	  	for (size_t k = i; k < bsize; k++)      // loop over jth constraint matrix
	  	{

	    	spda_out << 0. << " ";

	  	}
	  	}
	 }


  // F7 Flag - Q relations
	 if (flag_pass.F7_flag)
	 {
	 	  for (size_t i = 0;   i < bsize; i++)      // loop over ith constraint matrix
		  {
		  for (size_t j = i+1; j < bsize; j++)      // loop over jth constraint matrix
		  {
		  for (size_t k = i;   k < bsize; k++)    // loop over matrix row
		  {
		  for (size_t l = k+1; l < bsize; l++)    // loop over matrix column
		  {
		    if (j < l && k <= i)
		      continue;

			double val = kron_del(i,k)*kron_del(j,l) - kron_del(i,l)*kron_del(j,k);

		    spda_out << val << " ";

		  }
		  }
		  }
		  }

	 }


	 // F10 Flag - G relations
	 if (flag_pass.F10_flag)
	 {
		  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
		  {
		  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
		  {
		  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
		  {
		  for (size_t l = j; l < bsize; l++)    // loop over matrix column
		  {
		    
		    if (j == l && k < i)
		      continue;

		    spda_out << 0. << " ";

		  }
		  }
		  }
		  }

	 }


  spda_out << std::endl;

}


void init_F1_flag (const con_flags flag_pass, std::ofstream & spda_out, const size_t bsize, size_t & con_count)
{


    for (size_t ip = 0;  ip < bsize; ip++)
    {
    for (size_t kp = ip; kp < bsize; kp++)
    {


      double val1 = kron_del (ip, kp);

      size_t n = ip + 1;
      size_t m = kp + 1; 

      if (val1 != 0.)
      	spda_out << con_count << " " << 1 << " " << n << " " << m << " " << val1 << "\n";


    }
    }


    con_count++;
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

void init_F2_flag (const con_flags flag_pass, std::ofstream & spda_out, const size_t bsize, size_t & con_count)
{


//  std::stringstream ss1;
//  std::stringstream ss2;
//  std::stringstream ss3;

  std::ios_base::sync_with_stdio(false);

//  std::string ss1;
//  std::string ss2;
//  std::string ss3;

//  ss1.reserve(bsize*bsize*32);
//  ss2.reserve(bsize*bsize*32);
//  ss3.reserve(bsize*bsize*bsize*bsize*32);  

//  spda_out << std::setprecision(4);

//  std::cout << "BEFORE" << std::endl;

//  size_t extent = 2*bsize*bsize + (bsize*bsize-1)*(bsize*bsize-1);

//  extent *= bsize*bsize;

//  char * buffer = new char [extent*10];

//  char * buffera = new char [bsize*bsize*10];
//  char * bufferb = new char [bsize*bsize*10];
//  char * bufferc = new char [bsize*bsize*bsize*bsize*10];

//  char* ss1 = new char [10000];
//  char* ss2 = new char [10000];
//  char* ss3 = new char [10000];

//  std::cout << "AFTER" << std::endl;

//  int j1, j2, j3;

//  j1 = sprintf (ss1, "%f", 0.);
//  j2 = sprintf (ss2, "%f", 0.);
//  j3 = sprintf (ss3, "%f", 0.);

//  std::cout << "LOOP START" << std::endl;


//  int track = 0;
//  int num = 0;

  for (size_t i = 0; i < bsize; i++)
  {
  for (size_t j = i; j < bsize; j++)
  {

//	std::cout << "TEST";

//  	int testa = 0, testb = 0, testc = 0;

//  	int a = 0, b = 0, c = 0;

    for (size_t k = 0; k < bsize; k++)
    {
    for (size_t l = k; l < bsize; l++)
    {

//    	std::cout << j1 << " ";

      double val1 = (1./2.)*(kron_del(i,k)*kron_del(j,l) + kron_del(i,l)*kron_del(j,k));

      size_t n = k + 1;
      size_t m = l + 1; 

      if (val1 != 0.)
      	spda_out << con_count << " " << 1 << " " << n << " " << m << " " << val1 << "\n";



//      track += sprintf (buffer+track, "%f", val1);
//      track += sprintf (buffer+track, " ");

//      num += 9;

//      testa += sprintf (buffera+testa, "%f", val1);

//      testa += sprintf (buffera+testa, " ");
//      spda_out.write (buffera+testa-9,9);

//      a++;

//      ss1.append(val1);
//      spda_out.write("%f ", val1);
//      spda_out << val1 << " ";
//      ss1 << val1 << " ";
//      	j1 += sprintf (ss1+j1, "%f", val1);
//      	j2 += sprintf (ss2+j2, "%f", val1);

/*
      ss2 << val1 << " ";

      if (flag_pass.F3_flag and k < l)
      {
      	float val3 = 0.;

	    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
	    {
	    for (size_t lp = kp; lp < bsize; lp++)    // loop over matrix column
	    {
	    	if (lp > kp)
	    		ss3 << val3 << " ";

	    }
		}


      }

	*/
    }
    }



    for (size_t k = 0; k < bsize; k++)
    {
    for (size_t l = k; l < bsize; l++)
    {


      double val2 = (1./2.)*(kron_del(i,k)*kron_del(j,l) + kron_del(i,l)*kron_del(j,k));

      size_t n = k + 1;
      size_t m = l + 1; 

      if (val2 != 0.)
      	spda_out << con_count << " " << 2 << " " << n << " " << m << " " << val2 << "\n";

//      testb += sprintf (bufferb+testb, "%f", val2);

//      testb += sprintf (bufferb+testb, " ");

//      spda_out.write (bufferb+testb-9,9);

//      b++;

//      track += sprintf (buffer+track, "%f", val2);
//      track += sprintf (buffer+track, " ");

//      num += 9;

//      spda_out << val2 << " ";

//      spda_out.write("%f ", val2);
    }
    }



//	track += sprintf (buffer+track, "\n");


//    num++;
//	spda_out.write(buffera, 9*a);
//	spda_out.write(bufferb, 9*b);
//	spda_out.write(bufferc, 9*c);
//	spda_out << "\n";
//  spda_out.write(ss1.str().c_str(), ss1.str().length());
//  spda_out.write(ss2.str().c_str(), ss2.str().length());
//  spda_out.write(ss3.str().c_str(), ss3.str().length());
//	fprintf(spda_out, "%s", ss1);
//	fprintf(spda_out, "%s", ss2);
//	fprintf(spda_out, "%s", ss3);
//  spda_out << "\n";

//  ss1.str("");
//  ss1.clear();
//  ss2.str("");
//  ss2.clear();  
//  ss3.str("");
//  ss3.clear();

//  spda_out << ss;

    con_count++;

  }
  }

//  spda_out.write(buffer, num);

}





int F3_3_matrix (const int i, const int k, const int ip, const int jp, const int kp, const int lp)
{
  const int value = kron_del (i,ip) * kron_del (k,kp) * kron_del (jp, lp);

  return value;  
}

int F3_3_matrix_A (const int i, const int k, const int ip, const int jp, const int kp, const int lp)
{
	const int value = 
	    F3_3_matrix (i, k, ip, jp, kp, lp) +  F3_3_matrix (i, k, kp, lp, ip, jp)
      - F3_3_matrix (i, k, jp, ip, kp, lp) -  F3_3_matrix (i, k, kp, lp, jp, ip)
      - F3_3_matrix (i, k, ip, jp, lp, kp) -  F3_3_matrix (i, k, lp, kp, ip, jp)
      + F3_3_matrix (i, k, jp, ip, lp, kp) +  F3_3_matrix (i, k, lp, kp, jp, ip);

	return value;
}

/***************************************************************

Function to create our F3 constraint matrix to fix the relation
between the 2RDM partial trace and the 1RDM.

***************************************************************/


void init_F3_flag (const con_flags flag_pass, std::ofstream & spda_out, const size_t bsize, const size_t N, size_t & con_count)
{

  std::ios_base::sync_with_stdio(false);


  for (size_t i = 0; i < bsize; i++)
  {
  for (size_t k = i; k < bsize; k++)
  {



    for (size_t ip = 0; ip < bsize; ip++)
    {
    for (size_t kp = ip; kp < bsize; kp++)
    {


      double val1 = 
      -1.0 * (N - 1.0) / 4.0 * (
        kron_del(i,ip)*kron_del(k,kp) + kron_del(k,ip)*kron_del(i,kp)
        );

      size_t n = ip + 1;
      size_t m = kp + 1; 

      if (val1 != 0.)
      	spda_out << con_count << " " << 1 << " " << n << " " << m << " " << val1 << "\n";



    }
    }



    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
      if (ip >= jp)
        continue;

    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
      if (kp >= lp)
        continue;

      size_t ips = ip + 1;
      size_t jps = jp + 1;
      size_t kps = kp + 1;
      size_t lps = lp + 1;

      size_t n = jps - ips + (2*bsize - ips) * (ips - 1)/2;
      size_t m = lps - kps + (2*bsize - kps) * (kps - 1)/2;

      double val3 = 1./8. * F3_3_matrix_A (i, k, ip, jp, kp, lp);


      if (val3 != 0. and n <= m)
      	spda_out << con_count << " " << 3 << " " << n << " " << m << " " << val3 << "\n";


    }
    }
    }
    }

    con_count++;

  }
  }



}

/***************************************************************

Function to create our F4 constraint matrix to fix the 2RDM trace.

***************************************************************/

template <typename one_array, typename three_array>
void init_F4_flag (three_array & F4_con, one_array & F4_val, const size_t bsize, const size_t particles)
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

        size_t left_P  = left  + 2*bsize;
        size_t right_P = right + 2*bsize;

        F4_con [0][left_P][right_P] = kron_del(i,k) * kron_del(j,l);

    }
    }

  }
  }

  F4_val[0] = particles * (particles - 1)/2;

}





int F7_2_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp)
{
  const int value = 
        kron_del(j,l)*kron_del(i,ip)*kron_del(k,jp)
        -
        kron_del(j,k)*kron_del(i,ip)*kron_del(l,jp)
        -
        kron_del(i,l)*kron_del(j,ip)*kron_del(k,jp)
        +
        kron_del(i,k)*kron_del(j,ip)*kron_del(l,jp);

  return value;

}

int F7_3_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = kron_del (i,ip) * kron_del (j,jp) * kron_del (k,kp) * kron_del (l,lp);

  return value;  
}

int F7_4_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = kron_del (l,ip) * kron_del (k,jp) * kron_del (j,kp) * kron_del (i,lp);

  return value;  
}

int F7_3_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
    const int value =
	    F7_3_matrix (i, j, k, l, ip, jp, kp, lp) +  F7_3_matrix (i, j, k, l, kp, lp, ip, jp)
  	- F7_3_matrix (i, j, k, l, jp, ip, kp, lp) -  F7_3_matrix (i, j, k, l, kp, lp, jp, ip)
    - F7_3_matrix (i, j, k, l, ip, jp, lp, kp) -  F7_3_matrix (i, j, k, l, lp, kp, ip, jp)
    + F7_3_matrix (i, j, k, l, jp, ip, lp, kp) +  F7_3_matrix (i, j, k, l, lp, kp, jp, ip);

	return value;
}

int F7_4_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = 
        F7_4_matrix (i, j, k, l, ip, jp, kp, lp) +  F7_4_matrix (i, j, k, l, kp, lp, ip, jp)
      - F7_4_matrix (i, j, k, l, jp, ip, kp, lp) -  F7_4_matrix (i, j, k, l, kp, lp, jp, ip)
      - F7_4_matrix (i, j, k, l, ip, jp, lp, kp) -  F7_4_matrix (i, j, k, l, lp, kp, ip, jp)
      + F7_4_matrix (i, j, k, l, jp, ip, lp, kp) +  F7_4_matrix (i, j, k, l, lp, kp, jp, ip);

  return value;  
}

/***************************************************************

Function to create our F7 constraint matrix to fix the linear 
relations between the 1RDM, q, 2RDM, and Q. 

***************************************************************/

void init_F7_flag (const con_flags flag_pass, std::ofstream & spda_out, const size_t bsize)
{


  for (size_t i = 0;   i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = i+1; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = i;   k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = k+1; l < bsize; l++)    // loop over matrix column
  {
    if (j < l && k <= i)
      continue;

    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    	double val1 = 0.;
    	spda_out << val1 << " ";
    }
	}


    
    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    	double val2 = 
    	(1./2.) * (F7_2_matrix (i, j, k, l, ip, jp) + F7_2_matrix (i, j, k, l, jp, ip));
    	spda_out << val2 << " ";
    }
	}



    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
      if (ip >= jp)
        continue;

    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
      if (kp >= lp)
        continue;

    	double val3 = (1./4.) * F7_3_matrix_A (i, j, k, l, ip, jp, kp, lp);

    	spda_out << val3 << " ";
	}
	}
	}
	}


    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
      if (ip >= jp)
        continue;

    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
      if (kp >= lp)
        continue;

    	double val4 = (-1./8.) * F7_4_matrix_A (i, j, k, l, ip, jp, kp, lp);

    	spda_out << val4 << " ";
	}
	}
	}
	}


	if (flag_pass.F10_flag)
	{
		for (size_t ip = 0;  ip < bsize; ip++)      // loop over ith constraint matrix
		{
		for (size_t jp = 0;  jp < bsize; jp++)      // loop over jth constraint matrix
		{
		for (size_t kp = 0;  kp < bsize; kp++)    // loop over matrix row
		{
		for (size_t lp = jp; lp < bsize; lp++)    // loop over matrix column
		{
		    
		    if (jp == lp && kp < ip)
		      continue;

		    spda_out << 0. << " ";

		}
		}
		}
		}

	}


  	spda_out << std::endl;

  }	
  }
  }
  }



}






int F10_1_matrix (const int i, const int j, const int k, const int l, const int ip, const int kp)
{
  const int value = kron_del(j,l) * kron_del(i,ip)*kron_del(k,kp);

  return value;

}

int F10_3_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = 
        kron_del(i,ip) * kron_del(l,jp) * kron_del(k,kp) * kron_del(j,lp);

  return value;
}

int F10_5_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = 
        kron_del(i,ip) * kron_del(j,jp) * kron_del(k,kp) * kron_del(l,lp);

  return value;
}


int F10_3_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
    const int value =
      F10_3_matrix (i, j, k, l, ip, jp, kp, lp) +  F10_3_matrix (i, j, k, l, kp, lp, ip, jp)
    - F10_3_matrix (i, j, k, l, jp, ip, kp, lp) -  F10_3_matrix (i, j, k, l, kp, lp, jp, ip)
    - F10_3_matrix (i, j, k, l, ip, jp, lp, kp) -  F10_3_matrix (i, j, k, l, lp, kp, ip, jp)
    + F10_3_matrix (i, j, k, l, jp, ip, lp, kp) +  F10_3_matrix (i, j, k, l, lp, kp, jp, ip);

  return value;
}

/***************************************************************

Function to create our F10 constraint matrix to fix the G linear
relations.

***************************************************************/

template <typename one_array, typename three_array, typename six_array, typename eight_array>
void init_F10_flag (six_array & F10_build_1, eight_array & F10_build_3, eight_array & F10_build_5, three_array & F10_con, one_array & F10_val, const size_t bsize)
{



  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = j; l < bsize; l++)    // loop over matrix column
  {
    
    if (j == l && k < i)
      continue;

    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t kp = 0; kp < bsize; kp++)      // loop over jth constraint matrix
    {

      F10_build_1 [i][j][k][l][ip][kp] =  (-1./2.) * (
        F10_1_matrix (i, j, k, l, ip, kp) + F10_1_matrix (i, j, k, l, kp, ip)
        );
  
    }
    }
  }
  }
  }
  }


  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = j; l < bsize; l++)    // loop over matrix column
  {
    
    if (j == l && k < i)
      continue;

    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
      F10_build_3 [i][j][k][l][ip][jp][kp][lp] =  (1./4.) * F10_3_matrix_A (i, j, k, l, ip, jp, kp, lp);

      F10_build_5 [i][j][k][l][ip][jp][kp][lp] =  (1./2.) * (
        F10_5_matrix (i, j, k, l, ip, jp, kp, lp) + F10_5_matrix (i, j, k, l, kp, lp, ip, jp)
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
  for (size_t l = j; l < bsize; l++)    // loop over matrix column
  {

    if (j == l && k < i)
      continue;

    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {

      size_t b = counter;

      size_t left_p  = ip;
      size_t right_p = jp;

      F10_con[b][left_p][right_p] = F10_build_1 [i][j][k][l][ip][jp];

    }
    }

    counter++;
  }
  }
  }
  }




  counter = 0;


  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = j; l < bsize; l++)    // loop over matrix column
  {
    
    if (j == l && k < i)
      continue;

    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
      size_t b = counter;

      size_t ips = ip + 1;
      size_t jps = jp + 1;
      size_t kps = kp + 1;
      size_t lps = lp + 1;

      size_t left  = jps - ips + (2*bsize - ips) * (ips - 1)/2 - 1;
      size_t right = lps - kps + (2*bsize - kps) * (kps - 1)/2 - 1;

      size_t left_P  = left  + 2*bsize;
      size_t right_P = right + 2*bsize;

      size_t lG = ip * bsize + jp;
      size_t rG = kp * bsize + lp;

      size_t left_G  = lG + 2*bsize + 2*bsize*(bsize-1)/2;
      size_t right_G = rG + 2*bsize + 2*bsize*(bsize-1)/2;

      if (jp > ip and lp > kp)
        F10_con[b][left_P][right_P] = F10_build_3 [i][j][k][l][ip][jp][kp][lp];
      

      F10_con[b][left_G][right_G] = F10_build_5 [i][j][k][l][ip][jp][kp][lp];











    }
    }
    }
    }

    F10_val [counter] = 0.;

//    std::cout << "i j k l" << "\t" << i << " " << j << " " << k << " " << l << std::endl;
//    print (std::cout, F10_con[counter]);
//    std::cout << F10_val [counter] << std::endl;

    counter++;
  }
  }
  }
  }

  std::cout << "COUNTER = " << counter << std::endl;

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

  std::cout << ("Building system... ") << std::endl;



  size_t Q_num = 0;

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = i+1; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = i; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = k+1; l < bsize; l++)    // loop over matrix column
  {
    if (j < l && k <= i)
      continue;

    Q_num++;
  }
  }
  }
  }

  size_t G_num = 0;

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = j; l < bsize; l++)    // loop over matrix column
  {
    
    if (j == l && k < i)
      continue;

    G_num++;
  }
  }
  }
  }

  // define flags for all different combinations of conditions in RDM
  const bool two_body_toggle = true;
  const bool redundant_check = false;
  const bool diag_toggle = false;

  // START CONSTRAINT FLAG DEFINE

  const bool F1_flag = true;  // p START - TRACE CONDITION
  const bool F2_flag = true;  // q START - LINEAR RELATIONS
  const bool F3_flag = true;  // P START - TRACE CONDITION
  const bool F4_flag = false; // REDUNDANT - P TRACE CONDITION
  const bool F5_flag = false;  // P ANTI-SYMMETRY
  const bool F6_flag = false; // REDUNDANT - P ANTI-SYMMETRY
  const bool F7_flag = false;  // Q START - LINEAR RELATIONS
  const bool F8_flag = false; // Q ANTI-SYMMETRY
  const bool F9_flag = false; // REDUNDANT - Q ANTI-SYMMETRY
  const bool F10_flag = false; // G START - LINEAR REALTIONS

  const bool Q_flag = false;
  const bool G_flag = false;



  if (!two_body_toggle and (F4_flag or F5_flag or F6_flag or F7_flag))
  {
  	std::cerr << "ERROR: TWO-BODY TOGGLE AND F CONSTRAINT FLAGS INCOMPATIBLE - GIVING UP" << std::endl;
  	return EXIT_FAILURE;
  }	

  if (!two_body_toggle and (Q_flag or G_flag))
  {
    std::cerr << "ERROR: TWO-BODY TOGGLE AND Q/G FLAGS INCOMPATIBLE - GIVING UP" << std::endl;
      return EXIT_FAILURE;
  }

  if ((F7_flag and !Q_flag) or (!F7_flag and Q_flag))
  {
    std::cerr << "ERROR: INCOMPATIBLE SETTINGS ON F7 FLAG AND Q FLAG - GIVING UP" << std::endl;
      return EXIT_FAILURE;
  }

  if (G_flag and !Q_flag)
  {
    std::cerr << "ERROR: CANNOT INIT G FLAG WITHOUT Q FLAG - GIVING UP" << std::endl;
    return EXIT_FAILURE;
  }




  struct con_flags flag_pass;

  flag_pass.F1_flag = F1_flag;
  flag_pass.F2_flag = F2_flag;
  flag_pass.F3_flag = F3_flag;
  flag_pass.F4_flag = F4_flag;
  flag_pass.F5_flag = F5_flag;
  flag_pass.F6_flag = F6_flag;
  flag_pass.F7_flag = F7_flag;
  flag_pass.F8_flag = F8_flag;
  flag_pass.F9_flag = F9_flag;
  flag_pass.F10_flag = F10_flag;

  flag_pass.two_body_toggle = two_body_toggle;
  flag_pass.redundant_check = redundant_check;

  flag_pass.Q_flag = Q_flag;
  flag_pass.G_flag = G_flag;

  flag_pass.diag_toggle = diag_toggle;



  const size_t F1num  = 1;
//  const size_t F2num = bsize*bsize;
  const size_t F2num  = bsize * (bsize + 1)/2;
  const size_t F3num  = bsize * (bsize + 1)/2;
//  const size_t F3num = bsize*bsize;
  const size_t F4num  = 1;
  const size_t F5num  = 0;//anti_num;//bsize*bsize*bsize*bsize;
  const size_t F6num  = 0;//index_num;//bsize*bsize*bsize*bsize;
  const size_t F7num  = Q_num;//bsize*bsize*bsize*bsize;
  const size_t F8num  = 0;//anti_num;//bsize*bsize*bsize*bsize;
  const size_t F9num  = 0;//index_num;//bsize*bsize*bsize*bsize;
  const size_t F10num = G_num;//bsize*bsize*bsize*bsize;
  // turn off or on two-body potential




  const std::string diag_file = "diagnostic_out/test_diag.dat";
  const std::string spda_file = "sdp_files/test_spd.dat-s";

  std::ofstream diag_out (diag_file);
  std::ofstream spda_out (spda_file);

//  double weight[nodes], x[nodes];							// create arrays to hold gaussian weights/positions
//  gauss (mesh_size, 2, 0., 4., x, weight);						// creating weights/positions for s-wave G-L quad.
/*  for (int i = 0; i < nodes; i++)
  {
	printf ("%f", x[i]);
	printf ("\n");
  }*/

  size_t cons = 0;
  size_t blocks = 0;

  if (F1_flag)
  {
  	cons += F1num;
  	blocks++;
  }

  if (F2_flag)
  {
  	cons += F2num;
  	blocks++;
  }
  
  if (F3_flag)
  {
  	cons += F3num;
  	blocks++;
  }

  if (F7_flag)
  {
  	cons += F7num;
  	blocks++;
  }

  if (F10_flag)
  {
  	cons += F10num;
  	blocks++;
  }




  two_array ref_m (boost::extents[bsize][7]);
  two_array h1_mat(boost::extents[bsize][bsize]);
  five_array h2_mat(boost::extents[bsize][bsize][bsize][bsize][5]);

  fullm_populate_hamiltonian (ref_m, h1_mat, h2_mat, m_ref, m_mat, hw, diag_out, diag_toggle);

  std::cout << "HAMILTONIAN BUILT" << std::endl;

  two_array comp_h2 (boost::extents[bsize*(bsize-1)/2][bsize*(bsize-1)/2]);

  compactify_h2 (ref_m, comp_h2, h2_mat, diag_out, diag_toggle);

  std::cout << "POTENTIAL COMPACTIFIED" << std::endl;

  spda_out << cons << std::endl;
  spda_out << blocks << std::endl;

  if (F1_flag)
  {
  	spda_out << bsize << " ";
  }

  if (F2_flag)
  {
  	spda_out << bsize << " ";
  }
  
  if (F3_flag)
  {
  	spda_out << bsize * (bsize-1)/2 << " ";
  }

  if (F7_flag)
  {
  	spda_out << bsize * (bsize-1)/2 << " ";
  }

  if (F10_flag)
  {
  	spda_out << bsize * bsize << " ";
  }

  spda_out << std::endl;

  
  init_con_values (flag_pass, spda_out, bsize, particles);

  size_t con_count = 0;


  init_C_matrix (flag_pass, spda_out, h1_mat, comp_h2, con_count);

  std::cout << "C MATRIX DONE" << std::endl;


  if (F1_flag)
  {
    init_F1_flag (flag_pass, spda_out, bsize, con_count);
    std::cout << "FLAG 1 DONE" << std::endl;
  }

  if (F2_flag)
  {
    init_F2_flag (flag_pass, spda_out, bsize, con_count);
    std::cout << "FLAG 2 DONE" << std::endl;
  }

  if (F3_flag)
  {
    init_F3_flag (flag_pass, spda_out, bsize, particles, con_count);
    std::cout << "FLAG 3 DONE" << std::endl;
  }

  if (F7_flag)
  {
    init_F7_flag (flag_pass, spda_out, bsize);
    std::cout << "FLAG 7 DONE" << std::endl;
  }

  if (F10_flag)
  {
//    init_F10_flag (flag_pass, spda_out, bsize);
    std::cout << "FLAG 10 DONE" << std::endl;
  }


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


