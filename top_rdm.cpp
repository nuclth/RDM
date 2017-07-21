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
	size_t basis;
	size_t particles;
	double hw;
    std::string m_ref;
    std::string m_mat;
};

struct con_flags
{
	bool N_flag;
	bool O_flag;
	bool P_flag;
  	bool Q_flag;
  	bool G_flag;
	bool diag_toggle;
	bool two_body_toggle;
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
typedef boost::multi_array<double, 7> seven_array;
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


 Function to print out an arbitrary multiarray


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

//  ss >> orbit_dummy_1 >> orbit_dummy_2 >> ref_num >> n >> l >> j >> m_j >> tz;  // assign values of the line

  ss >> ref_num >> n >> l >> j >> m_j >> tz;  // assign values of the line

//  std::cout << ref_num << " " << n << " " << l << " " << j << " " << m_j << " " << tz << "\n";


  // only extract neutron-neutron states
  if (tz == -1.)
  {
    ref_m [ele_in][0] = ref_num;    // reference number of the line
    ref_m [ele_in][1] = n;          // principle quantum number
    ref_m [ele_in][2] = l;          // orbital angular mom.
    ref_m [ele_in][3] = j * 0.5;    // total angular mom.
    ref_m [ele_in][4] = m_j * 0.5;  // total angular mom. projection
    ref_m [ele_in][5] = tz;         // isospin projection (should all be -1.0)
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

Kronecker delta function. Returns 1 if i=j and 0 otherwise.

***************************************************************/

double kron_del(const size_t i, const size_t j)
{

  if (i == j)
    return 1.;

  return 0.;
}

/***************************************************************

Function to create the constraint matrix values for the SDP solver.

These are the real numbers on the right hand side of the SDP.

***************************************************************/



void init_con_values (const con_flags flag_pass, std::ofstream & spda_out, const size_t bsize, const size_t particles)
{

//  size_t cmat_extent = F1_con.size();


  // F1 Flag - N Trace condition
	if (flag_pass.N_flag)
	  spda_out << particles << " ";

  // F2 Flag - Linear p q relations
	if (flag_pass.O_flag)
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
	 if (flag_pass.P_flag)
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
	 if (flag_pass.Q_flag)
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
	 if (flag_pass.G_flag)
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

    spda_out << std::setprecision(16);

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





	if (flag_pass.P_flag)
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

Function to create our F1 constraint matrix to fix the total
particle number N (fixes trace of the 1RDM to be N). 

In keeping with other constraint functions, first populate the
F1_con_1 matrix (only one constraint) then populate the full F1_con
matrix. 

Note that F1_con_1 is a basis x basis size matrix while F1_con 
matches the dimensions of the F0 constraint matrix.

***************************************************************/

void init_N_flag (std::ofstream & spda_out, const size_t bsize, size_t & con_count)
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

Function to create our O constraint matrices to fix the relation
between the 1RDM and its twin q. 

First populate the F2_con_1 matrix which serves as the basis x basis
size constraint matrix for the 1RDM sector and the q sector (both 
have the same constraint matrix).

Then populate this matrix into the full F0 sized constraint term 
F2_con. 

***************************************************************/

void init_O_flag (std::ofstream & spda_out, const size_t bsize, size_t & con_count)
{

  std::ios_base::sync_with_stdio(false);


  for (size_t i = 0; i < bsize; i++)
  {
  for (size_t j = i; j < bsize; j++)
  {

    for (size_t k = 0; k < bsize; k++)
    {
    for (size_t l = k; l < bsize; l++)
    {


      double val1 = (1./2.)*(kron_del(i,k)*kron_del(j,l) + kron_del(i,l)*kron_del(j,k));

      size_t n = k + 1;
      size_t m = l + 1; 

      if (val1 != 0.)
      	spda_out << con_count << " " << 1 << " " << n << " " << m << " " << val1 << "\n";


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


    }
    }



    con_count++;

  }
  }


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

Function to create our P constraint matrix to fix the relation
between the 2RDM partial trace and the 1RDM.

***************************************************************/


void init_P_flag (std::ofstream & spda_out, const size_t bsize, const size_t N, size_t & con_count)
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

void init_Q_flag (std::ofstream & spda_out, const size_t bsize, size_t con_count)
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

    
    for (size_t ip = 0;  ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = ip; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    	double val2 = 
    	(1./2.) * (F7_2_matrix (i, j, k, l, ip, jp) + F7_2_matrix (i, j, k, l, jp, ip));

    	size_t n = ip + 1;
        size_t m = jp + 1; 

        if (val2 != 0.)
      		spda_out << con_count << " " << 2 << " " << n << " " << m << " " << val2 << "\n";
    }
	}

    for (size_t ip = 0;    ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = ip+1; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0;   kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = kp+1; lp < bsize; lp++)    // loop over matrix column
    {
  //  	if (jp < lp && kp <= ip)
  //    		continue;


    	double val3 = (1./4.) * F7_3_matrix_A (i, j, k, l, ip, jp, kp, lp);


        size_t ips = ip + 1;
        size_t jps = jp + 1;
        size_t kps = kp + 1;
        size_t lps = lp + 1;

        size_t n = jps - ips + (2*bsize - ips) * (ips - 1)/2;
        size_t m = lps - kps + (2*bsize - kps) * (kps - 1)/2;


        if (val3 != 0. and n <= m)
   		   	spda_out << con_count << " " << 3 << " " << n << " " << m << " " << val3 << "\n";


	}
	}
	}
	}


    for (size_t ip = 0;    ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = ip+1; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0;   kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = kp+1; lp < bsize; lp++)    // loop over matrix column
    {
//    	if (jp < lp && kp <= ip)
//      		continue;

    	double val4 = (-1./8.) * F7_4_matrix_A (i, j, k, l, ip, jp, kp, lp);

    	size_t ips = ip + 1;
        size_t jps = jp + 1;
        size_t kps = kp + 1;
        size_t lps = lp + 1;

        size_t n = jps - ips + (2*bsize - ips) * (ips - 1)/2;
        size_t m = lps - kps + (2*bsize - kps) * (kps - 1)/2;


        if (val4 != 0. and n <= m)
   		   	spda_out << con_count << " " << 4 << " " << n << " " << m << " " << val4 << "\n";


	}
	}
	}
	}



  con_count++;

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

void init_G_flag (std::ofstream & spda_out, const size_t bsize, size_t con_count)
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

    for (size_t ip = 0;  ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t kp = ip; kp < bsize; kp++)      // loop over jth constraint matrix
    {


      double val1 =  
      (-1./2.) * (F10_1_matrix (i, j, k, l, ip, kp) + F10_1_matrix (i, j, k, l, kp, ip));

      size_t n = ip + 1;
      size_t m = kp + 1; 

      if (val1 != 0.)
      	spda_out << con_count << " " << 1 << " " << n << " " << m << " " << val1 << "\n";

  
    }
    }


    for (size_t ip = 0;    ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = ip+1; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0;   kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = kp+1; lp < bsize; lp++)    // loop over matrix column
    {
//    	if (jp < lp && kp <= ip)
//      		continue;
      	
      	double val3 = (1./4.) * F10_3_matrix_A (i, j, k, l, ip, jp, kp, lp);

        size_t ips = ip + 1;
        size_t jps = jp + 1;
        size_t kps = kp + 1;
        size_t lps = lp + 1;

        size_t n = jps - ips + (2*bsize - ips) * (ips - 1)/2;
        size_t m = lps - kps + (2*bsize - kps) * (kps - 1)/2;


        if (val3 != 0. and n <= m)
   		   	spda_out << con_count << " " << 3 << " " << n << " " << m << " " << val3 << "\n";

    }
    }
    }
	}



    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
//        if (jp == lp && kp < ip)
//    	  continue;

      	double val5 =  
      	(1./2.) * (F10_5_matrix (i, j, k, l, ip, jp, kp, lp) + F10_5_matrix (i, j, k, l, kp, lp, ip, jp));

        size_t n = ip * bsize + jp + 1;
        size_t m = kp * bsize + lp + 1;

        if (val5 != 0. and n <= m)
   		   	spda_out << con_count << " " << 5 << " " << n << " " << m << " " << val5 << "\n";    

    }
    }
    }
    }

    con_count++;

  }
  }
  }
  }

}
/*

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
*/


/********************************************

BEGIN MAIN PROGRAM

********************************************/

int main ()
{
  
  parameters input_params;						             // create a struct to hold user specified parameters
  input_params = read_in_inputs (); 					     // read in the values from inputs.inp, store in input_params	

  // Now create/copy input parameters for ease of reading the code into newly defined local variables



  const size_t bsize = input_params.basis;	   // the size of the HO basis used
  const size_t particles = input_params.particles;		 // the number of neutrons (particles) in the trap
  const double hw = input_params.hw;					     // hbar * omega 
  const std::string m_ref = input_params.m_ref;    // single particle reference file - m scheme
  const std::string m_mat = input_params.m_mat;    // m scheme matrix elements



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
  const bool diag_toggle = false;

  // START CONSTRAINT FLAG DEFINE

  const bool N_flag = true; // p START - TRACE CONDITION
  const bool O_flag = true; // q START - LINEAR RELATIONS
  const bool P_flag = true; // P START - TRACE CONDITION
  const bool Q_flag = false; // Q START - LINEAR RELATIONS
  const bool G_flag = false; // G START - LINEAR REALTIONS



  if (!two_body_toggle and (P_flag or Q_flag or G_flag))
  {
  	std::cerr << "ERROR: TWO-BODY TOGGLE AND CONSTRAINT FLAGS INCOMPATIBLE - GIVING UP" << std::endl;
  	return EXIT_FAILURE;
  }	

  if (!P_flag and (Q_flag or G_flag))
  {
    std::cerr << "ERROR: CANNOT INIT G OR Q FLAG WITHOUT P FLAG - GIVING UP" << std::endl;
    return EXIT_FAILURE;
  }


  if (G_flag and !Q_flag)
  {
    std::cerr << "ERROR: CANNOT INIT G FLAG WITHOUT Q FLAG - GIVING UP" << std::endl;
    return EXIT_FAILURE;
  }




  struct con_flags flag_pass;

  flag_pass.N_flag = N_flag;
  flag_pass.O_flag = O_flag;
  flag_pass.P_flag = P_flag;
  flag_pass.Q_flag = Q_flag;
  flag_pass.G_flag = G_flag;

  flag_pass.two_body_toggle = two_body_toggle;
  flag_pass.diag_toggle = diag_toggle;



  const size_t F1num  = 1;
  const size_t F2num  = bsize * (bsize + 1)/2;
  const size_t F3num  = bsize * (bsize + 1)/2;
  const size_t F7num  = Q_num;//bsize*bsize*bsize*bsize;
  const size_t F10num = G_num;//bsize*bsize*bsize*bsize;





  const std::string diag_file = "diagnostic_out/test_diag.dat";
  const std::string spda_file = "sdp_files/test_spd.dat-s";

  std::ofstream diag_out (diag_file);
  std::ofstream spda_out (spda_file);



  size_t cons = 0;
  size_t blocks = 0;

  if (N_flag)
  {
  	cons += F1num;
  	blocks++;
  }

  if (O_flag)
  {
  	cons += F2num;
  	blocks++;
  }
  
  if (P_flag)
  {
  	cons += F3num;
  	blocks++;
  }

  if (Q_flag)
  {
  	cons += F7num;
  	blocks++;
  }

  if (G_flag)
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

  if (P_flag)
  {
  	spda_out << bsize << " ";
  }

  if (O_flag)
  {
  	spda_out << bsize << " ";
  }
  
  if (P_flag)
  {
  	spda_out << bsize * (bsize-1)/2 << " ";
  }

  if (Q_flag)
  {
  	spda_out << bsize * (bsize-1)/2 << " ";
  }

  if (G_flag)
  {
  	spda_out << bsize * bsize << " ";
  }

  spda_out << std::endl;

  
  init_con_values (flag_pass, spda_out, bsize, particles);

  size_t con_count = 0;


  init_C_matrix (flag_pass, spda_out, h1_mat, comp_h2, con_count);

  std::cout << "C MATRIX DONE" << std::endl;


  if (N_flag)
  {
    init_N_flag (spda_out, bsize, con_count);
    std::cout << "N FLAG DONE" << std::endl;
  }

  if (O_flag)
  {
    init_O_flag (spda_out, bsize, con_count);
    std::cout << "O FLAG DONE" << std::endl;
  }

  if (P_flag)
  {
    init_P_flag (spda_out, bsize, particles, con_count);
    std::cout << "P FLAG DONE" << std::endl;
  }

  if (Q_flag)
  {
    init_Q_flag (spda_out, bsize, con_count);
    std::cout << "Q FLAG DONE" << std::endl;
  }

  if (G_flag)
  {
    init_G_flag (spda_out, bsize, con_count);
    std::cout << "G FLAG DONE" << std::endl;
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
	  ss >> input.basis;

	if (counter == 1)
	  ss >> input.particles;

	if (counter == 2)
	  ss >> input.hw;

	if (counter == 3)
	  ss >> input.m_ref;

  	if (counter == 4)
      ss >> input.m_mat;


	counter++;					// increase the counter after reading in a parameter	
  }

  return input;						// return the struct 
}


