/********************************************

File to hold the necessary function definitions for our derived
fullm class. 

********************************************/


#include "phys_system.h"
#include "fullm.h"
#include <fstream>
#include <string>
#include <sstream>

/*******************************************

Constructor for the system. Read in the reference for the single particle
states then read in the actual matrix elements.

Throw exception if the number of particles not equal to 2, 6, or 8. 

*******************************************/

fullm::fullm (std::string choice, std::string reference_file, std::string matrix_file) : phys_system (choice, reference_file, matrix_file)
{

  try
  {
    nmax = read_in_reference_m_scheme (reference_file);
    read_in_matrix_m_scheme (nmax, matrix_file);
  }

  catch (const char * msg)
  {
    std::cerr << msg << std::endl;	
  }
 
}

/***************************************************

Function to access the relevant parts of the twobody potential. 
Called in the Hartree Fock solver. 

NOTE: FOR THE FULL CHOICE, THIS IS ONLY WELL-DEFINED FOR CLOSED-SHELL
SYSTEMS.

***************************************************/

double fullm::pot_access_modifier (int a, int b, int c, int d, int particles)
{

  double potential_factor = 0;

  // extract out the relevant single particle numbers for the 4 input numbers (v_{abcd})
  int refer_a = ref_m(a,0);
  double l_a = ref_m(a,2);
  double j_a = ref_m(a,3);
  double m_a = ref_m(a,4);
	
  int refer_b = ref_m(b,0);
  double l_b = ref_m(b,2);
  double j_b = ref_m(b,3);
  double m_b = ref_m(b,4);	

  int refer_c = ref_m(c,0);
  double l_c = ref_m(c,2);
  double j_c = ref_m(c,3);
  double m_c = ref_m(c,4);

  int refer_d = ref_m(d,0);
  double l_d = ref_m(d,2);
  double j_d = ref_m(d,3);
  double m_d = ref_m(d,4);


  // orthogonality condition check on l,j & l',j'
  // (should be superfluous check if matrix elements are right)
  if (l_a != l_c || j_a != j_c || l_b != l_d || j_b != j_d || m_a != m_c || m_b != m_d)
    return 0.;

  // perform check to see if matrix element is even non-zero 
  // due to l' and j' only being summed over occupied states

  double lmax = 0.;
  double jmax = 0.;

  if (particles == 2)
  {
    lmax = 0.;
    jmax = 0.5;
  }
	
  else if (particles == 8)
  {
    lmax = 1.;
    jmax = 1.5;
  }

  else if (particles == 20)
  {
    lmax = 2.;
    jmax = 2.5;
  }

  else
	throw "ERROR: Not a valid number of particles!";

  // check for l,j occupation in the system
  if (l_a > lmax || j_a > jmax || l_c > lmax || j_c > jmax)
    return 0.;	

  if (m_a == m_c && m_d == m_b)
    potential_factor =  m_scheme_matrix (refer_a, refer_b) (refer_c, refer_d);

  return potential_factor; 
}

/***********************************************

Function to output the kinetic energy and HO 1-body potential
energy contributions. E = (2n + l + 3/2) * hbar * omega where
n is the principle quantum number and l is the orbital angular 
momentum. This term is diagonal in the principle quantum number n
e.g. <n1|H_0|n2> = E*\delta_{n1,n2}.

Called in the hartree fock solver. 

***********************************************/

double fullm::kinetic_access (int a, int b, double hw)
{   

  double onebody_fact = 0.;

  // initialization of variables (garbage values here, values for the sake of definite declaration)
  double n1 = 100.; 
  double l1 = 100.;
  double j1 = 100.;
  double m1 = 100.;
  double iso1 = 100.;
  double n3 = 0.;
  double l3 = 0.;
  double j3 = 0.;
  double m3 = 0.;
  double iso3 = 0.;

  // loop over reference matrix: take out n, orbital ang. mom., and isospin projection for references
  // a & b --> check to see if they're all equivalent: if yes return KE, if not return 0
  for (int i = 0; i < nmax; i++)
  {

    if (ref_m(i,0) == (double)a)
    {
  	n1 = ref_m(i,1);
	l1 = ref_m(i,2);
	j1 = ref_m(i,3);
	m1 = ref_m(i,4);
	iso1 = ref_m(i,5);
    }

    if (ref_m(i,0) == (double)b)
    {
  	n3 = ref_m(i,1);
	l3 = ref_m(i,2);
	j3 = ref_m(i,3);
	m3 = ref_m(i,4);
	iso3 = ref_m(i,5);
    }
	  
  } // end for loop

  // check for equivalence of associated quantum numbers
  if (n1 == n3 && iso1 == iso3 && l1 == l3 && j1 == j3 && m1 == m3 && iso1 == 1.)
	onebody_fact = (2.*n1 + l1 + 1.5) * hw;

  return onebody_fact;
  

}

/**********************************************

Function to read in the reference for m scheme matrix elements 
calculated from Morten's code. 

If unexpected errors pop up, check the format on Morten's code.

THIS WILL BREAK IF THE INPUT FILE .dat CHANGES TITLE OR FORMAT.

**********************************************/

int fullm::read_in_reference_m_scheme (std::string m_ref_file)
{
  const char * m_reference_file = (m_ref_file).c_str();
  // input file stream for m_scheme
  std::ifstream m_ref_in (m_reference_file);
 
  if (m_ref_in.fail())
	throw "ERROR: Cannot open m-scheme reference file";
 
  std::size_t total_lines = 0;
  std::string dummy;
  double ref_num, n, l , j , m_j, tz;

  // find total number of defined reference lines
  while (std::getline (m_ref_in, dummy))
	++total_lines;

  // Matrix to hold reference values. Divide by 2 for size because we only want to store
  // and process elements for neutrons. Input file contains both proton and neutron terms.
  //ref_m = arma::zeros (total_lines/2, 6); 
  ref_m = arma::zeros (8, 6); 

  // clear the file stream, reset to read in the elements
  m_ref_in.clear();
  m_ref_in.seekg (0, std::ios::beg);

  std::size_t ele_in = 0;

  // read in and assign the references line by line
  // to the matrix ref_m
  for (std::size_t i = 0; i < total_lines; i++)
  {
	std::string orbit_dummy_1;
	std::string orbit_dummy_2;
	std::getline (m_ref_in, dummy);

	if (!dummy.length() || dummy[0] == '#')			// skip zero length lines and lines that start with #
	  continue;

	std::stringstream ss;

	ss << dummy;						// read in the line to stringstream ss

	ss >> orbit_dummy_1 >> orbit_dummy_2 >> ref_num >> n >> l >> j >> m_j >> tz; 	// assign values of the line

	// only extract neutron-neutron states
	if (tz == 1. && ele_in < 8)
	{
	  ref_m (ele_in,0) = ref_num;				// reference number of the line
    ref_m (ele_in,1) = n;					// principle quantum number
 	  ref_m (ele_in,2) = l;					// orbital angular mom.
	  ref_m (ele_in,3) = j * 0.5;				// total angular mom.
	  ref_m (ele_in,4) = m_j * 0.5;				// total angular mom. projection
	  ref_m (ele_in,5) = tz;				// isospin projection (should all be +1.0)
	  ele_in++;
	}
	
	
  }

  ref_m.print();						// print the resulting matrix

  // return final dimensionality of reference matrix
  return ele_in;
}

/**********************************************

Function to read in the actual m scheme matrix  elements 
calculated from Morten's code. 

THIS WILL BREAK IF THE INPUT FILE .dat CHANGES TITLE OR FORMAT.

**********************************************/


void fullm::read_in_matrix_m_scheme (int m_size, std::string m_mat_file)
{
  // input file stream
  const char * m_matrix_file = (m_mat_file).c_str();
  std::ifstream m_matrix_in (m_matrix_file);

  std::size_t total_lines = 0;
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

  // set the size of the m_scheme holder
  m_scheme_matrix.set_size (m_size, m_size);

  for (int a = 0; a < m_size; a++)
  {
  for (int b = 0; b < m_size; b++)
  {
	m_scheme_matrix (a,b).set_size (m_size, m_size);
  }
  }

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

 	for (int w = 0; w < m_size; w++)
	{
	  if (alpha == ref_m(w,0)) 
	    i = w;
	}
	
	for (int w = 0; w < m_size; w++)
	{
	  if (beta == ref_m(w,0)) 
	    j = w;
	}

	for (int w = 0; w < m_size; w++)
	{
	  if (gamma == ref_m(w,0)) 
	    k = w;
	}

	for (int w = 0; w < m_size; w++)
	{
	  if (delta == ref_m(w,0)) 
	    l = w;
	}

	if (i >= 0 && j >= 0 && k >= 0 && l >= 0)	// only activates if all 4 "if" statements above are true
	  m_scheme_matrix (i,j)(k,l) = value;

  }

  // reset index term on reference matrix to span 0 to m_size instead of e.g. 3, 4, 9, 10, etc... for neutrons

  for (int n = 0; n < m_size; n++)
	ref_m(n,0) = n;

  ref_m.print();


}

/*************************************************

Function to take care of degeneracies in basis states of
the matrix. Inactive for all but j-scheme calculation.

*************************************************/

double fullm::orb_degen (int, int, std::string, int)
{
  return 1.;
}

/************************************************

Function to return the size of our basis for the HF solver

************************************************/

int fullm::basis_extent ()
{
  return nmax;
}
