/*****************************************************************************

Programmer: Alex Dyhdalo

Main program to create the SDP output file. 

*****************************************************************************/

// necessary definitions

// standard libraries
#include <cstdio>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

// custom header files
#include "matrix_define.h"
#include "auxiliary.h"
#include "flags.h"
#include "hamiltonian.h"

using std::cout;
using std::endl;
using std::boolalpha;


/********************************************

BEGIN MAIN PROGRAM

********************************************/

int main ()
{
  

  parameters input_params = read_in_inputs (); 							// create a struct to hold user specified parameters
  
  
  const size_t nmax  = input_params.nmax;								// create local copies of input params
  const size_t basis_hw = input_params.basis_hw;   
  const size_t particles = input_params.particles;		 
  const bool N_flag = input_params.N_flag;
  const bool O_flag = input_params.O_flag;
  const bool P_flag = input_params.P_flag;
  const bool two_body_toggle = input_params.two_body_toggle;


  const std::string sdpa_file = "sdp_files/test_sdp.dat-s";				// set up output file for SDP
  FILE * sdpa_out;
  const char * sdpa_char = (sdpa_file).c_str();
  sdpa_out = fopen (sdpa_char, "w");


  string_holder input_strings = string_reader (nmax, basis_hw);			// create struct for string names of data files

 

  const size_t bsize = total_obme_states (input_strings.ref_obme);		// get total # of 1-body and 2-body states
  const size_t tbme_size = total_tbme_states (input_strings.ref_tbme);


  const size_t O_num = count_no_flags (input_strings.no_flag);			// get total # of constraint matrices in the O and P flag
  const size_t P_num = count_P_cons (input_strings.pflag_info);	




  cout << "N FLAG: " << boolalpha << N_flag << endl;
  cout << "O FLAG: " << boolalpha << O_flag << endl;
  cout << "P FLAG: " << boolalpha << P_flag << endl;
  cout << "NMAX: " << nmax << endl;
  cout << "TOTAL SP ORBITALS: " << bsize << endl;
  cout << "TOTAL TB STATES: " << tbme_size << endl;
  cout << "TOTAL O FLAG CONSTRAINTS: " << O_num << endl;
  cout << "TOTAL P FLAG CONSTRAINTS: " << P_num << endl << endl;
  cout << "Building system... " << endl;



  if (!two_body_toggle and P_flag)
  {
  	std::cerr << "ERROR: TWO-BODY TOGGLE AND CONSTRAINT FLAGS INCOMPATIBLE - GIVING UP" << std::endl;
  	return EXIT_FAILURE;
  }	




  two_array array_ref_obme (boost::extents[bsize][7]);					// 1- and 2-body reference arrays
  two_array array_ref_tbme (boost::extents[tbme_size][13]);

  two_array h1_mat(boost::extents[bsize][bsize]);						// 1- and 2-body hamiltonians
  two_array h2_mat (boost::extents[tbme_size][tbme_size]);



  populate_hamiltonian (array_ref_obme, array_ref_tbme, h1_mat, h2_mat, input_strings.ref_obme, input_strings.me_obme, input_strings.ref_tbme, input_strings.me_tbme, two_body_toggle, nmax, input_strings.h2_flag);


  std::cout << "HAMILTONIAN BUILT" << std::endl;



  size_t NO_blocks = count_NO_blocks (array_ref_obme);								// count number of blocks in 1-body matrix 
  one_array obme_blocks (boost::extents[NO_blocks]);								// create and populate array to hold 1-body block values
  populate_obme_blocks (obme_blocks, array_ref_obme);



  output_con_blocks (input_params, NO_blocks, O_num, P_num, sdpa_out);				// output total # of constraints and blocks to SDP file
  output_blocks (input_params, array_ref_obme, array_ref_tbme, bsize, sdpa_out); 	// output matrix block sizes to SDP file

  
  init_con_values (input_params, sdpa_out, bsize, tbme_size, particles, P_num, input_strings.no_flag); // output constraint values to SDP file

  size_t con_count = 0; // constraint matrix counter

  init_C_matrix (input_params, sdpa_out, h1_mat, h2_mat, con_count, obme_blocks, input_strings.no_flag, input_strings.h2_flag, 2 * NO_blocks); // output H constraint to SDP file

  std::cout << "C MATRIX DONE" << std::endl;


  if (N_flag)
  {
    init_N_flag (sdpa_out, bsize, con_count, input_strings.no_flag); // output N particle matrix constraint to SDP file
    std::cout << "N FLAG DONE" << std::endl;
  }

  if (O_flag)
  {
    init_O_flag (sdpa_out, bsize, con_count, input_strings.no_flag, NO_blocks); // output O flag matrix constraints to SDP file
    std::cout << "O FLAG DONE" << std::endl;
  }


  if (P_flag)
  {
    init_P_flag (sdpa_out, bsize, con_count, particles, tbme_size, array_ref_tbme, input_strings.pflag_info, 2 * NO_blocks); // output P flag matrix constraints to SDP file
    std::cout << "P FLAG DONE" << std::endl;
  }


  return EXIT_SUCCESS;

}
 

/************************************************

END MAIN PROGRAM

************************************************/