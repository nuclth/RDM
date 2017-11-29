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


size_t count_P_cons (const std::string);
size_t total_tbme_states (const std::string tbme_filename);
size_t total_obme_states (const std::string obme_filename);

size_t count_NO_blocks (const two_array & array_ref_obme);
void populate_obme_blocks (one_array & obme_blocks, const two_array & array_ref_obme);

size_t count_no_flags (const std::string);




/********************************************

BEGIN MAIN PROGRAM

********************************************/

int main ()
{
  

  // create a struct to hold user specified parameters
  parameters input_params = read_in_inputs ();						             
  
  // Now create/copy input parameters for ease of reading the code into newly defined local variables
  const size_t nmax  = input_params.nmax;	
  const size_t basis_hw = input_params.basis_hw;   
  const size_t particles = input_params.particles;		 


  const bool N_flag = input_params.N_flag;
  const bool O_flag = input_params.O_flag;
  const bool P_flag = input_params.P_flag;

  const bool two_body_toggle = input_params.two_body_toggle;

  cout << "N FLAG: " << boolalpha << N_flag << endl;
  cout << "O FLAG: " << boolalpha << O_flag << endl;
  cout << "P FLAG: " << boolalpha << P_flag << endl;


  string_holder input_strings = string_reader (nmax, basis_hw);

  cout << ("Building system... ") << endl;

 

  const size_t bsize = total_obme_states (input_strings.ref_obme);

  const size_t tbme_size = total_tbme_states (input_strings.ref_tbme);
  two_array array_ref_tbme (boost::extents[tbme_size][13]);


  size_t P_num = 0;

  if (two_body_toggle) P_num = count_P_cons (input_strings.pflag_info);





  cout << "BSIZE: " << bsize << endl;
  cout << "TBME SIZE: " << tbme_size << endl;
  cout << "PNUM: " << P_num << endl << endl;
  



  if (!two_body_toggle and P_flag)
  {
  	std::cerr << "ERROR: TWO-BODY TOGGLE AND CONSTRAINT FLAGS INCOMPATIBLE - GIVING UP" << std::endl;
  	return EXIT_FAILURE;
  }	






  struct con_flags flag_pass;

  flag_pass.N_flag  = N_flag;
  flag_pass.O_flag  = O_flag;
  flag_pass.P_flag  = P_flag;

  flag_pass.two_body_toggle = two_body_toggle;




  const std::string sdpa_file = "sdp_files/test_sdp.dat-s";



  FILE * sdpa_out;


  const char * sdpa_char = (sdpa_file).c_str();

  sdpa_out = fopen (sdpa_char, "w");


  two_array array_ref_obme (boost::extents[bsize][7]);
  two_array h1_mat(boost::extents[bsize][bsize]);

  two_array h2_mat (boost::extents[tbme_size][tbme_size]);



  populate_hamiltonian (array_ref_obme, array_ref_tbme, h1_mat, h2_mat, input_strings.ref_obme, input_strings.me_obme, input_strings.ref_tbme, input_strings.me_tbme, two_body_toggle, nmax, input_strings.h2_flag);

  std::cout << "HAMILTONIAN BUILT" << std::endl;



  size_t NO_blocks = count_NO_blocks (array_ref_obme);

  one_array obme_blocks (boost::extents[NO_blocks]);

  populate_obme_blocks (obme_blocks, array_ref_obme);

  size_t cons = 0;
  size_t blocks = 0;

  if (N_flag)
  {
  	cons += 1;
  	blocks+= NO_blocks;
  }

  if (O_flag)
  {
  	cons += count_no_flags (input_strings.no_flag);
  	blocks+= NO_blocks;
  }

  if (two_body_toggle)
  	blocks+= 2;

  if (P_flag)
  	cons += P_num;



  fprintf (sdpa_out, "%lu\n", cons);
  fprintf (sdpa_out, "%lu\n", blocks);


  if (N_flag)
  {
  	for (size_t loop = 0; loop < bsize; loop++)
  	{
  		size_t block_size = array_ref_obme[loop][5];

  		if (block_size > 0) fprintf (sdpa_out, "%lu ", block_size);
  	}
  }



  if (O_flag)
  {
  	for (size_t loop = 0; loop < bsize; loop++)
  	{
  		size_t block_size = array_ref_obme[loop][5];

  		if (block_size > 0) fprintf (sdpa_out, "%lu ", block_size);
  	}
  }
  

  if (two_body_toggle)
  {
  	fprintf (sdpa_out, "%lu ", (size_t)array_ref_tbme[0][11]);
  	fprintf (sdpa_out, "%lu ", (size_t)array_ref_tbme[1][11]);  	
  }





  fprintf (sdpa_out, "\n");

  
  init_con_values (flag_pass, sdpa_out, bsize, tbme_size, particles, P_num, input_strings.no_flag);

  size_t con_count = 0;


  init_C_matrix (flag_pass, sdpa_out, h1_mat, h2_mat, con_count, obme_blocks, input_strings.no_flag, input_strings.h2_flag, 2 * NO_blocks);

  std::cout << "C MATRIX DONE" << std::endl;


  if (N_flag)
  {
    init_N_flag (sdpa_out, bsize, con_count, input_strings.no_flag);
    std::cout << "N FLAG DONE" << std::endl;
  }

  if (O_flag)
  {
    init_O_flag (sdpa_out, bsize, con_count, input_strings.no_flag, NO_blocks);
    std::cout << "O FLAG DONE" << std::endl;
  }


  if (P_flag)
  {
    init_P_flag (sdpa_out, bsize, con_count, particles, tbme_size, array_ref_tbme, input_strings.pflag_info, 2 * NO_blocks);
    std::cout << "P FLAG DONE" << std::endl;
  }


  return EXIT_SUCCESS;

}
 

/************************************************

END MAIN PROGRAM

************************************************/

size_t count_P_cons (const std::string p_filename)
{
  const char * p_file = p_filename.c_str();
  // input file stream for m_scheme
  std::ifstream ref_in (p_file);
 
  std::string dummy;

  size_t P_num = 0;

  size_t ob_b1, ob_b2, ob_block;
  size_t tb_b1, tb_b2, tb_block;
  bool new_flag;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy)) 
  {
	if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
    	continue;

	  	std::stringstream ss;

	 	ss << dummy;            // read in the line to stringstream ss
	  	ss >> ob_b1 >> ob_b2 >> ob_block >> tb_b1 >> tb_b2 >> tb_block >> new_flag;

	  	if (new_flag)
	  		++P_num;
  }

  return P_num;
}


/***************************************************************



***************************************************************/

size_t total_tbme_states (const std::string tbme_filename)
{
  const char * ref_file = tbme_filename.c_str();
  std::ifstream ref_in (ref_file);
 
  size_t total_lines = 0;
  std::string dummy;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy)) 
  {
	if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
    	continue;

  	++total_lines;
  }

  return total_lines;
}


/***************************************************************



***************************************************************/

size_t total_obme_states (const std::string ref_obme)
{
  const char * ref_file = ref_obme.c_str();
  std::ifstream ref_in (ref_file);
 
  size_t total_lines = 0;
  std::string dummy;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy)) 
  {
	if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
    	continue;

  	++total_lines;
  }

  return total_lines;
}


size_t count_NO_blocks (const two_array & array_ref_obme)
{
	size_t count = 0;
	size_t states = array_ref_obme.size();

	for (size_t loop = 0; loop < states; loop++) 
	{
		if (array_ref_obme[loop][5] > 0) count += 1;
	}

	return count;
}


void populate_obme_blocks (one_array & obme_blocks, const two_array & array_ref_obme)
{
	size_t count = 0;
	size_t states = array_ref_obme.size();

	for (size_t loop = 0; loop < states; loop++) 
	{
		if (array_ref_obme[loop][5] > 0)
		{
			obme_blocks[count] = array_ref_obme[loop][5];
			count++;
		}
	}
}


size_t count_no_flags (const std::string no_flag)
{
	const char * no_file = no_flag.c_str();
	// input file stream for m_scheme
	std::ifstream ref_in (no_file);
	 
	std::string dummy;
	size_t total_lines = 0;

	// find total number of defined reference lines
	while (std::getline (ref_in, dummy))
	{

		if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
	    	continue;

	    total_lines++;
	}

	return total_lines;
}