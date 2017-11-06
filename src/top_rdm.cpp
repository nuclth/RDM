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

//#include "boost/cstdlib.hpp"

#include "matrix_define.h"
#include "auxiliary.h"
#include "flags.h"
#include "hamiltonian.h"

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

using std::cout;
using std::endl;


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
  
  // define flags for all different combinations of conditions in RDM
  const bool two_body_toggle = true;


  // START CONSTRAINT FLAG DEFINE

  const bool N_flag  = true; // p START - TRACE CONDITION
  const bool O_flag  = true; // q START - LINEAR RELATIONS
  const bool P_flag  = true; // P START - TRACE CONDITION
  const bool Q_flag  = false; // Q START - LINEAR RELATIONS
  const bool G_flag  = false; // G START - LINEAR REALTIONSi

  const bool T1_flag = false;
  const bool T2_flag = false;

  const bool NN_flag = false;




  parameters input_params;						             // create a struct to hold user specified parameters
  input_params = read_in_inputs (); 					     // read in the values from inputs.inp, store in input_params	

  // Now create/copy input parameters for ease of reading the code into newly defined local variables
  const size_t nmax  = input_params.nmax;	
  const size_t basis_hw = input_params.basis_hw;   
  const size_t particles = input_params.particles;		 

  std::cout << ("Building system... ") << std::endl;

  // Create filenames from inputs
  std::stringstream parser;
  
  parser << "me_files/ref_files/nmax" << nmax << "_spm.dat";
  const std::string morten_spm = parser.str();

  // clear parser and repeat for other files
  parser.str("");
  parser.clear();

  parser << "me_files/ref_files/nmax" << nmax << "_python_sp.dat";
  const std::string ref_obme = parser.str();

  parser.str("");
  parser.clear();

  parser << "me_files/obme/nmax" << nmax << "_obme_hw" << basis_hw << ".dat";
  const std::string me_obme = parser.str();

  parser.str("");
  parser.clear();

  parser << "me_files/ref_files/nmax" << nmax << "_python_tb.dat";
  const std::string ref_tbme = parser.str();

  parser.str("");
  parser.clear();

  parser << "me_files/tbme/nmax" << nmax << "_tbme_hw" << basis_hw << ".dat";
  const std::string me_tbme = parser.str();

  parser.str("");
  parser.clear();

  parser << "flag_files/nmax" << nmax << "_python_noflag.dat";
  const std::string no_flag = parser.str();

  parser.str("");
  parser.clear();

  parser << "flag_files/nmax" << nmax << "_python_h2flag.dat";
  const std::string h2_flag = parser.str();

  parser.str("");
  parser.clear();

  parser << "flag_files/nmax" << nmax << "_python_pflag.dat";
  const std::string pflag_info = parser.str();

  parser.str("");
  parser.clear();


  const size_t bsize = total_obme_states (ref_obme);

  const size_t tbme_size = total_tbme_states (ref_tbme);
  two_array array_ref_tbme (boost::extents[tbme_size][13]);



  if (two_body_toggle) readin_ref_tbme (array_ref_tbme, ref_tbme);


  size_t P_num = 0;

  if (two_body_toggle) P_num = count_P_cons (pflag_info);





  cout << "BSIZE: " << bsize << endl;
  cout << "TBME SIZE: " << tbme_size << endl;
  cout << "PNUM: " << P_num << endl << endl;
  

  size_t Q_num = 0;
/*
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
*/
  size_t G_num = 0;
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

    G_num++;
  }
  }
  }
  }
*/



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

  flag_pass.N_flag  = N_flag;
  flag_pass.O_flag  = O_flag;
  flag_pass.P_flag  = P_flag;
  flag_pass.Q_flag  = Q_flag;
  flag_pass.G_flag  = G_flag;
  flag_pass.T1_flag = T1_flag;
  flag_pass.T2_flag = T2_flag;

  flag_pass.NN_flag = NN_flag;

  flag_pass.two_body_toggle = two_body_toggle;



//  const size_t F1num  = 1;
//  const size_t F2num  = bsize * (bsize + 1)/2;
  const size_t NN_num = 1;
//  const size_t F3num  = P_num;//tbme_size * (tbme_size + 1)/2;
  const size_t F7num  = Q_num;//bsize*bsize*bsize*bsize;
  const size_t F10num = G_num;//bsize*bsize*bsize*bsize;

  const size_t T1_num = 1;
  const size_t T2_num = 1; //T2_count (bsize);


//  const std::string diag_file = "diagnostic_out/test_diag.dat";
  const std::string sdpa_file = "sdp_files/test_sdp.dat-s";

//  std::ofstream diag_out (diag_file);
//  std::ofstream sdpa_out (sdpa_file);


//  const char * diag_char = (diag_file).c_str();
  const char * sdpa_char = (sdpa_file).c_str();

//  FILE * diag_out;
  FILE * sdpa_out;

//  diag_out = fopen (diag_char, "w");
  sdpa_out = fopen (sdpa_char, "w");


  two_array array_ref_obme (boost::extents[bsize][7]);
  two_array h1_mat(boost::extents[bsize][bsize]);

//  std::cout << "TBME SIZE: " << tbme_size << std::endl;

  two_array h2_mat (boost::extents[tbme_size][tbme_size]);


//  five_array h2_mat(boost::extents[bsize][bsize][bsize][bsize][5]);

  fullm_populate_hamiltonian (array_ref_obme, array_ref_tbme, h1_mat, h2_mat, ref_obme, me_obme, ref_tbme, me_tbme, two_body_toggle, nmax, h2_flag);

  std::cout << "HAMILTONIAN BUILT" << std::endl;

//  two_array comp_h2 (boost::extents[bsize*(bsize-1)/2][bsize*(bsize-1)/2]);

//  if (two_body_toggle) 
//  {
//  	compactify_h2 (ref_m, comp_h2, h2_mat, diag_out, diag_toggle);
//
//  	std::cout << "POTENTIAL COMPACTIFIED" << std::endl;
//  }


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
  	cons += count_no_flags (no_flag);
  	blocks+= NO_blocks;
  }

  if (two_body_toggle)
  	blocks+= 2;
  
  if (NN_flag)
  	cons += NN_num;

  if (P_flag)
  	cons += P_num;

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

  if (T1_flag)
  {
  	cons += T1_num;
  	blocks++;
  }

  if (T2_flag)
  {
  	cons += T2_num;
  	blocks++;
  }

  fprintf (sdpa_out, "%lu\n", cons);
  fprintf (sdpa_out, "%lu\n", blocks);

//  sdpa_out << cons << std::endl;
//  sdpa_out << blocks << std::endl;

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


  if (Q_flag)
  	fprintf (sdpa_out, "%lu ", tbme_size);

  if (G_flag)
  	fprintf (sdpa_out, "%lu ", bsize*bsize);



  if (T2_flag)
  	fprintf (sdpa_out, "%lu ", T2_DIM_count (bsize));


  fprintf (sdpa_out, "\n");

  
  init_con_values (flag_pass, sdpa_out, bsize, tbme_size, particles, P_num, no_flag);

  size_t con_count = 0;


  init_C_matrix (flag_pass, sdpa_out, h1_mat, h2_mat, con_count, obme_blocks, no_flag, h2_flag, 2 * NO_blocks);

  std::cout << "C MATRIX DONE" << std::endl;


  if (N_flag)
  {
    init_N_flag (sdpa_out, bsize, con_count, no_flag);
    std::cout << "N FLAG DONE" << std::endl;
  }

  if (O_flag)
  {
    init_O_flag (sdpa_out, bsize, con_count, no_flag, NO_blocks);
    std::cout << "O FLAG DONE" << std::endl;
  }

  if (NN_flag)
  {
  	init_NN_flag (sdpa_out, bsize, con_count);
  	std::cout << "NN FLAG DONE" << std::endl;
  }

  if (P_flag)
  {
    init_P_flag (sdpa_out, bsize, con_count, particles, tbme_size, array_ref_tbme, pflag_info, 2 * NO_blocks);
    std::cout << "P FLAG DONE" << std::endl;
  }

  if (Q_flag)
  {
    init_Q_flag (sdpa_out, bsize, con_count);
    std::cout << "Q FLAG DONE" << std::endl;
  }

  if (G_flag)
  {
    init_G_flag (sdpa_out, bsize, con_count);
    std::cout << "G FLAG DONE" << std::endl;
  }

  if (T1_flag)
  {
  	std::cout << "T1 FLAG DONE" << std::endl;
  }

  if (T2_flag)
  {
  	init_T2_flag (sdpa_out, bsize, con_count);
  	std::cout << "T2 FLAG DONE" << std::endl;
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