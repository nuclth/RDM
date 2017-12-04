#include "auxiliary.h"
#include <sstream>
#include <fstream>

using std::boolalpha;
using std::string;

/********************************************************************************************************

Function which reads in the user specified values in inputs.inp to an instance of the struct parameters.

Returns:

	input - instance of the parameters struct

*********************************************************************************************************/

parameters read_in_inputs ()
{
	
  struct parameters input;						// local struct to hold user values
 
  std::ifstream input_in ("inputs.inp");		// the relevant input file stream
  string dummy;
  int counter = 0;								// counter lets us know which variable we're reading in

  while (std::getline (input_in, dummy)) 		// while runs as long as there are lines to read in the input file
  {
	if (!dummy.length() || dummy[0] == '#')		// skip zero length lines and lines that start with #
	  continue;

	std::stringstream ss;

	ss << dummy;								// read in the line to the stringstream ss

	if (counter == 0)							// read in parameters by the way they're ordered in "inputs.inp"
		ss >> input.nmax;

	if (counter == 1)
		ss >> input.basis_hw;

	if (counter == 2)
		ss >> input.particles;

	if (counter == 3)
		ss >> boolalpha >> input.N_flag;

	if (counter == 4)
		ss >> boolalpha >> input.O_flag;

	if (counter == 5)
		ss >> boolalpha >> input.two_body_toggle;

	if (counter == 6)
		ss >> boolalpha >> input.P_flag;

	counter++;									// increase the counter after reading in a parameter	
  }

  return input;									// return the struct 
}


/********************************************************************************************************

Function which reads in the filenames of various auxiliary files into a string_holder struct.

Arguments:
	
	nmax - Nmax truncation value

	basis_hw - energy parameter (hbar omega) of the basis expansion

Returns:
	
	input - instance of the string_holder struct

*********************************************************************************************************/


string_holder string_reader (const size_t nmax, const size_t basis_hw)
{

  struct string_holder input;

  
  std::stringstream parser;											// create stringstream object
  
  parser << "me_files/ref_files/nmax" << nmax << "_spm.dat";        // read in filename and save string from parser
  input.morten_spm = parser.str();


  parser.str("");											  		// clear parser and repeat for other files
  parser.clear();

  parser << "me_files/ref_files/nmax" << nmax << "_python_sp.dat";
  input.ref_obme = parser.str();

  parser.str("");
  parser.clear();

  parser << "me_files/obme/nmax" << nmax << "_obme_hw" << basis_hw << ".dat";
  input.me_obme = parser.str();

  parser.str("");
  parser.clear();

  parser << "me_files/ref_files/nmax" << nmax << "_python_tb.dat";
  input.ref_tbme = parser.str();

  parser.str("");
  parser.clear();

  parser << "me_files/tbme/nmax" << nmax << "_tbme_hw" << basis_hw << ".dat";
  input.me_tbme = parser.str();

  parser.str("");
  parser.clear();

  parser << "flag_files/nmax" << nmax << "_python_noflag.dat";
  input.no_flag = parser.str();

  parser.str("");
  parser.clear();

  parser << "flag_files/nmax" << nmax << "_python_h2flag.dat";
  input.h2_flag = parser.str();

  parser.str("");
  parser.clear();

  parser << "flag_files/nmax" << nmax << "_python_pflag.dat";
  input.pflag_info = parser.str();

  parser.str("");
  parser.clear();

  return input;
}


/***************************************************************

Function to count the total number of one-body basis states.

***************************************************************/

size_t total_obme_states (const string ref_obme)
{
  const char * ref_file = ref_obme.c_str();
  std::ifstream ref_in (ref_file);
 
  size_t total_lines = 0;
  string dummy;

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

Function to count the total number of two-body basis states.

***************************************************************/

size_t total_tbme_states (const string tbme_filename)
{
  const char * ref_file = tbme_filename.c_str();
  std::ifstream ref_in (ref_file);
 
  size_t total_lines = 0;
  string dummy;

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

Function to count the total number of O flag constraint matrices. 

***************************************************************/

size_t count_no_flags (const string no_flag)
{
	const char * no_file = no_flag.c_str();
	// input file stream for m_scheme
	std::ifstream ref_in (no_file);
	 
	string dummy;
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


/***************************************************************

Function to count the total number of P flag constraint matrices.

***************************************************************/

size_t count_P_cons (const string p_filename)
{
  const char * p_file = p_filename.c_str();
  // input file stream for m_scheme
  std::ifstream ref_in (p_file);
 
  string dummy;

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



/*****************************************************************************

Function to output the total number of constraints and blocks in our SDP file

*****************************************************************************/

void init_con_blocks (const parameters flag_pass, const size_t NO_blocks, const size_t O_num, const size_t P_num, FILE * sdpa_out)
{

  size_t cons = 0;
  size_t blocks = 0;

  if (flag_pass.N_flag)
  {
  	cons += 1;
  	blocks+= NO_blocks;
  }

  if (flag_pass.O_flag)
  {
  	cons += O_num;
  	blocks+= NO_blocks;
  }

  if (flag_pass.two_body_toggle)
  	blocks+= 2;

  if (flag_pass.P_flag)
  	cons += P_num;



  fprintf (sdpa_out, "%lu\n", cons);
  fprintf (sdpa_out, "%lu\n", blocks);

  return;
}