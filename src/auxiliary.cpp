#include "auxiliary.h"
#include <sstream>
#include <fstream>

using std::boolalpha;

/********************************************************************************************************

Function which reads in the user specified values in inputs.inp to an instance of the struct parameters.

Returns:

	input - instance of the parameters struct

*********************************************************************************************************/

parameters read_in_inputs ()
{
	
  struct parameters input;						// local struct to hold user values
 
  std::ifstream input_in ("inputs.inp");		// the relevant input file stream
  std::string dummy;
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