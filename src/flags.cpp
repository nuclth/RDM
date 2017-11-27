#include "flags.h"
#include "auxiliary.h"
#include <iomanip>



/***************************************************************

Function to create the constraint matrix values for the SDP solver.

These are the real numbers on the right hand side of the SDP.

***************************************************************/



void init_con_values (const con_flags flag_pass, FILE * sdpa_out, const size_t bsize, const size_t tbme_size, const size_t particles, const size_t P_num, const std::string no_flag)
{



  // N Trace condition
  if (flag_pass.N_flag)
  	fprintf(sdpa_out, "%lu ", particles);

  // Linear p q relations
  if (flag_pass.O_flag)
  {

	  const char * no_file = no_flag.c_str();
	  // input file stream for m_scheme
	  std::ifstream ref_in (no_file);
	 
	  std::string dummy;
	  size_t b1, b2, m1, m2;

	  // find total number of defined reference lines
	  while (std::getline (ref_in, dummy))
	  {

		if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
	    	continue;

	  	std::stringstream ss;

	 	ss << dummy;            // read in the line to stringstream ss
	  	ss >> b1 >> b2 >> m1 >> m2;

	    if (b1 == b2)
      		fprintf(sdpa_out, "%f ", 1.0);
	    else
      		fprintf(sdpa_out, "%f ", 0.0);
	  }

   }


  // P and p trace relation
   if (flag_pass.P_flag)
   {
      for (size_t i = 0; i < P_num; i++)      
      {
      	  fprintf(sdpa_out, "%f ", 0.0);
      }
   }



  fprintf(sdpa_out, "\n");

  return;
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


void init_C_matrix (const con_flags flag_pass, FILE * sdpa_out, const two_array & h1_mat, const two_array & h2_mat, size_t & con_count, const one_array & obme_blocks, const std::string no_flag, const std::string h2_flag, const size_t two_body_block_bias)
{


  const char * no_file = no_flag.c_str();
  // input file stream for m_scheme
  std::ifstream ref_in (no_file);
 
  std::string dummy;
  size_t b1a, b2a, m1a, m2a;

  size_t bcount = 0;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy))
  {

	if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
    	continue;

  	std::stringstream ss;

 	ss << dummy;            // read in the line to stringstream ss
  	ss >> b1a >> b2a >> m1a >> m2a;

  	double val1 = h1_mat [m1a][m2a] * -1.0;

  	if (b1a == 1 and b2a == 1) bcount++;

    if (val1 != 0. and b1a <= b2a)
    	fprintf(sdpa_out, "%lu %lu %lu %lu %f\n", con_count, bcount, b1a, b2a, val1);
  }


  if (flag_pass.two_body_toggle)
  {
	  const char * h2_file = h2_flag.c_str();
	  // input file stream for m_scheme
	  std::ifstream h2_in (h2_file);
	 
//	  std::string dummy;
	  size_t b1, b2, m1, m2, sp1, sp2, sp3, sp4, msp1, msp2, msp3, msp4, block;

	  // find total number of defined reference lines
	  while (std::getline (h2_in, dummy))
	  {

		if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
	    	continue;

	  	std::stringstream ss;

	 	ss << dummy;            // read in the line to stringstream ss
	  	ss >> b1 >> b2 >> m1 >> m2 >> sp1 >> sp2 >> sp3 >> sp4 >> msp1 >> msp2 >> msp3 >> msp4 >> block;

	  	double val3 = h2_mat [m1][m2] * -1./2.;

	  	block += two_body_block_bias;

	    if (val3 != 0. and b1 <= b2)
	    	fprintf(sdpa_out, "%lu %lu %lu %lu %f\n", con_count, block, b1, b2, val3);
	  }
  }


  con_count++;

  return;
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

void init_N_flag (FILE * sdpa_out, const size_t bsize, size_t & con_count, const std::string no_flag)
{


  const char * no_file = no_flag.c_str();
  // input file stream for m_scheme
  std::ifstream ref_in (no_file);
 
  std::string dummy;
  size_t b1, b2, m1, m2;

  size_t bcount = 0;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy))
  {

	if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
    	continue;

  	std::stringstream ss;

 	ss << dummy;            // read in the line to stringstream ss
  	ss >> b1 >> b2 >> m1 >> m2;

  	double val1 = kron_del (m1, m2);

  	if (b1 == 1 and b2 == 1) bcount++;

    if (val1 != 0. and b1 <= b2)
    	fprintf(sdpa_out, "%lu %lu %lu %lu %f\n", con_count, bcount, b1, b2, val1);
  }


    con_count++;


    return;
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

void init_O_flag (FILE * sdpa_out, const size_t bsize, size_t & con_count, const std::string no_flag, const size_t obme_block_count)
{

  const char * no_file = no_flag.c_str();
  // input file stream for m_scheme
  std::ifstream ref_in (no_file);
 
  std::string dummy;
  size_t b1, b2, m1, m2;

  size_t bcount1 = 0;
  size_t bcount2 = obme_block_count;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy))
  {

	if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
    	continue;

  	std::stringstream ss;

 	ss << dummy;            // read in the line to stringstream ss
  	ss >> b1 >> b2 >> m1 >> m2;


  	if (b1 == 1 and b2 == 1)
  	{
  		bcount1++;
  		bcount2++;
  	}

  	if (b1 <= b2)
  	{
    	fprintf(sdpa_out, "%lu %lu %lu %lu %f\n", con_count, bcount1, b1, b2, 1.0);
    	fprintf(sdpa_out, "%lu %lu %lu %lu %f\n", con_count, bcount2, b1, b2, 1.0);
  	}

    con_count++;
  }


  return;
}


/***************************************************************

Function to create our P constraint matrix to fix the relation
between the 2RDM partial trace and the 1RDM.

***************************************************************/


void init_P_flag (FILE * sdpa_out, const size_t bsize, size_t & con_count, const size_t N, const size_t tbme_size, 
	const two_array & array_ref_tbme, const std::string pflag_info, const size_t two_body_block_bias)
{

  const char * p_file = pflag_info.c_str();
  // input file stream for m_scheme
  std::ifstream ref_in (p_file);
 
  std::string dummy;

  size_t ob_b1, ob_b2, ob_block;
  size_t tb_b1, tb_b2, tb_block;
  bool new_flag;

  double prefac = -1.0 * (N - 1.0) / 2.0;

  bool start = true;

  // find total number of defined reference lines
  while (std::getline (ref_in, dummy))
  {

	if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
    	continue;

  	std::stringstream ss;

 	ss << dummy;            // read in the line to stringstream ss
  	ss >> ob_b1 >> ob_b2 >> ob_block >> tb_b1 >> tb_b2 >> tb_block >> new_flag;

  	tb_block += two_body_block_bias;

  	if (new_flag and !start)
    	con_count++;

    start = false;

  	if (new_flag)
    	fprintf(sdpa_out, "%lu %lu %lu %lu %f\n", con_count, ob_block, ob_b1, ob_b2, prefac);


    fprintf(sdpa_out, "%lu %lu %lu %lu %f\n", con_count, tb_block, tb_b1, tb_b2, 1.0);


  }



  return;
}
