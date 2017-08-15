/*****************************************************************************

Programmer: Alex Dyhdalo
Last Mod: 8/2017

This is the main wrapper file that calls all the other machinery to generate the SDP file.

The struct below holds all the user specified parameters which are read in via the 
read_in_inputs function. User inputs are read in from the file "inputs.inp". This is 
done to avoid having to recompile when the user changes input parameters. 

*****************************************************************************/

// necessary definitions

#include "matrix_define.h"
#include "auxiliary.h"
#include "flags.h"
#include "hamiltonian.h"

#include <cstdio>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>




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

  const bool two_body_toggle = input_params.two_body_toggle;
  const bool diag_toggle = input_params.diag_toggle;

  const bool N_flag = input_params.N_flag; // p START - TRACE CONDITION
  const bool O_flag = input_params.O_flag; // q START - LINEAR RELATIONS
  const bool P_flag = input_params.P_flag; // P START - TRACE CONDITION
  const bool Q_flag = input_params.Q_flag; // Q START - LINEAR RELATIONS
  const bool G_flag = input_params.G_flag; // G START - LINEAR REALTIONS

  const std::string spda_file = input_params.sdp_file;
  const std::string diag_file = input_params.diag_file;

  const bool mscheme_toggle = input_params.mscheme_toggle;

  std::cout << ("Building system... ") << std::endl;


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



  const size_t N_num  = 1;
  const size_t O_num  = bsize * (bsize + 1)/2;
  const size_t P_num  = bsize * (bsize + 1)/2;
  const size_t Q_num = Q_loop (bsize);
  const size_t G_num = G_loop (bsize);



  std::ofstream diag_out (diag_file);
  std::ofstream spda_out (spda_file);


  header_sdp_file (flag_pass, bsize, N_num, O_num, P_num, Q_num, G_num, spda_out);

  two_array ref_m  (boost::extents[bsize][7]);
  two_array h1_mat (boost::extents[bsize][bsize]);

  five_array h2_mat(boost::extents[bsize][bsize][bsize][bsize][5]);
  two_array comp_h2 (boost::extents[bsize*bsize][bsize*bsize]);



  if (mscheme_toggle)
  {
    fullm_populate_hamiltonian (ref_m, h1_mat, h2_mat, m_ref, m_mat, hw, diag_out, diag_toggle);
    std::cout << "HAMILTONIAN BUILT" << std::endl;
  	compactify_h2 (ref_m, comp_h2, h2_mat, diag_out, diag_toggle);
  	std::cout << "POTENTIAL COMPACTIFIED" << std::endl;
  }

  else if (!mscheme_toggle)
  {
    h2_mat.resize(boost::extents[0][0][0][0][0]);
    fullj_populate_hamiltonian (ref_m, h1_mat, comp_h2, m_ref, m_mat, hw, diag_out, diag_toggle);
    std::cout << "HAMILTONIAN BUILT" << std::endl;
  }


  
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
    init_P_flag (spda_out, bsize, con_count, particles);
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

