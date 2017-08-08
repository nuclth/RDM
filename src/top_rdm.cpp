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
#include "aux.h"
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

  if (two_body_toggle)
  	blocks++;
  
  if (P_flag)
  	cons += F3num;

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

  if (N_flag)
  {
  	spda_out << bsize << " ";
  }

  if (O_flag)
  {
  	spda_out << bsize << " ";
  }
  
  if (two_body_toggle)
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

