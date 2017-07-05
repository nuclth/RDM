/*****************************************************************************

Programmer: Alex Dyhdalo


*****************************************************************************/

// necessary definitions
#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "boost/cstdlib.hpp"
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <cstdio>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "flags.h"
#include "auxiliary.h"
#include "hamiltonian.h"
#include "spd.h"





/********************************************

BEGIN MAIN PROGRAM

********************************************/

int main ()
{
  
  parameters input_params;						             // create a struct to hold user specified parameters
  input_params = read_in_inputs (); 					     // read in the values from inputs.inp, store in input_params	

  // Now create/copy input parameters for ease of reading the code into newly defined local variables


  const size_t bsize = input_params.basis;	   // the size of the HO basis used, ONLY ACTIVE FOR S-WAVE
  const size_t particles = input_params.particles;		 // the number of neutrons (particles) in the trap
  const double hw = input_params.hw;					     // hbar * omega 
  const std::string m_ref = input_params.m_ref;    // single particle reference file - m scheme
  const std::string m_mat = input_params.m_mat;    // m scheme matrix elements


  std::cout << ("Building system... ") << std::endl;



  size_t Q_num_loop = 0;

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

    Q_num_loop++;
  }
  }
  }
  }

  size_t G_num_loop = 0;

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

    G_num_loop++;
  }
  }
  }
  }



  // define flags for all different combinations of conditions in RDM


  // START CONSTRAINT FLAG DEFINE

  const bool N_flag = true;  // p START - TRACE CONDITION
  const bool O_flag = true;  // q START - LINEAR RELATIONS
  const bool P_flag = true;  // P START - TRACE CONDITION
  const bool Q_flag = false; // Q START - LINEAR RELATIONS
  const bool G_flag = false; // G START - LINEAR REALTIONS

  const bool two_body_toggle = true;
  const bool diag_toggle 	 = true;
  const bool redundant_check = false;


  if (!two_body_toggle and (P_flag or Q_flag or G_flag))
  {
  	std::cerr << "ERROR: TWO-BODY TOGGLE AND F CONSTRAINT FLAGS INCOMPATIBLE - GIVING UP" << std::endl;
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
  flag_pass.diag_toggle 	= diag_toggle;
  flag_pass.redundant_check = redundant_check;



  const size_t N_num = 1;
  const size_t O_num = bsize * (bsize + 1)/2;
  const size_t P_num = bsize * (bsize + 1)/2;
  const size_t Q_num = Q_num_loop;	//bsize*bsize*bsize*bsize;
  const size_t G_num = G_num_loop;	//bsize*bsize*bsize*bsize;




  size_t cmat_extent = 2 * bsize;

  if (two_body_toggle)
  {
      cmat_extent = 2*bsize + bsize * (bsize-1)/2;

    if (Q_flag)
      cmat_extent = 2*bsize + 2 * bsize * (bsize-1)/2;

    if (G_flag)
      cmat_extent = 2*bsize + 2 * bsize * (bsize-1)/2 + bsize * bsize;
  }


  const std::string diag_file = "diagnostic_out/test_diag.dat";
  const std::string spda_file = "sdp_files/test_spd.dat";

  std::ofstream diag_out (diag_file);
  std::ofstream spda_out (spda_file);

 


  two_array ref_m     (boost::extents[bsize][7]);
  two_array h1_mat    (boost::extents[bsize][bsize]);
  five_array h2_mat   (boost::extents[bsize][bsize][bsize][bsize][5]);


  two_array block_mat (boost::extents[15][4]);
  two_array bcopy     (boost::extents[15][4]);

  two_array comp_h2  (boost::extents[bsize*(bsize-1)/2][bsize*(bsize-1)/2]);
  two_array block_h2 (boost::extents[bsize*(bsize-1)/2][bsize*(bsize-1)/2]);
  two_array trans_h2 (boost::extents[bsize*(bsize-1)/2][bsize*(bsize-1)/2]);

  two_array comp_basis_ref  (boost::extents[bsize*(bsize-1)/2][2]);
  two_array block_basis_ref (boost::extents[bsize*(bsize-1)/2][2]);




  fullm_populate_hamiltonian (ref_m, h1_mat, h2_mat, m_ref, m_mat, hw, diag_out, diag_toggle, block_mat);

  size_t p_subblocks = block_list (ref_m, block_mat, bcopy);

  block_mat.resize(boost::extents[p_subblocks][4]);
  bcopy.resize(boost::extents[0][0]);

//  print (std::cout, block_mat);

  one_array oned_blocks (boost::extents[p_subblocks]);

  size_t block_num = fill_oned_blocks (p_subblocks, oned_blocks, block_mat);

  oned_blocks.resize(boost::extents[block_num]);

//  std::cout << "HAMILTONIAN POPULATED" << std::endl;



  compactify_h2 (ref_m, comp_h2,  h2_mat, diag_out, diag_toggle, comp_basis_ref);

  blockdiag_h2  (ref_m, block_h2, h2_mat, diag_out, diag_toggle, block_mat, block_basis_ref);


  create_transformation (comp_basis_ref, block_basis_ref, trans_h2, bsize*(bsize-1)/2);

  verify_tranformation  (comp_h2, block_h2, trans_h2);

//  print(std::cout, comp_h2);





  two_array c_matrix (boost::extents[cmat_extent][cmat_extent]);

  three_array N_con (boost::extents[0][0][0]);
  three_array O_con (boost::extents[0][0][0]);
  three_array P_con (boost::extents[0][0][0]);
  three_array Q_con (boost::extents[0][0][0]);
  three_array G_con (boost::extents[0][0][0]);

  one_array   N_val (boost::extents[0]);
  one_array   O_val (boost::extents[0]);
  one_array   P_val (boost::extents[0]);
  one_array   Q_val (boost::extents[0]);
  one_array   G_val (boost::extents[0]);


  if (N_flag)
  {
    N_con.resize(boost::extents[N_num][cmat_extent][cmat_extent]);
    N_val.resize(boost::extents[N_num]);
  }

  if (O_flag)
  {
    O_con.resize(boost::extents[O_num][cmat_extent][cmat_extent]);
    O_val.resize(boost::extents[O_num]);
  }

  if (P_flag)
  {
    P_con.resize(boost::extents[P_num][cmat_extent][cmat_extent]);
    P_val.resize(boost::extents[P_num]);
  }


  if (Q_flag)
  {
    Q_con.resize(boost::extents[Q_num][cmat_extent][cmat_extent]);
    Q_val.resize(boost::extents[Q_num]);
  }

 
  if (G_flag)
  {
    G_con.resize(boost::extents[G_num][cmat_extent][cmat_extent]);
    G_val.resize(boost::extents[G_num]);
  }





  create_c_matrix (h1_mat, block_h2, c_matrix, two_body_toggle);

  if (N_flag)
  {
    two_array  F1_build_1 (boost::extents[bsize][bsize]);
    init_N_flag (N_con, N_val, bsize, particles);
    F1_build_1.resize(boost::extents[0][0]);
    std::cout << "N FLAG DONE" << std::endl;
  }



  if (O_flag)
  {
    four_array F2_build_1 (boost::extents[bsize][bsize][bsize][bsize]);
    init_O_flag (F2_build_1, O_con, O_val, bsize);
    F2_build_1.resize(boost::extents[0][0][0][0]);
    std::cout << "O FLAG DONE" << std::endl;
  }

  if (P_flag)
  {
    four_array F3_build_1 (boost::extents[bsize][bsize][bsize][bsize]);
    six_array  F3_build_3 (boost::extents[bsize][bsize][bsize][bsize][bsize][bsize]);
    init_P_flag (F3_build_1, F3_build_3, P_con, P_val, bsize, particles, trans_h2);
    F3_build_1.resize(boost::extents[0][0][0][0]);
    F3_build_3.resize(boost::extents[0][0][0][0][0][0]);
    std::cout << "P FLAG DONE" << std::endl;  
  }



  if (Q_flag)
  {
    size_t skip = 0;
    six_array   F7_build_1 (boost::extents[bsize][bsize][bsize][bsize][bsize][bsize]);
    eight_array F7_build_3 (boost::extents[bsize][bsize][bsize][bsize][bsize][bsize][bsize][bsize]);
    eight_array F7_build_5 (boost::extents[bsize][bsize][bsize][bsize][bsize][bsize][bsize][bsize]);
    init_Q_flag (F7_build_1, F7_build_3, F7_build_5, Q_con, Q_val, bsize, cmat_extent, skip);
    F7_build_1.resize(boost::extents[0][0][0][0][0][0]);
    F7_build_3.resize(boost::extents[0][0][0][0][0][0][0][0]);
    F7_build_5.resize(boost::extents[0][0][0][0][0][0][0][0]);
    std::cout << "Q FLAG DONE" << std::endl;
  }


  if (G_flag)
  {
    six_array F10_build_1 (boost::extents[bsize][bsize][bsize][bsize][bsize][bsize]);
    eight_array F10_build_3 (boost::extents[bsize][bsize][bsize][bsize][bsize][bsize][bsize][bsize]);
    eight_array F10_build_5 (boost::extents[bsize][bsize][bsize][bsize][bsize][bsize][bsize][bsize]);
    init_G_flag (F10_build_1, F10_build_3, F10_build_5, G_con, G_val, bsize);
    F10_build_1.resize(boost::extents[0][0][0][0][0][0]);
    F10_build_3.resize(boost::extents[0][0][0][0][0][0][0][0]);
    F10_build_5.resize(boost::extents[0][0][0][0][0][0][0][0]);
    std::cout << "G FLAG DONE" << std::endl;
  }



  std::cout << "FLAG INITIALIZATION DONE" << std::endl;


  std::cout << "CREATING SPDA FILE" << std::endl;

  create_spda_file (block_mat, oned_blocks, c_matrix, flag_pass, N_con, N_val, O_con, O_val, P_con, P_val, Q_con, Q_val, G_con, G_val, bsize, cmat_extent, spda_out);



  if (diag_toggle)
  {
  	diag_format (diag_out, h1_mat, comp_h2, block_h2, comp_basis_ref, block_basis_ref, trans_h2, c_matrix, N_num, N_con, N_val, O_num, O_con, O_val, P_num, P_con, P_val, Q_num, Q_con, Q_val, G_num, G_con, G_val, flag_pass);
  }
 
  return EXIT_SUCCESS;
}

/************************************************

END MAIN PROGRAM

************************************************/

