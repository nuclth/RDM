/*****************************************************************************

Programmer: Alex Dyhdalo


*****************************************************************************/

// necessary inclusion files
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include "flags.h"
#include "auxiliary.h"
#include "hamiltonian.h"
#include "output.h"



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
  const bool Q_flag = true; // Q START - LINEAR RELATIONS
  const bool G_flag = false; // G START - LINEAR REALTIONS

  const bool two_body_toggle = true;
  const bool diag_toggle 	 = false;
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




  size_t cmat_extent = full_matrix_extent (bsize, flag_pass);





  const std::string diag_file = "diagnostic_out/test_diag.dat";
  const std::string spda_file = "sdp_files/test_spd.dat";

  std::ofstream diag_out (diag_file);
  std::ofstream spda_out (spda_file);

 



  two_array h1_mat      (boost::extents[bsize][bsize]);
  two_array block_h2    (boost::extents[bsize*(bsize-1)/2][bsize*(bsize-1)/2]);
  two_array trans_h2    (boost::extents[bsize*(bsize-1)/2][bsize*(bsize-1)/2]);
  two_array block_mat   (boost::extents[2*bsize][4]);
  one_array oned_blocks (boost::extents[2*bsize]);

  fullm_populate_hamiltonian (bsize, h1_mat, block_h2, trans_h2, block_mat, oned_blocks, m_ref, m_mat, hw, diag_out, diag_toggle);







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
    init_Q_flag (Q_con, Q_val, bsize, cmat_extent, skip, trans_h2);
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
  	diag_format (diag_out, c_matrix, N_num, N_con, N_val, O_num, O_con, O_val, P_num, P_con, P_val, Q_num, Q_con, Q_val, G_num, G_con, G_val, flag_pass);
  }
 
  return EXIT_SUCCESS;
}

/************************************************

END MAIN PROGRAM

************************************************/

