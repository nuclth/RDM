/*****************************************************************************
Programmer: Alex Dyhdalo
Last Mod: 8/2017
*****************************************************************************/

#include "hamiltonian.h"
#include "mscheme.h"
#include "jscheme.h"
#include <iostream>


/***************************************************************

Wrapper function to read in single-particle m scheme reference 
file, populate 1-body Hamiltonian matrix elements, and then 
populate 2-body matrix elements.

***************************************************************/


void fullm_populate_hamiltonian (two_array & ref_m, two_array & h1_mat, five_array & h2_mat, const std::string reference_file, const std::string matrix_file, const double hw, std::ofstream & diag_out, const bool diag_toggle) 
{

  try
  {
    read_in_reference_m_scheme (ref_m, reference_file, diag_out, diag_toggle);
    std::cout << "REFERENCE READ" << std::endl;
    mpopulate_1body (ref_m, h1_mat, hw);
    std::cout << "1 BODY POPULATED" << std::endl;
    read_in_matrix_m_scheme (ref_m, h2_mat, matrix_file);
    std::cout << "2 BODY POPULATED" << std::endl;
  }

  catch (const char * msg)
  {
    std::cerr << msg << std::endl;  
  }
 
}


/***************************************************************


 Function to populate the 1-body part of the Hamiltonian. 
 The 1-body part here is T + U for kinetic energy T and external potential U


***************************************************************/


void mpopulate_1body (const two_array & ref_m, two_array & h1_mat, const double hw)
{

  size_t mat_length = h1_mat.size();


  for (size_t i = 0; i < mat_length; ++i)
  {
    double n = ref_m [i][1];
    double l = ref_m [i][2];

    h1_mat[i][i] = (2.*n + l + 1.5) * hw;
  }

}


/***************************************************************

Wrapper function to read in single-particle j scheme reference 
file, populate 1-body Hamiltonian matrix elements, and then 
populate 2-body matrix elements.

***************************************************************/


void fullj_populate_hamiltonian (two_array & ref_j, two_array & twop_basis, two_array & h1_mat, two_array & h2_mat, const std::string reference_file, const std::string matrix_file, const double hw, std::ofstream & diag_out, const bool diag_toggle) 
{

  try
  {
    read_in_reference_j_scheme (ref_j, reference_file, diag_out, diag_toggle);
    std::cout << "REFERENCE READ" << std::endl;
    jpopulate_1body (ref_j, h1_mat, hw);
    std::cout << "1 BODY POPULATED" << std::endl;

    size_t total_twop_basis = count_twopart_jscheme (ref_j);
    std::cout << "J SCHEME TWO-PARTICLE BASIS EXTENT: " << total_twop_basis << std::endl;

    twop_basis.resize(boost::extents[total_twop_basis][3]);
    h2_mat.resize(boost::extents[total_twop_basis][total_twop_basis]);

    create_2pbasis_jscheme (twop_basis, ref_j);

    read_in_matrix_j_scheme (h2_mat, twop_basis, matrix_file);
    std::cout << "2 BODY POPULATED" << std::endl;
  }

  catch (const char * msg)
  {
    std::cerr << msg << std::endl;  
  }
 
}

/***************************************************************


 Function to populate the 1-body part of the Hamiltonian. 
 The 1-body part here is T + U for kinetic energy T and external potential U


***************************************************************/

void jpopulate_1body (const two_array & ref_j, two_array & h1_mat, const double hw)
{

  size_t mat_length = h1_mat.size();


  for (size_t i = 0; i < mat_length; ++i)
  {
 // 	double j     = ref_j[i][3];
    double degen = 1.;//(2.*j + 1.);
    h1_mat[i][i] = ref_j[i][6] * degen;
  }

}

/***************************************************************



***************************************************************/


void compactify_h2 (const two_array & ref_m, two_array & comp_h2, five_array & h2_mat, std::ofstream & diag_out, const bool diag_toggle)
{
  const size_t bsize = h2_mat.size();

  double value;

  size_t alpha, beta, gamma, delta;

  if(diag_toggle)
  {
    diag_out << "Output of nonzero 2-body ME" << std::endl << std::endl;    
  }

  for (size_t i = 0; i < bsize; ++i)
  {
  for (size_t j = 0; j < bsize; ++j)
  {
//    if (i >= j)
//      continue;

  for (size_t k = 0; k < bsize; ++k)
  {
  for (size_t l = 0; l < bsize; ++l)
  {
//    if (k >= l)
//      continue;

            value = h2_mat [i][j][k][l][0];

            alpha = h2_mat [i][j][k][l][1];
            beta  = h2_mat [i][j][k][l][2];
            gamma = h2_mat [i][j][k][l][3];
            delta = h2_mat [i][j][k][l][4];

            if (diag_toggle)
            {
              if (value != 0.)
              {
              	diag_out << alpha << " " << beta << " " << gamma << " " << delta << " " << value << " \t";

                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (alpha == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\t";
                }


                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (beta == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\t";
                }


                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (gamma == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\t";
                }


                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (delta == ref_m[loop][0])
                    diag_out << ref_m[loop][6] << "\n";
                }

              }
            }

//            size_t ips = i + 1;
//            size_t jps = j + 1;
//            size_t kps = k + 1;
//            size_t lps = l + 1;

//            size_t left  = jps - ips + (2*bsize - ips) * (ips - 1)/2 - 1;
//            size_t right = lps - kps + (2*bsize - kps) * (kps - 1)/2 - 1;

            size_t left  = i * bsize + j;
            size_t right = k * bsize + l;

  /*         if (i == 0 && j == 4 && k == 0)
            {
            	std::cout << ref_m[l][1] << "\t" << ref_m[l][2] << "\n";
            }

           if (i == 0 && j == 1)
            {
            	std::cout << right << " " << l << " " << k << "\t" << value << "\n"; 
            }
*/

//            std::cout << value << std::endl;

            comp_h2 [left][right] = value;

  }
  }
  }
  }

  if (diag_toggle)
    diag_out << std::endl << std::endl;
}


