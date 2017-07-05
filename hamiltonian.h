#ifndef hamiltonian_h
#define hamiltonian_h

#include <cstdio>
#include <string>
#include <fstream>
#include "arraydefs.h"

void populate_1body (const two_array & ref_m, two_array & h1_mat, const double hw, std::ofstream & diag_out, const bool diag_toggle);


void read_in_reference_m_scheme (two_array & ref_m, const std::string m_ref_file, std::ofstream & diag_out, const bool diag_toggle);


void read_in_matrix_m_scheme (const two_array & ref_m, five_array & h2_mat, const std::string m_mat_file);

void block_list (const size_t bsize, const two_array & ref_m, two_array & block_mat, std::ofstream & diag_out, const bool diag_toggle);

void fill_oned_blocks (one_array & oned_blocks, const two_array & block_mat);

void fullm_populate_hamiltonian (const size_t bsize, two_array & h1_mat, two_array & block_h2, two_array & trans_h2, two_array & block_mat, one_array & oned_blocks, const std::string reference_file, const std::string matrix_file, const double hw, std::ofstream & diag_out, const bool diag_toggle);

void compactify_h2 (const two_array & ref_m, two_array & comp_h2, five_array & h2_mat, two_array & basis_ref, std::ofstream & diag_out, const bool diag_toggle);

void blockdiag_h2 (const two_array & ref_m, two_array & block_h2, five_array & h2_mat, const two_array & block_mat, two_array & basis_ref, std::ofstream & diag_out, const bool diag_toggle);

void create_transformation (const two_array & compact_ref, const two_array & block_ref, two_array & trans_h2, size_t extent, std::ofstream & diag_out, const bool diag_toggle);

void verify_tranformation  (const two_array & comp_h2, const two_array & block_h2, const two_array & trans_h2);

#endif