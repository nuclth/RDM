#ifndef hamiltonian_h
#define hamiltonian_h

#include <cstdio>
#include <string>
#include <fstream>
#include "arraydefs.h"

void populate_1body (const two_array & ref_m, two_array & h1_mat, const double hw);

void read_in_reference_m_scheme (two_array & ref_m, const std::string m_ref_file, std::ofstream & diag_out, const bool diag_toggle);

void read_in_matrix_m_scheme (const two_array & ref_m, five_array & h2_mat, const std::string m_mat_file);

size_t block_list (const two_array & ref_m, two_array & block_mat, two_array & bcopy);

size_t fill_oned_blocks (const size_t sub_blocks, one_array & oned_blocks, const two_array & block_mat);

void fullm_populate_hamiltonian (two_array & ref_m, two_array & h1_mat, five_array & h2_mat, const std::string reference_file, const std::string matrix_file, const double hw, std::ofstream & diag_out, const bool diag_toggle, two_array & block_mat);

void compactify_h2 (const two_array & ref_m, two_array & comp_h2, five_array & h2_mat, std::ofstream & diag_out, const bool diag_toggle, two_array & basis_ref);

void blockdiag_h2 (const two_array & ref_m, two_array & block_h2, five_array & h2_mat, std::ofstream & diag_out, const bool diag_toggle, const two_array & block_mat, two_array & basis_ref);

void create_transformation (const two_array & compact_ref, const two_array & block_ref, two_array & trans_h2, size_t extent);

void verify_tranformation  (const two_array & comp_h2, const two_array & block_h2, const two_array & trans_h2);

#endif