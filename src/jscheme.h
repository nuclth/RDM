/*****************************************************************************
Programmer: Alex Dyhdalo
Last Mod: 8/2017
*****************************************************************************/

#ifndef JSCHEME_H
#define JSCHEME_H

#include "matrix_define.h"
#include <string>
#include <fstream>

void read_in_reference_j_scheme (two_array & ref_j, const std::string j_ref_file, std::ofstream & diag_out, const bool diag_toggle);
void read_in_matrix_j_scheme (two_array & h2_mat, const two_array & twop_basis, const std::string j_mat_file);
size_t count_twopart_jscheme (const two_array & ref_j);
void create_2pbasis_jscheme (two_array & twop_basis, const two_array & ref_j);

#endif