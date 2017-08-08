#ifndef MSCHEME_H
#define MSCHEME_H

#include <string>
#include <fstream>
#include "matrix_define.h"


void read_in_reference_m_scheme (two_array & ref_m, const std::string m_ref_file, std::ofstream & diag_out, const bool diag_toggle);
void read_in_matrix_m_scheme (const two_array & ref_m, five_array & h2_mat, const std::string m_mat_file);

#endif