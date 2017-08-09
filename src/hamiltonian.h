/*****************************************************************************
Programmer: Alex Dyhdalo
Last Mod: 8/2017
*****************************************************************************/

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <string>
#include <fstream>
#include "matrix_define.h"




// function prototypes
void mpopulate_1body (const two_array & ref_m, two_array & h1_mat, const double hw);
void jpopulate_1body (const two_array & ref_j, two_array & h1_mat, const double hw);
void compactify_h2 (const two_array & ref_m, two_array & comp_h2, five_array & h2_mat, std::ofstream & diag_out, const bool diag_toggle);
void fullm_populate_hamiltonian (two_array & ref_m, two_array & h1_mat, five_array & h2_mat, const std::string reference_file, const std::string matrix_file, const double hw, std::ofstream & diag_out, const bool diag_toggle);

void fullj_populate_hamiltonian (two_array & ref_j, two_array & h1_mat, five_array & h2_mat, const std::string reference_file, const std::string matrix_file, const double hw, std::ofstream & diag_out, const bool diag_toggle);


#endif