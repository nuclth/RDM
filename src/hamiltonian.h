#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <string>
#include <iostream>
#include "matrix_define.h"




// function prototypes
void populate_1body (const two_array & ref_m, two_array & h1_mat, const double hw);
void fullm_populate_hamiltonian (two_array & ref_m, two_array & h1_mat, five_array & h2_mat, const std::string reference_file, const std::string matrix_file, const double hw, std::ofstream & diag_out, const bool diag_toggle);
void compactify_h2 (const two_array & ref_m, two_array & comp_h2, five_array & h2_mat, std::ofstream & diag_out, const bool diag_toggle);
void read_in_reference_m_scheme (two_array & ref_m, const std::string m_ref_file, std::ofstream & diag_out, const bool diag_toggle);
void read_in_matrix_m_scheme (const two_array & ref_m, five_array & h2_mat, const std::string m_mat_file);

double find_obme_me (const double n1, const double n2, const double l, const double j, const double mj, two_array & obme);
int get_obme_lines ();
void read_in_obme (two_array & obme);

//template <> void print (std::ostream& os, const double & x);

#endif