#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <string>
#include <iostream>
#include "matrix_define.h"




// function prototypes
void populate_1body (const two_array & ref_m, two_array & h1_mat, const std::string, const int nmax);
void fullm_populate_hamiltonian (two_array & array_ref_obme, two_array & array_ref_tbme, two_array & h1_mat, two_array & h2_mat, const std::string ref_obme, const std::string me_obme, const std::string ref_tbme, const std::string me_tbme, const bool two_body_toggle, const int nmax, const std::string h2_flag);
void compactify_h2 (const two_array & ref_m, two_array & comp_h2, five_array & h2_mat, std::ofstream & diag_out, const bool diag_toggle);
void readin_ref_obme (two_array & ref_m, const std::string m_ref_file);
void read_in_matrix_m_scheme (const two_array & ref_m, five_array & h2_mat, const std::string m_mat_file);

double find_obme_me (const double n1, const double n2, const double l, const double j, const double mj, two_array & obme);
int get_obme_lines (const std::string);
void read_in_obme (two_array & obme, const std::string);


void readin_ref_tbme (two_array & ref_tbme, const std::string tbme_filename);
void populate_2body (const two_array & ref_m, const two_array & ref_tbme, two_array & h2_mat, const std::string me_tbme, const int nmax, const std::string h2_flag);

//template <> void print (std::ostream& os, const double & x);

template <typename Array> void print(std::ostream& os, const Array & A);

template <> void print (std::ostream& os, const double & x);

int cantor (size_t a, size_t b);

#endif