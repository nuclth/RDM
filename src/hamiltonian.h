#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <string>
#include <iostream>
#include "matrix_define.h"




// function prototypes

void readin_ref_obme (two_array & ref_obme, const std::string obme_filename);
void readin_ref_tbme (two_array & ref_tbme, const std::string tbme_filename);


void populate_1body (const two_array & ref_m, two_array & h1_mat, const std::string, const int nmax);
void populate_2body (const two_array & ref_m, const two_array & ref_tbme, two_array & h2_mat, const std::string me_tbme, const int nmax, const std::string h2_flag);

void populate_hamiltonian (two_array & array_ref_obme, two_array & array_ref_tbme, two_array & h1_mat, two_array & h2_mat, const std::string ref_obme, const std::string me_obme, const std::string ref_tbme, const std::string me_tbme, const bool two_body_toggle, const int nmax, const std::string h2_flag);

void read_in_obme (two_array & obme, const std::string);


int cantor (size_t a, size_t b);

#endif