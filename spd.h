#ifndef spd_h
#define spd_h

#include <cstdio>
#include <fstream>
#include <string>
#include "arraydefs.h"
#include "auxiliary.h"

void con_matrix_out (const two_array & m_pass, size_t con_count, size_t bsize, struct con_flags flag_pass, std::ofstream & spda_out, const one_array & block_mat);

void check_constraints (const three_array & Con, const std::string & label, const bool output, const size_t num, const size_t cmat_extent);

void blockdiag_spdaout (const size_t sub_blocks, const two_array & block_mat, std::ofstream & spda_out);

size_t num_sub_blocks (size_t sub_blocks, const two_array & block_mat);

void create_spda_file (const two_array & block_mat, const one_array & oned_blocks, const two_array & c_matrix, const struct con_flags flag_pass, const three_array & F1_con, const one_array & F1_val, const three_array & F2_con, const one_array & F2_val, const three_array & F3_con, const one_array & F3_val, const three_array & F7_con, const one_array & F7_val, const three_array & F10_con, const one_array & F10_val, const size_t bsize, const size_t cmat_extent, std::ofstream & spda_out);

#endif