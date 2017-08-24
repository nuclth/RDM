#ifndef FLAGS_H
#define FLAGS_H

#include <fstream>
#include "matrix_define.h"
#include "auxiliary.h"

void init_con_values (const con_flags flag_pass, FILE * spda_out, const size_t bsize, const size_t particles);

void init_C_matrix (const con_flags flag_pass, FILE * spda_out, const two_array & h1_mat, const two_array & h2_mat, size_t & con_count);

void init_N_flag (FILE * spda_out, const size_t bsize, size_t & con_count);
void init_O_flag (FILE * spda_out, const size_t bsize, size_t & con_count);
void init_P_flag (FILE * spda_out, const size_t bsize, size_t & con_count, const size_t N);
void init_Q_flag (FILE * spda_out, const size_t bsize, size_t & con_count);
void init_G_flag (FILE * spda_out, const size_t bsize, size_t & con_count);

#endif