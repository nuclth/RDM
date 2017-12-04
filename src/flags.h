#ifndef FLAGS_H
#define FLAGS_H

#include <fstream>
#include "matrix_define.h"
#include "auxiliary.h"

void init_con_values (const parameters flag_pass, FILE * spda_out, const size_t bsize, const size_t tbme_size, const size_t particles, const size_t P_num, const std::string no_flag);

void init_C_matrix (const parameters flag_pass, FILE * spda_out, const two_array & h1_mat, const two_array & h2_mat, size_t & con_count, const one_array & obme_blocks, const std::string no_flag, const std::string h2_flag, const size_t);

void init_N_flag (FILE * spda_out, const size_t bsize, size_t & con_count, const std::string no_flag);
void init_O_flag (FILE * spda_out, const size_t bsize, size_t & con_count, const std::string no_flag, const size_t obme_block_count);

void init_P_flag (FILE * sdpa_out, const size_t bsize, size_t & con_count, const size_t N, const size_t tbme_size, 
	const two_array & array_ref_tbme, const std::string, const size_t);




#endif