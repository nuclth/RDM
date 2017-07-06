 // flag.h file

#ifndef flags_h
#define flags_h

#include <cstdio>
#include "arraydefs.h"

double kron_del(const size_t i, const size_t j);

//template <typename two_array>
void create_c_matrix (const two_array & h1_mat, const two_array & h2_mat, two_array & c_matrix, const bool two_body_toggle);

//template <typename one_array, typename three_array>
void init_N_flag (three_array & F1_con, one_array & F1_val, const size_t bsize, const size_t particles);

//template <typename three_array, typename one_array, typename four_array>
void init_O_flag (three_array & F2_con, one_array & F2_val, const size_t bsize);

//template <typename one_array, typename two_array, typename three_array, typename four_array, typename six_array>
void init_P_flag (three_array & F3_con, one_array & F3_val, const size_t bsize, const size_t N, const two_array & trans_h2);

//template <typename one_array, typename three_array, typename six_array, typename eight_array>
void init_Q_flag (three_array & F7_con, one_array & F7_val, const size_t bsize, const size_t cmat_extent, size_t & skip, const two_array & trans_h2);

//template <typename one_array, typename three_array, typename six_array, typename eight_array>
void init_G_flag (six_array & F10_build_1, eight_array & F10_build_3, eight_array & F10_build_5, three_array & F10_con, one_array & F10_val, const size_t bsize);


int F3_3_matrix (const int i, const int k, const int ip, const int jp, const int kp, const int lp);

int F3_3_matrix_A (const int i, const int k, const int ip, const int jp, const int kp, const int lp);
 
int F7_2_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp);

int F7_3_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);

int F7_5_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);

int F7_3_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);

int F7_5_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);

int F10_1_matrix (const int i, const int j, const int k, const int l, const int ip, const int kp);

int F10_3_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);

int F10_5_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);

int F10_3_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);


 #endif