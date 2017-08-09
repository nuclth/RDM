/*****************************************************************************
Programmer: Alex Dyhdalo
Last Mod: 8/2017
*****************************************************************************/

#ifndef SYM_ROUTINES_H
#define SYM_ROUTINES_H


int F3_3_matrix   (const int i, const int k, const int ip, const int jp, const int kp, const int lp);
int F3_3_matrix_S (const int i, const int k, const int ip, const int jp, const int kp, const int lp);
int F3_3_matrix_A (const int i, const int k, const int ip, const int jp, const int kp, const int lp);


int F7_2_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp);

int F7_3_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);

int F7_4_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);

int F7_3_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);

int F7_4_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);


int F10_1_matrix (const int i, const int j, const int k, const int l, const int ip, const int kp);

int F10_3_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);

int F10_5_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);

int F10_3_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp);

#endif