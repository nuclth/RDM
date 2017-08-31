#ifndef SYM_ROUTINES_H
#define SYM_ROUTINES_H

#include "auxiliary.h"


inline double F_NN_matrix (const int i, const int j, const int k, const int l)
{
	const double value = kron_del (i,k) * kron_del (j,l);

	return value;
}

inline double F_NN_matrix_A (const int i, const int j, const int k, const int l)
{
	const double value = (1./8.) * 
		( F_NN_matrix (i, j, k, l) + F_NN_matrix (k, l, i, j)
		- F_NN_matrix (j, i, k, l) - F_NN_matrix (k, l, j, i)
		- F_NN_matrix (i, j, l, k) - F_NN_matrix (l, k, i, j)
		+ F_NN_matrix (j, i, l, k) + F_NN_matrix (l, k, j, i) );

	return value;
}

inline int F3_3_matrix (const int i, const int k, const int ip, const int jp, const int kp, const int lp)
{
  const int value = kron_del (i,ip) * kron_del (k,kp) * kron_del (jp, lp);

  return value;  
}

inline int F3_3_matrix_A (const int i, const int k, const int ip, const int jp, const int kp, const int lp)
{
	const int value = 
	    F3_3_matrix (i, k, ip, jp, kp, lp) +  F3_3_matrix (i, k, kp, lp, ip, jp)
      - F3_3_matrix (i, k, jp, ip, kp, lp) -  F3_3_matrix (i, k, kp, lp, jp, ip)
      - F3_3_matrix (i, k, ip, jp, lp, kp) -  F3_3_matrix (i, k, lp, kp, ip, jp)
      + F3_3_matrix (i, k, jp, ip, lp, kp) +  F3_3_matrix (i, k, lp, kp, jp, ip);

	return value;
}

inline int F7_2_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp)
{
  const int value = 
        kron_del(j,l)*kron_del(i,ip)*kron_del(k,jp)
        -
        kron_del(j,k)*kron_del(i,ip)*kron_del(l,jp)
        -
        kron_del(i,l)*kron_del(j,ip)*kron_del(k,jp)
        +
        kron_del(i,k)*kron_del(j,ip)*kron_del(l,jp);

  return value;
}

inline int F7_3_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = kron_del (i,ip) * kron_del (j,jp) * kron_del (k,kp) * kron_del (l,lp);

  return value; 
}

inline int F7_4_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = kron_del (l,ip) * kron_del (k,jp) * kron_del (j,kp) * kron_del (i,lp);

  return value; 
}

inline int F7_3_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
    const int value =
	  F7_3_matrix (i, j, k, l, ip, jp, kp, lp) +  F7_3_matrix (i, j, k, l, kp, lp, ip, jp)
  	- F7_3_matrix (i, j, k, l, jp, ip, kp, lp) -  F7_3_matrix (i, j, k, l, kp, lp, jp, ip)
    - F7_3_matrix (i, j, k, l, ip, jp, lp, kp) -  F7_3_matrix (i, j, k, l, lp, kp, ip, jp)
    + F7_3_matrix (i, j, k, l, jp, ip, lp, kp) +  F7_3_matrix (i, j, k, l, lp, kp, jp, ip);

	return value;
}

inline int F7_4_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = 
        F7_4_matrix (i, j, k, l, ip, jp, kp, lp) +  F7_4_matrix (i, j, k, l, kp, lp, ip, jp)
      - F7_4_matrix (i, j, k, l, jp, ip, kp, lp) -  F7_4_matrix (i, j, k, l, kp, lp, jp, ip)
      - F7_4_matrix (i, j, k, l, ip, jp, lp, kp) -  F7_4_matrix (i, j, k, l, lp, kp, ip, jp)
      + F7_4_matrix (i, j, k, l, jp, ip, lp, kp) +  F7_4_matrix (i, j, k, l, lp, kp, jp, ip);

  return value;  
}


inline int F10_1_matrix (const int i, const int j, const int k, const int l, const int ip, const int kp)
{
  const int value = kron_del(j,l) * kron_del(i,ip)*kron_del(k,kp);

  return value;
}

inline int F10_3_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = 
        kron_del(i,ip) * kron_del(l,jp) * kron_del(k,kp) * kron_del(j,lp);

  return value;
}

inline int F10_5_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = 
        kron_del(i,ip) * kron_del(j,jp) * kron_del(k,kp) * kron_del(l,lp);

  return value;
}

inline int F10_3_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
    const int value =
      F10_3_matrix (i, j, k, l, ip, jp, kp, lp) +  F10_3_matrix (i, j, k, l, kp, lp, ip, jp)
    - F10_3_matrix (i, j, k, l, jp, ip, kp, lp) -  F10_3_matrix (i, j, k, l, kp, lp, jp, ip)
    - F10_3_matrix (i, j, k, l, ip, jp, lp, kp) -  F10_3_matrix (i, j, k, l, lp, kp, ip, jp)
    + F10_3_matrix (i, j, k, l, jp, ip, lp, kp) +  F10_3_matrix (i, j, k, l, lp, kp, jp, ip);

  return value;
}


inline double T2_1_matrix (const int i, const int j, const int k, const int l, const int m, const int n, const int ip, const int lp)
{
	const double value = 
		(1./2.) * kron_del (i, ip) * kron_del (l, lp) *
		( kron_del (j, m) * kron_del (k, n) - kron_del (k, m) * kron_del (j, n) 
		- kron_del (j, n) * kron_del (k, m)	+ kron_del (k, n) * kron_del (j, m) );

	return value;
}


inline double T2_3_matrix (const int i, const int j, const int k, const int l, const int m, const int n, const int ip, const int jp, const int kp, const int lp)
{
	const double value = 
		(1./2.) * kron_del (i, l) *
		( kron_del (ip, m) * kron_del (jp, n) * kron_del (kp, j) * kron_del (lp, k)
		- kron_del (ip, m) * kron_del (jp, n) * kron_del (kp, k) * kron_del (lp, j)
		- kron_del (ip, n) * kron_del (jp, m) * kron_del (kp, j) * kron_del (lp, k)
		+ kron_del (ip, n) * kron_del (jp, m) * kron_del (kp, k) * kron_del (lp, j) ) 
		-
		2. *
		( kron_del (j, m) * kron_del (ip, i) * kron_del (jp, n) * kron_del (kp, l) * kron_del (lp, k)
		- kron_del (k, m) * kron_del (ip, i) * kron_del (jp, n) * kron_del (kp, l) * kron_del (lp, j)
		- kron_del (j, n) * kron_del (ip, i) * kron_del (jp, m) * kron_del (kp, l) * kron_del (lp, k)
		+ kron_del (k, n) * kron_del (ip, i) * kron_del (jp, m) * kron_del (kp, l) * kron_del (lp, j) );

	return value;
}

inline double T2_1_matrix_S (const int i, const int j, const int k, const int l, const int m, const int n, const int ip, const int lp)
{
	const double value = (1./2.) * 
		( T2_1_matrix (i, j, k, l, m, n, ip, lp) + T2_1_matrix (i, j, k, l, m, n, lp, ip) ); 

	return value;
}


inline double T2_3_matrix_A (const int i, const int j, const int k, const int l, const int m, const int n, const int ip, const int jp, const int kp, const int lp)
{
	const double value = (1./8.) *
		( T2_3_matrix (i, j, k, l, m, n, ip, jp, kp, lp) + T2_3_matrix (i, j, k, l, m, n, kp, lp, ip, jp) 
		- T2_3_matrix (i, j, k, l, m, n, jp, ip, kp, lp) - T2_3_matrix (i, j, k, l, m, n, kp, lp, jp, ip) 
		- T2_3_matrix (i, j, k, l, m, n, ip, jp, lp, kp) - T2_3_matrix (i, j, k, l, m, n, lp, kp, ip, jp)
		+ T2_3_matrix (i, j, k, l, m, n, jp, ip, lp, kp) + T2_3_matrix (i, j, k, l, m, n, lp, kp, jp, ip) );

	return value;
}


inline double T2_6_matrix (const int i, const int j, const int k, const int l, const int m, const int n, const int ip, const int jp, const int kp, const int lp, const int mp, const int np)
{
	const double value = (-1.) *
		( kron_del (i, ip) * kron_del (j, jp) * kron_del (k, kp) * kron_del (l, lp) * kron_del (m, mp) * kron_del (n ,np) );

	return value;
}

inline double T2_6_matrix_A (const int i, const int j, const int k, const int l, const int m, const int n, const int ip, const int jp, const int kp, const int lp, const int mp, const int np)
{
	const double value = (1./8.) * 
		( T2_6_matrix (i, j, k, l, m, n, ip, jp, kp, lp, mp, np) + T2_6_matrix (i, j, k, l, m, n, lp, mp, np, ip, jp, kp)
		- T2_6_matrix (i, j, k, l, m, n, ip, kp, jp, lp, mp, np) - T2_6_matrix (i, j, k, l, m, n, lp, mp, np, ip, kp, jp) 
		- T2_6_matrix (i, j, k, l, m, n, ip, jp, kp, lp, np, mp) - T2_6_matrix (i, j, k, l, m, n, lp, np, mp, ip, jp, kp)
		+ T2_6_matrix (i, j, k, l, m, n, ip, kp, jp, lp, np, mp) + T2_6_matrix (i, j, k, l, m, n, lp, np, mp, ip, kp, jp));

		return value;
}

#endif
