/*****************************************************************************
Programmer: Alex Dyhdalo
Last Mod: 8/2017
*****************************************************************************/

#include "sym_routines.h"
#include "aux.h"


int F3_3_matrix (const int i, const int k, const int ip, const int jp, const int kp, const int lp)
{
  const int value = kron_del (i,ip) * kron_del (k,kp) * kron_del (jp, lp);

  return value;  
}

int F3_3_matrix_A (const int i, const int k, const int ip, const int jp, const int kp, const int lp)
{
	const int value = 
	    F3_3_matrix (i, k, ip, jp, kp, lp) +  F3_3_matrix (i, k, kp, lp, ip, jp)
      - F3_3_matrix (i, k, jp, ip, kp, lp) -  F3_3_matrix (i, k, kp, lp, jp, ip)
      - F3_3_matrix (i, k, ip, jp, lp, kp) -  F3_3_matrix (i, k, lp, kp, ip, jp)
      + F3_3_matrix (i, k, jp, ip, lp, kp) +  F3_3_matrix (i, k, lp, kp, jp, ip);

	return value;
}




int F7_2_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp)
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

int F7_3_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = kron_del (i,ip) * kron_del (j,jp) * kron_del (k,kp) * kron_del (l,lp);

  return value;  
}

int F7_4_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = kron_del (l,ip) * kron_del (k,jp) * kron_del (j,kp) * kron_del (i,lp);

  return value;  
}

int F7_3_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
    const int value =
	    F7_3_matrix (i, j, k, l, ip, jp, kp, lp) +  F7_3_matrix (i, j, k, l, kp, lp, ip, jp)
  	- F7_3_matrix (i, j, k, l, jp, ip, kp, lp) -  F7_3_matrix (i, j, k, l, kp, lp, jp, ip)
    - F7_3_matrix (i, j, k, l, ip, jp, lp, kp) -  F7_3_matrix (i, j, k, l, lp, kp, ip, jp)
    + F7_3_matrix (i, j, k, l, jp, ip, lp, kp) +  F7_3_matrix (i, j, k, l, lp, kp, jp, ip);

	return value;
}

int F7_4_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = 
        F7_4_matrix (i, j, k, l, ip, jp, kp, lp) +  F7_4_matrix (i, j, k, l, kp, lp, ip, jp)
      - F7_4_matrix (i, j, k, l, jp, ip, kp, lp) -  F7_4_matrix (i, j, k, l, kp, lp, jp, ip)
      - F7_4_matrix (i, j, k, l, ip, jp, lp, kp) -  F7_4_matrix (i, j, k, l, lp, kp, ip, jp)
      + F7_4_matrix (i, j, k, l, jp, ip, lp, kp) +  F7_4_matrix (i, j, k, l, lp, kp, jp, ip);

  return value;  
}







int F10_1_matrix (const int i, const int j, const int k, const int l, const int ip, const int kp)
{
  const int value = kron_del(j,l) * kron_del(i,ip)*kron_del(k,kp);

  return value;

}

int F10_3_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = 
        kron_del(i,ip) * kron_del(l,jp) * kron_del(k,kp) * kron_del(j,lp);

  return value;
}

int F10_5_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = 
        kron_del(i,ip) * kron_del(j,jp) * kron_del(k,kp) * kron_del(l,lp);

  return value;
}


int F10_3_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
    const int value =
      F10_3_matrix (i, j, k, l, ip, jp, kp, lp) +  F10_3_matrix (i, j, k, l, kp, lp, ip, jp)
    - F10_3_matrix (i, j, k, l, jp, ip, kp, lp) -  F10_3_matrix (i, j, k, l, kp, lp, jp, ip)
    - F10_3_matrix (i, j, k, l, ip, jp, lp, kp) -  F10_3_matrix (i, j, k, l, lp, kp, ip, jp)
    + F10_3_matrix (i, j, k, l, jp, ip, lp, kp) +  F10_3_matrix (i, j, k, l, lp, kp, jp, ip);

  return value;
}