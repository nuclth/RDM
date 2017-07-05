
#include "flags.h"
#include <gsl/gsl_blas.h>
#include <iostream>

/***************************************************************

Function to create the F0 constraint matrix for the SDP solver.
This matrix is just the Hamiltonian matrix elements for the 1RDM
and 2RDM parts (where applicable) with zeros in the other parts
of the matrix.

Note that since the SDP solver maximizes the objective function,
here we multiply all elements by -1 so that we are maximizing
the negative energy.

***************************************************************/

//template <typename two_array>
void create_c_matrix (const two_array & h1_mat, const two_array & h2_mat, two_array & c_matrix, const bool two_body_toggle)
{
  size_t h1_len = h1_mat.size();
  size_t h2_len = h2_mat.size();

  for (size_t i = 0; i < h1_len; i++) 
  {
    for (size_t j = 0; j < h1_len; j++)
    {
      c_matrix[i][j] = h1_mat [i][j] * -1.0;
    }

  }

  if (two_body_toggle)
  {
    for (size_t i = 0; i < h2_len; i++) 
    {
    for (size_t j = 0; j < h2_len; j++)
    {
      size_t ic = i + 2*h1_len;
      size_t jc = j + 2*h1_len;

// CHANGE
//      c_matrix[ic][jc] = h2_mat [i][j] * -1.0;
        c_matrix[ic][jc] = h2_mat [i][j] * -1./2.;
    }
    }
  }


//  print(std::cout, c_matrix);
}

/***************************************************************

Kronecker delta function. Returns 1 if i=j and 0 otherwise.

***************************************************************/

double kron_del(const size_t i, const size_t j)
{

  if (i == j)
    return 1.;

  return 0.;
}


/***************************************************************

Function to create our F1 constraint matrix to fix the total
particle number N (fixes trace of the 1RDM to be N). 

In keeping with other constraint functions, first populate the
F1_con_1 matrix (only one constraint) then populate the full F1_con
matrix. 

Note that F1_con_1 is a basis x basis size matrix while F1_con 
matches the dimensions of the F0 constraint matrix.

***************************************************************/

//template <typename one_array, typename three_array>
void init_N_flag (three_array & F1_con, one_array & F1_val, const size_t bsize, const size_t particles)
{

//  size_t cmat_extent = F1_con.size();


  for (size_t i = 0; i < bsize; i++)
  {
    if (i < bsize)
      F1_con[0][i][i] = 1.;
  }

  F1_val[0] = particles;

}


/***************************************************************

Function to create our F2 constraint matrices to fix the relation
between the 1RDM and its twin q. 

First populate the F2_con_1 matrix which serves as the basis x basis
size constraint matrix for the 1RDM sector and the q sector (both 
have the same constraint matrix).

Then populate this matrix into the full F0 sized constraint term 
F2_con. 

***************************************************************/


//template <typename three_array, typename one_array, typename four_array>
void init_O_flag (four_array & F2_con_1, three_array & F2_con, one_array & F2_val, const size_t bsize)
{


//  size_t F2num = F2_con.size();           // number of F2 constraint matrices
//  size_t bsize = F2_con_1.size();         // basis size

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
    for (size_t k = 0; k < bsize; k++)    // loop over matrix row
    {
    for (size_t l = 0; l < bsize; l++)    // loop over matrix column
    {

        F2_con_1[i][j][k][l] = (1./2.)*(kron_del(i,k)*kron_del(j,l) + kron_del(i,l)*kron_del(j,k));
    }
    }

  }
  }

  size_t counter = 0;

  for (size_t i1 = 0;  i1 < bsize; i1++)
  {
  for (size_t i2 = i1; i2 < bsize; i2++)
  {
    for (size_t j = 0; j < bsize; j++)
    {
    for (size_t k = 0; k < bsize; k++)
    {
      size_t jq = j + bsize;
      size_t kq = k + bsize;

      size_t i = counter;

      F2_con[i][j][k]   = F2_con_1 [i1][i2][j][k];
      F2_con[i][jq][kq] = F2_con_1 [i1][i2][j][k];
    }
    }

    if (i1 == i2)
      F2_val[counter] = 1.;
    else
      F2_val[counter] = 0.;

    counter++;
  }
  }

//  printf("%i", F2num);
//  printf("\n");
//  printf("%i", bsize);
//  printf("\n");
}

int F3_3_matrix (const int i, const int k, const int ip, const int jp, const int kp, const int lp)
{
  const int value = kron_del (i,ip) * kron_del (k,kp) * kron_del (jp, lp);

  return value;  
}

/***************************************************************

Function to create our F3 constraint matrix to fix the relation
between the 2RDM partial trace and the 1RDM.

***************************************************************/

//template <typename one_array, typename two_array, typename three_array, typename four_array, typename six_array>
void init_P_flag (four_array & F3_build_1, six_array & F3_build_3, three_array & F3_con, one_array & F3_val, const size_t bsize, const size_t N, const two_array & trans_h2)
{

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t k = i; k < bsize; k++)      // loop over jth constraint matrix
  {
    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
      F3_build_1 [i][k][ip][kp] = -1.0 * (N - 1.0) / 4.0 * (
        kron_del(i,ip)*kron_del(k,kp) + kron_del(k,ip)*kron_del(i,kp)
        );
    }
    }
  }
  }

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t k = i; k < bsize; k++)      // loop over jth constraint matrix
  {
    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
     F3_build_3 [i][k][ip][jp][kp][lp] = 1./8. * (
        F3_3_matrix (i, k, ip, jp, kp, lp) +  F3_3_matrix (i, k, kp, lp, ip, jp)

      - F3_3_matrix (i, k, jp, ip, kp, lp) -  F3_3_matrix (i, k, kp, lp, jp, ip)

      - F3_3_matrix (i, k, ip, jp, lp, kp) -  F3_3_matrix (i, k, lp, kp, ip, jp)

      + F3_3_matrix (i, k, jp, ip, lp, kp) +  F3_3_matrix (i, k, lp, kp, jp, ip)        
      );

//     	std::cout << "i k; ip jp kp lp: ";
//     	std::cout << i << " " << k << "\t" << ip << " " << jp << " " << kp << " " << lp;
//     	std::cout << "\t" << F3_build_3 [i][k][ip][jp][kp][lp] << std::endl;
    }
    }
    }
    }
  }
  }


  size_t counter = 0;

  for (size_t i1 = 0; i1 < bsize; i1++)
  {
  for (size_t i2 = i1; i2 < bsize; i2++)
  {
    for (size_t j = 0; j < bsize; j++)
    {
    for (size_t k = 0; k < bsize; k++)
    {

      size_t b = counter;

      F3_con[b][j][k] = F3_build_1 [i1][i2][j][k];


    }
    }

    counter++;
  }
  }


  counter = 0;

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t k = i; k < bsize; k++)      // loop over jth constraint matrix
  {
    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
      if (ip >= jp)
        continue;

    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
      if (kp >= lp)
        continue;


      size_t b = counter;

      size_t ips = ip + 1;
      size_t jps = jp + 1;
      size_t kps = kp + 1;
      size_t lps = lp + 1;

      size_t left  = jps - ips + (2*bsize - ips) * (ips - 1)/2 - 1;
      size_t right = lps - kps + (2*bsize - kps) * (kps - 1)/2 - 1;

//      std::cout << ips << " " << jps << " " << kps << " " << lps << std::endl;

//      std::cout << left << "\t" << right << std::endl;

      size_t left_P  = left  + 2*bsize;
      size_t right_P = right + 2*bsize;


      F3_con[b][left_P][right_P] = F3_build_3 [i][k][ip][jp][kp][lp];


//      std::cout << "SUCCESS" << std::endl;
    }
    }
    }
    }

    F3_val[counter] = 0.;

    counter++;
  }
  }



  for (size_t a = 0; a < counter; a++)
  {
	  const size_t extent = bsize * (bsize-1) / 2;
	  double perm_data   [extent*extent];
	  double flag_data   [extent*extent];
	  double inter_data  [extent*extent];
	  double answer_data [extent*extent];

	  size_t q = 0;

	  for (size_t i = 0; i < extent; i++)
	  {
	  for (size_t j = 0; j < extent; j++)
	  {

	  		size_t l = i + 2 * bsize;
	  		size_t r = j + 2 * bsize;

			perm_data   [q] = trans_h2 [i][j];
			flag_data   [q] = F3_con[a][l][r];
			inter_data  [q] = 0.0;
			answer_data [q] = 0.0;
			q++;
	  }
	  }

	  gsl_matrix_view p  = gsl_matrix_view_array (perm_data, extent, extent);
	  gsl_matrix_view f  = gsl_matrix_view_array (flag_data, extent, extent);

	  gsl_matrix_view inter  = gsl_matrix_view_array (inter_data,  extent, extent);
	  gsl_matrix_view answer = gsl_matrix_view_array (answer_data, extent, extent);


	  gsl_blas_dgemm (CblasTrans,   CblasNoTrans, 1.0, &p.matrix, &f.matrix, 0.0, &inter.matrix);
	  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &inter.matrix, &p.matrix, 0.0, &answer.matrix);




	  q = 0;

	  for (size_t i = 0; i < extent; i++)
	  {
	  for (size_t j = 0; j < extent; j++)
	  {
			double value = answer_data [q];
			q++;

	  		size_t l = i + 2 * bsize;
	  		size_t r = j + 2 * bsize;

			F3_con[a][l][r] = value;
	  }
	  } 
    
  }

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

int F7_5_matrix (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
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

int F7_5_matrix_A (const int i, const int j, const int k, const int l, const int ip, const int jp, const int kp, const int lp)
{
  const int value = 
        F7_5_matrix (i, j, k, l, ip, jp, kp, lp) +  F7_5_matrix (i, j, k, l, kp, lp, ip, jp)
      - F7_5_matrix (i, j, k, l, jp, ip, kp, lp) -  F7_5_matrix (i, j, k, l, kp, lp, jp, ip)
      - F7_5_matrix (i, j, k, l, ip, jp, lp, kp) -  F7_5_matrix (i, j, k, l, lp, kp, ip, jp)
      + F7_5_matrix (i, j, k, l, jp, ip, lp, kp) +  F7_5_matrix (i, j, k, l, lp, kp, jp, ip);

  return value;  
}

/***************************************************************

Function to create our F7 constraint matrix to fix the linear 
relations between the 1RDM, q, 2RDM, and Q. 

***************************************************************/

//template <typename one_array, typename three_array, typename six_array, typename eight_array>
void init_Q_flag (six_array & F7_build_2, eight_array & F7_build_3, eight_array & F7_build_5, three_array & F7_con, one_array & F7_val, const size_t bsize, const size_t cmat_extent, size_t & skip, const two_array & trans_h2)
{


  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = i+1; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = i; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = k+1; l < bsize; l++)    // loop over matrix column
  {
    if (j < l && k <= i)
      continue;

    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
        F7_build_2 [i][j][k][l][ip][jp] = (1./2.) * (
          F7_2_matrix (i, j, k, l, ip, jp) + F7_2_matrix (i, j, k, l, jp, ip)
        );


/*        std::cout << "i j k l; ip jp: ";
     	std::cout << i << " " << j << " " << k << " " << l;
     	std::cout << "\t";
     	std::cout << ip << " " << jp;
     	std::cout << "\t" << F7_build_2 [i][j][k][l][ip][jp] << std::endl;
*/
    }
    }
   }
   }
   }
   } 


//  std::cout << std::endl << std::endl;

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = i+1; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = i; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = k+1; l < bsize; l++)    // loop over matrix column
  {
    if (j < l && k <= i)
      continue;

    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
//    	if (ip >= jp)
//    		continue;

    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
//    	if (kp >= lp)
//    		continue;


    	F7_build_3 [i][j][k][l][ip][jp][kp][lp] = (1./4.) *
    	F7_3_matrix_A (i, j, k, l, ip, jp, kp, lp);
    	F7_build_5 [i][j][k][l][ip][jp][kp][lp] = (-1./8.) *
    	F7_5_matrix_A (i, j, k, l, ip, jp, kp, lp);


/*     	std::cout << "i j k l; ip jp kp lp: ";
     	std::cout << i << " " << j << " " << k << " " << l;
     	std::cout << "\t";
     	std::cout << ip << " " << jp << " " << kp << " " << lp;
     	std::cout << "\t" << F7_build_3 [i][j][k][l][ip][jp][kp][lp] << std::endl;
*/
//    	if (F7_build_3 [i][j][k][l][ip][jp][kp][lp] != 0)
//    			std::cout << "ip jp kp lp" << "\t" << ip << " " << jp << " " << kp << " " << lp << std::endl;

//    	F7_build_3 [i][j][k][l][ip][jp][kp][lp] = 
//    	F7_3_matrix (i, j, k, l, ip, jp, kp, lp) + F7_3_matrix (i, j, k, l, kp, lp, ip, jp);

//    	F7_build_5 [i][j][k][l][ip][jp][kp][lp] = F7_5_matrix (i, j, k, l, ip, jp, kp, lp);
/*      F7_build_3 [i][j][k][l][ip][jp][kp][lp] = (1./32.) * (
      	  F7_3_matrix_A (i, j, k, l, ip, jp, kp, lp) + F7_3_matrix_A (k, l, i, j, ip, jp, kp, lp)
      	- F7_3_matrix_A (j, i, k, l, ip, jp, kp, lp) - F7_3_matrix_A (k, l, j, i, ip, jp, kp, lp)
      	- F7_3_matrix_A (i, j, l, k, ip, jp, kp, lp) - F7_3_matrix_A (l, k, i, j, ip, jp, kp, lp)
      	+ F7_3_matrix_A (j, i, l, k, ip, jp, kp, lp) + F7_3_matrix_A (l, k, j, i, ip, jp, kp, lp)
      );

      F7_build_5 [i][j][k][l][ip][jp][kp][lp] = (1./64.) * (
      	  F7_5_matrix_A (i, j, k, l, ip, jp, kp, lp) + F7_5_matrix_A (k, l, i, j, ip, jp, kp, lp)
      	- F7_5_matrix_A (j, i, k, l, ip, jp, kp, lp) - F7_5_matrix_A (k, l, j, i, ip, jp, kp, lp)
      	- F7_5_matrix_A (i, j, l, k, ip, jp, kp, lp) - F7_5_matrix_A (l, k, i, j, ip, jp, kp, lp)
      	+ F7_5_matrix_A (j, i, l, k, ip, jp, kp, lp) + F7_5_matrix_A (l, k, j, i, ip, jp, kp, lp)
      );
*/
    }
    }
    }
    }

//    std::cout << "i j k l" << "\t" << i << " " << j << " " << k << " " << l << std::endl;
//    print (std::cout, F7_build_3[i][j][k][l]);
//    std::cout << std::endl;

  }
  }
  }
  }




  size_t counter = 0;

  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = i+1; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = i; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = k+1; l < bsize; l++)    // loop over matrix column
  {
    if (j < l && k <= i)
      continue;
/*
    if (j == l && k < i)
      continue;

    if (i == k && l < j)
      continue;

    if (i > k and j > l)
      continue;
*/
    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
      if (ip >= jp)
        continue;

    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
      if (kp >= lp)
        continue;

      size_t b = counter;

      size_t ips = ip + 1;
      size_t jps = jp + 1;
      size_t kps = kp + 1;
      size_t lps = lp + 1;

      size_t left  = jps - ips + (2*bsize - ips) * (ips - 1)/2 - 1;
      size_t right = lps - kps + (2*bsize - kps) * (kps - 1)/2 - 1;

      size_t left_P  = left  + 2*bsize;
      size_t right_P = right + 2*bsize;

      size_t left_Q  = left  + 2*bsize + bsize*(bsize-1)/2;
      size_t right_Q = right + 2*bsize + bsize*(bsize-1)/2;

      F7_con[b][left_P][right_P] = F7_build_3 [i][j][k][l][ip][jp][kp][lp];
      F7_con[b][left_Q][right_Q] = F7_build_5 [i][j][k][l][ip][jp][kp][lp];


    }
    }
    }
    }

    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
      size_t b = counter;

      size_t left_q  = ip + bsize;
      size_t right_q = jp + bsize;

      F7_con[b][left_q][right_q] = F7_build_2 [i][j][k][l][ip][jp];
    }
    }

    F7_val [counter] = 
    kron_del(i,k)*kron_del(j,l) - kron_del(i,l)*kron_del(j,k);

    bool mat_used = false;
    bool val_used = false;

    for (size_t l1 = 0;  l1 < cmat_extent; l1++)
    {
    for (size_t l2 = l1; l2 < cmat_extent; l2++)
    {
      if (F7_con[counter][l1][l2] != 0.)
        mat_used = true;
    }
    }

    if (F7_val [counter] != 0.)
      val_used = true;


    if (!mat_used and !val_used)
    {
      skip++;
      continue;
    }


    if (!mat_used and val_used)
    {
      std::cout << "ERROR: CONSTRAINT MATRIX ZERO AND VALUE NON-ZERO" << std::endl;
//      print (std::cout, F7_con[counter]);
      std::cout << F7_val [counter] << std::endl;
      std::cout << "END ERROR" << std::endl;
    }

//    std::cout << "i j k l" << "\t" << i << " " << j << " " << k << " " << l << std::endl;
//    print (std::cout, F7_con[counter]);
//    std::cout << F7_val [counter] << std::endl;


    counter++;
  }
  }
  }
  }
  std::cout << "COUNTER = " << counter << std::endl;
  std::cout << "SKIP " << skip << " MATRICES" << std::endl;


  for (size_t a = 0; a < counter; a++)
  {
    const size_t extent = bsize * (bsize-1) / 2;

    double perm_data     [extent*extent];

    double flag_P_data   [extent*extent];
    double flag_Q_data   [extent*extent];

    double inter_P_data  [extent*extent];
    double answer_P_data [extent*extent];

    double inter_Q_data  [extent*extent];
    double answer_Q_data [extent*extent];

    size_t q = 0;

    for (size_t i = 0; i < extent; i++)
    {
    for (size_t j = 0; j < extent; j++)
    {

      size_t lP = i + 2 * bsize;
      size_t rP = j + 2 * bsize;

      size_t lQ = i + 2 * bsize + bsize*(bsize-1)/2;
      size_t rQ = j + 2 * bsize + bsize*(bsize-1)/2;

      perm_data     [q] = trans_h2 [i][j];

      flag_P_data   [q] = F7_con[a][lP][rP];
      flag_Q_data   [q] = F7_con[a][lQ][rQ];

      inter_P_data  [q] = 0.0;
      answer_P_data [q] = 0.0;

      inter_Q_data  [q] = 0.0;
      answer_Q_data [q] = 0.0;
      q++;
    }
    }

    gsl_matrix_view p  = gsl_matrix_view_array (perm_data, extent, extent);
    gsl_matrix_view f1 = gsl_matrix_view_array (flag_P_data, extent, extent);
    gsl_matrix_view f2 = gsl_matrix_view_array (flag_Q_data, extent, extent);

    gsl_matrix_view inter_P  = gsl_matrix_view_array (inter_P_data,  extent, extent);
    gsl_matrix_view answer_P = gsl_matrix_view_array (answer_P_data, extent, extent);

    gsl_matrix_view inter_Q  = gsl_matrix_view_array (inter_Q_data,  extent, extent);
    gsl_matrix_view answer_Q = gsl_matrix_view_array (answer_Q_data, extent, extent);


    gsl_blas_dgemm (CblasTrans,   CblasNoTrans, 1.0, &p.matrix, &f1.matrix, 0.0, &inter_P.matrix);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &inter_P.matrix, &p.matrix, 0.0, &answer_P.matrix);

    gsl_blas_dgemm (CblasTrans,   CblasNoTrans, 1.0, &p.matrix, &f2.matrix, 0.0, &inter_Q.matrix);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &inter_Q.matrix, &p.matrix, 0.0, &answer_Q.matrix);


    q = 0;

    for (size_t i = 0; i < extent; i++)
    {
    for (size_t j = 0; j < extent; j++)
    {
      double value_P = answer_P_data [q];
      double value_Q = answer_Q_data [q];
      q++;

      size_t lP = i + 2 * bsize;
      size_t rP = j + 2 * bsize;

      size_t lQ = i + 2 * bsize + bsize*(bsize-1)/2;
      size_t rQ = j + 2 * bsize + bsize*(bsize-1)/2;

      F7_con[a][lP][rP] = value_P;
      F7_con[a][lQ][rQ] = value_Q;
    }
    } 
    
  }

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

/***************************************************************

Function to create our F10 constraint matrix to fix the G linear
relations.

***************************************************************/

//template <typename one_array, typename three_array, typename six_array, typename eight_array>
void init_G_flag (six_array & F10_build_1, eight_array & F10_build_3, eight_array & F10_build_5, three_array & F10_con, one_array & F10_val, const size_t bsize)
{



  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = j; l < bsize; l++)    // loop over matrix column
  {
    
    if (j == l && k < i)
      continue;

    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t kp = 0; kp < bsize; kp++)      // loop over jth constraint matrix
    {

      F10_build_1 [i][j][k][l][ip][kp] =  (-1./2.) * (
        F10_1_matrix (i, j, k, l, ip, kp) + F10_1_matrix (i, j, k, l, kp, ip)
        );
  
    }
    }
  }
  }
  }
  }


  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = j; l < bsize; l++)    // loop over matrix column
  {
    
    if (j == l && k < i)
      continue;

    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
      F10_build_3 [i][j][k][l][ip][jp][kp][lp] =  (1./4.) * F10_3_matrix_A (i, j, k, l, ip, jp, kp, lp);

      F10_build_5 [i][j][k][l][ip][jp][kp][lp] =  (1./2.) * (
        F10_5_matrix (i, j, k, l, ip, jp, kp, lp) + F10_5_matrix (i, j, k, l, kp, lp, ip, jp)
        );
    }
    }
    }
    }
  }
  }
  }
  }



  size_t counter = 0;


  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = j; l < bsize; l++)    // loop over matrix column
  {

    if (j == l && k < i)
      continue;

    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {

      size_t b = counter;

      size_t left_p  = ip;
      size_t right_p = jp;

      F10_con[b][left_p][right_p] = F10_build_1 [i][j][k][l][ip][jp];

    }
    }

    counter++;
  }
  }
  }
  }




  counter = 0;


  for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = 0; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = 0; k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = j; l < bsize; l++)    // loop over matrix column
  {
    
    if (j == l && k < i)
      continue;

    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
      size_t b = counter;

      size_t ips = ip + 1;
      size_t jps = jp + 1;
      size_t kps = kp + 1;
      size_t lps = lp + 1;

      size_t left  = jps - ips + (2*bsize - ips) * (ips - 1)/2 - 1;
      size_t right = lps - kps + (2*bsize - kps) * (kps - 1)/2 - 1;

      size_t left_P  = left  + 2*bsize;
      size_t right_P = right + 2*bsize;

      size_t lG = ip * bsize + jp;
      size_t rG = kp * bsize + lp;

      size_t left_G  = lG + 2*bsize + 2*bsize*(bsize-1)/2;
      size_t right_G = rG + 2*bsize + 2*bsize*(bsize-1)/2;

      if (jp > ip and lp > kp)
        F10_con[b][left_P][right_P] = F10_build_3 [i][j][k][l][ip][jp][kp][lp];
      

      F10_con[b][left_G][right_G] = F10_build_5 [i][j][k][l][ip][jp][kp][lp];











    }
    }
    }
    }

    F10_val [counter] = 0.;

//    std::cout << "i j k l" << "\t" << i << " " << j << " " << k << " " << l << std::endl;
//    print (std::cout, F10_con[counter]);
//    std::cout << F10_val [counter] << std::endl;

    counter++;
  }
  }
  }
  }

  std::cout << "COUNTER = " << counter << std::endl;

}

