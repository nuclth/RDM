/*****************************************************************************
Programmer: Alex Dyhdalo
Last Mod: 8/2017
*****************************************************************************/

#include "flags.h"
#include "sym_routines.h"
#include "aux.h"
#include <iomanip>
#include <cmath>



/***************************************************************

Function to create the constraint matrix values for the SDP solver.

These are the real numbers on the right hand side of the SDP.

***************************************************************/



void init_con_values (const con_flags flag_pass, std::ofstream & spda_out, const size_t bsize, const size_t particles)
{

//  size_t cmat_extent = F1_con.size();


  // F1 Flag - N Trace condition
  if (flag_pass.N_flag)
    spda_out << particles << " ";

  // F2 Flag - Linear p q relations
  if (flag_pass.O_flag)
  {
    for (size_t i1 = 0;  i1 < bsize; i1++)
    {
    for (size_t i2 = i1; i2 < bsize; i2++)
    {

      double degen1 = sqrt(degen(i1+1));
      double degen2 = sqrt(degen(i2+1));

      if (i1 == i2)
        spda_out << degen1*degen2 << " ";
      else
        spda_out << 0. << " ";


    }
    }
   }


  // F3 Flag - P and p trace relation
   if (flag_pass.P_flag)
   {
      for (size_t i = 0; i < bsize; i++)      // loop over ith constraint matrix
      {
      for (size_t k = i; k < bsize; k++)      // loop over jth constraint matrix
      {

        spda_out << 0. << " ";

      }
      }
   }


  // F7 Flag - Q relations
   if (flag_pass.Q_flag)
   {
      for (size_t i = 0;   i < bsize; i++)      // loop over ith constraint matrix
      {
      for (size_t j = i+1; j < bsize; j++)      // loop over jth constraint matrix
      {
      for (size_t k = i;   k < bsize; k++)    // loop over matrix row
      {
      for (size_t l = k+1; l < bsize; l++)    // loop over matrix column
      {
        if (j < l && k <= i)
          continue;

      double val = kron_del(i,k)*kron_del(j,l) - kron_del(i,l)*kron_del(j,k);

        spda_out << val << " ";

      }
      }
      }
      }

   }


   // F10 Flag - G relations
   if (flag_pass.G_flag)
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

        spda_out << 0. << " ";

      }
      }
      }
      }

   }


  spda_out << std::endl;

}



/***************************************************************

Function to create the F0 constraint matrix for the SDP solver.
This matrix is just the Hamiltonian matrix elements for the 1RDM
and 2RDM parts (where applicable) with zeros in the other parts
of the matrix.

Note that since the SDP solver maximizes the objective function,
here we multiply all elements by -1 so that we are maximizing
the negative energy.

***************************************************************/


void init_C_matrix (const con_flags flag_pass, std::ofstream & spda_out, const two_array & h1_mat, const two_array & h2_mat, size_t & con_count)
{
  size_t h1_len = h1_mat.size();
  size_t h2_len = h2_mat.size();

    spda_out << std::setprecision(16);

    for (size_t ip = 0;  ip < h1_len; ip++)
    {
    for (size_t jp = ip; jp < h1_len; jp++)
    {



      size_t n = ip + 1;
      size_t m = jp + 1; 

      double degen1 = sqrt(degen(n));
      double degen2 = sqrt(degen(m));

      double val1 = h1_mat [ip][jp] * -1.0 * degen1 * degen2;

      if (val1 != 0.)
        spda_out << con_count << " " << 1 << " " << n << " " << m << " " << val1 << "\n";
    }
    }





  if (flag_pass.two_body_toggle)
  {
      for (size_t ip = 0; ip < h2_len; ip++)      // loop over ith constraint matrix
      {
      for (size_t jp = ip; jp < h2_len; jp++)      // loop over jth constraint matrix
      {

        double val3 = h2_mat [ip][jp] * -1./2.;

        size_t n = ip + 1;
        size_t m = jp + 1; 

        if (val3 != 0.)
          spda_out << con_count << " " << 3 << " " << n << " " << m << " " << val3 << "\n";

      }
      }
  }

  con_count++;

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

void init_N_flag (std::ofstream & spda_out, const size_t bsize, size_t & con_count)
{


    for (size_t ip = 0;  ip < bsize; ip++)
    {
    for (size_t kp = ip; kp < bsize; kp++)
    {

      size_t n = ip + 1;
      size_t m = kp + 1; 

      double degen1 = sqrt(degen(n));
      double degen2 = sqrt(degen(m));

      double val1 = kron_del (ip, kp) * degen1 * degen2;



      if (val1 != 0.)
      	spda_out << con_count << " " << 1 << " " << n << " " << m << " " << val1 << "\n";


    }
    }


    con_count++;
}



/***************************************************************

Function to create our O constraint matrices to fix the relation
between the 1RDM and its twin q. 

First populate the F2_con_1 matrix which serves as the basis x basis
size constraint matrix for the 1RDM sector and the q sector (both 
have the same constraint matrix).

Then populate this matrix into the full F0 sized constraint term 
F2_con. 

***************************************************************/

void init_O_flag (std::ofstream & spda_out, const size_t bsize, size_t & con_count)
{

  std::ios_base::sync_with_stdio(false);


  for (size_t i = 0; i < bsize; i++)
  {
  for (size_t j = i; j < bsize; j++)
  {

    for (size_t k = 0; k < bsize; k++)
    {
    for (size_t l = k; l < bsize; l++)
    {

      size_t n = k + 1;
      size_t m = l + 1; 

      double degen1 = sqrt(degen(n));
      double degen2 = sqrt(degen(m));

      double val1 = (1./2.)*(kron_del(i,k)*kron_del(j,l) + kron_del(i,l)*kron_del(j,k)) * degen1 * degen2;

      if (val1 != 0.)
        spda_out << con_count << " " << 1 << " " << n << " " << m << " " << val1 << "\n";


    }
    }



    for (size_t k = 0; k < bsize; k++)
    {
    for (size_t l = k; l < bsize; l++)
    {

      size_t n = k + 1;
      size_t m = l + 1; 

      double degen1 = sqrt(degen(n));
      double degen2 = sqrt(degen(m));

      double val2 = (1./2.)*(kron_del(i,k)*kron_del(j,l) + kron_del(i,l)*kron_del(j,k)) * degen1 * degen2;

      if (val2 != 0.)
        spda_out << con_count << " " << 2 << " " << n << " " << m << " " << val2 << "\n";


    }
    }



    con_count++;

  }
  }


}



/***************************************************************

Function to create our P constraint matrix to fix the relation
between the 2RDM partial trace and the 1RDM.

***************************************************************/


void init_P_flag (std::ofstream & spda_out, const size_t bsize, size_t & con_count, const size_t N)
{

  std::ios_base::sync_with_stdio(false);


  for (size_t i = 0; i < bsize; i++)
  {
  for (size_t k = i; k < bsize; k++)
  {



    for (size_t ip = 0; ip < bsize; ip++)
    {
    for (size_t kp = ip; kp < bsize; kp++)
    {

      size_t n = ip + 1;
      size_t m = kp + 1; 

      double degen1 = sqrt (degen(n));
      double degen2 = sqrt (degen(m));

      double val1 = 
      -1.0 * (N - 1.0) / 4.0 * (
        kron_del(i,ip)*kron_del(k,kp) + kron_del(k,ip)*kron_del(i,kp)
        ) * degen1 * degen2;

      if (val1 != 0.)
        spda_out << con_count << " " << 1 << " " << n << " " << m << " " << val1 << "\n";



    }
    }



    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
//      if (ip >= jp)
//        continue;

    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
//      if (kp >= lp)
//        continue;

      size_t ips = ip + 1;
      size_t jps = jp + 1;
      size_t kps = kp + 1;
      size_t lps = lp + 1;

      double degen1 = sqrt(degen (ips));
      double degen2 = sqrt(degen (jps));
      double degen3 = sqrt(degen (kps));
      double degen4 = sqrt(degen (lps));

//      size_t n = jps - ips + (2*bsize - ips) * (ips - 1)/2;
//      size_t m = lps - kps + (2*bsize - kps) * (kps - 1)/2;

    size_t n = ip * bsize + jp + 1;
    size_t m = kp * bsize + lp + 1;

//      double val3 = 1./8. * F3_3_matrix_A (i, k, ip, jp, kp, lp);

      double val3 = 1./2. * F3_3_matrix_S (i, k, ip, jp, kp, lp) * degen1 * degen2 * degen3 * degen4;

      if (val3 != 0. and n <= m)
        spda_out << con_count << " " << 3 << " " << n << " " << m << " " << val3 << "\n";


    }
    }
    }
    }

    con_count++;

  }
  }



}


/***************************************************************

Function to create our F7 constraint matrix to fix the linear 
relations between the 1RDM, q, 2RDM, and Q. 

***************************************************************/

void init_Q_flag (std::ofstream & spda_out, const size_t bsize, size_t & con_count)
{


  for (size_t i = 0;   i < bsize; i++)      // loop over ith constraint matrix
  {
  for (size_t j = i+1; j < bsize; j++)      // loop over jth constraint matrix
  {
  for (size_t k = i;   k < bsize; k++)    // loop over matrix row
  {
  for (size_t l = k+1; l < bsize; l++)    // loop over matrix column
  {
    if (j < l && k <= i)
      continue;

    
    for (size_t ip = 0;  ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = ip; jp < bsize; jp++)      // loop over jth constraint matrix
    {
      double val2 = 
      (1./2.) * (F7_2_matrix (i, j, k, l, ip, jp) + F7_2_matrix (i, j, k, l, jp, ip));

      size_t n = ip + 1;
        size_t m = jp + 1; 

        if (val2 != 0.)
          spda_out << con_count << " " << 2 << " " << n << " " << m << " " << val2 << "\n";
    }
  }

    for (size_t ip = 0;    ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = ip+1; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0;   kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = kp+1; lp < bsize; lp++)    // loop over matrix column
    {
  //    if (jp < lp && kp <= ip)
  //        continue;


      double val3 = (1./4.) * F7_3_matrix_A (i, j, k, l, ip, jp, kp, lp);


        size_t ips = ip + 1;
        size_t jps = jp + 1;
        size_t kps = kp + 1;
        size_t lps = lp + 1;

        size_t n = jps - ips + (2*bsize - ips) * (ips - 1)/2;
        size_t m = lps - kps + (2*bsize - kps) * (kps - 1)/2;


        if (val3 != 0. and n <= m)
          spda_out << con_count << " " << 3 << " " << n << " " << m << " " << val3 << "\n";


  }
  }
  }
  }


    for (size_t ip = 0;    ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = ip+1; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0;   kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = kp+1; lp < bsize; lp++)    // loop over matrix column
    {
//      if (jp < lp && kp <= ip)
//          continue;

      double val4 = (-1./8.) * F7_4_matrix_A (i, j, k, l, ip, jp, kp, lp);

      size_t ips = ip + 1;
        size_t jps = jp + 1;
        size_t kps = kp + 1;
        size_t lps = lp + 1;

        size_t n = jps - ips + (2*bsize - ips) * (ips - 1)/2;
        size_t m = lps - kps + (2*bsize - kps) * (kps - 1)/2;


        if (val4 != 0. and n <= m)
          spda_out << con_count << " " << 4 << " " << n << " " << m << " " << val4 << "\n";


  }
  }
  }
  }



  con_count++;

  } 
  }
  }
  }



}



/***************************************************************

Function to create our F10 constraint matrix to fix the G linear
relations.

***************************************************************/

void init_G_flag (std::ofstream & spda_out, const size_t bsize, size_t & con_count)
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

    for (size_t ip = 0;  ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t kp = ip; kp < bsize; kp++)      // loop over jth constraint matrix
    {


      double val1 =  
      (-1./2.) * (F10_1_matrix (i, j, k, l, ip, kp) + F10_1_matrix (i, j, k, l, kp, ip));

      size_t n = ip + 1;
      size_t m = kp + 1; 

      if (val1 != 0.)
        spda_out << con_count << " " << 1 << " " << n << " " << m << " " << val1 << "\n";

  
    }
    }


    for (size_t ip = 0;    ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = ip+1; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0;   kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = kp+1; lp < bsize; lp++)    // loop over matrix column
    {
//      if (jp < lp && kp <= ip)
//          continue;
        
        double val3 = (1./4.) * F10_3_matrix_A (i, j, k, l, ip, jp, kp, lp);

        size_t ips = ip + 1;
        size_t jps = jp + 1;
        size_t kps = kp + 1;
        size_t lps = lp + 1;

        size_t n = jps - ips + (2*bsize - ips) * (ips - 1)/2;
        size_t m = lps - kps + (2*bsize - kps) * (kps - 1)/2;


        if (val3 != 0. and n <= m)
          spda_out << con_count << " " << 3 << " " << n << " " << m << " " << val3 << "\n";

    }
    }
    }
  }



    for (size_t ip = 0; ip < bsize; ip++)      // loop over ith constraint matrix
    {
    for (size_t jp = 0; jp < bsize; jp++)      // loop over jth constraint matrix
    {
    for (size_t kp = 0; kp < bsize; kp++)    // loop over matrix row
    {
    for (size_t lp = 0; lp < bsize; lp++)    // loop over matrix column
    {
//        if (jp == lp && kp < ip)
//        continue;

        double val5 =  
        (1./2.) * (F10_5_matrix (i, j, k, l, ip, jp, kp, lp) + F10_5_matrix (i, j, k, l, kp, lp, ip, jp));

        size_t n = ip * bsize + jp + 1;
        size_t m = kp * bsize + lp + 1;

        if (val5 != 0. and n <= m)
          spda_out << con_count << " " << 5 << " " << n << " " << m << " " << val5 << "\n";    

    }
    }
    }
    }

    con_count++;

  }
  }
  }
  }

}
/*

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
*/
