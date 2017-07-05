#include <iostream>
#include "output.h"

/***************************************************************

 
 Function to automate output of constraint matrices.


***************************************************************/


void con_matrix_out (const two_array & m_pass, size_t con_count, size_t bsize, struct con_flags flag_pass, std::ofstream & spda_out, const one_array & block_mat)
{

  const size_t PQsize = bsize*(bsize-1)/2;

  size_t block = 1;
  size_t lower = 0;
  size_t upper = bsize;

  for (size_t i = lower; i < upper; i++)
  {
  for (size_t j = i;     j < upper; j++)
  {
    size_t k = i + 1 - lower;
    size_t l = j + 1 - lower;

    if (m_pass[i][j] != 0.0)
      spda_out << con_count << " " << block << " " << k << " " << l << " " << m_pass[i][j] << std::endl;
  }
  }

  block++;

  lower += bsize;
  upper += bsize;

  for (size_t i = lower; i < upper; i++)
  {
  for (size_t j = i;     j < upper; j++)
  {
    size_t k = i + 1 - lower;
    size_t l = j + 1 - lower;

    if (m_pass[i][j] != 0.0)
      spda_out << con_count << " " << block << " " << k << " " << l << " " << m_pass[i][j] << std::endl;
  }
  }

  block++;

  lower += bsize;
  upper += PQsize;

  if (flag_pass.two_body_toggle)
  {

  	size_t a = 0;
  	size_t offset = 0;

    for (size_t i = lower; i < upper; i++)
    {
    for (size_t j = i;     j < upper; j++)
    {
//      size_t k = i + 1 - lower;
//      size_t l = j + 1 - lower;
    
      size_t k = i - lower;
      size_t l = j - lower;

      int m = k - offset + 1;
      int n = l - offset + 1;

      size_t sub = (size_t) abs(block_mat[a]);


      if (abs(m) > abs(sub) and a < block_mat.size())
      {
        offset += abs(sub);
        block++;
        a++;

        m -= abs(sub);
        n -= abs(sub);
      }

      if (m_pass[i][j] != 0.0 and abs(m) <= abs(sub) and abs (n) <= abs(sub))
        spda_out << con_count << " " << block << " " << m << " " << n << " " << m_pass[i][j] << std::endl;

    }
    }




    if (flag_pass.Q_flag)
    {

      block++;
      a = 0;
      offset = 0;

      lower += PQsize; 
      upper += PQsize;


      for (size_t i = lower; i < upper; i++)
      {
      for (size_t j = i;     j < upper; j++)
      {
  //      size_t k = i + 1 - lower;
  //      size_t l = j + 1 - lower;
      
        size_t k = i - lower;
        size_t l = j - lower;

        int m = k - offset + 1;
        int n = l - offset + 1;

        size_t sub = (size_t) abs(block_mat[a]);


        if (abs(m) > abs(sub) and a < block_mat.size())
        {
          offset += abs(sub);
          block++;
          a++;

          m -= abs(sub);
          n -= abs(sub);
        }

        if (m_pass[i][j] != 0.0 and abs(m) <= abs(sub) and abs (n) <= abs(sub))
          spda_out << con_count << " " << block << " " << m << " " << n << " " << m_pass[i][j] << std::endl;

      }
      }

    }


    

    if (flag_pass.G_flag)
    {

      lower += PQsize; 
      upper += bsize*bsize;

      for (size_t i = lower; i < upper; i++)
      {
      for (size_t j = i;     j < upper; j++)
      {
        size_t k = i + 1 - lower;
        size_t l = j + 1 - lower;

        if (m_pass[i][j] != 0.0)
          spda_out << con_count << " " << 5 << " " << k << " " << l << " " << m_pass[i][j] << std::endl;
      }
      }

    }

  }



}

/***************************************************************

 
 Function to look for redundancies in constraint terms.


***************************************************************/


void check_constraints (const three_array & Con, const std::string & label, const bool output, const size_t num, const size_t cmat_extent)
{

  std::cout << label << " REDUNDANCY CHECK" << std::endl;

  size_t total = 0;


  for (size_t b = 0; b < num; b++)
  {
  for (size_t q = b; q < num; q++)
  {
    if (q == b)
      continue;

    bool mat_used = false;

    for (size_t l1 = 0;  l1 < cmat_extent; l1++)
    {
    for (size_t l2 = l1; l2 < cmat_extent; l2++)
    {
      if (Con[q][l1][l2] != 0.)
        mat_used = true;
    }
    }

    if (!mat_used)
    	continue;

    if (Con[b] == Con[q])
    {

      total++;
      if (output)
      {
        std::cout << "b " << b << "\t" << "q " << q << std::endl;
        if (false)
        {
          print (std::cout, Con[b]);
          std::cout << "\n";
          print (std::cout, Con[q]);
        }
      }
    }
  }
  }

  std::cout << "Total redundancies: " << total << std::endl;
}

/***************************************************************

 
 Function to output the two-body M scheme sub blocks in the correct
 SPDA format. 


***************************************************************/



void blockdiag_spdaout (const size_t sub_blocks, const two_array & block_mat, std::ofstream & spda_out)
{


	for (size_t a = 0; a < sub_blocks; a++)
	{
		int sub = block_mat [a][3];

		if (sub != 0)
			spda_out << sub << " "; 

	}

/*
	size_t t = 0;

	for (size_t a = 0; a < sub_blocks; a++)
		{

	  		const size_t sub_term = (size_t) block_mat[a][2];


	  		if (sub_term == 1)
	  		{
	  			bool diag_check = true;

	  			double num = 1;


	  			while (diag_check and (a+1) < sub_blocks)
	  			{
					const size_t b = a + 1;

	  				if((size_t)block_mat[b][2] == 1)
	  				{
	  					num++;
	  					a++;
	  				}

	  				else 
	  					diag_check = false;

	  			}


	  			spda_out << (int)(-1 * num) << " "; 

	  		}

	  		else 
	  			spda_out << sub_term << " "; 
	  		
		}
*/
}

/***************************************************************
 
 Function to output the number of two-body sub-blocks. 

***************************************************************/



size_t num_sub_blocks (size_t sub_blocks, const two_array & block_mat)
{

  size_t number = 0;

  for (size_t a = 0; a < sub_blocks; a++)
  {
    double sub = block_mat [a][3];

    if (sub != 0.)
    {
      number++;
    }

  }

  return number;
}

/***************************************************************

 
 Function to output the SPD file in correct format.


***************************************************************/


void create_spda_file (const two_array & block_mat, const one_array & oned_blocks, const two_array & c_matrix, const struct con_flags flag_pass, const three_array & F1_con, const one_array & F1_val, const three_array & F2_con, const one_array & F2_val, const three_array & F3_con, const one_array & F3_val, const three_array & F7_con, const one_array & F7_val, const three_array & F10_con, const one_array & F10_val, const size_t bsize, const size_t cmat_extent, std::ofstream & spda_out)
{
  const size_t F1num = F1_con.size();
  const size_t F2num = F2_con.size();
  const size_t F3num = F3_con.size();

  const size_t F7num = F7_con.size();


  const size_t F10num = F10_con.size();
  //  size_t cmat_extent = F1_con.size();

//  const size_t PQsize = bsize*(bsize-1)/2;

  const size_t sub_blocks = block_mat.size();

  if(flag_pass.redundant_check)
  { 
    check_constraints (F1_con, "F1", true, F1num, cmat_extent);
    check_constraints (F2_con, "F2", true, F2num, cmat_extent);
    check_constraints (F3_con, "F3", true, F3num, cmat_extent);
    check_constraints (F7_con, "F7", false, F7num, cmat_extent);
    check_constraints (F10_con, "F10", true, F10num, cmat_extent);         
  }


  // number of constraints start
  size_t num_cons = 0;  
  
  if (flag_pass.N_flag)
  	num_cons += F1num;	// N trace constraint 

  if (flag_pass.O_flag)			
  	num_cons += F2num; //  p + q constraints  

  if (flag_pass.two_body_toggle)
  {
    if (flag_pass.P_flag)
      num_cons += F3num;

    if (flag_pass.Q_flag)
      num_cons += F7num;

    if (flag_pass.G_flag)
      num_cons += F10num;
  }

  spda_out << num_cons << std::endl;         // output number of constraints

print (std::cout, oned_blocks);


  // number of blocks start
  size_t blocks = 2;

  if (flag_pass.two_body_toggle)
  {
    blocks += num_sub_blocks (sub_blocks, block_mat);

    if (flag_pass.Q_flag)
      blocks += num_sub_blocks (sub_blocks, block_mat);

    if (flag_pass.G_flag)
      blocks += num_sub_blocks (sub_blocks, block_mat);
  }

  spda_out << blocks << std::endl; 						 // output number of blocks in X
 


  // block sizes start
  if (flag_pass.N_flag)
    spda_out << bsize << " ";

  if (flag_pass.O_flag)
    spda_out << bsize << " ";     // output block sizes
  
  if (flag_pass.two_body_toggle)
  {
  	
  	blockdiag_spdaout (sub_blocks, block_mat, spda_out);


    if (flag_pass.Q_flag)
       	blockdiag_spdaout (sub_blocks, block_mat, spda_out);


    if (flag_pass.G_flag)
    	blockdiag_spdaout (sub_blocks, block_mat, spda_out);

  }
  
  spda_out << std::endl;	             


  // constraint values start

  if (flag_pass.N_flag)
  {
    for (size_t i = 0; i < F1num; i++)
      spda_out << F1_val[i] << " ";
  }
						 // output N trace constraint

  if (flag_pass.O_flag)
  {
    for (size_t i = 0; i < F2num; i++)
      spda_out << F2_val[i] << " ";
  }

  if (flag_pass.P_flag)
  {
    for (size_t i = 0; i < F3num; i++)
      spda_out << F3_val[i] << " ";
  }

  if (flag_pass.Q_flag)
  {
    for (size_t i = 0; i < F7num; i++)
      spda_out << F7_val[i] << " ";
  }


  if (flag_pass.G_flag)
  {
    for (size_t i = 0; i < F10num; i++)
      spda_out << F10_val[i] << " ";
  }

  spda_out << std::endl;






  // constraint matrix start

  size_t con_count = 0;


  con_matrix_out (c_matrix, con_count, bsize, flag_pass, spda_out, oned_blocks);
  con_count++;


  if (flag_pass.N_flag)
  {
  for (size_t cnum = 0; cnum < F1num; cnum++)
  {
     con_matrix_out (F1_con[cnum], con_count, bsize, flag_pass, spda_out, oned_blocks);
     con_count++;
  }
  }

  if (flag_pass.O_flag)
  {
	for (size_t cnum = 0; cnum < F2num; cnum++)
	{
	   con_matrix_out (F2_con[cnum], con_count, bsize, flag_pass, spda_out, oned_blocks);
	   con_count++;
	}
  }

  if (flag_pass.P_flag)
  {
	for (size_t cnum = 0; cnum < F3num; cnum++)
	{
	   con_matrix_out (F3_con[cnum], con_count, bsize, flag_pass, spda_out, oned_blocks);
	   con_count++;
	}
  }


 

  if (flag_pass.Q_flag)
  {
  for (size_t cnum = 0; cnum < F7num; cnum++)
  {
     con_matrix_out (F7_con[cnum], con_count, bsize, flag_pass, spda_out, oned_blocks);
     con_count++;
  }
  }

 


  if (flag_pass.G_flag)
  {
  for (size_t cnum = 0; cnum < F10num; cnum++)
  {
     con_matrix_out (F10_con[cnum], con_count, bsize, flag_pass, spda_out, oned_blocks);
     con_count++;
  }
  }


}
