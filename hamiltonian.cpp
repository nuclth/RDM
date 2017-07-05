#include "hamiltonian.h"
#include <gsl/gsl_blas.h>
#include <iostream>
#include "auxiliary.h"





/***************************************************************


 Function to populate the 1-body part of the Hamiltonian. 
 The 1-body part here is T + U for kinetic energy T and external potential U


***************************************************************/


void populate_1body (const two_array & ref_m, two_array & h1_mat, const double hw)
{

  size_t mat_length = h1_mat.size();


  for (size_t i = 0; i < mat_length; ++i)
  {
    double n = ref_m [i][1];
    double l = ref_m [i][2];

    h1_mat[i][i] = (2.*n + l + 1.5) * hw;
  }

}


/***************************************************************

 
 Function to read in the reference file for m scheme.


***************************************************************/



void read_in_reference_m_scheme (two_array & ref_m, const std::string m_ref_file, std::ofstream & diag_out, const bool diag_toggle)
{
  const char * m_reference_file = (m_ref_file).c_str();
  // input file stream for m_scheme
  std::ifstream m_ref_in (m_reference_file);
 
  if (m_ref_in.fail())
  throw "ERROR: Cannot open m-scheme reference file";
 
  size_t total_lines = 0;
  std::string dummy;
  double ref_num, n, l, j , m_j, tz;

  // find total number of defined reference lines
  while (std::getline (m_ref_in, dummy))
  ++total_lines;


  // clear the file stream, reset to read in the elements
  m_ref_in.clear();
  m_ref_in.seekg (0, std::ios::beg);

  size_t ref_size = ref_m.size();
  size_t ele_in = 0;

  // read in and assign the references line by line
  // to the matrix ref_m
  for (size_t i = 0; i < total_lines; i++)
  {
  std::string orbit_dummy_1;
  std::string orbit_dummy_2;
  std::getline (m_ref_in, dummy);

  if (!dummy.length() || dummy[0] == '#')     // skip zero length lines and lines that start with #
    continue;

  std::stringstream ss;

  ss << dummy;            // read in the line to stringstream ss

  ss >> orbit_dummy_1 >> orbit_dummy_2 >> ref_num >> n >> l >> j >> m_j >> tz;  // assign values of the line

//  std::cout << ref_num << " " << n << " " << l << " " << j << " " << m_j << " " << tz << "\n";


  // only extract neutron-neutron states
  if (tz == 1.)
  {
    ref_m [ele_in][0] = ref_num;    // reference number of the line
    ref_m [ele_in][1] = n;          // principle quantum number
    ref_m [ele_in][2] = l;          // orbital angular mom.
    ref_m [ele_in][3] = j * 0.5;    // total angular mom.
    ref_m [ele_in][4] = m_j * 0.5;  // total angular mom. projection
    ref_m [ele_in][5] = tz;         // isospin projection (should all be +1.0)
    ref_m [ele_in][6] = ele_in;     // new index for sp orbital
    ele_in++;
  }

  if (ele_in >= ref_size)
    break;

  }
  
  
  
  if(diag_toggle)
  {
    diag_out << "Single particle orbitals pulled from REF file" << std::endl << std::endl;
    print(diag_out, ref_m);            // print the resulting matrix
    diag_out << std::endl << std::endl;
  }


  
}

/***************************************************************

Function to read in the actual m scheme matrix  elements 
calculated from Morten's code. 

THIS WILL BREAK IF THE INPUT FILE .dat CHANGES TITLE OR FORMAT.

***************************************************************/


void read_in_matrix_m_scheme (const two_array & ref_m, five_array & h2_mat, const std::string m_mat_file)
{
  // input file stream
  const char * m_matrix_file = (m_mat_file).c_str();
  std::ifstream m_matrix_in (m_matrix_file);

  size_t total_lines = 0;
  std::string dummy;
  int alpha, beta, gamma, delta;
  double value;

  if (m_matrix_in.fail())
  throw "ERROR: Cannot read m-scheme matrix file";

  // get the total number of lines in the matrix elements file
  while(std::getline(m_matrix_in, dummy))
  ++total_lines;
 
  // clear and reset the file stream
  m_matrix_in.clear();
  m_matrix_in.seekg(0, std::ios::beg);

  const size_t m_size = h2_mat.size();

  // set the size of the m_scheme holder
  //m_scheme_matrix.set_size (m_size, m_size);


  int i,j,k,l;

  for (size_t n = 0; n < total_lines; n++)
  {
  std::getline (m_matrix_in, dummy);

  if (!dummy.length() || dummy[0] == '#')
    continue;

  i = -1;
  j = -1;
  k = -1;
  l = -1;

  std::stringstream ss;

  ss << dummy;

  ss >> alpha >> beta >> gamma >> delta >> value;

  // check to see if alpha, beta, gamma, delta we're reading in
  // matches the reference numbers in our reference matrix
  // if yes --> read it in, otherwise discard it

  for (size_t w = 0; w < m_size; w++)
  {
    if (alpha == ref_m[w][0]) 
      i = w;
  }
  
  for (size_t w = 0; w < m_size; w++)
  {
    if (beta == ref_m[w][0]) 
      j = w;
  }

  for (size_t w = 0; w < m_size; w++)
  {
    if (gamma == ref_m[w][0]) 
      k = w;
  }

  for (size_t w = 0; w < m_size; w++)
  {
    if (delta == ref_m[w][0]) 
      l = w;
  }

  if (i >= 0 && j >= 0 && k >= 0 && l >= 0) // only activates if all 4 "if" statements above are true
    {
      h2_mat [i][j][k][l][0] = value;
      h2_mat [i][j][k][l][1] = alpha;
      h2_mat [i][j][k][l][2] = beta;
      h2_mat [i][j][k][l][3] = gamma;
      h2_mat [i][j][k][l][4] = delta;
    }
  }

  // reset index term on reference matrix to span 0 to m_size instead of e.g. 3, 4, 9, 10, etc... for neutrons

//  for (size_t n = 0; n < m_size; n++)
//  ref_m[n][0] = n;

//  print(std::cout, ref_m);


}



/***************************************************************

Function for M projection and parity blocks.

***************************************************************/


size_t block_list (const two_array & ref_m, two_array & block_mat, two_array & bcopy)
{
	size_t num 	= 0;


	const size_t orbital_number = ref_m.size();

	for (size_t i = 0;   i < orbital_number; i++)
	{
	for (size_t j = i+1; j < orbital_number; j++)
	{
		bool unique = true;


		double P = -1.;

		const size_t l1 = ref_m[i][2];
		const double m1 = ref_m[i][4];

		const size_t l2 = ref_m[j][2];
		const double m2 = ref_m[j][4];

		if ((l1 + l2) % 2 == 0)
           	P = 0.;

        else if ((l1 + l2) % 2 == 1)
           	P = 1.;

		const double mtot = m1 + m2;

		for (size_t k = 0; k < num; k++)
		{
			const double taken_M = block_mat [k][0];
			const double taken_P = block_mat [k][1];

			if (taken_M == mtot and taken_P == P)
				unique = false;

		}

		if (unique)
		{
			block_mat [num][0] = mtot;
			block_mat [num][1] = P;

			bcopy     [num][0] = mtot;
			bcopy     [num][1] = P;
			num++;

		}
	}
	}

//	std::cout << num << std::endl;

	//print (std::cout, bcopy);

	size_t even = 0;

	for (size_t i = 0; i < num; i++)
	{
		bool check = false;

		double M_start = bcopy [i][0];
		double P_start = bcopy [i][1];


//		std::cout << M_start << " " << P_start << std::endl;

		size_t mark = i;

//		std::cout << M_start << " " << P_start << std::endl;

		for (size_t j = 0; j < num; j++)
		{
			double M_loop = bcopy [j][0];
			double P_loop = bcopy [j][1];

			if (M_loop <= M_start and P_loop == 0.)
			{
				M_start = M_loop;
				P_start = P_loop;
				mark = j;
				check = true;
			}
		}

//		std::cout << M_start << " " << P_start << std::endl << std::endl;
		if (check)
		{
			block_mat [even][0] = M_start;
			block_mat [even][1] = P_start;

			bcopy [mark][1] = 2.;

//		print (std::cout, bcopy);

			even++;
		}
	}

	size_t odd = even;

//	print (std::cout, bcopy);
//	print (std::cout, block_mat);


	for (size_t i = 0; i < num; i++)
	{
		bool check = false;

		double M_start = bcopy [i][0];
		double P_start = bcopy [i][1];

		size_t mark = i;

		for (size_t j = 0; j < num; j++)
		{
			double M_loop = bcopy [j][0];
			double P_loop = bcopy [j][1];


			if (M_loop <= M_start and P_loop == 1.)
			{
				M_start = M_loop;
				P_start = P_loop;
				mark = j;
				check = true;
			}
		}

		if (check)
		{
			block_mat [odd][0] = M_start;
			block_mat [odd][1] = P_start;

			bcopy [mark][1] = 2.;

			odd++;
		}
	}





	for (size_t a = 0; a < odd; a++)
	{

	double N = 0.;

	for (size_t i = 0;   i < orbital_number; i++)
	{
	for (size_t j = i+1; j < orbital_number; j++)
	{

	for (size_t k = i; 	 k < orbital_number; k++)
	{
	for (size_t l = k+1; l < orbital_number; l++)
	{

		if (j < l && k <= i)
      		continue;

		double P1 = -1.;

		const size_t l1 = ref_m[i][2];
		const double m1 = ref_m[i][4];

		const size_t l2 = ref_m[j][2];
		const double m2 = ref_m[j][4];

		if ((l1 + l2) % 2 == 0)
           	P1 = 0.;

        else if ((l1 + l2) % 2 == 1)
           	P1 = 1.;

		const double mleft = m1 + m2;



		double P2 = -1.;

		const size_t l3 = ref_m[k][2];
		const double m3 = ref_m[k][4];

		const size_t l4 = ref_m[l][2];
		const double m4 = ref_m[l][4];

		if ((l3 + l4) % 2 == 0)
           	P2 = 0.;

        else if ((l3 + l4) % 2 == 1)
           	P2 = 1.;

		const double mright = m3 + m4;

		if (mleft == mright and P1 == P2)
		{

			if (block_mat[a][0] == mleft and block_mat[a][1] == P1)
				N++;
		}

	}
	}

	}
	}

	double r = (sqrt(1. + 8. * N) - 1.) / 2;

	block_mat [a][2] = r;

	}



	size_t t = 0;


	for (size_t a = 0; a < odd; a++)
	{
		const double sub_term = block_mat[a][2];


	  	if (sub_term == 1.)
	  	{
	  		bool diag_check = true;

  			double sub_loop = 1.;


			while (diag_check and (a+1) < odd)
			{
				const size_t b = a + 1;

	  			if(block_mat[b][2] == 1.)
	  			{
	  				sub_loop++;
	  				a++;
	  			}

	  			else 
	  				diag_check = false;

	  		}

	  		block_mat[t][3] = (-1. * sub_loop);
	  		t++;

	  	}

	  	else 
	  	{
  			block_mat[t][3] = sub_term;
  			t++;
  		}
	}

	return num;
}



size_t fill_oned_blocks (const size_t sub_blocks, one_array & oned_blocks, const two_array & block_mat)
{
  size_t num = 0; 

	for (size_t a = 0; a < sub_blocks; a++)
	{
  	oned_blocks [a] = block_mat[a][3];

    if (block_mat[a][3] != 0.)
      num++;
  }

  return num;
}


/***************************************************************

Wrapper function to read in single-particle m scheme reference 
file, populate 1-body Hamiltonian matrix elements, and then 
populate 2-body matrix elements.

***************************************************************/



void fullm_populate_hamiltonian (two_array & ref_m, two_array & h1_mat, five_array & h2_mat, const std::string reference_file, const std::string matrix_file, const double hw, std::ofstream & diag_out, const bool diag_toggle, two_array & block_mat) 
{

  try
  {
    read_in_reference_m_scheme (ref_m, reference_file, diag_out, diag_toggle);
    populate_1body (ref_m, h1_mat, hw);
    read_in_matrix_m_scheme (ref_m, h2_mat, matrix_file);
  }

  catch (const char * msg)
  {
    std::cerr << msg << std::endl;  
  }
 

}


/***************************************************************

Function to put the potential in two index form. 

***************************************************************/


void compactify_h2 (const two_array & ref_m, two_array & comp_h2, five_array & h2_mat, std::ofstream & diag_out, const bool diag_toggle, two_array & basis_ref)
{
  const size_t bsize = h2_mat.size();

  size_t z = 0;

  double value;

  size_t alpha, beta, gamma, delta;

  size_t l1, l2, l3, l4;
  double m1, m2, m3, m4;
  double mleft, mright;

  if(diag_toggle)
  {
    diag_out << "Output of nonzero 2-body ME" << std::endl << std::endl;    
  }

  for (size_t i = 0; i < bsize; ++i)
  {
  for (size_t j = i+1; j < bsize; ++j)
  {

  for (size_t k = 0; k < bsize; ++k)
  {
  for (size_t l = k+1; l < bsize; ++l)
  {

//  	if (j < l && k <= i)
//      continue;

            value = h2_mat [i][j][k][l][0];

            alpha = h2_mat [i][j][k][l][1];
            beta  = h2_mat [i][j][k][l][2];
            gamma = h2_mat [i][j][k][l][3];
            delta = h2_mat [i][j][k][l][4];

//            std::cout << alpha << " " << beta << "\t" << gamma << " " << delta << " " << value << "\n";

            if (diag_toggle)
            {
              if (value != 0.)
              {
              	diag_out << alpha << " " << beta << "\t" << gamma << " " << delta << " " << value << " \t\t";

                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (alpha == ref_m[loop][0])
                  {
                    diag_out << ref_m[loop][6] << "\t";
                    l1 = ref_m[loop][2];
                    m1 = ref_m[loop][4];
                  }
                }


                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (beta == ref_m[loop][0])
                  {
                    diag_out << ref_m[loop][6] << "\t";
                    l2 = ref_m[loop][2];
                    m2 = ref_m[loop][4];
                  }
                }


                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (gamma == ref_m[loop][0])
                  {
                    diag_out << ref_m[loop][6] << "\t";
                    l3 = ref_m[loop][2];
                    m3 = ref_m[loop][4];
                  }
                }


                for (size_t loop = 0; loop < bsize; loop++)
                {
                  if (delta == ref_m[loop][0])
                  {
                    diag_out << ref_m[loop][6] << "\t";
                    l4 = ref_m[loop][2];
                    m4 = ref_m[loop][4];
                  }
                }

                char P1 = 'W';

                if ((l1 + l2) % 2 == 0)
                	P1 = 'E';

                else if ((l1 + l2) % 2 == 1)
                	P1 = 'O';


                char P2 = 'W';

                if ((l3 + l4) % 2 == 0)
                	P2 = 'E';

                else if ((l3 + l4) % 2 == 1)
                	P2 = 'O';

                mleft  = m1 + m2;

                mright = m3 + m4;

                diag_out << mleft << " " << mright;
                diag_out << "\t" << P1 << " " << P2;

                diag_out << "\n";
              }


            }

            size_t ips = i + 1;
            size_t jps = j + 1;
            size_t kps = k + 1;
            size_t lps = l + 1;

            size_t left  = jps - ips + (2*bsize - ips) * (ips - 1)/2 - 1;
            size_t right = lps - kps + (2*bsize - kps) * (kps - 1)/2 - 1;


  /*         if (i == 0 && j == 4 && k == 0)
            {
            	std::cout << ref_m[l][1] << "\t" << ref_m[l][2] << "\n";
            }

           if (i == 0 && j == 1)
            {
            	std::cout << right << " " << l << " " << k << "\t" << value << "\n"; 
            }
*/

//            std::cout << value << std::endl;

            comp_h2 [left][right] = value;

            if (alpha != 0 and gamma != 0 and beta != 0 and delta != 0)
            {

            	bool new_val = true;

            	for (size_t q = 0; q < basis_ref.size(); q++)
            	{
            		if (basis_ref [q][0] == alpha and basis_ref [q][1] == beta)
            			new_val = false;

            		if (basis_ref [q][1] == alpha and basis_ref [q][0] == beta)
            			new_val = false;            		
            		
            	}


            	if (new_val)
	            {
	            	basis_ref [z][0] = alpha;
    	        	basis_ref [z][1] = beta;
        	    	z++;
//        	    	std::cout << "HIT" << std::endl;
        	    }
            }

  }
  }
  }
  }


//  print (std::cout, basis_ref);

  if (diag_toggle)
    diag_out << std::endl << std::endl;
}





/***************************************************************

Function to put the potential in two index form with block 
diagonlization. 

***************************************************************/

void blockdiag_h2 (const two_array & ref_m, two_array & block_h2, five_array & h2_mat, std::ofstream & diag_out, const bool diag_toggle, const two_array & block_mat, two_array & basis_ref)
{
  const size_t bsize = h2_mat.size();

  const size_t sub_blocks = block_mat.size();

  double value;

  size_t alpha, beta, gamma, delta;

  size_t l1, l2, l3, l4;
  double m1, m2, m3, m4;
  double mleft, mright;

//  if(diag_toggle)
//  {
//    diag_out << "Output of nonzero 2-body ME" << std::endl << std::endl;    
//  }

  size_t offset = 0;

  size_t z = 0;

  for (size_t a = 0; a < sub_blocks; a++)
  {

  	size_t left  = 0;
  	size_t right = 0;


  for (size_t i = 0;   i < bsize; ++i)
  {
  for (size_t j = i+1; j < bsize; ++j)
  {

  for (size_t k = 0;   k < bsize; ++k)
  {
  for (size_t l = k+1; l < bsize; ++l)
  {

//  	if (j < l && k <= i)
//      continue;

            value = h2_mat [i][j][k][l][0];

            alpha = h2_mat [i][j][k][l][1];
            beta  = h2_mat [i][j][k][l][2];
            gamma = h2_mat [i][j][k][l][3];
            delta = h2_mat [i][j][k][l][4];

//              	diag_out << alpha << " " << beta << "\t" << gamma << " " << delta << " " << value << " \t\t";

            for (size_t loop = 0; loop < bsize; loop++)
            {
                if (alpha == ref_m[loop][0])
                {
//                    diag_out << ref_m[loop][6] << "\t";
                    l1 = ref_m[loop][2];
                    m1 = ref_m[loop][4];
                }
            }


            for (size_t loop = 0; loop < bsize; loop++)
            {
                if (beta == ref_m[loop][0])
                {
//                    diag_out << ref_m[loop][6] << "\t";
                    l2 = ref_m[loop][2];
                    m2 = ref_m[loop][4];
                }
            }


            for (size_t loop = 0; loop < bsize; loop++)
            {
                if (gamma == ref_m[loop][0])
                {
//                    diag_out << ref_m[loop][6] << "\t";
                    l3 = ref_m[loop][2];
                    m3 = ref_m[loop][4];
                }
            }


            for (size_t loop = 0; loop < bsize; loop++)
            {
                if (delta == ref_m[loop][0])
                {
//                    diag_out << ref_m[loop][6] << "\t";
                    l4 = ref_m[loop][2];
                    m4 = ref_m[loop][4];
                }
            }


            double P1 = 2;

            if ((l1 + l2) % 2 == 0)
                P1 = 0;

            else if ((l1 + l2) % 2 == 1)
               	P1 = 1;


            double P2 = 2;

            if ((l3 + l4) % 2 == 0)
                P2 = 0;

            else if ((l3 + l4) % 2 == 1)
               	P2 = 1;


            mleft  = m1 + m2;

            mright = m3 + m4;

                //diag_out << mleft << " " << mright;
                //diag_out << "\t" << P1 << " " << P2;

                //diag_out << "\n";



            if (mleft == mright and mleft == block_mat[a][0] and P1 == P2 and P1 == block_mat[a][1] and alpha != 0)
            {

            	block_h2 [left + offset][right + offset] = value;

            	right++;

            	if (left == 0)
            	{
            	    basis_ref [z][0] = gamma;
            		basis_ref [z][1] = delta;
            		z++;
            	}

            	if (right >= block_mat [a][2])
            	{
            		left++;
            		right = 0;
            	}
//            	std::cout << mleft << " " << P1 << std::endl;
//            	std::cout << alpha << " " << beta << "\t" << gamma << " " << delta << "\t";
//            	std::cout << value << std::endl << std::endl;

            }

/*            size_t ips = i + 1;
            size_t jps = j + 1;
            size_t kps = k + 1;
            size_t lps = l + 1;

            size_t left  = jps - ips + (2*bsize - ips) * (ips - 1)/2 - 1;
            size_t right = lps - kps + (2*bsize - kps) * (kps - 1)/2 - 1;
*/

  /*         if (i == 0 && j == 4 && k == 0)
            {
            	std::cout << ref_m[l][1] << "\t" << ref_m[l][2] << "\n";
            }

           if (i == 0 && j == 1)
            {
            	std::cout << right << " " << l << " " << k << "\t" << value << "\n"; 
            }
*/

//            std::cout << value << std::endl;

//            comp_h2 [left][right] = value;

  }
  }
  }
  }

  offset += block_mat [a][2];

  }

//  print (std::cout, basis_ref);

//  print (std::cout, block_h2);
//  if (diag_toggle)
//    diag_out << std::endl << std::endl;
}


/***************************************************************

***************************************************************/

void create_transformation (const two_array & compact_ref, const two_array & block_ref, two_array & trans_h2, size_t extent)
{
	for (size_t i = 0; i < extent; i++)
	{
		size_t com_a = compact_ref [i][0];
		size_t com_b = compact_ref [i][1];

	for (size_t j = 0; j < extent; j++)
	{
		size_t blo_a = block_ref [j][0];
		size_t blo_b = block_ref [j][1];

		if (blo_a == com_a and blo_b == com_b)
			trans_h2 [i][j] = 1.0;
	}
	}

//	print (std::cout, trans_h2);
}


void verify_tranformation  (const two_array & comp_h2, const two_array & block_h2, const two_array & trans_h2)
{
	const size_t extent = comp_h2.size();

	double perm_data  [extent*extent];
	double comp_data  [extent*extent];
	double zero1_data [extent*extent];
	double zero2_data [extent*extent];
/*

	gsl_matrix * p  = gsl_matrix_alloc (extent, extent);
	gsl_matrix * pt = gsl_matrix_alloc (extent, extent);
	gsl_matrix * c  = gsl_matrix_alloc (extent, extent);

	gsl_matrix * inter  = gsl_matrix_alloc (extent, extent);
	gsl_matrix * answer = gsl_matrix_alloc (extent, extent);

*/

	size_t q = 0;

	for (size_t i = 0; i < extent; i++)
	{
	for (size_t j = 0; j < extent; j++)
	{

		perm_data  [q] = trans_h2 [i][j];
		comp_data  [q] = comp_h2  [i][j];
		zero1_data [q] = 0.0;
		zero2_data [q] = 0.0;
		q++;
//		gsl_matrix_set (p,  i, j, trans_h2 [i][j]);
//		gsl_matrix_set (pt, i, j, trans_h2 [i][j]);
//		gsl_matrix_set (c,  i, j, comp_h2  [i][j]);

	}
	}

	gsl_matrix_view p  = gsl_matrix_view_array (perm_data, extent, extent);
	gsl_matrix_view c  = gsl_matrix_view_array (comp_data, extent, extent);

	gsl_matrix_view inter  = gsl_matrix_view_array (zero1_data, extent, extent);
	gsl_matrix_view answer = gsl_matrix_view_array (zero2_data, extent, extent);


  gsl_blas_dgemm (CblasTrans,   CblasNoTrans, 1.0, &p.matrix, &c.matrix, 0.0, &inter.matrix);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &inter.matrix, &p.matrix, 0.0, &answer.matrix);

/*
	for (size_t i = 0; i < extent; i++)
	{
	gsl_vector_view col_i = gsl_matrix_row (pt.matrix, i);

		for (size_t j = 0; j < extent; j++)
		{
			double value = gsl_vector_get(&col_i.vector, j);

			std::cout << value << " ";

		}

	std::cout << std::endl;
	} 

	for (size_t i = 0; i < extent; i++)
	{
	gsl_vector_view col_i = gsl_matrix_row (c.matrix, i);

		for (size_t j = 0; j < extent; j++)
		{
			double value = gsl_vector_get(&col_i.vector, j);

			std::cout << value << " ";

		}

	std::cout << std::endl;
	} 


	for (size_t i = 0; i < extent; i++)
	{
	gsl_vector_view col_i = gsl_matrix_row (pt.matrix, i);

		for (size_t j = 0; j < extent; j++)
		{
			double value = gsl_vector_get(&col_i.vector, j);

			std::cout << value << " ";

		}

	std::cout << std::endl;
	} 




	for (size_t i = 0; i < extent; i++)
	{
	gsl_vector_view col_i = gsl_matrix_row (pt.matrix, i);

		for (size_t j = 0; j < extent; j++)
		{
			double value = gsl_vector_get(&col_i.vector, j);

			std::cout << value << " ";

		}

	std::cout << std::endl;
	} 

*/


  	q = 0;

	for (size_t i = 0; i < extent; i++)
	{
	for (size_t j = 0; j < extent; j++)
	{
		double value = zero2_data [q];
		q++;

		if (value != block_h2 [i][j])
		{
//			std::cout << value << " " << block_h2 [i][j] << std::endl;
			std::cerr << "ERROR IN VERIFICATION: TRANSFORMATION DOES NOT REPRODUCE BLOCK" << std::endl;
		}
	}
	} 

}
