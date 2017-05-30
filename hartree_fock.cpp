/******************************************************

File to hold the function definitions for our hartree_fock
class.

******************************************************/

#include "hartree_fock.h"
#include <armadillo>
#include <iomanip>
#include <string>

/********************************************************

Constructor for the Hartree-Fock class. Constructor takes the 
arguments of total basis states ntotal, total HF iterations, the 
type of system to solve in (e.g. swave or fullm or fullj) and the number
of particles (part_number).

*********************************************************/

hartree_fock::hartree_fock (int ntotal, int total_loops, std::string type, int part_number)
{
  nmax = ntotal;
  total_iter = total_loops;
  density_mat = arma::zeros (nmax, nmax);
  eps = 1.0E-8;					// number to determine when we've reached convergence
  choice = type;
  particles = part_number; 
}

/********************************************************

Function to update the density matrix after diagonalizing the
HF hamiltonian. The new density matrix is taken from the the lowest
N eigenvectors corresponding to the existing N particles.  
Take N lowest eigenvectors from matrix unitary into temp1 and then form the new
density matrix using temp1 * transpose(temp1). The final resulting density matrix
is a linear combination of this new temp1 * transpose(temp1) and the previous 
density matrix. This is done to avoid oscillatory behavior. 

Note that this is slightly more complicated for j-scheme calculations due to
the degeneracy of the basis states. 

********************************************************/

void hartree_fock::update_density_mat (arma::mat & unitary, phys_system & instance)
{
  arma::mat temp1 = arma::zeros (nmax, particles);

  arma::mat dens_prev = density_mat;		// density matrix before updating

  density_mat = arma::zeros (nmax, nmax);

  for (int i = 0; i < particles; i++)
  {
	temp1.col(i) = unitary.col(i);		// pull out the lowest N eigenvectors
  }

  arma::mat full_temp = temp1 * temp1.t();  	// part of the new density matrix

  for (int i = 0; i < nmax; i++)		// j-scheme degeneracy complication
  {
  for (int j = 0; j < nmax; j++)
  {
	full_temp (i,j) = full_temp (i, j) * instance.orb_degen (i, j, "density", particles); 
  }
  }

/*  double trace = 0;

  for (int i = 0; i < nmax; i++)
  {
	trace += full_temp (i,i);
  }*/

  density_mat = 0.5*(dens_prev + full_temp); 	// final resulting density matrix


}

/************************************************************

Function that runs the HF self-consistent procedure. 

************************************************************/


void hartree_fock::run (phys_system & instance, double hw)
{
  double sp_pot;

  arma::mat h_mat;
  arma::vec energy = arma::zeros (nmax, 1);
  arma::vec e_prev = arma::zeros (nmax, 1);
  arma::vec diff = arma::zeros (nmax, 1);
  arma::mat eig_vec = arma::zeros (nmax, nmax);

  // populate the density matrix with an initial guess (here just 1's along the diagonal) 
  for (int i = 0; i < particles; i++)
	  density_mat (i,i) = (1. * instance.orb_degen(i, i, "initial", particles));

  int iter = 0;

  while (iter < total_iter)
  {

    h_mat = arma::zeros (nmax, nmax);			// hamiltonian matrix 
    energy = arma::zeros (nmax, 1);			// vector to hold the eigenvalues
    eig_vec = arma::zeros (nmax, nmax);			// matrix to hold the eigenvectors

    for (int alpha = 0; alpha < nmax; alpha++)
    {
	for (int beta = 0; beta < nmax; beta++)
	{
	
	  sp_pot = 0.;
	  
	  for (int gamma = 0; gamma < nmax; gamma++)
	  {
		for (int delta = 0; delta < nmax; delta++)
		{
		  sp_pot += (density_mat (delta, gamma) * instance.pot_access_modifier (gamma, alpha, delta, beta, particles)); 
		
		} 

	  } // END INTERNAL SINGLE PARTICLE POTENTIAL LOOPS

	// populate hamiltonian matrix (Redundacy in alpha, beta loops? should it go over only beta(0,alpha) alpha(0, nmax)?)
	  h_mat (alpha, beta) = instance.kinetic_access (alpha, beta, hw) + sp_pot;

	}
    } // end top two for loops


  if (iter == 0)
  {
    std::cout << std::endl << "Initial Hamiltonian" << std::endl;
    h_mat.print();
    std::cout << std::endl << "Initial Density Matrix" << std::endl;
    density_mat.print();
  }

  arma::eig_sym (energy, eig_vec, h_mat);

  if (iter == 0)
  {
    std::cout << std::endl << "Initial Eigenvectors" << std::endl;
    eig_vec.print();
    std::cout << std::endl << "Initial Energies" << std::endl; 
    energy.print();
  }

  iter++; // add 1 to iterations


  // Record the difference between this iteration's energy eigenvalues and the previous ones
  for (int i = 0; i < nmax; i++)
  {
	diff (i)  = fabs(energy (i) - e_prev (i));
  }
  
  // break while loop if largest element of difference smaller than eps (tolerance number)
  if (fabs(diff.max()) < eps && iter > 1)
	break;

  // update the density matrix with the new eigenvectors
  update_density_mat (eig_vec, instance); 

  
  std::cout << "Iter-" << iter << std::setw(40) << " Intermed. G.S. Energy = " << calc_energy(instance, hw) << std::endl;

  e_prev = energy;

  } // end WHILE loop

  std::cout << "Final Eigenvalues" << std::endl;

  energy.print();

  std::cout << "Total loops = " << iter << std::endl;

  std::cout << "Final Density Matrix" << std::endl;
  density_mat.print();

  std::cout << std::endl << "Final Ground State Energy = " << std::setprecision(10) << calc_energy(instance, hw) << std::endl;

}

/******************************************************

Function to calculate the final ground state energy of the
neutron drop system. Result is to take the trace of the 1-body terms 
(kinetic term and harmonic oscillator trap single particle potential) 
with the final density matrix. Then, add 0.5 times the trace of the 
antisymmetric two body potential with two density matrices. 

e.g. with repeated indices summed over

E_HF = (T_ij * p_ji) + (1/2 * V_{ij,kl} * p_ki * p_lj)

where:
  
  T = 1-body terms
  p = final density matrix
  V = antisymmetric two-body matrix elements

*******************************************************/

double hartree_fock::calc_energy (phys_system & instance, double hw) 
{
  double one_body = 0.;
  double two_body = 0.;
  double total_E = 0.;

  for (int i = 0; i < nmax; i++)
  {
  for (int j = 0; j < nmax; j++)
  {
	one_body += instance.kinetic_access (i, j, hw) * density_mat (j, i);
  }
  }

  double gamma;

  for (int a = 0; a < nmax; a++)
  {
  for (int b = 0; b < nmax; b++)
  {
  
  gamma = 0.;

  for (int c = 0; c < nmax; c++)
  {
  for (int d = 0; d < nmax; d++)
  {
//	gamma += density_mat (d,c) * instance.pot_access_modifier (b, c, a, d, particles);
//	gamma += density_mat (d,c) * instance.pot_access_modifier (a, c, b, d, particles);
	gamma += density_mat (d,c) * instance.pot_access_modifier (c, b, d, a, particles);
  }
  }
	two_body += 0.5*density_mat (a,b) * gamma;
  }
  }

 total_E = one_body + two_body;

  return total_E;
}

