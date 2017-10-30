#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 11:03:52 2017

@author: alex
"""

import pandas as pd
import numpy as np

def nmax_count (nmax, sp_list, tb_list):
    """Wrapper function to call the single-particle basis
    creating function and then the two-particle basis creating function."""
    sp_list = create_sp_list (nmax, sp_list)
    tb_list = create_tb_list (nmax, sp_list, tb_list)
    
    return (sp_list, tb_list)
    
def create_sp_list (nmax, sp_list):
    """Function to create a list of single-particle orbitals for a given Nmax
    truncation in m-scheme harmonic oscillator basis."""
    
    # initialize values    
    ncount = 0
    num= 1
    n = 0
    l = 0
    j = abs (l - 1/2)
    mj = -j


    while nmax >= ncount:    
        # add new single-particle orbital to dataframe
        sp_list = sp_list.append({'number': num, 'n': n, 'l': l, 
                                  'j': j, 'mj' : mj}, ignore_index=True)
    
        num += 1
    
        # cycle through mj and j values
        if mj < j:
            mj += 1
            continue
    
        if j < abs (l + 1/2):
            j += 1
            mj = -j
            continue
      
        if 2*n + (l+1) <= nmax:
            l += 1
            j = abs (l-1/2)
            mj = -j
            continue
        
        n += 1
        l = 0
   
        ncount = 2*n
    
        # reset j,mj for new n,l values
        j = abs (l - 1/2)
        mj = -j
    
        # break if we get above nmax
        if ncount > nmax:
            break
    
    return sp_list

def create_tb_list (nmax, sp_list, tb_list):
    """Function to create the two-particle basis set. Two loops are run over 
    the single-particle basis set taking two elements in the set taking care
    not to do repeats, e.g., 1,2 and 2,1 and not to take the same element
    twice e.g., 1,1"""
    
    total_orbs = len (sp_list)
    num = 1
    
    for orb1 in range (0, total_orbs):
        for orb2 in range (orb1+1, total_orbs):
            
            n1  = sp_list.iloc[orb1].n
            l1  = sp_list.iloc[orb1].l
            j1  = sp_list.iloc[orb1].j
            mj1 = sp_list.iloc[orb1].mj

            n2  = sp_list.iloc[orb2].n
            l2  = sp_list.iloc[orb2].l
            j2  = sp_list.iloc[orb2].j
            mj2 = sp_list.iloc[orb2].mj

            # The hamiltonian is rotationally invariant so            
            # only take elements from total M=0 block
            if mj1 + mj2 != 0:
                continue
            
            ncount = 2*n1 + l1 + 2*n2 + l2
            
            if ncount > nmax:
                continue
            
            tb_list = tb_list.append({'number': num, 
                                      'n1': n1, 'l1': l1, 'j1': j1, 'mj1': mj1,
                                      'n2': n2, 'l2': l2, 'j2': j2, 'mj2': mj2
                                      }, ignore_index=True) 
    
            num +=1
    
    return tb_list
    

def increment_nl (n, l, ncount):
    """Function to increment either n or l in the single-particle space
    depending on which is lower energy. Incrementing n resets l back to 
    zero."""
    
    n_incr = 2 * (n + 1)
    l_incr = 2 * n + (l + 1)
    
    if n_incr < l_incr:
        n += 1
        l = 0
        ncount = n_incr
        
    else:
        l += 1
        ncount = l_incr
    
    return (n, l, ncount)


if __name__ == '__main__':
    
    sp_list = pd.DataFrame(columns= ['number', 'n', 'l', 'j', 'mj'])
    tb_list = pd.DataFrame(columns= ['number', 'n1', 'l1', 'j1', 'mj1', 
                                     'n2', 'l2', 'j2', 'mj2'])


    (sp_list, tb_list) = nmax_count(4, sp_list, tb_list)

    sp_list.number = sp_list.number.astype(int)
    sp_list.n = sp_list.n.astype(int)
    sp_list.l = sp_list.l.astype(int)
    
    tb_list.number = tb_list.number.astype(int)
    tb_list.n1 = tb_list.n1.astype(int)
    tb_list.l1 = tb_list.l1.astype(int)
    tb_list.n2 = tb_list.n2.astype(int)
    tb_list.l2 = tb_list.l2.astype(int)
    
    np.savetxt('../me_files/nmax4_python_sp.dat', sp_list.values, fmt='%6.2f')
    np.savetxt('../me_files/nmax4_python_tb.dat', tb_list.values, fmt='%6.2f')

    print(sp_list) ; print ('\n')
    print(tb_list)