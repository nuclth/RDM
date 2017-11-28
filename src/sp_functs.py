#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 11:14:38 2017

@author: alex
"""

import sys

# import third party libraries
import pandas as pd
import numpy as np



def create_sp_list (nmax, sp_list):
    """Function to create a list of single-particle orbitals for a given Nmax
    truncation in m-scheme harmonic oscillator basis. The sp loop first goes over the
    angular projection and then the total angular momentum quantum number. Afterwards, 
    the loop increments to a greater total N = 2n + l (e.g. N = 1 with n,l = 0,1 to 
    N = 2 with n,l = 0,2 and 1,0).
    
    Arguments:
        nmax    - nmax truncation size
        sp_list - dataframe to hold single-particle orbital quantum numbers 
        
    Returns:
        sp_list
            number - sp orbital number
            n - principal quantum number
            l - orbital quantum  number
            2j - twice the angular momentum quantum number 
            2mj - twice ang. mom. quantum num. projection 
    """
    
    # initialize values
    new_Ncount = True    
    ncount = 0
    num= 1
    n = 0
    l = 0
    twoj = 2 * abs (l + 1/2)
    twomj = -twoj
    icount = 0
    l_list = []

    while nmax >= ncount:    
        
        # add new single-particle orbital to dataframe
        sp_list = sp_list.append({'number': num, 'n': n, 'l': l, 
                                  '2j': twoj, '2mj' : twomj}, ignore_index=True)
    
        # increment sp orbital counter
        num += 1
    
        # cycle through mj values for a given n,l,j
        if twomj < twoj:
            twomj += 2
            continue
    
        # cycle through j values for a given n,l
        if twoj > 2 * abs (l - 1/2):
            twoj = twoj - 2
            twomj = -twoj
            continue

        # conditions to increment N count 
        if icount >= len(l_list) and not new_Ncount:
            new_Ncount = True

        # increment N count
        if new_Ncount:
            new_Ncount = False
            ncount += 1
            icount = 0
            
            # create matching lists of total N for n,l 
            # e.g., l_list[a] + n_list[a] = N for any a
            l_list = [a for a in range (ncount, -1, -2)]
            n_list = [a for a in range (0, int(np.floor(ncount/2))+1, 1)]
      
        if ncount > nmax:
            break
        
        # set new n and l
        l = l_list[icount]
        n = n_list[icount]
        
        # increment n,l counter for a given total N
        icount += 1

        # reset j,mj for new n,l values
        twoj = 2 * abs (l + 1/2)
        twomj = -twoj
    
    
    return sp_list




def sort_sp (sp_list):
    """Sort sp_list dataframe by 2j, 2mj values in ascending order. Then sort by
    parity of the given sp orbital, even first then odd.
    
    Argument:
        sp_list - dataframe to hold sp orbital quantum numbers
        
    Returns:
        sp_list sorted by parity and parity blocks sorted by 2j, 2mj
    """

  
    sp_list = sp_list.sort_values(['2j','2mj'])
        
    even_parity = sp_list['l'] % 2 == 0
    odd_parity = sp_list['l'] % 2 == 1
    
    shifted_list = sp_list[even_parity]
    shifted_list = shifted_list.append(sp_list[odd_parity])

    # reset index for dataframe
    shifted_list = shifted_list.reset_index(drop=True)

    
    return shifted_list
    



def add_blocks (sp_list):
    """Add new column called 'blocks' to sp dataframe. Value is the block
    size of a set of given sp orbitals. The block size is only given for the
    last orbital in a given block; all other values are set to 0. 
    
    Argument:
        sp_list - dataframe to hold sp orbital quantum numbers
        
    Returns:
        sp_list with additional blocks column
    """
  
    # dummy list, all values initialized to -2
    n_list = pd.DataFrame(sp_list[['n']])
    n_list.loc[:,'blocks'] = -2

    # counter for # of orbitals in a block
    ncount = 1
    
    # loop over n quantum number for sp orbital and the next sp orbital
    for val1, val2 in zip(n_list.iterrows(), n_list.iloc[1:].iterrows()):
        
        # n values
        n_set = val1[1][0]
        n_nxt = val2[1][0]
        
        # if both in same block
        if n_set < n_nxt:
            ncount += 1
            val1[1][1] = 0
        
        # if new block case 1
        if n_set == n_nxt:
            ncount = 1
            val1[1][1] = 1

        # if new block case 2
        if n_set > n_nxt:
            val1[1][1] = ncount
            ncount = 1
  
    # set last block value (always 1)
    n_list.iloc[len(n_list)-1][1] = 1

    sp_list.loc[:,'blocks'] = n_list.loc[:,'blocks']
    
    return sp_list



def sp_relational_db_morten (nmax, sp_list):
    """Function to create a relational table for the sp orbital numbering scheme
    used by Morten and myself. This is needed to determine which sp state is which in
    our two schemes and convert seamlessly. This function is necessary for two reasons:
        1. Morten's sp states include proton states (not needed by me) 
        2. I reorder my sp states for block diagonalization purposes
    
    Arguments:
        nmax    - nmax truncation size
        sp_list - dataframe to hold single-particle orbital quantum numbers
        
    Returns:
        relational_db
            morten_sp - morten's sp orbital number for a given state
            alex_sp - my sp orbital number for the same state
    """
    
    
    morten_splist = pd.DataFrame(columns = ['number', 'n', 'l', '2j', '2mj'])
    
    # open and read in Morten sp data file
    with open('../me_files/ref_files/nmax' + str(nmax) + "_spm.dat", 'r') as f:
        content = f.readlines()
        for x in content:
            row = x.split()

            if row[0] == '#': continue 
            num   = int(row[2])
            n     = int(row[3])
            l     = int(row[4])
            twoj  = int(row[5])
            twomj = int(row[6])
            
            morten_splist = morten_splist.append({'number':num, 'n':n, 'l':l,'2j':twoj,
                                                  '2mj':twomj}, ignore_index=True)
    
    # create relational dataframe
    relational_db = pd.DataFrame(columns = ['morten_sp', 'alex_sp'])
    

    # fill in relational dataframe with my sp numbers and Morten's
    for index, number in zip(morten_splist.index, morten_splist['number']):
        relational_db = relational_db.append({'morten_sp':number, 
                                              'alex_sp': index+1}, ignore_index=True)
    
    return relational_db

