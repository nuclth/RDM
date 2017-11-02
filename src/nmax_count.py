#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 11:03:52 2017

@author: alex
"""

import sys

import pandas as pd
import numpy as np

    
def create_sp_list (nmax, sp_list):
    """Function to create a list of single-particle orbitals for a given Nmax
    truncation in m-scheme harmonic oscillator basis."""
    
    # initialize values
    new_ncount = True    
    ncount = 0
    num= 1
    n = 0
    l = 0
    j = 2 * abs (l + 1/2)
    mj = -j
    icount = 0
    l_list = []

    while nmax >= ncount:    
        # add new single-particle orbital to dataframe
        sp_list = sp_list.append({'number': num, 'n': n, 'l': l, 
                                  'j': j, 'mj' : mj}, ignore_index=True)
    
        num += 1
    
        # cycle through mj and j values
        if mj < j:
            mj += 2
            continue
    
        if j > 2 * abs (l - 1/2):
            j = j - 2
            mj = -j
            continue

 #       (n, l, ncount) = increment_nl (n, l, ncount)

        if icount >= len(l_list) and not new_ncount:
            new_ncount = True

        if new_ncount:
            new_ncount = False
            ncount += 1
            icount = 0
            l_list = [a for a in range (ncount, -1, -2)]
            n_list = [a for a in range (0, int(np.floor(ncount/2))+1, 1)]
      
        if ncount > nmax:
            break
        
#        print (icount)
#        print (len(l_list))
        
        l = l_list[icount]
        n = n_list[icount]
        
        icount += 1

        # reset j,mj for new n,l values
        j = 2 * abs (l + 1/2)
        mj = -j
    
        # break if we get above nmax
#        if ncount > nmax:
#            break
    
    return sp_list



def increment_nl (n, l, ncount):
    """Function to increment either n or l in the single-particle space
    depending on which is lower energy. Incrementing n resets l back to 
    zero."""
    
    ncount = 5
    
  
    l_list = [a for a in range (ncount, -1, -2)]
    n_list = [a for a in range (0, int(np.floor(ncount/2))+1, 1)]

    print (len(l_list))
    sys.exit()

    for n1 in n_list:
        for l1 in l_list:
            print (l1)

    print (n_list)
    print (l_list)
    sys.exit()
    
    if (2*(n+1) + l > ncount and 2*n + l + 1 > ncount): ncount += 1
    
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




def create_tb_list (nmax, sp_list, tb_list):
    """Function to create the two-particle basis set. Two loops are run over 
    the single-particle basis set taking two elements in the set taking care
    not to do repeats, e.g., 1,2 and 2,1 and not to take the same element
    twice e.g., 1,1"""
    
    total_orbs = len (sp_list)
    num = 1
    
    for orb1 in range (0, total_orbs):
        for orb2 in range (orb1+1, total_orbs):
            
            sp1 = sp_list.iloc[orb1].number
            
            n1  = sp_list.iloc[orb1].n
            l1  = sp_list.iloc[orb1].l
            j1  = sp_list.iloc[orb1].j
            mj1 = sp_list.iloc[orb1].mj

            sp2 = sp_list.iloc[orb2].number

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
            
            mj1 *= 0.5
            mj2 *= 0.5
            
            j1 *= 0.5
            j2 *= 0.5
            
            tb_list = tb_list.append({'number': num, 
                                      'n1': n1, 'l1': l1, 'j1': j1, 'mj1': mj1,
                                      'n2': n2, 'l2': l2, 'j2': j2, 'mj2': mj2,
                                      'sp1' : sp1, 'sp2' : sp2
                                      }, ignore_index=True) 
    
            num +=1
    
    return tb_list


def P_flag_num (tb_list):
    
    sp_vals = pd.DataFrame({'sp1': tb_list.sp1, 'sp2':tb_list.sp2})

    repeats = []
    count = 0
    

    for i1, sp1, sp2 in zip(sp_vals.index, sp_vals.sp1, sp_vals.sp2):
        for i2, sp3, sp4 in zip (sp_vals.index, sp_vals.sp1, sp_vals.sp2):
            if sp2 == sp4 and i1 <= i2 and (sp1,sp3) not in repeats:
                repeats.append((sp1, sp3))
#                print (sp1, sp2, sp3, sp4)
                count+=1

    return (count, repeats)

def sort_sp (sp_list):
    
#    shifted_list = pd.DataFrame(columns= ['number', 'n', 'l', 'j', 'mj'])

    sp_list = sp_list.sort_values(['j','mj'])
    sp_list = sort_sp_l (sp_list)
    
    return sp_list


def sort_sp_l (sp_list):
    
    even_parity = sp_list['l'] % 2 == 0
    odd_parity = sp_list['l'] % 2 == 1
    
    shifted_list = sp_list[even_parity]
    shifted_list = shifted_list.append(sp_list[odd_parity])

    return shifted_list


def add_blocks (sp_list):
    
  
    n_list = sp_list[['n']]
    n_list['blocks'] = -2

    ncount = 1
    
    for val1, val2 in zip(n_list.iterrows(), n_list.iloc[1:].iterrows()):
        
        n_set = val1[1][0]
        n_nxt = val2[1][0]
        
        if n_set < n_nxt:
            ncount += 1
            val1[1][1] = 0
        
        if n_set == n_nxt:
            ncount = 1
            val1[1][1] = 1

        if n_set > n_nxt:
            val1[1][1] = ncount
            ncount = 1
  
          
    n_list.iloc[len(n_list)-1][1] = 1

    sp_list['blocks'] = n_list['blocks']
    

    
    return sp_list

if __name__ == '__main__':
    
    sp_list = pd.DataFrame(columns= ['number', 'n', 'l', 'j', 'mj'])
    tb_list = pd.DataFrame(columns= ['number', 'n1', 'l1', 'j1', 'mj1', 
                                     'n2', 'l2', 'j2', 'mj2', 'sp1', 'sp2'])

    nmax = 4

    sp_list = create_sp_list (nmax, sp_list)


    sp_list.number = sp_list.number.astype(int)
    sp_list.n = sp_list.n.astype(int)
    sp_list.l = sp_list.l.astype(int)
    sp_list.j = sp_list.j.astype(int)
    sp_list.mj = sp_list.mj.astype(int)

    
    sp_list = sort_sp (sp_list)
    
    sp_list = add_blocks (sp_list)
    
    print (sp_list)

    
    np.savetxt('../me_files/ref_files/nmax' + str(nmax) + 
               '_python_sp.dat', sp_list.values, fmt='%6i')

    sys.exit()
    
    tb_list = create_tb_list (nmax, sp_list, tb_list)
    
    tb_list.number = tb_list.number.astype(int)
    tb_list.n1 = tb_list.n1.astype(int)
    tb_list.l1 = tb_list.l1.astype(int)
    tb_list.n2 = tb_list.n2.astype(int)
    tb_list.l2 = tb_list.l2.astype(int)
    
    tb_list.sp1 = tb_list.sp1.astype(int)
    tb_list.sp2 = tb_list.sp2.astype(int)
    
    np.savetxt('../me_files/ref_files/nmax' + str(nmax) + 
               '_python_tb.dat', tb_list.values, fmt='%6.2f')
    
    print ('Memory size sp list {} bytes'.format(sys.getsizeof(sp_list)))
    print ('Memory size tb list {} bytes'.format(sys.getsizeof(tb_list)))

    (Pnum, P_spterms) = P_flag_num (tb_list)

    np.savetxt('../flag_files/nmax' + str (nmax) +
               '_python_pflag.dat', P_spterms, fmt='%6i')
    
    print(sp_list) ; print ('\n')
    print(tb_list)