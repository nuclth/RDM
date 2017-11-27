#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 11:03:52 2017

@author: alex
"""

# standard libraries import
import sys

# import third party libraries
import pandas as pd
import numpy as np

    
def create_sp_list (nmax, sp_list):
    """Function to create a list of single-particle orbitals for a given Nmax
    truncation in m-scheme harmonic oscillator basis. 
    
    Arguments:
    nmax    - nmax truncation size
    sp_list - dataframe to hold single-particle orbital quantum numbers 
    """
    
    # initialize values
    new_ncount = True    
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
                                  'j': twoj, 'mj' : twomj}, ignore_index=True)
    
        # increment sp orbital counter
        num += 1
    
        # cycle through mj values for a given n,l,j
        if twomj < twoj:
            twomj += 2
            continue
    
        # cylce through j values for a given n,l
        if twoj > 2 * abs (l - 1/2):
            twoj = twoj - 2
            twomj = -twoj
            continue

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
        
      
        l = l_list[icount]
        n = n_list[icount]
        
        icount += 1

        # reset j,mj for new n,l values
        twoj = 2 * abs (l + 1/2)
        twomj = -twoj
    

    
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



def sort_sp (sp_list):
    """Sort sp_list dataframe by j, mj values in ascending order. Then 
    
    Argument:
        sp_list - dataframe to hold sp orbital quantum numbers
    """


    print (sp_list)
    
    sp_list = sp_list.sort_values(['j','mj'])
    
    print (sp_list)
        
    sp_list = sort_sp_l (sp_list)
    
    print (sp_list)

    sys.exit()

    
    return sp_list


def sort_sp_l (sp_list):
    
    even_parity = sp_list['l'] % 2 == 0
    odd_parity = sp_list['l'] % 2 == 1
    
    shifted_list = sp_list[even_parity]
    shifted_list = shifted_list.append(sp_list[odd_parity])

    # drop index
    shifted_list = shifted_list.reset_index(drop=True)

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
    
#    print (sp_list)
    
    return sp_list

def sort_tb (tb_list):

    even_parity = (tb_list['l1'] + tb_list['l2']) % 2 == 0
    odd_parity  = (tb_list['l1'] + tb_list['l2']) % 2 == 1

    num_even = sum (a for a in even_parity)
    num_odd = sum (a for a in odd_parity)
    
#    print (num_even, num_odd)
    
    shifted_list = tb_list[even_parity]
    shifted_list = shifted_list.append(tb_list[odd_parity])

    shifted_list['blocks'] = 0
    
    shifted_list = shifted_list.reset_index(drop=True)
    
    shifted_list['blocks'][0] = num_even
    shifted_list['blocks'][1] = num_odd


    return shifted_list

def no_flag (sp_list):
      
    sp_vals = pd.DataFrame(columns=['b1','b2','m1','m2','sp1','sp2', 'block'])

    block_list = sp_list[['blocks']]
    
    m1 = -1
    m2 = 0 
  
    block = 0
  
    for num in block_list['blocks']:
        
        if num < 1:
            continue

        block += 1

        for val1 in range (1, num+1, 1):
            m1 += 1
            m2 = m1
            
            sp_num1 = sp_list.iloc[m1].number
            
            
            for val2 in range (val1, num+1, 1):
                sp_num2 = sp_list.iloc[m2].number
                sp_vals = sp_vals.append({'b1': val1, 'b2': val2, 
                                          'm1': m1, 'm2': m2, 'sp1' : sp_num1,
                                          'sp2' : sp_num2, 'block': block
                                          }, ignore_index=True)
                m2 += 1
                
                

    return sp_vals




def h2_flag (tb_list, sp_relational):
    
    sp_vals = pd.DataFrame(columns=['b1','b2','m1','m2',
                                    'sp1','sp2','sp3','sp4', 
                                    'm_sp1', 'm_sp2', 'm_sp3', 'm_sp4',
                                    'block'])

    block_list = tb_list[['blocks']]
    
    m1 = -1
    m2 = 0 
    block = 0
    
    for num in block_list['blocks']:
        
        if num < 1:
            continue

        block += 1

        for val1 in range (1, num+1, 1):
            m1 += 1
            m2 = m1
            
            sp1 = tb_list.iloc[m1].sp1
            sp2 = tb_list.iloc[m1].sp2
            
            m_sp1 = sp_relational[sp_relational['alex_sp'] == sp1].morten_sp.iloc[0]
            m_sp2 = sp_relational[sp_relational['alex_sp'] == sp2].morten_sp.iloc[0]
            
            for val2 in range (val1, num+1, 1):
                sp3 = tb_list.iloc[m2].sp1
                sp4 = tb_list.iloc[m2].sp2
                
                m_sp3 = sp_relational[sp_relational['alex_sp'] == sp3].morten_sp.iloc[0]
                m_sp4 = sp_relational[sp_relational['alex_sp'] == sp4].morten_sp.iloc[0]
                
                sp_vals = sp_vals.append({'b1': val1, 'b2': val2, 
                                          'm1': m1, 'm2': m2, 'sp1': sp1,
                                          'sp2': sp2, 'sp3':sp3, 
                                          'sp4' : sp4, 'm_sp1' : m_sp1, 'm_sp2': m_sp2,
                                          'm_sp3' : m_sp3, 'm_sp4': m_sp4, 'block': block
                                          }, ignore_index=True)
                m2 += 1
         
#    print (sp_vals)
#    sys.exit()
         
    return sp_vals


def p_flag (no_flag, h2_flag):
    
    
#    sp_vals = pd.DataFrame({'sp1': tb_list.sp1, 'sp2':tb_list.sp2})

    sp_vals = pd.DataFrame(columns=['ob_b1','ob_b2','ob_block',
                                    'tb_b1','tb_b2','tb_block','new'])
    repeats = []
    

    for sp1, sp3 in zip(h2_flag.sp1, h2_flag.sp3):
        
        
        if (sp1, sp3) in repeats:
            continue
        
        repeats.append((sp1,sp3))
        
#        print (repeats)
        
        new = True

        

#        print ('SP1, SP3 {} {}'.format(sp1, sp3))
        
        ob_row = no_flag[(no_flag['sp1'] == sp1) & (no_flag['sp2'] == sp3)]
        
        if ob_row.empty:
            continue
        
#        print ('OB ROW')
#        print (ob_row)
        
        ob_b1 = ob_row.b1.iloc[0]
        ob_b2 = ob_row.b2.iloc[0]
        ob_block = ob_row.block.iloc[0]
        
        loop_repeats = []
        
        for sp2, sp4 in zip(h2_flag.sp2, h2_flag.sp4):
            
            if (sp2 != sp4):
                continue


            if (sp2, sp4) in loop_repeats:
                continue
            
            
            loop_repeats.append((sp2,sp4))
            
            tb_row = h2_flag[(h2_flag['sp1'] == sp1) & (h2_flag['sp2'] == sp2)
            & (h2_flag['sp3'] == sp3) & (h2_flag['sp4'] == sp4)]

            if tb_row.empty:
                continue
            
#            print ('TB_row')
#            print (tb_row)
            
            tb_b1 = tb_row.b1.iloc[0]
            tb_b2 = tb_row.b2.iloc[0]
            tb_block = tb_row.block.iloc[0]

#            print ('ENTRIES {} {} {} {} {} {} {}'.format(
#                    ob_b1, ob_b2, ob_block, tb_b1, tb_b2, tb_block, new))

            sp_vals = sp_vals.append({'ob_b1':ob_b1, 'ob_b2':ob_b2, 
                                      'ob_block':ob_block, 'tb_b1':tb_b1, 'tb_b2':tb_b2,
                                      'tb_block':tb_block, 'new':new}, ignore_index=True)
            
            new = False
        
#        print (no_flag)

        
#        for i2, sp3, sp4 in zip (sp_vals.index, sp_vals.sp1, sp_vals.sp2):
#            if sp2 == sp4 and i1 <= i2 and (sp1,sp3) not in repeats:
#                repeats.append((sp1, sp3))
#                print (sp1, sp2, sp3, sp4)
#                count+=1


    return sp_vals


def sort_pflag (p_terms):
    
    sp_vals = pd.DataFrame(columns=['ob_b1','ob_b2','ob_block',
                                    'tb_b1','tb_b2','tb_block','new'])

    hold_vals = pd.DataFrame(columns=['ob_b1','ob_b2','ob_block',
                                    'tb_b1','tb_b2','tb_block','new'])
    
    start = True   
    
#    p_terms['new'].iloc[len(p_terms)-1] = False

    for index, term in zip(p_terms.index, p_terms['new']):
        
        if term == False or start:
            hold_vals = hold_vals.append(p_terms.iloc[index], ignore_index=True)
        
#        print (index, term)
        
        if term == True and not start:
            hold_vals = hold_vals.sort_values(['tb_block', 'tb_b1','tb_b2'])
            hold_vals['new'] = False
            hold_vals['new'].iloc[0] = True
            sp_vals = sp_vals.append(hold_vals, ignore_index=True)
            hold_vals.drop(hold_vals.index, inplace=True)
            hold_vals = hold_vals.append(p_terms.iloc[index], ignore_index=True)

        
        if index+1 == len(p_terms):
            hold_vals = hold_vals.sort_values(['tb_block', 'tb_b1','tb_b2'])
            hold_vals['new'] = False
            hold_vals['new'].iloc[0] = True
            sp_vals = sp_vals.append(hold_vals, ignore_index=True)
            
        start = False
     
    
    return sp_vals


def sp_relational_db_morten (nmax, sp_list):
    
    morten_splist = pd.DataFrame(columns = ['number', 'n', 'l', '2j', '2mj'])
    
    with open('../me_files/ref_files/nmax' + str(nmax) + "_spm.dat", 'r') as f:
        content = f.readlines()
        for x in content:
            row = x.split()
#            print (row)
            if row[0] == '#': continue 
            num   = int(row[2])
            n     = int(row[3])
            l     = int(row[4])
            twoj  = int(row[5])
            twomj = int(row[6])
            
            morten_splist = morten_splist.append({'number':num, 'n':n, 'l':l,'2j':twoj,
                                                  '2mj':twomj}, ignore_index=True)
#    print (morten_splist)
     
    relational_db = pd.DataFrame(columns = ['morten_sp', 'alex_sp'])
    

    for index, number in zip(morten_splist.index, morten_splist['number']):
        relational_db = relational_db.append({'morten_sp':number, 
                                              'alex_sp': index+1}, ignore_index=True)
    
    return relational_db




if __name__ == '__main__':
    
    sp_list = pd.DataFrame(columns= ['number', 'n', 'l', 'j', 'mj'])
    tb_list = pd.DataFrame(columns= ['number', 'n1', 'l1', 'j1', 'mj1', 
                                     'n2', 'l2', 'j2', 'mj2', 'sp1', 'sp2'])

    
    nmax = 2

    sp_list = create_sp_list (nmax, sp_list)


    sp_list.number = sp_list.number.astype(int)
    sp_list.n = sp_list.n.astype(int)
    sp_list.l = sp_list.l.astype(int)
    sp_list.j = sp_list.j.astype(int)
    sp_list.mj = sp_list.mj.astype(int)

    
    sp_list = sort_sp (sp_list)
    
    sp_list = add_blocks (sp_list)
    
    print (sp_list)


    sp_relational = sp_relational_db_morten (nmax, sp_list)

    
    np.savetxt('../me_files/ref_files/nmax' + str(nmax) + 
               '_python_sp.dat', sp_list.values, fmt='%6i')

    no_terms = no_flag (sp_list)

    np.savetxt('../flag_files/nmax' + str(nmax) + 
               '_python_noflag.dat', no_terms, fmt='%6i')


    print (no_terms)

    print ("OBME SUCCESS!")
    
    tb_list = create_tb_list (nmax, sp_list, tb_list)

    tb_list = sort_tb (tb_list)

    tb_list.number = tb_list.number.astype(int)
    tb_list.n1 = tb_list.n1.astype(int)
    tb_list.l1 = tb_list.l1.astype(int)
    tb_list.n2 = tb_list.n2.astype(int)
    tb_list.l2 = tb_list.l2.astype(int)
    
    tb_list.sp1 = tb_list.sp1.astype(int)
    tb_list.sp2 = tb_list.sp2.astype(int)

    print (tb_list)
    
    np.savetxt('../me_files/ref_files/nmax' + str(nmax) + 
               '_python_tb.dat', tb_list.values, fmt='%6.2f')
    
    print ('Memory size sp list {} bytes'.format(sys.getsizeof(sp_list)))
    print ('Memory size tb list {} bytes'.format(sys.getsizeof(tb_list)))

    h2_terms = h2_flag (tb_list, sp_relational)

    h2_terms = h2_terms.astype(int)

    print (h2_terms)

    np.savetxt('../flag_files/nmax' + str (nmax) +
               '_python_h2flag.dat', h2_terms, fmt='%6i')

    p_terms = p_flag (no_terms, h2_terms)

#    print (p_terms)

    p_terms = sort_pflag (p_terms)

    print (p_terms)

    np.savetxt('../flag_files/nmax' + str (nmax) +
               '_python_pflag.dat', p_terms, fmt='%6i')
    
    
    print ("FINISHED")
#    sys.exit()
    
#    print(sp_list) ; print ('\n')
#    print(tb_list)