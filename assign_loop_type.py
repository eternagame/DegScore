import numpy as np
import re

def convert_structure_to_bps(secstruct):

    bps = []

    left_delimiters = ['(','{','[']
    right_delimiters = [')','}',']']

    for (left_delim, right_delim) in list(zip(left_delimiters, right_delimiters)):

        left_list = []
        for i, char in enumerate(secstruct):
            if char == left_delim:
                left_list.append(i)

            elif char == right_delim:
                bps.append([left_list[-1],i])
                left_list = left_list[:-1]

        assert len(left_list)==0

    return bps

def secstruct_to_partner(secstruct):
    '''Convert secondary structure string to partner array.
    I.E. ((.)) -> [4,3,-1,1,0]
    '''
    bps = convert_structure_to_bps(secstruct)
    partner_vec = -1*np.ones([len(secstruct)]) 

    for (i,j) in bps:
        partner_vec[i] = j
        partner_vec[j] = i

    return partner_vec

def write_loop_assignments(dbn_string):
    '''Input: dot-parenthesis string
    Output: bpRNA-style loop type assignments'''
    
    pair_partners = secstruct_to_partner(dbn_string)
    
    #print(pair_partners)
    bprna_string=['u']*len(dbn_string)

    # assign stems
    for s_ind, s in enumerate(dbn_string):
        if s != '.':
            bprna_string[s_ind] = 'S'
                
    # get loop regions
    
    while 'u' in ''.join(bprna_string):
        #print(''.join(bprna_string))

        obj = re.search(r"uu*", ''.join(bprna_string))
        start_ind, end_ind = obj.start(), obj.end()
        
        n_open_hps = dbn_string[:start_ind].count(')') - dbn_string[:start_ind].count('(')
        
        if n_open_hps == 0:
            bprna_string[start_ind:end_ind] = 'E'*(end_ind-start_ind)

        else:

            last_stem_pairing = int(pair_partners[start_ind - 1])
            next_stem_pairing = int(pair_partners[end_ind ])
            
            if last_stem_pairing == end_ind:
                bprna_string[start_ind:end_ind] = 'H'*(end_ind-start_ind)

            elif (last_stem_pairing - 1 == next_stem_pairing):
                bprna_string[start_ind:end_ind] = 'B'*(end_ind-start_ind)
                
            elif dbn_string[start_ind-1]==')' and dbn_string[end_ind]=='(':
                bprna_string[start_ind:end_ind] = 'M'*(end_ind-start_ind)
                
            else:
                if dbn_string[next_stem_pairing+1:last_stem_pairing] == '.'*(last_stem_pairing - next_stem_pairing-1):
                    bprna_string[start_ind:end_ind] = 'I'*(end_ind-start_ind)
                    bprna_string[next_stem_pairing+1:last_stem_pairing] = 'I'*(last_stem_pairing - next_stem_pairing-1)

                else:
                    bprna_string[start_ind:end_ind] = 'M'*(end_ind - start_ind)
    return ''.join(bprna_string)