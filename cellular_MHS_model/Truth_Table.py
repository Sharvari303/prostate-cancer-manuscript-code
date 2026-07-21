#!/usr/bin/env python
import itertools, json, math

def GenerateTable(p53_Network, writefile='no', filename='Truth_Table.txt'):

    Combined_Rule = []
    Init_Cond = ''
    for key, val in p53_Network.items():
        #if key[-1]: Init_Cond = Init_Cond + key[0] + '='
        single_node_rule = []
        for i in itertools.product([False, True], repeat=len(val)):
            species, weights = zip(*val)
            combelements = zip(species,i)
            sumelem = sum(map(lambda xi,yi : xi*yi, [int(j) for j in i],weights)) if species != "None" else 0 + key[1]
            for self_node in [0, 1]:
                node_status_f = 1 if sumelem > 0 else (0 if sumelem < 0 else self_node)
                if node_status_f:
                    species_mod_1 =[(species[k] if list(i)[k] else 'not ' + species[k]) for k in range(len(species))] if species[0] != "None" else []
                    species_mod_2 = [key[0]] if self_node else []
                    species_mod =species_mod_1 + species_mod_2
                    single_node_rule.append(' and '.join(species_mod))

        Combined_Rule.append('1: '+ key[0] + ' *= (' + ') or ('.join(single_node_rule) + ')' + '\r\n')

    #Init_Cond = Init_Cond + 'True'
    #text = Init_Cond + '\r\n' + ''.join(Combined_Rule)
    text = ''.join(Combined_Rule)

    if writefile is not 'no':
        with open(filename,'w') as fp: fp.write(text)

    return text

def Calculate_State(p53_Network, data):
    """
    Calculates the future state of the node based
    on the network parameters

    """

    Init_Cond = ''
    data_next = {}
    state_trn_prob = {}
    data_next.update(dict((key,val) for key,val in data.items()))
    state_trn_prob.update(dict((key,val) for key,val in data.items()))
    c = 0.001
    for key, val in p53_Network.items():
        species, weights = zip(*val)
        sumelem = sum(map(lambda xi,yi : xi*yi, [data[j] for j in species],weights)) if species[0] != "None" else 0 + key[1]
        data_next[key[0]] = 1 if sumelem > 0 else (0 if sumelem < 0 else data[key[0]])
        state_trn_prob[key[0]] = 0.5 + 0.5*math.tanh(sumelem) if sumelem > 0 else (0.5 - 0.5*math.tanh(sumelem) if sumelem < 0 else 1 - c)

    return data_next, state_trn_prob
