##############################################################################
# Copyright (c) 2023, Ruben Becker, Sajjad Ghobadi
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * The name of the author may not be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
# EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##############################################################################

"""Module proxy for the cpp implementation."""

import os
import uuid
import influence_max as im
import operator


def grdy_im(G, graph_file, k, file_cpp):
    """Execute TIM implementation by Tang et al (2014) from cpp implementations.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the seed set)    
    file_cpp : str
        string specifying the name of output file written in cpp implementation.
        
    Returns
    -------
    S : set of seed nodes
    running_time : The execution time of TIM to return the computed solution

    """
    file_cpp = file_cpp + '_' + str(uuid.uuid1()) + '.txt'
    tim_command = 'cd' + ' ' + './TIM_final' + '&& ./adding_links -model IC -dataset' + ' ' +'.'+ graph_file + ' ' + '-epsilon 0.1 -k' + ' ' + str(k) + ' ' + '-B' + ' ' + str(1) + ' ' + '-setting' + ' '+ 'tim' + ' ' + '-file_name' + ' ' + str(file_cpp)# + ' ' + '-folder'# + ' ' + str(folder_name)                  
    file = './TIM_final/' + file_cpp
    
    process = os.popen(tim_command)
    output = process.read()
    process.close()
    print(output) 
    
    S = set()
    with open(file) as f:
        line1 = f.readline()
        running_time = float(line1)
            
        line2 = f.readline()
        nodes = map(str, line2.split())
        for v in nodes:
            S.add(int(v))

    os.remove(file)
    return S, running_time
    

def mult_weight(G, graph_file, k, file_cpp):
    """Execute multiplicative weight routine for the set-based problem of Becker et al (2022) from cpp implementation.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the ssed set)    
    file_cpp : str
        string specifying the name of output file written in cpp implementation.
        
    Returns
    -------
    p : dict
        dictionary with keys integers encoding sets and values probabilities.
    running_time : The execution time of mult_weight to return the computed solution

    """
    file_cpp = file_cpp + '_' + str(uuid.uuid1()) + '.txt'    
    tim_command = 'cd' + ' ' + './TIM_final' + '&& ./adding_links -model IC -dataset' + ' ' +'.' + graph_file + ' ' + '-epsilon 0.1 -k' + ' ' + str(k) + ' ' + '-B' + ' ' + str(1) + ' ' + '-setting' + ' '+ 'set' + ' ' + '-file_name' + ' ' + str(file_cpp)
    file = './TIM_final/'+ file_cpp
    process = os.popen(tim_command)
    output = process.read()
    process.close()
    print(output)   
        
    p = {}
    with open(file) as f:
        line1 = f.readline()
        running_time = float(line1)
            
    with open(file) as f:
        lines = f.readlines()[1:]
        for line in lines:
            *sets , pro = map(float, line.split())
            if len(sets) == 0.0:
                p[0.0] = pro
            else:
                p[im.set_to_number(G, sets)] = pro
    os.remove(file)
    sum_p = sum(p.values())
    if sum_p != 1:
        p = dict(sorted(p.items(), key=lambda x: x[1]))
        if sum_p > 1:
            x = sum_p - 1
            for S_pro in p:
                if p[S_pro] >= x:
                      p[S_pro] -= x
                      break
        if sum_p < 1:
            x = 1 - sum_p
            max_S = max(p.items(), key=operator.itemgetter(1))[0]
            p[max_S] += x;                       
    return p, running_time


def to_minC_infl(graph_file, k, B, file_cpp):
    """ Execute the implementation of to_minC_infl proposed in the paper from cpp implementations.
    """
    return read_file_ae(graph_file, k, B, file_cpp, 'to_minC_infl')
    
def to_minC_min(graph_file, k, B, file_cpp): 
    """ Execute the implementation of to_minC_min proposed in the paper from cpp implementations
    """
    return read_file_ae(graph_file, k, B, file_cpp, 'to_minC_min')   

def grdy_al(graph_file, k, B, file_cpp):
    """ Execute the implementation of grdy_al proposed in the paper from cpp implementations
    """
    return read_file_ae(graph_file, k, B, file_cpp, 'grdy_al')

def random(graph_file, k, B, file_cpp):
    """ Execute the implementation of random from cpp implementations
    """
    return read_file_ae(graph_file, k, B, file_cpp, 'random')

def max_weight(graph_file, k, B, file_cpp): 
    """ Execute the implementation of max_weight from cpp implementations
    """
    return read_file_ae(graph_file, k, B, file_cpp, 'max_weight')


def read_file_ae(graph_file, k, B, file_cpp, setting): 
    """ Read the output file written by cpp implementations.

    Parameters
    ----------
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the seed set)
    B : int
        Budget (size of the non-edges to be added to the graph) 
    file_cpp : str
        string specifying the name of output file written by cpp implementation.
    setting : str
            string specifying the name of function to be executed in cpp implementations.
        
    Returns
    -------
    graphname : str
        name of the graph
    fcn_name : str
        name of the fuction to be executed
    n : int
        number of nodes of the graph
    m : int
        number of edges of the graph
    k : int
        Budget (size of the seed set)
    B : int
        Budget (size of the non-edges to be added to the graph)
    running_time : double
        The execution time of fcn_name to return the computed solution
    spread : int
        expected number of reached nodes from seed set
    min_pro : double
        mimimum probability of being reached of communities from seed set
    max_pro : double
        maximum probability of being reached of communities from seed set
    """
    file_cpp = file_cpp + '_' + str(uuid.uuid1()) + '.txt'
    tim_command = 'cd' + ' ' + './TIM_final' + '&& ./adding_links -model IC -dataset' + ' ' +'.'+ graph_file + ' ' + '-epsilon 0.1 -k' + ' ' + str(k) + ' ' + '-B' + ' ' + str(B) + ' ' + '-setting' + ' '+ setting + ' ' + '-file_name' + ' ' + str(file_cpp)# + ' ' + '-folder'# + ' ' + str(folder_name)                  
    file = './TIM_final/' + file_cpp

    process = os.popen(tim_command)
    output = process.read()
    process.close()
    print(output) 
    
    with open(file) as f:
        line1 = f.readline()
        graphname = str.strip(line1)
        
        line2 = f.readline()
        fcn_name = str.strip(line2)
        
        line3 = f.readline()
        n = int(line3)

        line4 = f.readline()
        m = int(line4)  
        
        line5 = f.readline()
        k = int(line5)
        
        line6 = f.readline()
        B = int(line6)
        
        line7 = f.readline()
        running_time = float(line7)

        line8 = f.readline()
        spread = float(line8)
        
        line9 = f.readline()
        min_pro = float(line9)
        
        line10 = f.readline()
        max_pro = float(line10)
        
    os.remove(file)
    return graphname, fcn_name, n, m, k, B, running_time, spread, min_pro, max_pro

 

