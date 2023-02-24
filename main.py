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


"""Main module that executes the experiments."""

import os
import sys
import threading
import time
import uuid
import random
import networkx as nx
import numpy as np
from numpy.random import choice, seed
import influence_max as im
import cpp_proxy as cpp
import print_functions as pf
import generation as gen
import maximin_fish as mf
from tsang import moso as tm
from tsang.algorithms import rounding

def sample_sets(G, vector, times, type):
    """Sample times many sets from vector depending on type.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    vector : depends on type
        If type is deterministic, vector is a set.
        If type is node_based or swap_rounding, vector is a dictionary with
            keys being the nodes and values the probability of the respective
            node.
        If type is set_based, vector is a dictionary with keys being ints
            representing the sets (binary) and values are the probability
            of the respective set.
    times : int
        How many times to sample a set.
    type : str
        The type that specifies how to sample the sets, this depends on whether
        this is a deterministic, node_based, or set_based problem.

    Returns
    -------
    list of lists
        The list of sampled sets.

    """
    if type == 'deterministic':
        sets = [vector for _ in range(times)]
    elif type == 'swap_rounding':
        x_items_list = sorted(vector.items())
        x = np.array([x_items[1] for x_items in x_items_list])
        rounded_xs = [rounding(x) for _ in range(times)]
        sets = [[v for v in G.nodes() if rounded_xs[i][v]]
                for i in range(times)]
    elif type == 'set_based':
        sets = [im.number_to_set(G, choice(list(vector.keys()),
                                           p=list(vector.values())))
                for _ in range(times)]
    else:
        print("Error: Unknown option.", type)
        assert(False)
    return sets


def comp_ex_post(G, solution, fct_name):
    """Compute fairness abd spread values by sampling one or several sets.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    solution : depends on fct_name
        The computed solution.
    fct_name : str
        The name of the function that computed solution.

    Returns
    -------
    int
        The ex_post fairness value and spread obtained for one and several sampled sets.

    """
    sets = sample_sets(G, solution, 1,
                       ex_post_sampling_types[fct_name])
    
    if fct_name in deterministic_algorithms + ['grdy_im']:
        sets_pof = sets
    else:
        sets_pof = sample_sets(G, solution, 100,
                                ex_post_sampling_types[fct_name])
    sigma_sets_pof = [im.sigma(G, S, weight='weight')
                              for S in sets_pof]
    sigma_pof = np.mean(sigma_sets_pof)
    sigma_sets_ep = [im.sigma(G, S, weight='weight')
                             for S in sets]
    sigma_ep = np.mean(sigma_sets_ep)
    min_probs_sets = [min(im.sigma_C(G, S).values())
                      for S in sets]
    ex_post_apx_val = np.mean(min_probs_sets)

    return ex_post_apx_val, sigma_pof, sigma_ep


def comp_ex_ante(G, solution, fct_name):
    """Compute ex_ante values ((0.1, 0.1)-approximately).

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    solution : depends on fct_name
        The computed solution.
    fct_name : str
        The name of the function that computed solution.

    Returns
    -------
    int
        The ex_ante value obtained by (0.1, 0.1)-approximating in case of
        the node_based problem and exact computation in case of the set_based
        problem.

    """
    if fct_name in deterministic_algorithms + ['grdy_im']:
        comm_probs = im.sigma_C(G, solution)
        return min(comm_probs.values()), max(comm_probs.values())
    else:
        assert(fct_name in probabilistic_algorithms + cpp_algos + ['mult_weight'])
        if fct_name in ['mult_weight']:
            comm_probs = im.sigma_C_p(G, solution)
        elif fct_name in ['moso']:
            comm_probs = im.sigma_C_x(G, solution, 0.1, 0.1)
        else:
            print("Error: Unknown option:", fct_name)
            assert(False)
    return min(comm_probs.values()) , max(comm_probs.values())


def read_graph(graph_file):
    """Read graph from files.

    Parameters
    ----------
    graph_file : The path of the folder containing graph files

    Returns
    -------
    G : networkx.DiGraph
        The underlying considered instance.
    """
    G = nx.read_edgelist(graph_file + '/graph_ic.txt', nodetype = int,
                          data=(('p',float),), create_using = nx.DiGraph())
      
    gen.init_communities(G)    
    G.graph['nodes_in_comms'] = set()    
    
    comm_file = open(graph_file + '/community.txt','r')
    lines = comm_file.readlines()[1:]    
    for line in lines:
        comm , *nodes = map(str, line.split())
        G.graph['communities'][comm] = [int(v) for v in nodes]
        
        for u in nodes:
            G.nodes[int(u)]['communities'].append(comm)
            G.graph['nodes_in_comms'].add(int(u))        
    
    comm_file.close()    
    gen.set_unit_node_weights(G)
             
    return G


def execute(function, graph_file, k, rep_per_graph_i, num, B): 
    """Execute function on G and writes results to out_file (global).

    Parameters
    ----------
    function : function that takes a networkx graph and returns a solution.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the seed set)
    rep_per_graph_i : int
        The number showing repetition per graph.
    num : int
        Number specifying the name of output file (used for cpp implementations).
    B : int
            Budget (size of the non-edges to be added to the graph)
    """
    
    fct_name = function.__name__
        
    G = read_graph(graph_file)
    G.graph['T'] = no_le_graphs
    G.graph['k'] = k
    
    start = time.time()
    G.graph['leg'], G.graph['leg_tsang'] = im.gen_leg(G) 
    
    if fct_name in python_algos:        
        solution = function(G)
        running_time = time.time() - start
    elif fct_name in ['mult_weight', 'grdy_im']:
        solution, ex_time = function(G, graph_file, k, str(num))
        running_time = time.time() - start
    else:
        graphname, fct_name, n, m, k, B, running_time, exp_spread, ex_post_min, max_pro = function(graph_file, k, B, fct_name + '_' + str(num)) # folder_cpp 
        spread = exp_spread
        ex_ante_min = ex_post_min
    
    if fct_name in python_algos + ['mult_weight', 'grdy_im']:
        G.graph['leg'], G.graph['leg_tsang'] = im.gen_leg(G) 
    
        n, m = len(G), len(G.edges())
        ex_ante_min, ex_ante_max = comp_ex_ante(G, solution, fct_name)  
        ex_post_min, exp_spread, spread = comp_ex_post(G, solution, fct_name)
        B = 0
      
    x = graph_file.split("/")
    graphname = x[3]
            
    res = (graphname,
           fct_name,
           n, 
           m,
           k,
           B,
           running_time,
           spread,
           exp_spread,
           ex_post_min,
           ex_ante_min)

    with open('./results_ae/'+out_file, 'a') as f:
        f.write('\t\t\t'.join('%s' % x for x in res) + '\n')
    

def generate_executables(
        functions,
        instance_folder,
        k_range,
        rep_per_graph,
        no_le_graphs,
        B_range
        ):     
    """Generate executables by generating graphs and combining them.

    Parameters
    ----------
    functions : list
        List of all the functions that are supposed to be executed.
    instance_folder : str
        Name of the folder containing graph files.
    k_range : list
        Range of k-values to be tested.
    rep_per_graph : ing
        Number of times the execution should be repeated per graph.
    no_le_graphs : int
        Number of live-edge graphs to be generated and used in all following
        computations.
    B_range : list
        Range of B-values to be tested.

    Returns
    -------
    list
        List of executables, i.e., (function, graph)-pairs.

    """
    num = 0.0    
    executables = []
    instance_names = sorted([f for f in os.listdir('./data_set/' + instance_folder) #  + '/' 
                            if not f.startswith(".")])
                  
    for graph_folder_name in instance_names:
        for rep_per_graph_i in range(rep_per_graph):
            for k in k_range:
                for function in functions:
                    fct_name = function.__name__
                    num += 1
                    if fct_name in cpp_algos:
                        for B in B_range:
                            executables.append((function,
                                                './data_set/' + instance_folder + '/' + graph_folder_name,
                                                k, rep_per_graph_i, num, B))
                    else:
                        executables.append((function,
                                            './data_set/' + instance_folder + '/' + graph_folder_name,
                                            k, rep_per_graph_i, num, 0))
                        
            print('collected', len(executables), 'executables')
    return executables


#############
# main
#############

# forbid python 2 usage
version = sys.version_info[0]
if version == 2:
    sys.exit("This script shouldn't be run by python 2 ")

# do not set seed specifically
s = None
seed(s)

# the following lists are used to specify the way in which sets are sampled
# for the respective algorithms

deterministic_algorithms = [
    'myopic',
    'grdy_maximin',
    ]

probabilistic_algorithms = [
    'moso']

python_algos = deterministic_algorithms + probabilistic_algorithms
cpp_algos = [
    'to_minC_infl',
    'to_minC_min',
    'grdy_al',
    'max_weight',
    'random',
    ]

sampling_types = [
    'deterministic',
    'set_based',
    'swap_rounding'
]

ex_ante_sampling_types = {
    'mult_weight': 'set_based',
    'moso': 'node_based'}
for alg in deterministic_algorithms + ['grdy_im']:
    ex_ante_sampling_types[alg] = 'deterministic'

ex_post_sampling_types = ex_ante_sampling_types
ex_post_sampling_types['moso'] = 'swap_rounding'

print('++++++++++++++++++++++++++++++++++++++++++++++++++++')
print('++++++ Expecting experiment_type as argument. ++++++')
print('++++++++++++++++++++++++++++++++++++++++++++++++++++')


# read number of desired processes from the shell
experiment_type = sys.argv[1]
if len(sys.argv) == 3:
    number_of_processes = int(sys.argv[2])
else:
    number_of_processes = 1

# default values for experiments     
functions = [cpp.grdy_im,
             cpp.mult_weight,
             mf.myopic,
             im.grdy_maximin,
             tm.moso,
             cpp.to_minC_infl,
             cpp.to_minC_min,
             cpp.random,
             cpp.max_weight,
             cpp.grdy_al
             ]
no_le_graphs = 100
rep_per_graph = 5

# specify experiment dependent parameters
if experiment_type == 'ba-singletons-0_0_4':
    instance_folder = 'ba-singletons-0_0_4'
    k_range = [25]
    B_range = [10, 20, 50]    
         
elif experiment_type == 'tsang-region-gender-0-0_4':           
    instance_folder = 'tsang-region-gender-0-0_4'
    k_range = [25]
    B_range = [10, 20, 50]
   
elif experiment_type == 'arena_0_0_2_bfs_comm_10':           
    instance_folder = 'arena_0_0_2_bfs_comm_10'
    k_range = [20]
    B_range = [10, 20, 50]
     
elif experiment_type == 'email-Eu-core_0_0_2':           
    instance_folder = 'email-Eu-core_0_0_2'
    k_range = [20]
    B_range = [10, 20, 50]
    
elif experiment_type == 'irvine_0_0_2_bfs_comm_10':           
    instance_folder = 'irvine_0_0_2_bfs_comm_10'
    k_range = [20]
    B_range = [10, 20, 50]

elif experiment_type == 'ca-GrQc_0_0.2_bfs_comm_10':           
    instance_folder = 'ca-GrQc_0_0.2_bfs_comm_10'
    k_range = [20]
    B_range = [10, 20, 50]
    
elif experiment_type == 'ca-HepTh_0_0.2_bfs_comm_10':           
    instance_folder = 'ca-HepTh_0_0.2_bfs_comm_10'
    k_range = [20]
    B_range = [10, 20, 50]
        
elif experiment_type == 'youtube_0_0_1_n_3000_no_singl':           
    instance_folder = 'youtube_0_0_1_n_3000_no_singl'
    k_range = [20]
    B_range = [10, 20, 50]
       
else:
    print("Error: Unknown option.")
    assert(0)
    
    
# create output file with header
folder = './results_ae/'
out_file = experiment_type + '.txt'
if os.path.exists(folder + out_file):
    out_file = out_file[:-4] + '_' + str(int(time.time())) + '.txt'
print('Output is written to:', out_file, '\n')
header = [
    'graphname',
    'algorithm',
    'n',
    'm',
    'k',
    'B',
    'running_time',
    'spread',
    'exp_spread',
    'ex_post_min',
    'ex_ante_min'
    ]

with open(folder+out_file, 'a') as f:
    f.write('\t\t\t'.join('%s' % x for x in header) + '\n')
    
    
path_make = 'cd' + ' ' + './TIM_final' + '&& make'
process_make = os.popen(path_make)
output_make = process_make.read()
process_make.close()
print(output_make)

#generate the various experiments
executables = generate_executables(
    functions,
    instance_folder,
    k_range,
    rep_per_graph,
    no_le_graphs,
    B_range
    )

#run experiments in parallel (if number_of_processes > 1)

thread_list = []
count = 0.0
for executable in executables:
    thread_list.insert(
        0,
        threading.Thread(
            target=execute,
            args=(executable[0], executable[1], executable[2],
                  executable[3], executable[4] , executable[5])))  

while thread_list:
    if threading.active_count() < number_of_processes + 1:
        thread_list.pop().start()





