This is the source code used for the experimental part of the paper 
			Improving Fairness in Information Exposure by Adding Links published at AAAI-23

For execution of the experiments use
    python main.py experiment_type N

where experiment_type is one of
        [   'ba-singletons-0_0_4',
            'tsang-region-gender-0-0_4',
            'arena_0_0_2_bfs_comm_10',
            'email-Eu-core_0_0_2',
            'irvine_0_0_2_bfs_comm_10',
            'ca-GrQc_0_0.2_bfs_comm_10',
            'ca-HepTh_0_0.2_bfs_comm_10',
            'youtube_0_0_1_n_3000_no_singl'
            ]
and N is the number of experiments that are supposed to be run in parallel.

For creating the instances use
    python main_generate.py experiment_type


- For creating the networks used by Fish et al. [1] download the data sets into a subfolder 'fish' and put it into the folder 'data_set'.
- For running tests on the youtube network download the data sets into a subfolder 'youtube' and put it into the folder 'data_set'.
- For running tests on the networks used by Tsang et al. [2], download the instances into the subfolder 'networks' and put it into the folder 'tsang'.

- The files have the following content:
    - main.py: contains the infrastructure for the execution of the experiments.
    - generation.py: contains the functions that generate the different instances.
    - influence_max.py: contains the necessary implementations of influence maximization functions: computation of the influence sigma of a set (and related functions sigma_v, sigma_C, see the paper), computation of nodes reachable from a given set, generation of live edge graphs, greedy algorithms for influence maximization and grdy_maximin
    - maximin_fish.py: contains the implementations of the myopic routine of Fish et al.
    - main_generate.py: contains the infrastructure for generating the instances used in the paper.
    - print_functions.py: contains some functions for prettier printing of graphs etc.
    - moso.py: contains the function used to call the code of Tsang et al. [2]
    - cpp_proxy.py: contains the functions for executing the c++ implementations, the algorithms proposed in the paper and the multiplicative weight routine of Young [3] for the set-based problem of Becker et al.[4].

- Folder 'TIM_final' contains the c++ implementation of the proposed algorithms in the paper. It uses TIM implementation for influence maximization [5].
- The execution generates an output file within the folder 'results_ae' with the name being the experiment_type.


- The experiment_types have the following meaning, see the paper for the interpretation of the community structure types:
    - 'ba-singletons-0_0_4':                + Barabasi albert graph (parameter m=2),
                                            + singleton community structure
                                            + edge weights uniformly at random in [0,0.4]
                                            + k = 25
                                            + n = 200
                                            + B = [10, 20, 50]

    - 'tsang-region-gender-0-0_4':          + instances of Tsang et al. (2019),
                                            + community structure induced by attributes gender and region
                                            + edge weights uniformly at random in [0,0.4]
                                            + k = 25
                                            + n = 500
                                            + B = [10, 20, 50]

    - 'arena_0_0_2_bfs_comm_10':            + instances used by Fish et al. (2019)
                                            + BFS community structure with 10 communities
                                            + edge weights uniformly at random in [0,0.2]
                                            + k = 20
                                            + n = 1133
                                            + B = [10, 20, 50]

    - 'email-Eu-core_0_0_2':                + instances used by Fish et al. (2019)      
                                            + real community structure
                                            + edge weights uniformly at random in [0,0.2]
                                            + k = 20
                                            + n = 1005
                                            + B = [10, 20, 50]
                                          
    - 'irvine_0_0_2_bfs_comm_10':           + instances used by Fish et al. (2019)
                                            + BFS community structure with 10 communities
                                            + edge weights uniformly at random in [0,0.2]
                                            + k = 20
                                            + n = 1899 
                                            + B = [10, 20, 50]                               

    - 'ca-GrQc_0_0.2_bfs_comm_10':          + instances used by Fish et al. (2019)
                                            + BFS community structure with 10 communities
                                            + edge weights uniformly at random in [0,0.2]
                                            + k = 20
                                            + n = 5242
                                            + B = [10, 20, 50]
                                            
    - 'ca-HepTh_0_0.2_bfs_comm_10':         + instances used by Fish et al. (2019)
                                            + community structure induced by departments
                                            + edge weights uniformly at random in [0,0.2]
                                            + k = 20
                                            + n = 9877
                                            + B = [10, 20, 50]

    - 'youtube_0_0_1_n_3000_no_singl':      + real community structure
                                            + edge weights uniformly at random in [0,0.1]
                                            + k = 20
                                            + n = 3000
                                            + B = [10, 20, 50]
                                            

[1] Fish, Bashardoust, Boyd, Friedler, Scheidegger, Venkatasubramanian. Gaps in Information Access in Social Networks. WWW 2019.
    http://sorelle.friedler.net/papers/access_networks_www19.pdf
[2] Tsang, Wilder, Rice, Tambe, Zick. Group-Fairness in Influence Maximization. IJCAI 2019.
    https://www.ijcai.org/Proceedings/2019/0831.pdf
[3] Young, Neal E. "Randomized Rounding Without Solving the Linear Program." SODA. Vol. 95. 1995.
[4] Becker, D’Angelo, Ghobadi, and Gilbert. Fairness in Influence Maximization through Randomization. J. Artif. Intell. Res., 73: 1251–1283. 2022.
[5] Tang, Youze, Xiaokui Xiao, and Yanchen Shi. "Influence maximization: Near-optimal time complexity meets practical efficiency." Proceedings of the 2014 ACM SIGMOD international conference on Management of data. 2014.
