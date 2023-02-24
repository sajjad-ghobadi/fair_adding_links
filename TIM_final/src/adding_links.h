//typedef double (*pf)(int,int);
#include <math.h> 
#include "head.h"
#include <map>
#include <random>
#include <numeric>
#include <algorithm>
#include <set>

class adding_links: public TimGraph
{
    public:
        int i_multi = 1;
        pair<vector<pair<int, int>>, double> set_max_min_pro;

        adding_links(string folder, string graph_file, string setting):TimGraph(folder, graph_file, setting){
                }

        vector<pair<int, int>> grdy_al(){
            // Implement the grdy_al algorithm proposed in the paper.
            //
            // return :
            //      a set of non-edges to be added to the graph

            vector<pair<int, int>> E_prime;
            set_of_hyperGT_hyperG_rootR_p();
            map<int, vector<int>> gT_temp;
            map<int, vector<double>> probT_temp;
            vector<vector<vector<int>>> set_of_hyperGT_temp;
            vector<map<int, vector<int>>> set_of_hyperG_temp;
            vector<map<int, vector<int>>> set_of_rootR_temp;
            vector<map<vector<int>, double>> p_temp;

            gT_temp = gT;
            probT_temp = probT;
            set_of_hyperGT_temp = set_of_hyperGT;
            set_of_hyperG_temp = set_of_hyperG;
            set_of_rootR_temp = set_of_rootR; 
            p_temp = p;
           
            while((int)E_prime.size() < B){
                pair<pair<int, int>, double> best_e_bar;
                best_e_bar.second = -1.0;
                map<int, double> pro_v_is_seed;
                set<int> seedNodes;
                for(auto item:p){
                    for(auto seed_p:item){
                        for(int v:seed_p.first){
                            seedNodes.insert(v);
                            pro_v_is_seed[v] += seed_p.second;
                        }
                    } 
                }

                vector<pair<int, int>> new_E_bar;            
                // consider non-edges that are incident to the seed nodes and belong to different communities
                for(auto node_in_nghbr:gT_non_edges){
                    int max_index;
                    double max = 0.0;
                    string C1 = node_communities[node_in_nghbr.first][0];
                    int seed_count = 0.0;
                    for(int i = 0; i < (int)node_in_nghbr.second.size(); i++){
                        string C2 = node_communities[node_in_nghbr.second[i]][0];
                        if(C1 != C2){
                            if(find(seedNodes.begin(), seedNodes.end(), node_in_nghbr.second[i]) != seedNodes.end()){
                                seed_count += 1;
                                if((double)probT_non_edges[node_in_nghbr.first][i] * pro_v_is_seed[node_in_nghbr.second[i]] > max){
                                    max = (double)probT_non_edges[node_in_nghbr.first][i] * pro_v_is_seed[node_in_nghbr.second[i]];
                                    max_index = i;
                                }
                                if(seed_count == (int)seedNodes.size()){
                                    break;
                                }
                            }
                        }
                    }
                    if(max > 0.0){
                        new_E_bar.push_back(make_pair(node_in_nghbr.second[max_index], node_in_nghbr.first));
                    }
                }

                E_bar = new_E_bar;
                for(auto e:E_prime){
                    if(find(E_bar.begin(), E_bar.end(), e) != E_bar.end()){
                        E_bar.erase(find(E_bar.begin(), E_bar.end(), e));
                    }
                }
  
                map<string, double> comm_probs_without_v;
                comm_probs_without_v = sigma_C_p();

                map<int, double>v_minTau;
                map<int, int>v_minTau_bool;
                for(int i = 0; i < n; i++){
                    v_minTau_bool[i] = 0.0;
                }
                
                for(int j=0; j < (int)E_bar.size(); j++){
                    if(v_minTau_bool[E_bar[j].second] == 0.0){
                        vector<map<vector<int>, double>> p_temp;
                        for(int i = 0; i < (int)p.size(); i++){
                            for(auto S_p:p[i]){
                                map<vector<int>, double> temp;
                                temp[{E_bar[j].second}] = p[i][S_p.first];
                                p_temp.push_back(temp);
                            }
                        }

                        map<string, double> comm_probs_v;
                        comm_probs_v = sigma_C_p();

                        map<string, double> tau_C;
                        for(auto C:communities){
                            tau_C[C] = comm_probs_without_v[C] + comm_probs_v[C];
                        }

                        auto min_pro_tau = *min_element(tau_C.begin(), tau_C.end(),
                            [](decltype(tau_C)::value_type& l, decltype(tau_C)::value_type& r) -> bool { return l.second < r.second; });
                   
                        v_minTau[E_bar[j].second] = min_pro_tau.second;
                        v_minTau_bool[E_bar[j].second] = 1.0;
                    }
                
                    if(v_minTau[E_bar[j].second] > best_e_bar.second){
                        probT[E_bar[j].second].push_back(E_bar_random_weight[E_bar[j]]);
                        gT[E_bar[j].second].push_back(E_bar[j].first);
                        update_hyperGT_hyperG_rootR_p(E_bar[j], E_bar_random_weight[E_bar[j]]);
                        map<string, double> C_pros;
                        C_pros = sigma_C_p();

                        auto min_pro = *min_element(C_pros.begin(), C_pros.end(),
                            [](decltype(C_pros)::value_type& l, decltype(C_pros)::value_type& r) -> bool { return l.second < r.second; });
                        if(min_pro.second > best_e_bar.second){
                            best_e_bar = make_pair(E_bar[j], min_pro.second);
                        } 
                        gT = gT_temp;
                        probT = probT_temp;
                        set_of_hyperGT = set_of_hyperGT_temp;
                        set_of_hyperG = set_of_hyperG_temp;
                        set_of_rootR = set_of_rootR_temp;
                        p = p_temp;

                    }
                }
                if(find(E_prime.begin(), E_prime.end(), best_e_bar.first) == E_prime.end()){
                    E_prime.push_back(best_e_bar.first);
                    probT[best_e_bar.first.second].push_back(E_bar_random_weight[best_e_bar.first]);
                    gT[best_e_bar.first.second].push_back(best_e_bar.first.first);
                    update_hyperGT_hyperG_rootR_p(best_e_bar.first, E_bar_random_weight[best_e_bar.first]);
                    gT_temp = gT;
                    probT_temp = probT;
                    set_of_hyperGT_temp = set_of_hyperGT;
                    set_of_hyperG_temp = set_of_hyperG;
                    set_of_rootR_temp = set_of_rootR;
                    p_temp = p;
                }

                best_e_bar.second = -1.0;
                if((int)E_bar.size() == 0.0){
                    break;
                }  
            }
            
            return E_prime;
        }


        vector<pair<int, int>> to_minC_infl(){
            // Implement the algorithm to_minC_infl proposed in the paper.
            // return :
            //      a set of non-edges to be added to the graph        
            vector<pair<int, int>> E_prime;
            set_of_hyperGT_hyperG_rootR_p();
            while((int)E_prime.size() < B){
                map<int, vector<int>> gT_non_node_seeds;
                map<int, vector<double>> probT_non_node_seeds;
                map<int, double> pro_v_is_seed;
                set<int> seedNodes;
                for(auto item:p){
                    for(auto seed_p:item){
                        for(int v:seed_p.first){
                            seedNodes.insert(v);
                            pro_v_is_seed[v] += seed_p.second;
                        }
                    } 
                }

                for(int v=0; v<n; v++){
                    if(find(seedNodes.begin(), seedNodes.end(), v) == seedNodes.end()){
                        int count_seeds = 0.0;
                        for(int i = 0; i < (int)gT_non_edges[v].size(); i++){
                            if(find(seedNodes.begin(), seedNodes.end(), gT_non_edges[v][i]) != seedNodes.end()){
                                count_seeds += 1;
                                gT_non_node_seeds[v].push_back(gT_non_edges[v][i]);
                                probT_non_node_seeds[v].push_back(probT_non_edges[v][i]);
                                if(count_seeds == seedNodes.size()){
                                    break;
                                }
                            }
                        }
                    }
                }

                map<string, double> comm_pros;
                comm_pros = sigma_C_p();

                auto comm_min_pro = *min_element(comm_pros.begin(), comm_pros.end(),
                    [](decltype(comm_pros)::value_type& l, decltype(comm_pros)::value_type& r) -> bool {return l.second < r.second;});
                
                map<int, long long> coverage_node_incr;
                for(int h = 0; h < set_of_hyperGT.size(); h++){
                    map<int, vector<int>> hyperG_C;  
                    vector<vector<int>> hyperGT_C;

                    for(int v:communities_node[comm_min_pro.first]){
                        for(int indes_of_RRsets:set_of_rootR[h][v]){
                            hyperGT_C.push_back(set_of_hyperGT[h][indes_of_RRsets]);
                        }   
                    }

                    for(int i = 0; i < hyperGT_C.size(); i++){
                        for(int t:hyperGT_C[i]){
                            hyperG_C[t].push_back(i);
                        }
                    }
                    hyperG = hyperG_C;  
                    hyperGT = hyperGT_C;
                    rootR = set_of_rootR[h];
                    // update the coverage by removing the RR sets covered by seedSets, then compute increment for each node
                    priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
                    map<int, int>coverage;

                    for(int i=0; i<n; i++){
                        pair<int, int>tep(make_pair(i, (int)hyperG[i].size()));
                        heap.push(tep);
                        coverage[i] = (int)hyperG[i].size();
                    }
                    
                    long long influence = 0;
                    long long numEdge = hyperGT.size();

                    vector<bool> edgeMark(numEdge, false);
                    map<int, bool> nodeMark;
                    for(int v=0; v<n; v++){
                        nodeMark[v] = true;
                    }

                    for(auto item:p[h]){
                        seedSet = item.first;
                    }
                    
                    for(int maxInd:seedSet){
                        vector<int>e = hyperG[maxInd];
                        influence += coverage[maxInd];
                        nodeMark[maxInd] = false;
                        for(unsigned int j = 0; j < e.size(); ++j){
                            if(edgeMark[e[j]])continue;
                            vector<int>nList = hyperGT[e[j]];
                            for(unsigned int l = 0; l < nList.size(); ++l){
                                if(nodeMark[nList[l]]){
                                    coverage[nList[l]]--;
                                }
                            }
                            edgeMark[e[j]] = true;
                        }
                    }
                    for(int i = 0; i < n; i++){
                        coverage_node_incr[i] += coverage[i];
                    }
                }

                pair<int, int> max_edge;
                int max_index;
                double max_increment = 0.0;
                for(auto v_in_seeds:gT_non_node_seeds){
                    for(int i=0; i<v_in_seeds.second.size(); i++){
                        if((double)coverage_node_incr[v_in_seeds.first] * probT_non_node_seeds[v_in_seeds.first][i] * pro_v_is_seed[v_in_seeds.second[i]] > max_increment){
                            if(find(E_prime.begin(), E_prime.end(), make_pair(gT_non_node_seeds[v_in_seeds.first][i], v_in_seeds.first)) == E_prime.end()){
                                max_increment = (double)coverage_node_incr[v_in_seeds.first] * probT_non_node_seeds[v_in_seeds.first][i] * pro_v_is_seed[v_in_seeds.second[i]];
                                max_edge = make_pair(gT_non_node_seeds[v_in_seeds.first][i], v_in_seeds.first);
                                max_index = i;
                            }
                        }
                    }
                }
        
                E_prime.push_back(max_edge);
                probT[max_edge.second].push_back(probT_non_node_seeds[max_edge.second][max_index]);
                gT[max_edge.second].push_back(max_edge.first);
                update_hyperGT_hyperG_rootR_p(max_edge, probT_non_node_seeds[max_edge.second][max_index]);
            }
            
            return E_prime;
        }

        vector<pair<int, int>> to_minC_min(){
            // Implement the algorithm to_minC_min proposed in the paper.
            // return :
            //      a set of non-edges to be added to the graph 
            vector<pair<int, int>> E_prime;
            set_of_hyperGT_hyperG_rootR_p();
            while((int)E_prime.size() < B){
                map<int, double> pro_v_is_seed;
                set<int> seedNodes;
                for(auto item:p){
                    for(auto seed_p:item){
                        for(int v:seed_p.first){
                            seedNodes.insert(v);
                            pro_v_is_seed[v] += seed_p.second;
                        }
                    } 
                } 

                map<string, double> C_probs;
                C_probs = sigma_C_p();
                map<double, string> pro_CC;
                for(auto item:C_probs){
                    pro_CC[item.second] = item.first;
                }
                
                for(auto pro_C:pro_CC){
                    map<double, int> pro_nodes_in_C;
                    for(int v:communities_node[pro_C.second]){
                        pro_nodes_in_C[node_pros_p[v]] = v;
                    }
                    for(auto p_v:pro_nodes_in_C){
                        int max_index;
                        double max = 0.0;
                        int count_seeds = 0.0;
                        for(int i=0; i < (int)gT_non_edges[p_v.second].size(); i++){
                            if(find(seedNodes.begin(), seedNodes.end(), gT_non_edges[p_v.second][i]) != seedNodes.end()){
                                if((double)probT_non_edges[p_v.second][i] * pro_v_is_seed[gT_non_edges[p_v.second][i]] > max){
                                    max = (double)probT_non_edges[p_v.second][i] * pro_v_is_seed[gT_non_edges[p_v.second][i]];
                                    max_index = i;
                                }
                                if(count_seeds == seedNodes.size()){
                                    break;
                                }
                            }
                        }

                        if(max == 0.0){
                            continue;
                        }
                    
                        if(max > 0.0){
                            E_prime.push_back(make_pair(gT_non_edges[p_v.second][max_index], p_v.second));
                            pair<int, int> e_bar;    
                            e_bar = make_pair(gT_non_edges[p_v.second][max_index], p_v.second);
                            probT[e_bar.second].push_back(probT_non_edges[p_v.second][max_index]);
                            gT[e_bar.second].push_back(e_bar.first);
                            update_hyperGT_hyperG_rootR_p(e_bar, probT_non_edges[p_v.second][max_index]);
                            gT_non_edges[p_v.second].erase(gT_non_edges[p_v.second].begin() + max_index);
                            probT_non_edges[p_v.second].erase(probT_non_edges[p_v.second].begin() + max_index);
                        }
                        break;   
                    }
                    break; 
                }
            } 
                    
            return E_prime;
        }

        
        vector<pair<int, int>> max_weight(){
            // Implement the baseline max_weight.
            // return :
            //      a set of non-edges to be added to the graph 
            vector<pair<int, int>> E_prime;
            map<pair<int, int>, double> E_bar_random_weight_temp;
            E_bar_random_weight_temp = E_bar_random_weight;
            for(int i = 0; i < B; i++){
                auto e_bar_maxW = *max_element(E_bar_random_weight_temp.begin(), E_bar_random_weight_temp.end(),
                    [](decltype(E_bar_random_weight_temp)::value_type& l, decltype(E_bar_random_weight_temp)::value_type& r) -> bool { return l.second < r.second; });
                E_prime.push_back(e_bar_maxW.first);
                E_bar_random_weight_temp.erase(find(E_bar_random_weight_temp.begin(), E_bar_random_weight_temp.end(), e_bar_maxW));
            }
            
            return E_prime;
        }


        vector<pair<int, int>> random(){
            // Implement the baseline random.
            // return :
            //      a set of non-edges to be added to the graph
            vector<pair<int, int>> E_prime;
            vector<int> index_weight((int)E_bar.size(), 1.0);
            vector<int> rand_numbers;
            random_device rd;
            mt19937 generator{rd()}; 
            discrete_distribution<> distribution(index_weight.begin(), index_weight.end());
            while(true){
                int rnd = distribution(generator);
                if(find(rand_numbers.begin(), rand_numbers.end(), rnd) == rand_numbers.end()){
                    rand_numbers.push_back(rnd);
                    E_prime.push_back(E_bar[rnd]);
                    if((int)E_prime.size() == B){
                        break;
                    }
                }
            }
            
            return E_prime;
        }


        void updating_RRsets(pair<int, int> e_bar, double p){
            for(int index_of_RR:hyperG[e_bar.second]){ // update RR sets using e_bar
                double randDouble=double(sfmt_genrand_uint32(&sfmtSeed))/double(RAND_MAX)/2; 
                if(randDouble <= p){
                    BuildHypergraphNode(e_bar.first, index_of_RR, false);
                    // compute intrsection of hyperGT[index_of_RR] and RRset
                    sort(RRset.begin(), RRset.end());
                    sort(hyperGT[index_of_RR].begin(), hyperGT[index_of_RR].end());  
                    vector<int> diff;
                    set_difference(RRset.begin(), RRset.end(), hyperGT[index_of_RR].begin(), hyperGT[index_of_RR].end(), back_inserter(diff));
                    for(int v:diff){
                        hyperGT[index_of_RR].push_back(v);
                        hyperG[v].push_back(index_of_RR);
                    }
                }
            }
        }

        
        void set_of_hyperGT_hyperG_rootR_p(){            
            // compute the probability distribution over seed sets p(F, k).
            for(int i = 0; i < 10; i++){    
                EstimateOPT(epsilon, true);
                sort(seedSet.begin(), seedSet.end());
                add_to_set_of_hyperGT_G_rootR_p();
            }
        }


        void update_hyperGT_hyperG_rootR_p(pair<int, int> e, double w_e){
            // update sets of hyperGT, hyperG, rootR and p after adding non-edge e to the graph 
            vector<vector<vector<int>>> set_of_hyperGT_update;
            vector<map<int, vector<int>>> set_of_hyperG_update;  
            vector<map<int, vector<int>>> set_of_rootR_update;
            vector<map<vector<int>, double>> p_update;

            p_update = p;
            set_of_hyperGT_update = set_of_hyperGT;
            set_of_hyperG_update = set_of_hyperG;
            set_of_rootR_update = set_of_rootR;

            p.clear();
            set_of_hyperGT.clear();
            set_of_hyperG.clear();
            set_of_rootR.clear();

            for(int i=0; i<(int)set_of_hyperGT_update.size(); i++){
                hyperGT = set_of_hyperGT_update[i];
                hyperG = set_of_hyperG_update[i];  
                rootR = set_of_rootR_update[i];

                double countR = 0.0;
                for(auto S_p:p_update[i]){
                    countR = (double)p_update[i][S_p.first] * 10.0;                                       
                }

                vector<vector<int>> hyperGT_per_RR;
                map<int, vector<int>> hyperG_per_RR; 
                map<int, vector<int>> rootR_per_RR;
                hyperGT_per_RR = hyperGT;
                hyperG_per_RR = hyperG; 
                rootR_per_RR = rootR;
                for(int per_RR=0; per_RR < countR; per_RR++){
                    hyperGT = hyperGT_per_RR;
                    hyperG = hyperG_per_RR;
                    rootR = rootR_per_RR;
                    updating_RRsets(e, w_e);
                        
                    BuildSeedSet();
                    sort(seedSet.begin(), seedSet.end());
                    add_to_set_of_hyperGT_G_rootR_p();
                }
            }
        }


        void add_to_set_of_hyperGT_G_rootR_p(){
            if((int)p.size() == 0.0){
                set_of_hyperGT.push_back(hyperGT);
                set_of_hyperG.push_back(hyperG);
                set_of_rootR.push_back(rootR);
                map<vector<int>, double> temp;
                temp[seedSet] = 0.1;
                p.push_back(temp);
            }else{
                for(int j = 0; j < (int)p.size(); j++){
                    int jj = 0.0;
                    for(auto seed_p:p[j]){                                                
                        if(seed_p.first == seedSet){
                            p[j][seed_p.first] += 0.1;
                            jj = 1;
                            break;
                        }else{
                            if(j == (int)p.size()-1){
                                set_of_hyperGT.push_back(hyperGT);
                                set_of_hyperG.push_back(hyperG);
                                set_of_rootR.push_back(rootR);
                                map<vector<int>, double> temp;
                                temp[seedSet] = 0.1;
                                p.push_back(temp);
                                jj = 1;
                                break;
                            }
                        }
                    }
                    if (jj == 1.0){
                        break;
                    }           
                }
            }
        }

        map<int, double> sigma_v_p(){
            // compute sigma_v(p)
            // Returns
            // -------
            // node_pros_p : dict
            //       Dictionary with keys nodes and values probability of being reached of nodes using p. 
            map<vector<int>, map<int, double>> sets_node_pros;
            for(int i=0; i<set_of_hyperGT.size(); i++){
                hyperGT = set_of_hyperGT[i];
                hyperG = set_of_hyperG[i];  
                rootR = set_of_rootR[i];
                for(auto set_pro:p[i]){
                    seedSet = set_pro.first;
                }  
                sets_node_pros[seedSet] = sigma_v(seedSet);
            }
            node_pros_p.clear();
            for(int i=0; i<p.size(); i++){
                for(auto set_p:p[i]){
                    for(auto v_pro:sets_node_pros[set_p.first]){
                        node_pros_p[v_pro.first] += (double)set_p.second * v_pro.second;
                    }        
                }
            }

            return node_pros_p;
        }


        map<string, double> sigma_C_p(){
            // compute sigma_C(p)
            // Returns
            // -------
            // comm_probs : dict
            //       Dictionary with keys nodes and values the average probability of being reached of nodes in coomunities using p. 

            map<int, double> node_pro;
            node_pro = sigma_v_p();
            map<string, double> comm_probs;
            for(int v:nodes){
                for(auto C:node_communities[v]){
                    comm_probs[C] += (double)(node_pro[v] / communities_node[C].size());
                }
            }
            return comm_probs;
        }

        map<vector<int>, double> set_based_maximin(double epsilon, double eps=0.1, string set_or_node = "set"){
            // Call multi_weight for the set-based problem of Becker et al (2022).
            // Parameters
            // ----------
            // epsilon : float
            //      epsilon of (eps, delta)-approximation (for tim implementation).
            // eps : float
            //      eps of multi_weight.
            // set_or_node: str
            //      string showing that multi_weight should be executed for the set-based problem of Becker et al.

            return multi_weight(epsilon, eps, set_or_node = "set");
        }

        map<vector<int>, double> multi_weight(double epsilon, double eps, string set_or_node = "set"){
            // Implement multi_weight for the set-based algorithm of Becker et al (2022).
            // Parameters
            // ----------
            // epsilon : float
            //      epsilon of (eps, delta)-approximation (for tim implementation)     
            // eps : float
            //      eps of multi_weight.
            // set_or_node: str
            //      string showing that multi_weight should be executed for the set-based problem of Becker et al.
            // Returns
            // -------
            // p : dict
            //      Dictionary with keys sets and values probabilities.
            
            map<string, double>z;
            map<string, double>s;

            for(auto C:communities){
                z[C] = 1.0;
                s[C] = 0.0;
            }

            map<vector<int>, double>p;
            i_multi = 1;
            double primal = -INFINITY;
            double dual = INFINITY;

            while(true)
            {
                double sum_z = 0;
                for(auto item: z){
                    sum_z += item.second;
                }

                double sum_weights = 0.0;
                for(int index=0; index<nodes.size(); index++){
                    double temp = 0.0;
                    for(auto c: node_communities[nodes[index]]){
                        temp += (double) (z[c] / communities_node[c].size());
                    }
                    weights[index] = (double) temp;
                    sum_weights += (double)weights[index];
                }
 
                for(int v:nodes){
                    weights[v] *= (double) (n/sum_weights);
                }
                EstimateOPT(epsilon, true);
                vector<int>oracle_solution = seedSet;
                double oracle_value = spread;

                if((double)(oracle_value * (sum_weights / nodes.size()) / sum_z) > 1.0000001){ 
                    cout << "_____________________________"<<endl;
                    for(auto item: z){
                        cout<< item.second<<endl;
                    }
                    cout<< "oracle_value, sum(z.values())" <<endl;
                    cout << oracle_value<< "," << sum_z<< endl;
                    cout << (double)(oracle_value * (sum_weights / nodes.size()) / sum_z) <<endl;
                    cout << "_____________________________"<<endl;
                }
                ASSERT((double)(oracle_value * (sum_weights/nodes.size()) / sum_z) <= 1.0000001); 

                if(set_or_node == "set"){
                    vector<int>difference;
                    sort(oracle_solution.begin(), oracle_solution.end());
                    int ii = 0;
                    for(auto item:p){
                        if(oracle_solution == item.first)
                        {
                            p[oracle_solution] += 1;
                            ii += 1;
                            break;
                        }
                    }
                    if (ii == 0)
                    {
                        p[oracle_solution] = 1;
                    }
                }
                else if(set_or_node == "node"){
                    int jj = 0;
                    for(int v:oracle_solution){
                        jj = 0;
                        for(auto item:p){
                            if(find(item.first.begin(), item.first.end(), v) != item.first.end()){
                                p[item.first] += 1;
                                jj += 1;
                                break;
                            }
                        }
                        if(jj == 0){
                            p[{v}] = 1;
                        }
                    }
                }
                else{
                    cout <<"Error: Unknown option."<< set_or_node <<endl;
                }

                dual = min((double)dual, (double)oracle_value * (sum_weights / nodes.size()) / sum_z);

                map<string, double> community_prob;
                community_prob = sigma_C(oracle_solution);

                for (auto C:communities){
                    z[C] *= (double)(1 - (double)(eps * community_prob[C]));
                }

                for (auto C:communities){
                    s[C] = (double)(i_multi - 1) / i_multi * s[C] + (double)1 / i_multi * community_prob[C];
                }

                auto min_s = *min_element(s.begin(), s.end(),
                    [](decltype(s)::value_type& l, decltype(s)::value_type& r) -> bool { return l.second < r.second; });
                primal = min_s.second;
                if(primal >= (1 - eps) * dual){
                    break;
                }
                i_multi += 1;
            }

            for(auto item:p){
                p[item.first] = (double)item.second/i_multi;
            }
            return p; 
        }

};
