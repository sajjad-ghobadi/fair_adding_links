class VV
{
    public:
        vector<int> head;
        vector<int> next;
        vector<int> data;
        vector<int> vsize;
        void clear()
        {
            head.clear();
            next.clear();
            data.clear();
            vsize.clear();
        }
        // trick for not change code
        void push_back(vector<int> x){
            ASSERT(x.size()==0);
            addVector();
        }
        void addVector()
        {
            head.push_back(-1);
            vsize.push_back(0);
        }
        int size(int t){
            return vsize[t];
        }
        //vv[a].push_back(b)
        void addElement( int a, int b)
        {
            //a.push_back(b);
            vsize[a]++;
            data.push_back(b);
            next.push_back(head[a]);
            head[a]=next.size()-1;

        }
};

#include <queue>	
#include <utility> 
#include <random>
#include <ctime>

struct CompareBySecond {
	bool operator()(pair<int, int> a, pair<int, int> b)
	{
		return a.second < b.second;
	}
};

class InfGraph:public Graph
{
    public:
        map<int, vector<int>> hyperG_tempp;  
        vector<vector<int>> hyperGT_tempp; 

        int64 R = 0.0;
        double spread = 0.0;
        vector<double> RRset;
        map<int, double> node_probs;
        vector<vector<int>> hyperGT;  
        map<int, vector<int>> hyperG;  
        map<int, vector<int>> rootR; 
        
        InfGraph(string folder, string graph_file, string setting):Graph(folder, graph_file, setting){
            hyperG.clear();

            for(int i=0; i<n; i++){
                hyperG.insert(pair<int,vector<int>>(i, vector<int>())); 
            }
            for(int i=0; i<12; i++)
                sfmt_init_gen_rand(&sfmtSeed, i+1234); 
        }

        enum ProbModel {TR, WC, TR001};
        ProbModel probModel;


        void hyper_GT_G(int64 R, vector<int>nodes_or_comm_nodes){
            hyperGT.clear();
            hyperG.clear();
            rootR.clear();   
            for(int i:nodes){
                hyperG.insert(pair<int,vector<int>>(i, vector<int>()));    
                rootR.insert(pair<int,vector<int>>(i, vector<int>()));  
            }
                
            int R_per_node = 0.0;
            R_per_node = (int) (ceil((double) R/nodes_or_comm_nodes.size()));  
            int up_R = R_per_node * nodes_or_comm_nodes.size();
            hyperId = up_R;

            while((int)hyperGT.size() < up_R)
                hyperGT.push_back(vector<int>());

            int index = 0.0;
            for (int v:nodes_or_comm_nodes){
                for (int t = 0; t < R_per_node; t++){
                    rootR[v].push_back(index);        
                    BuildHypergraphNode(v, index, true);
                    index += 1;
                }
            }

            int totAddedElement=0; 
            for(int i=0; i<up_R; i++){     
                for(int t:hyperGT[i]){
                    hyperG[t].push_back(i);
                    totAddedElement++;
                }
            }
        }

        void BuildHypergraphR(int64 R){
            // create RR sets 
            if (setting == "set"){
                hyperG.clear();
                rootR.clear();   
                for(int i=0; i<n; i++){
                    hyperG.insert(pair<int,vector<int>>(i, vector<int>()));    
                    rootR.insert(pair<int,vector<int>>(i, vector<int>()));  
                }
                random_device rd;
                mt19937 generator{rd()}; 
                discrete_distribution<> distribution(weights.begin(),weights.end());

                hyperGT.clear();
                    while((int)hyperGT.size() < R)
                        hyperGT.push_back(vector<int>());

                hyperId = R; 
                    
                for(int i=0; i<R; i++){
                    int number = distribution(generator);
                    rootR[index_node[number]].push_back(i); 
                    BuildHypergraphNode(index_node[number], i, true); 
                }
                int totAddedElement=0; 
                for(int i=0; i<R; i++){     
                    for(int t:hyperGT[i]){
                        hyperG[t].push_back(i);
                        totAddedElement++;
                    }
                }
                ASSERT(hyperId == R);
            }
            else{
                hyper_GT_G(R, nodes);
            }            
        }

        int BuildHypergraphNode(int uStart, int hyperiiid, bool addHyperEdge){
            RRset.clear();
            int n_visit_edge=1;
            if(addHyperEdge)
            {
                ASSERT((int)hyperGT.size() > hyperiiid);
                hyperGT[hyperiiid].push_back(uStart);
            }else{
                RRset.push_back(uStart);
            }

            int n_visit_mark=0;
            q.clear();
            q.push_back(uStart);
            ASSERT(n_visit_mark < n);
            visit_mark[n_visit_mark++]=uStart;
            visit[uStart]=true;
            while(!q.empty()) {
                int expand=q.front();
                q.pop_front();
                if(influModel==IC){
                    int i=expand;
                    for(int j=0; j<(int)gT[i].size(); j++){
                        int v=gT[i][j];
                        n_visit_edge++;
                        double randDouble=double(sfmt_genrand_uint32(&sfmtSeed))/double(RAND_MAX)/2; // 
                        if(randDouble > probT[i][j])
                            continue;
                        if(visit[v])
                            continue;
                        if(!visit[v])
                        {
                            ASSERT(n_visit_mark < n);
                            visit_mark[n_visit_mark++]=v;
                            visit[v]=true;
                        }
                        q.push_back(v);
                        if(addHyperEdge)
                        {
                            ASSERT((int)hyperGT.size() > hyperiiid);
                            hyperGT[hyperiiid].push_back(v);
                        }else{
                            RRset.push_back(v);
                        }
                    }
                }
                else if(influModel==LT){
                    if(gT[expand].size()==0)
                        continue;
                    ASSERT(gT[expand].size()>0);
                    n_visit_edge+=gT[expand].size();
                    double randDouble=double(sfmt_genrand_uint32(&sfmtSeed))/double(RAND_MAX)/2;
                    for(int i=0; i<(int)gT[expand].size(); i++){
                        ASSERT( i< (int)probT[expand].size());
                        randDouble -= probT[expand][i];
                        if(randDouble>0)
                            continue;
                        int v=gT[expand][i];

                        if(visit[v])
                            break;
                        if(!visit[v])
                        {
                            visit_mark[n_visit_mark++]=v;
                            visit[v]=true;
                        }
                        q.push_back(v);
                        if(addHyperEdge)
                        {
                            ASSERT((int)hyperGT.size() > hyperiiid);
                            hyperGT[hyperiiid].push_back(v);
                        }
                        break;
                    }
                }
                else
                    ASSERT(false);
            }
            for(int i=0; i<n_visit_mark; i++)
                visit[visit_mark[i]]=false;
            return n_visit_edge;
        }

        //return the number of edges visited
        int64 hyperId = 0;
        deque<int> q;
        sfmt_t sfmtSeed;
		vector<int> seedSet;

		double BuildSeedSet(){
            // Compute the seed set and the expected spread from the seed set
			priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
            map<int, int>coverage;
        
            for (int i=0; i<n; i++){
				pair<int, int>tep(make_pair(i, (int)hyperG[i].size()));
				heap.push(tep);
				coverage[i] = (int)hyperG[i].size();
			}
			int maxInd;

			long long influence = 0;
			long long numEdge = hyperGT.size();

			// check if an edge is removed
			vector<bool> edgeMark(numEdge, false);
			// check if an node is remained in the heap
            map<int, bool> nodeMark;
            for (int v=0; v<n; v++){
                nodeMark[v] = true;
            }

			seedSet.clear();            
			while ((int)seedSet.size()<k){
				pair<int, int>ele = heap.top();
				heap.pop();
				if (ele.second > coverage[ele.first])
				{
					ele.second = coverage[ele.first];
					heap.push(ele);
					continue;
				}

				maxInd = ele.first;
				vector<int>e = hyperG[maxInd];  //the edge influence
				influence += coverage[maxInd];
				seedSet.push_back(maxInd);
				nodeMark[maxInd] = false;
            
				for (unsigned int j = 0; j < e.size(); ++j){
					if (edgeMark[e[j]])continue;
					vector<int>nList = hyperGT[e[j]];
					for (unsigned int l = 0; l < nList.size(); ++l){
						if (nodeMark[nList[l]]){
                            coverage[nList[l]]--;
                        }
					}
					edgeMark[e[j]] = true;
				}
			}

            spread = 1.0*influence / hyperGT.size()*nodes.size();
			return  spread;
		}

        map<int, double>sigma_v(vector<int> S){
            // Compute sigma_v(S).
            // Parameters
            // ----------
            // S : list of seed nodes
            // Returns
            // -------
            // pro_v : dict
            //       Dictionary with keys node and values probabilities of being reached of nodes from S. 

            sort(S.begin(), S.end());
            map<int, double> pro_v;
            for(int index=0; index<(int)nodes.size(); index++){
                pro_v[nodes[index]] = 0.0;
                for(int r:rootR[nodes[index]]){
                    sort(hyperGT[r].begin(), hyperGT[r].end());
                    vector<int> intersect;   
                    set_intersection(hyperGT[r].begin(), hyperGT[r].end(), S.begin(),
                                                            S.end(), back_inserter(intersect));
                    if((int)intersect.size()>0){
                       pro_v[nodes[index]] += 1.0;
                    }
                }
               pro_v[nodes[index]] = pro_v[nodes[index]] * (double) n/(hyperGT.size() * weights[index]);
            }
           return pro_v;
        }
       
        map<string, double>sigma_C(vector<int> S){
            // Compute sigma_C(S).
            // Parameters
            // ----------
            // S : list of seed nodes
            // Returns
            // -------
            // comm_probs : dict
            //       Dictionary with keys communities and values average probability of being reached of nodes in the community from S. 

            node_probs.clear();
            node_probs = sigma_v(S);
            map<string, double>comm_probs;
            for(auto C:communities){
                comm_probs[C] = 0.0;
            }
            for(int v:nodes){
                for(auto C:node_communities[v]){
                        comm_probs[C] += (double)(node_probs[v] / communities_node[C].size());
                }
            }
            return comm_probs;
        }

        double sigma_S(vector<int> S){
            // Compute sigma(S).
            // Parameters
            // ----------
            // p : list of seed nodes
            // Returns
            // -------
            // infl : float
            //       number showing the expected number of nodes reached from S. 

           map<int, double>node_pro = sigma_v(S);           
           double infl = 0.0;
           for (int v:nodes){
               infl += node_pro[v];
           }
           return infl;
        }

};