#define HEAD_INFO
#include "sfmt/SFMT.h"
using namespace std;
typedef double (*pf)(int,int);
#include <algorithm>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <set>

class Graph
{
    public:
        vector<vector<vector<int>>> set_of_hyperGT;
        vector<map<int, vector<int>>> set_of_hyperG;  
        vector<map<int, vector<int>>> set_of_rootR;
        vector<map<vector<int>, double>> p;
    
        map<pair<int, int>, double> all_edges;
        map<int, vector<int>> out_neighbors_v;
        int n, m, k, B;
        double epsilon;
        map<int, vector<int>> gT;
        map<int, vector<double>> probT;
        map<int, vector<int>> gT_non_edges;
        map<int, vector<double>> probT_non_edges;
        map<string, vector<int>> communities_node;
        map<int, vector<string>> node_communities;
        map<pair<int, int>, double> E_bar_random_weight;
        map<int, double> node_pros_p;
        map<int, int> nodeW;
        int sum_w = 0;
        vector<int> nodes;
        vector<pair<int, int>> E_bar;
        vector<double> weights;
        map<int, int> index_node;
        vector<string> communities;
        vector<bool> visit;
        vector<int> visit_mark;
        enum InfluModel {IC, LT};
        InfluModel influModel;
        void setInfuModel(InfluModel p){
            influModel=p;
        }

        string folder;
        string graph_file;
        string setting;
        void readNM(){
         // Read attribute.txt containing the number of nodes (n) and edges (m) of the graph
            ifstream cin((folder+"attribute.txt").c_str());
            ASSERT(!cin == false);
            string s;
            while(cin >> s){
                if(s.substr(0,2)=="n="){
                    n=atoi(s.substr(2).c_str());
                    continue;
                }
                if(s.substr(0,2)=="m="){
                    m=atoi(s.substr(2).c_str());
                    continue;
                }
                ASSERT(false);
            }
            visit_mark=vector<int>(n);     
            visit=vector<bool>(n);    
            cin.close();
        }
        void add_edge(int a, int b, double p){
            probT[b].push_back(p);
            gT[b].push_back(a);
        }
        
        void node_weight(int a, int b){
        // Assign weight 1 to each node
            if(nodeW.count(a)<1){
                nodeW[a] = 1;
                sum_w += nodeW[a];
            }
            if(nodeW.count(b)<1){
                nodeW[b] = 1;
                sum_w += nodeW[b];
            }  
        }

        void init_communities(){
            for(int j: nodes){      
                node_communities.insert(pair<int,vector<string>>(j, vector<string>()));
            }
        }

        void community(){
            // Read community.txt showing the community structure of the graph
            ifstream file(folder+"community.txt");
            string delimiter = " ";

            size_t pos = 0;
            string token;
            string line;
            init_communities();
            getline(file, line);
            while (getline(file, line)){
                string key; 
                int l = 0;
                while((pos = line.find(delimiter)) != string::npos){
                    l += 1;
                    token = line.substr(0, pos);
                    if(l == 1){
                        key = token;
                    }
                    if (l > 1){
                        communities_node[key].push_back(stoi(token));
                        node_communities[stoi(token)].push_back(key);
                    }
                    line.erase(0, pos + delimiter.length());
                }
                communities_node[key].push_back(stoi(line));
                node_communities[stoi(line)].push_back(key);
            }
            for(auto item: communities_node)
                communities.push_back(item.first);    
        }
        
        void readGraph(){
            FILE * fin= fopen((graph_file).c_str(), "r");
            long int readCnt=0;
            for(int i=0; i<m; i++){
                readCnt ++;
                int a, b;
                double p;
                int c=fscanf(fin, "%d%d%lf", &a, &b, &p);
                ASSERT(c==3);
                ASSERTT(c==3, a, b, p, c);
                
                add_edge(a, b, p);
                node_weight(a, b);
            }
            
            int h = 0.0; 
            for(auto item: nodeW){
                nodes.push_back(item.first);
                weights.push_back((double) item.second * n/sum_w); 
                index_node[h] = item.first;
                h += 1;
            }
            nodeW.clear();

            if(readCnt !=m)
                ExitMessage("m not equal to the number of edges in file "+graph_file);
            fclose(fin);
        }
        void read_non_edges(){
            // read all non-edges
            ifstream file(folder+"non_edges.txt");
            string line;
            while(getline(file, line)){
                stringstream linestream(line);
                string a;
                int b;
                double p;
              
                getline(linestream, a, ' ');
                linestream >> b >> p;

                if(setting == "grdy_al" or setting == "max_weight" or setting == "random"){
                    E_bar_random_weight[make_pair(stoi(a), b)] = p;
                    E_bar.push_back(make_pair(stoi(a), b));
                }
                if(setting != "random" and setting != "max_weight"){
                    probT_non_edges[b].push_back(p);
                    gT_non_edges[b].push_back(stoi(a));
                }
            }
        }

        sfmt_t sfmtSeed;

        Graph(string folder, string graph_file, string setting):folder(folder), graph_file(graph_file), setting(setting){
            readNM();
            community();
            readGraph();

            if(setting != "set"){
                read_non_edges();
            }        
        }

};
double sqr(double t)
{
    return t*t;
}

#include "infgraph.h"
#include "timgraph.h"


