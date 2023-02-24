//#define HEAD_TRACE
#define HEAD_INFO

#define HEAD_INFO
//#define HEAD_TRACE
#include "sfmt/SFMT.h"
#include "head.h"
#include "memoryusage.h"
#include "graph.h"
#include "adding_links.h"
#include <set>
#include <chrono>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>

void run(adding_links & x, string dataset, int k, int B, double epsilon, string model, string setting, string file_name){ // , string folder 
    //Infrastructure for the execution of the experiments
    // Parameters
    // ----------
    // dataset : str
    //    path to the dataset to be executed
    // k : int
    //    budget (seed set size)
    // B : int
    //    budget (number of non-edges to be added to the graph)
    // epsilon : float
    //    epsilon of (eps, delta)-approximation (for tim implementation).
    // model : str
    //    underlying diffusion model (IC or LT)
    // setting : str
    //    string showing which function should be executed
    // file_name : str 
    //    name of the file the results should be written


    cout << "dataset:" << dataset << " k:" << k << " -B:" << B << " epsilon:" << epsilon << " model:" << model << " setting:" << setting << " file_name:" << file_name << endl; //  << " folder: " << folder
    x.k=k;
    x.B=B;
    x.setting = setting;
    x.epsilon = epsilon;

    if(model=="IC")
        x.setInfuModel(InfGraph::IC);
    else if(model=="LT")
        x.setInfuModel(InfGraph::LT);
    else
        ASSERT(false);

    vector<int> S;
    map<vector<int>, double> p_set_based;

    auto start = chrono::high_resolution_clock::now();
    ios_base::sync_with_stdio(false);

    vector<pair<int, int>> E_prime;

    if(setting == "tim"){
        x.EstimateOPT(epsilon, true);
        S = x.seedSet;
    }
    else if (setting == "set"){
        p_set_based = x.set_based_maximin(epsilon, 0.1, "set");
    }
    else if (setting == "to_minC_infl"){
        E_prime = x.to_minC_infl();
    }
    else if (setting == "to_minC_min"){
        E_prime = x.to_minC_min();
    }
    else if (setting == "grdy_al"){
        E_prime = x.grdy_al();
    }
    else if (setting == "random"){
        E_prime = x.random();
    }
    else if (setting == "max_weight"){
        E_prime = x.max_weight();
    }
    else{
        cout <<"Error: Unknown option (setting)."<< endl;
        ASSERT(0); 
    }

    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;

    if(setting == "random" or setting == "max_weight"){ 
        for(auto e:E_prime){
            x.probT[e.second].push_back(x.E_bar_random_weight[e]);
            x.gT[e.second].push_back(e.first);
        }
    }

    if(setting == "tim"){
        ofstream out_file;
        out_file.open(file_name);
        out_file << time_taken << endl;
        int i = 0;
        for(int s:S){
            i += 1;
            if(i < (int)S.size()){
                out_file << s << " ";
            }
            else{
                out_file << s << endl;
            }
        }
    }
    else if (setting == "set"){
        ofstream out_file;
        out_file.open(file_name);
        out_file << time_taken << endl;
        for(auto S_pro:p_set_based){
            for (int s:S_pro.first){
                out_file << s << " ";
            }
            out_file << fixed << setprecision(20) << S_pro.second << endl;
        }
    }
    else{
        x.EstimateOPT(epsilon, true);
        map<string, double> C_pros;
        C_pros = x.sigma_C(x.seedSet);

        auto min_s = *min_element(C_pros.begin(), C_pros.end(),
            [](decltype(C_pros)::value_type& l, decltype(C_pros)::value_type& r) -> bool { return l.second < r.second; });
        double min_C = min_s.second; 

        auto max_s = *max_element(C_pros.begin(), C_pros.end(),
            [](decltype(C_pros)::value_type& l, decltype(C_pros)::value_type& r) -> bool { return l.second < r.second; });
        double max_C  = max_s.second;

        stringstream test(dataset);
        string segment;
        vector<string> seglist;
        while(getline(test, segment, '/')){
            seglist.push_back(segment);
        }
        string graphname = seglist[seglist.size()-1];

        ofstream out_file;
        out_file.open(file_name);

        out_file << graphname << endl;
        out_file << x.setting << endl;
        out_file << x.n << endl;
        out_file << x.m << endl;
        out_file << x.k << endl;
        out_file << x.B << endl;        
        out_file << time_taken << endl;
        out_file << x.spread << endl;
        out_file << min_C << endl;
        out_file << max_C << endl;

        out_file.close();

        // cout << "function " << ": " << setting << endl;
        // cout << "solution " << ": ";
        // int j = 0.0;
        // for(auto e:E_prime){
        //     j += 1;
        //     if (j != (int)E_prime.size()){
        //         cout << "(" << e.first << ", " << e.second << ")" << ", ";
        //     }else{
        //         cout << "(" << e.first << ", " << e.second << ")" << endl;
        //     }
        // }
        // cout << "spread " << ": " << x.spread << endl;  
        // cout << "min_comm_prob " << ": " << min_C << endl;
        // cout << "max_comm_prob " << ": " << max_C << endl; 
        // cout << "ex-time " << ": " << time_taken << endl;
        
    }

    Counter::show();
}



void parseArg(int argn, char ** argv)
{
    string file_name=""; 
    string setting=""; 
    string dataset="";

    double epsilon=0;
    string model="";
    int k=0.0;
    int B=0.0;

    for(int i=0; i<argn; i++)
    {
        if(argv[i]==string("-file_name")) file_name=string(argv[i+1]);
        if(argv[i]==string("-setting")) setting=string(argv[i+1]);
        if(argv[i]==string("-dataset")) dataset=string(argv[i+1])+"/";
        if(argv[i]==string("-epsilon")) epsilon=atof(argv[i+1]);
        if(argv[i]==string("-k")) k=atoi(argv[i+1]);
        if(argv[i]==string("-B")) B=atoi(argv[i+1]);
        if(argv[i]==string("-model")) {
            if(argv[i+1]==string("LT"))
            {
                model=argv[i+1];
            }
            else if(argv[i+1]==string("IC"))
            {
                model=argv[i+1];
            }
            else
                ExitMessage("model should be IC or LT");
        }
    }
    if (dataset=="")
        ExitMessage("argument dataset missing");
    if (k==0)
        ExitMessage("argument k missing");
    if (epsilon==0)
        ExitMessage("argument epsilon missing");
    if (model=="")
        ExitMessage("argument model missing");
    if (setting=="")                             
        ExitMessage("argument setting missing");
    if (file_name=="")                             
        ExitMessage("argument file_name missing");  

    string graph_file;
    if(model=="IC")
        graph_file=dataset + "graph_ic.txt"; 
    else if(model=="LT")
        graph_file=dataset + "graph_lt.inf";

    adding_links x(dataset, graph_file, setting);  
   
    run(x, dataset, k, B, epsilon, model, setting, file_name); 
}

int main(int argn, char ** argv)
{
    OutputInfo info(argn, argv);
    parseArg( argn, argv );
}
