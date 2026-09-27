// Bounded, reproducible adapter to the user-supplied OILL routines.
#include "reduce.h"
#include <random>
#include <fstream>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <stdexcept>
#include <cstdio>
extern mt19937 rand_generator;
static vector<ROW> original;
vector<ROW> get_matrix() { return original; }

static vector<int> validate(const vector<xpair>& seq) {
    auto a=original;
    for(auto p:seq) {
        if(p.src<0 || p.dst<0 || p.src>=SIZE || p.dst>=SIZE || p.src==p.dst)
            throw runtime_error("Invalid XOR operation");
        a[p.dst]^=a[p.src];
    }
    vector<int> perm(SIZE), seen(SIZE,0);
    for(int i=0;i<SIZE;i++) {
        if(a[i].count()!=1) throw runtime_error("Reduction is not a permutation");
        for(int j=0;j<SIZE;j++) if(a[i][SIZE-1-j]) {perm[i]=j;break;}
        if(seen[perm[i]]++) throw runtime_error("Duplicate permutation column");
    }
    return perm;
}

int main(int argc,char**argv) {
  try {
    if(argc!=7) throw runtime_error("Usage: evaluator input.txt output.json seed iterations windows window_max");
    ifstream f(argv[1]); if(!f) throw runtime_error("Cannot read matrix");
    int size; f>>size; if(size!=SIZE) throw runtime_error("Wrong matrix dimension");
    string s;
    for(int i=0;i<SIZE;i++) {
        if(!(f>>s) || s.size()!=SIZE || s.find_first_not_of("01")!=string::npos)
            throw runtime_error("Bad matrix row");
        original.emplace_back(s);
    }
    if(f>>s) throw runtime_error("Trailing input");
    auto rank_rows=original; int rank=0;
    for(int col=0;col<SIZE;col++) {
        int r=rank; while(r<SIZE && !rank_rows[r][col]) ++r;
        if(r==SIZE) continue;
        swap(rank_rows[r],rank_rows[rank]);
        for(int k=rank+1;k<SIZE;k++) if(rank_rows[k][col]) rank_rows[k]^=rank_rows[rank];
        ++rank;
    }
    if(rank!=SIZE) throw runtime_error("Singular input matrix");
    unsigned seed=stoul(argv[3]); int iters=stoi(argv[4]), windows=stoi(argv[5]), window_max=stoi(argv[6]);
    if(iters<0 || windows<0 || window_max<3) throw runtime_error("Invalid optimization limits");
    rand_generator.seed(seed);
    auto start=chrono::steady_clock::now();
    vector<xpair> best; bool have=false;
    auto save=[&](const vector<xpair>&seq,const string&stage,int iteration) {
        auto perm=validate(seq);
        if(have && seq.size()>=best.size()) return;
        best=seq; have=true;
        string out=argv[2], tmp=out+".tmp";
        ofstream o(tmp); if(!o) throw runtime_error("Cannot write result");
        o<<"{\"size\":"<<SIZE<<",\"xor_count\":"<<best.size()<<",\"seed\":"<<seed
         <<",\"stage\":\""<<stage<<"\",\"iteration\":"<<iteration
         <<",\"elapsed_to_best_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()
         <<",\"input_permutation\":[";
        for(int i=0;i<SIZE;i++) {if(i)o<<",";o<<perm[i];}
        o<<"],\"gates\":[";
        bool first=true;
        for(auto it=best.rbegin();it!=best.rend();++it) {if(!first)o<<",";first=false;o<<"["<<it->dst<<","<<it->src<<"]";}
        o<<"]}\n";o.close(); if(!o)throw runtime_error("Result write failed");
        if(rename(tmp.c_str(),out.c_str())!=0) throw runtime_error("Result rename failed");
    };
    // A valid checkpoint exists before any expensive heuristic work.
    auto a=original; save(strgy1(a),"gaussian",-1);
    for(int i=0;i<iters;i++) {
        a=original; auto seq=strgy3(a); save(seq,"greedy",i);
        reduce_step(seq); save(seq,"local_reduction",i);
        for(int w=0;w<windows && seq.size()>=3;w++) {
            int cap=min<int>(window_max,seq.size());
            int gap=3+rand_generator()%(cap-2);
            int at=rand_generator()%(seq.size()-gap+1);
            auto candidate=seq;
            get_equivalent_seq(candidate,gap,at);
            reduce_step(candidate);
            save(candidate,"window_reduction",i);
            if(candidate.size()<seq.size()) seq=candidate;
        }
    }
    cout<<"XOR="<<best.size()<<"\n";
    return 0;
  } catch(const exception&e) {cerr<<e.what()<<"\n";return 2;}
}
