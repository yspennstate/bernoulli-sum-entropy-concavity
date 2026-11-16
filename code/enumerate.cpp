// Exhaustive enumeration of entropy profiles. Each fibre is encoded in base 64;
// the digit at position k counts subsets of cardinality k. No digit can carry
// because there are at most 32 subsets. Equal sums are grouped after sorting.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

using Profile = std::vector<std::uint64_t>;
std::map<Profile, std::vector<int>> found;
std::uint64_t total = 0;
int m, bound;
std::vector<int> weights;

void visit() {
    std::array<std::pair<int,int>,32> terms;
    for (int mask=0; mask<(1<<m); ++mask) {
        int s=0, k=0;
        for (int i=0;i<m;++i) if (mask & (1<<i)) {s+=weights[i]; ++k;}
        terms[mask]={s,k};
    }
    std::sort(terms.begin(), terms.begin()+(1<<m));
    Profile key;
    for(int i=0; i<(1<<m);) {
        int j=i; std::uint64_t row=0;
        while(j<(1<<m) && terms[j].first==terms[i].first) {
            row += std::uint64_t(1) << (6*terms[j].second); ++j;
        }
        key.push_back(row); i=j;
    }
    std::sort(key.begin(),key.end());
    found.emplace(key,weights); ++total;
}
void walk(int i,int lower) {
    if(i==m) {visit();return;}
    for(int x=lower;x<=bound;++x) {weights[i]=x;walk(i+1,x);}
}
int main(int argc,char**argv) {
    if(argc!=4) return 2;
    m=std::stoi(argv[1]); bound=std::stoi(argv[2]);
    if(m<1 || m>5 || bound<1 || bound>55) return 2;
    weights.resize(m);walk(0,1);
    std::ofstream out(argv[3]);
    out<<m<<' '<<bound<<' '<<total<<' '<<found.size()<<'\n';
    for(auto const& item: found) {
        for(int x:item.second) out<<x<<' ';
        out<<item.first.size();
        for(auto x:item.first) out<<' '<<x;
        out<<'\n';
    }
    std::cout<<"m="<<m<<" bound="<<bound<<" vectors="<<total<<" profiles="<<found.size()<<'\n';
}
