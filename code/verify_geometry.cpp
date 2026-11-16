// Check closure completeness using nullspace projections and signed ray tests.
// Orbit keys are sets of contained ternary normals, not exterior coordinates.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>
using Row=std::array<int,6>;
using Key=std::array<uint64_t,5>;
std::vector<Row> planes;
std::vector<std::array<int,301>> actions;
std::vector<int> masks[7];int pos[7][64];
std::map<Key,Key> cache;
Key canonical(const Key& k){
 auto old=cache.find(k);if(old!=cache.end())return old->second;
 std::vector<int> on;
 for(int i=0;i<301;++i)if((k[i/64]>>(i%64))&1)on.push_back(i);
 Key best=k;
 for(auto& action:actions){Key q{};for(int i:on){int j=action[i];q[j/64]|=uint64_t(1)<<(j%64);}if(q<best)best=q;}
 cache.emplace(k,best);return best;
}
int sign_of(const Row& a){for(int x:a)if(x)return x>0?1:-1;return 0;}
struct Flat{int r;std::vector<int> p;std::vector<Row> null,rays,sig;Key closed{},canon{};};
Flat decode(const std::string& line){
 Flat f;std::istringstream in(line);in>>f.r;int x;while(in>>x)f.p.push_back(x);
 assert(f.p.size()==masks[f.r].size());int pivot=0;while(!f.p[pivot])++pivot;int I=masks[f.r][pivot],D=f.p[pivot];assert(D>0);
 std::vector<int> inds;for(int i=0;i<6;++i)if((I>>i)&1)inds.push_back(i);
 for(int j=0;j<6;++j)if(!((I>>j)&1)){
  Row v{};v[j]=D;
  for(int h=0;h<f.r;++h){int J=(I^(1<<inds[h]))|(1<<j);int t=__builtin_popcount((unsigned)(J&((1<<j)-1)));v[inds[h]]=-((h+t)%2?-1:1)*f.p[pos[f.r][J]];}
  f.null.push_back(v);
 }
 // All cofactor vectors of this already independently validated flat.
 int covered=0;
 for(int J:masks[f.r+1]){
  Row v{};bool plus=false,minus=false;int t=0;
  for(int j=0;j<6;++j)if((J>>j)&1){v[j]=(t%2?-1:1)*f.p[pos[f.r][J^(1<<j)]];plus|=v[j]>0;minus|=v[j]<0;++t;}
  if(plus==minus)continue;if(minus)for(int& z:v)z=-z;
  for(int j=0;j<6;++j)if(v[j])covered|=1<<j;
  f.rays.push_back(v);
 }
 assert(covered==63);
 for(int i=0;i<301;++i){Row s{};bool zero=true;for(int j=0;j<(int)f.null.size();++j){for(int h=0;h<6;++h)s[j]+=planes[i][h]*f.null[j][h];zero&=s[j]==0;}f.sig.push_back(s);if(zero)f.closed[i/64]|=uint64_t(1)<<(i%64);}
 f.canon=canonical(f.closed);return f;
}
int main(){
 auto start=std::chrono::steady_clock::now();
 for(int r=0;r<7;++r)for(int s=0;s<64;++s)if(__builtin_popcount((unsigned)s)==r){pos[r][s]=masks[r].size();masks[r].push_back(s);}
 std::map<Row,int> ids;
 for(int code=0;code<729;++code){int c=code;Row a{};bool p=false,n=false;for(int& x:a){x=c%3-1;c/=3;p|=x>0;n|=x<0;}if(p&&n&&sign_of(a)==1){ids[a]=planes.size();planes.push_back(a);}}
 assert(planes.size()==301);Row perm{0,1,2,3,4,5};
 do{std::array<int,301> action{};for(int i=0;i<301;++i){Row b{};for(int j=0;j<6;++j)b[perm[j]]=planes[i][j];if(sign_of(b)<0)for(int& x:b)x=-x;action[i]=ids.at(b);}actions.push_back(action);}while(std::next_permutation(perm.begin(),perm.end()));
 std::ifstream in("flats6_flats.txt");assert(in);std::vector<Flat> fs;std::set<Key> listed[6];std::string line;
 while(std::getline(in,line)){auto f=decode(line);assert(listed[f.r].insert(f.canon).second);fs.push_back(std::move(f));}
 uint64_t extensions[5]{},infeasible[5]{},dependent[5]{};
 for(const auto& f:fs){if(f.r==5)continue;int dim=f.null.size();
  for(int i=0;i<301;++i){const auto& a=f.sig[i];int pivot=0;while(pivot<dim&&a[pivot]==0)++pivot;if(pivot==dim){++dependent[f.r];continue;}
   // A new hyperplane intersects the strictly positive part iff the linear
   // form has both signs on the old cone's rays. Zero on all is handled above.
   bool plus=false,minus=false;
   for(const auto& ray:f.rays){int d=0;for(int h=0;h<6;++h)d+=planes[i][h]*ray[h];plus|=d>0;minus|=d<0;}
   if(!(plus&&minus)){++infeasible[f.r];continue;}
   Key k{};
   for(int j=0;j<301;++j){bool proportional=true;for(int h=0;h<dim;++h)if(int64_t(f.sig[j][h])*a[pivot]!=int64_t(a[h])*f.sig[j][pivot]){proportional=false;break;}if(proportional)k[j/64]|=uint64_t(1)<<(j%64);}
   assert(listed[f.r+1].count(canonical(k))==1);++extensions[f.r];
  }
  if((extensions[f.r]+infeasible[f.r]+dependent[f.r])%10000==0)std::cout<<"rank "<<f.r<<" checked "<<extensions[f.r]+infeasible[f.r]+dependent[f.r]<<std::endl;
 }
 std::ofstream out("geometry_results.json");out<<"{\"pass\":true,\"flats\":"<<fs.size()<<",\"orbits_by_rank\":[";
 for(int r=0;r<6;++r){if(r)out<<',';out<<listed[r].size();}out<<"],\"feasible_extensions\":[";
 for(int r=0;r<5;++r){if(r)out<<',';out<<extensions[r];}out<<"],\"infeasible_extensions\":[";
 for(int r=0;r<5;++r){if(r)out<<',';out<<infeasible[r];}out<<"],\"already_contained\":[";
 for(int r=0;r<5;++r){if(r)out<<',';out<<dependent[r];}out<<"],\"canonical_cache_size\":"<<cache.size()<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
 std::cout<<"PASS: no duplicate orbits; every feasible successor has a listed orbit"<<std::endl;
}
