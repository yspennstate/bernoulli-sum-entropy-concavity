// Exact collision flats modulo coordinate permutations. Single threaded.
#include <algorithm>
#include <cassert>
#include <cassert>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>
using V=std::vector<int>;
int n;
std::vector<V> masks, pos;
struct Action { V index,sign; };
std::vector<std::vector<Action>> actions;
void normalize(V& v){int g=0;for(int x:v)g=std::gcd(g,std::abs(x));if(!g)return;for(int x:v)if(x){if(x<0)g=-g;break;}for(int&x:v)x/=g;}
V extend(const V&p,int r,const V&a){
 V q(masks[r+1].size());
 for(int j=0;j<(int)q.size();++j){int s=masks[r+1][j],i=0;for(int k=0;k<n;++k)if(s>>k&1){q[j]+=(i%2?-1:1)*a[k]*p[pos[r][s^(1<<k)]];++i;}}
 normalize(q);return q;
}
V canonical(const V&p,int r){
 V best=p,q(p.size());
 for(auto const&a:actions[r]){
  int orientation=1;for(int j=0;j<(int)p.size();++j)if(p[a.index[j]]){orientation=a.sign[j]*(p[a.index[j]]>0?1:-1);break;}
  bool less=false,greater=false;
  for(int j=0;j<(int)p.size();++j){int x=orientation*a.sign[j]*p[a.index[j]];q[j]=x;if(!less&&!greater){less=x<best[j];greater=x>best[j];if(greater)break;}}
  if(less)best=q;
 }return best;
}
std::vector<V> rays(const V&p,int r){
 std::set<V> out;
 for(int s:masks[r+1]){
  V v(n);int i=0;bool plus=false,minus=false;
  for(int k=0;k<n;++k)if(s>>k&1){v[k]=(i%2?-1:1)*p[pos[r][s^(1<<k)]];plus|=v[k]>0;minus|=v[k]<0;++i;}
  if(plus==minus)continue;
  normalize(v);out.insert(v);
 }
 int covered=0;for(auto const&v:out)for(int i=0;i<n;++i)if(v[i])covered|=1<<i;
 if(covered!=(1<<n)-1)return {};
 return {out.begin(),out.end()};
}
using Profile=std::vector<std::uint64_t>;
Profile profile(const std::vector<V>&vs){
 std::map<V,std::uint64_t> atoms;
 for(int s=0;s<(1<<n);++s){V sig(vs.size());for(int j=0;j<(int)vs.size();++j)for(int i=0;i<n;++i)if(s>>i&1)sig[j]+=vs[j][i];atoms[sig]+=std::uint64_t(1)<<(6*__builtin_popcount((unsigned)s));}
 Profile p;for(auto const&kv:atoms)p.push_back(kv.second);std::sort(p.begin(),p.end());return p;
}
int main(int argc,char**argv){
 if(argc!=3)return 2;n=std::stoi(argv[1]);if(n<2||n>6)return 2;
 auto start=std::chrono::steady_clock::now();
 masks.resize(n+1);pos.resize(n+1,V(1<<n,-1));actions.resize(n);
 for(int s=0;s<(1<<n);++s){int r=__builtin_popcount((unsigned)s);pos[r][s]=masks[r].size();masks[r].push_back(s);}
 V perm(n);std::iota(perm.begin(),perm.end(),0);
 do{for(int r=0;r<n;++r){Action a;for(int s:masks[r]){V inds;int t=0,inv=0;for(int i=0;i<n;++i)if(s>>i&1){inds.push_back(perm[i]);t|=1<<perm[i];}for(int i=0;i<r;++i)for(int j=i+1;j<r;++j)inv+=inds[i]>inds[j];a.index.push_back(pos[r][t]);a.sign.push_back(inv%2?-1:1);}actions[r].push_back(a);}}while(std::next_permutation(perm.begin(),perm.end()));
 std::vector<V> planes;
 int limit=1;for(int i=0;i<n;++i)limit*=3;
 for(int code=0;code<limit;++code){int t=code;V a(n);bool plus=false,minus=false;for(int&i:a){i=t%3-1;t/=3;plus|=i>0;minus|=i<0;}if(!plus||!minus)continue;V b=a;normalize(b);if(a==b)planes.push_back(a);}
 std::map<Profile,std::vector<V>> profiles;
 std::set<V> current;current.insert(V{1});
 std::ofstream flats(std::string(argv[2])+"_flats.txt");
 std::vector<std::size_t> counts;
 for(int r=0;r<n;++r){
  counts.push_back(current.size());std::set<V> next;
  std::size_t visited=0;
  for(auto const&p:current){
   for(int x:p)assert(x>=-30&&x<=30); // Every subsequent extension coordinate has absolute value at most 180.
   // This invariant also bounds each unnormalized extension by 6*30.
   for(int x:p)assert(std::abs(x)<=30);
   auto vs=rays(p,r);if(vs.empty()){std::cerr<<"BAD FEASIBILITY\n";return 3;}
   auto key=profile(vs);profiles.emplace(key,vs);
   flats<<r;for(int x:p)flats<<' '<<x;flats<<'\n';
   if(r<n-1)for(auto const&a:planes){auto q=extend(p,r,a);if(std::all_of(q.begin(),q.end(),[](int x){return !x;}))continue;if(rays(q,r+1).empty())continue;next.insert(canonical(q,r+1));}
   ++visited;if(visited%1000==0)std::cout<<"progress rank="<<r<<" done="<<visited<<" / "<<current.size()<<" next="<<next.size()<<std::endl;
  }
  std::cout<<"rank="<<r<<" orbits="<<current.size()<<" profiles="<<profiles.size()<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
  current=std::move(next);
 }
 std::ofstream out(std::string(argv[2])+"_profiles.json");
 out<<"{\"m\":"<<n<<",\"planes\":"<<planes.size()<<",\"orbits_by_rank\":[";
 for(int r=0;r<n;++r){if(r)out<<',';out<<counts[r];}out<<"],\"profiles\":[";
 bool first=true;for(auto const&kv:profiles){if(!first)out<<',';first=false;out<<"{\"counts\":[";bool f=true;for(auto x:kv.first){if(!f)out<<',';f=false;out<<'[';for(int k=0;k<=n;++k){if(k)out<<',';out<<((x>>(6*k))&63);}out<<']';}out<<"],\"rays\":[";f=true;for(auto const&v:kv.second){if(!f)out<<',';f=false;out<<'[';for(int i=0;i<n;++i){if(i)out<<',';out<<v[i];}out<<']';}out<<"]}";}
 out<<"],\"complete\":true,\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
}
