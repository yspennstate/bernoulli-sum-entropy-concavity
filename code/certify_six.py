"""Exact rational certification of the six-weight profiles, with checkpoints."""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
from math import comb
from functools import lru_cache
import argparse, json, time

import bernstein_bounds as b

@lru_cache(maxsize=131072)
def atom_fixed(c,a,z):
    # Round EACH atom outward before summation. This avoids multiplying many
    # unrelated logarithm denominators; the direction of the proof is unchanged.
    v=b.atom_upper(c,a,z)*(1<<40)
    return -((-v.numerator)//v.denominator)

def upper(cs,a,z):
    return F(sum(atom_fixed(c,a,z) for c in cs),1<<40)

def counts(w):
    n=len(w);atoms={}
    for mask in range(1<<n):
        s=sum(w[i] for i in range(n) if mask>>i&1)
        atoms.setdefault(s,[0]*(n+1))[bin(mask).count('1')]+=1
    return tuple(sorted(tuple(x) for x in atoms.values()))

def witness(p):
    cs=tuple(sorted(tuple(x) for x in p['counts']));vs=p['rays']
    # A nonzero univariate polynomial of degree v-1 has at most v-1 roots.
    # Testing (number of pairs)*(v-1)+1 positive integers must suffice.
    n=len(vs[0]);limit=comb(1<<n,2)*(len(vs)-1)+1
    for t in range(1,limit+1):
        w=[sum(v[i]*t**j for j,v in enumerate(vs)) for i in range(n)]
        assert min(w)>0
        if counts(w)==cs:return w,t
    raise AssertionError('root bound failed')

def main():
    if not __debug__:raise RuntimeError('Run with assertions enabled')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',default='generated_six')
    output=Path(parser.parse_args().output_dir)
    output.mkdir(parents=True,exist_ok=True)
    if (output/'certificate_six.jsonl').exists():
        raise FileExistsError('Preserve the previous certificate; generate in a fresh output folder')
    start=time.monotonic();data=json.loads((Path(__file__).resolve().parent/'flats6_profiles.json').read_text())
    assert data['complete'] and data['m']==6
    ps=data['profiles'];assert len({tuple(map(tuple,p['counts'])) for p in ps})==len(ps)
    # A single common endpoint cutoff, independent of profile or weight size.
    delta=F(1,32768);m=6
    ep=4*m*m*(m-1)*delta*15
    en=(1-delta)**(m-3)*(1-m*delta)**2
    assert ep<F(1,2) and en>F(3,4)
    certificates=[];ncells=0;worst=None;max_t=0
    for idx,p in enumerate(ps,1):
        cs=tuple(sorted(tuple(x) for x in p['counts']))
        assert [sum(c[k] for c in cs) for k in range(m+1)]==[comb(m,k) for k in range(m+1)]
        assert Counter(cs)==Counter(tuple(reversed(c)) for c in cs)
        w,t=witness(p);max_t=max(t,max_t)
        stack=[(delta,F(1,2),0)];leaves=[]
        while stack:
            a,z,depth=stack.pop();u=upper(cs,a,z)
            if u<0:
                scaled=u*10**6
                rounded=F(-((-scaled.numerator)//scaled.denominator),10**6)
                if rounded<0:
                    assert u<=rounded
                    leaves.append([str(a),str(z),str(rounded)])
                    ncells+=1;worst=rounded if worst is None else max(worst,rounded)
                    continue
            assert depth<26,(idx,w,a,z,str(u))
            mid=(a+z)/2;stack.extend([(mid,z,depth+1),(a,mid,depth+1)])
        assert F(leaves[0][0])==delta and F(leaves[-1][1])==F(1,2)
        assert all(x[1]==y[0] for x,y in zip(leaves,leaves[1:]))
        certificates.append({'id':idx,'weights':w,'counts':cs,'cells':leaves})
        if idx%100==0 or idx==len(ps):
            summary={'m':m,'complete':idx==len(ps),'profiles_processed':idx,'profiles_total':len(ps),'cells':ncells,'largest_upper_bound':str(worst),'delta':str(delta),'endpoint_positive_upper':str(ep),'endpoint_negative_lower':str(en),'max_witness_parameter':max_t,'seconds':time.monotonic()-start}
            print(json.dumps(summary),flush=True)
            (output/'certification_progress.json').write_text(json.dumps(summary)+'\n')
            # Append-only records survive interruption without large rewrites.
            with (output/'certificate_six.jsonl').open('a') as out:
                for item in certificates:out.write(json.dumps(item,separators=(',',':'))+'\n')
            certificates=[]
    (output/'certificate_summary.json').write_text(json.dumps(summary,indent=2)+'\n')

if __name__=='__main__':main()
