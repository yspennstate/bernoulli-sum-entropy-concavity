"""Exact cell verification, using Bernstein finite differences and de Casteljau.
Uses normalized Bernstein differences and rational subdivision.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from pathlib import Path
from collections import Counter
import json,time
ROOT=Path(__file__).resolve().parent
S=ROOT

def profile(weights):
    # Encode support sum and cardinality as one polynomial exponent.
    n=len(weights);D=n+1;coeff={0:1}
    for w in weights:
        shift=D*w+1;new=coeff.copy()
        for exponent,c in coeff.items():new[exponent+shift]=new.get(exponent+shift,0)+c
        coeff=new
    atoms={}
    for exponent,c in coeff.items():
        s,k=divmod(exponent,D);atoms.setdefault(s,[0]*(n+1))[k]=c
    return tuple(sorted(map(tuple,atoms.values())))

def split(c,t):
    left=[c[0]];right=[c[-1]]
    while len(c)>1:
        c=tuple((1-t)*x+t*y for x,y in zip(c,c[1:]));left.append(c[0]);right.append(c[-1])
    return tuple(left),tuple(reversed(right))

@lru_cache(2048)
def bernstein(row):
    n=len(row)-1;c=tuple(Q(v,comb(n,k)) for k,v in enumerate(row))
    d=tuple(n*(y-x) for x,y in zip(c,c[1:]))
    e=tuple((n-1)*(y-x) for x,y in zip(d,d[1:]))
    return c,d,e

def enclose(c,a,b):
    _,right=split(c,a)
    left,_=split(right,(b-a)/(1-a))
    return min(left),max(left)

@lru_cache(16000)
def logq(x):
    assert x>0
    k=0;y=x
    while y<1:y*=2;k-=1
    while y>2:y/=2;k+=1
    def series(y):
        z=(y-1)/(y+1);power=z;total=Q(0)
        for j in range(22):total+=2*power/(2*j+1);power*=z*z
        return total,total+2*power/(45*(1-z*z))
    l,h=series(y);l2,h2=series(Q(2))
    if k>=0:return l+k*l2,h+k*h2
    return l+k*h2,h+k*l2

@lru_cache(131072)
def atom_upper(row,a,b):
    cs=bernstein(row)
    (ql,qh),(dl,dh),(el,eh)=[enclose(c,a,b) for c in cs]
    assert ql>0
    ll=logq(ql)[0];lh=logq(qh)[1]
    upper=-min(el*ll,el*lh,eh*ll,eh*lh)
    low_square=0 if dl<=0<=dh else min(dl*dl,dh*dh)
    v=(upper-low_square/qh)*(1<<42)
    return -((-v.numerator)//v.denominator)

def upper(profile,a,b):return Q(sum(atom_upper(row,a,b) for row in profile),1<<42)

def key(p):return tuple(sorted(map(tuple,p['counts'])))

def main():
    if not __debug__:raise RuntimeError('Run with assertions enabled')
    started=time.monotonic()
    ps=json.loads((S/'flats6_profiles.json').read_text())['profiles'];expected={key(p) for p in ps}
    compact=json.loads((S/'profiles_six_compact.json').read_text())['profiles']
    assert len(expected)==6767==len(compact)
    leading=[]
    for p in compact:
        w=p['weights'];cs=key(p)
        assert len(w)==6 and all(isinstance(x,int) and 1<=x<=36 for x in w)
        assert profile(w)==cs
        assert [sum(row[k] for row in cs) for k in range(7)]==[comb(6,k) for k in range(7)]
        assert Counter(cs)==Counter(tuple(reversed(c)) for c in cs)
        assert sum(c[1] for c in cs if c[0]==0)==6
        assert sum(c[0] for c in cs)==1
        leading.append(sum(c[1] for c in cs if c[0]==0))
    assert {key(p) for p in compact}==expected
    seen=set();ncells=0;storedmax=None;maxactual=None
    for line in (S/'certificate_six.jsonl').open():
        p=json.loads(line);cs=key(p)
        assert p['id']==len(seen)+1 and cs not in seen and cs in expected
        assert profile(p['weights'])==cs
        seen.add(cs);last=Q(1,32768)
        for sa,sb,su in p['cells']:
            a,b,u=map(Q,(sa,sb,su));assert a==last and a<b<=Q(1,2) and u<0
            actual=upper(cs,a,b);assert actual<=u,(p['id'],a,b,actual,u)
            last=b;ncells+=1
            storedmax=u if storedmax is None else max(u,storedmax)
            maxactual=actual if maxactual is None else max(actual,maxactual)
        assert last==Q(1,2)
        if len(seen)%1000==0:print('profiles',len(seen),'cells',ncells,'seconds',round(time.monotonic()-started,2),flush=True)
    assert seen==expected and ncells==130172 and storedmax==Q(-1,20000)
    report={'pass':True,'profiles':len(seen),'cells':ncells,'largest_stored_upper':str(storedmax),'largest_recomputed_upper':str(maxactual),'compact_max_coordinate':max(max(p['weights']) for p in compact),'asymptotic_singleton_coefficient_range':[min(leading),max(leading)],'seconds':time.monotonic()-started}
    (ROOT/'verification_six_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
if __name__=='__main__':main()
