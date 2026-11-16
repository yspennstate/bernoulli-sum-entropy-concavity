"""Independent rational cell replayer: power shifts, Bernstein conversion.

Uses monomial shifts rather than generator subdivision.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from pathlib import Path
import json,time

def derivative(c):return tuple((i+1)*c[i+1] for i in range(len(c)-1))

@lru_cache(maxsize=1024)
def polynomials(counts):
    n=len(counts)-1
    c=tuple(sum(counts[k]*comb(n-k,j-k)*(-1)**(j-k) for k in range(j+1)) for j in range(n+1))
    d=derivative(c)
    return c,d,derivative(d)

def enclosure(c,a,b):
    # Expand f(a+(b-a)t), then convert powers t^j into degree-d Bernstein.
    d=len(c)-1
    powers_a=[a**i for i in range(d+1)]
    powers_h=[(b-a)**i for i in range(d+1)]
    shifted=[powers_h[j]*sum(c[i]*comb(i,j)*powers_a[i-j] for i in range(j,d+1)) for j in range(d+1)]
    beta=[sum((shifted[j]*F(comb(k,j),comb(d,j)) for j in range(k+1)),F(0)) for k in range(d+1)]
    return min(beta),max(beta)

@lru_cache(maxsize=8192)
def log_interval(x):
    assert x>0
    # Bit lengths choose scaling; exact comparisons correct its two boundaries.
    k=x.numerator.bit_length()-x.denominator.bit_length()
    y=x/F(2)**k
    while y<1:y*=2;k-=1
    while y>=2:y/=2;k+=1
    def series(v):
        z=(v-1)/(v+1);z2=z*z;power=z;s=F(0)
        for j in range(18):s+=2*power/F(2*j+1);power*=z2
        return s,s+2*power/(37*(1-z2))
    lo,hi=series(y);l2,u2=series(F(2))
    return (lo+k*l2,hi+k*u2) if k>=0 else (lo+k*u2,hi+k*l2)

@lru_cache(maxsize=131072)
def atom_bound(counts,a,b):
    c,d,e=polynomials(counts)
    ql,qu=enclosure(c,a,b);dl,du=enclosure(d,a,b);el,eu=enclosure(e,a,b)
    assert ql>0 and qu>0
    ll=log_interval(ql)[0];lu=log_interval(qu)[1]
    product=min(el*ll,el*lu,eu*ll,eu*lu)
    square=F(0) if dl<=0<=du else min(dl*dl,du*du)
    v=(-product-square/qu)*(1<<40)
    return -((-v.numerator)//v.denominator)

def bound(cs,a,b):return F(sum(atom_bound(tuple(c),a,b) for c in cs),1<<40)

def reconstruct(w):
    fibres={0:[1,0,0,0,0,0,0]}
    # Dynamic polynomial product, independent of the bit-mask generator.
    for x in w:
        old={s:c[:] for s,c in fibres.items()}
        for s,c in old.items():
            dest=fibres.setdefault(s+x,[0]*7)
            for k in range(6):dest[k+1]+=c[k]
    return tuple(sorted(tuple(v) for v in fibres.values()))

def main():
    if not __debug__:raise RuntimeError('Run with assertions enabled')
    start=time.monotonic();summary=json.loads(Path('certificate_summary.json').read_text())
    assert summary['complete'] and summary['profiles_total']==6767
    listed=json.loads(Path('flats6_profiles.json').read_text())
    expected={tuple(sorted(tuple(c) for c in p['counts'])) for p in listed['profiles']}
    seen=set();cells=0;worst=None
    for idx,line in enumerate(Path('certificate_six.jsonl').open(),1):
        p=json.loads(line);w=p['weights'];cs=tuple(sorted(tuple(c) for c in p['counts']))
        assert p['id']==idx and min(w)>0 and cs==reconstruct(w)
        assert cs in expected and cs not in seen;seen.add(cs)
        end=F(1,32768)
        for aa,bb,uu in p['cells']:
            a,b,u=F(aa),F(bb),F(uu)
            assert a==end and a<b<=F(1,2) and u<0
            actual=bound(cs,a,b)
            assert actual<=u,(idx,a,b,actual,u)
            end=b;cells+=1;worst=u if worst is None else max(worst,u)
        assert end==F(1,2)
        if idx%500==0:print(json.dumps({'profiles':idx,'cells':cells,'seconds':time.monotonic()-start}),flush=True)
    assert seen==expected and cells==summary['cells'] and str(worst)==summary['largest_upper_bound']
    delta=F(1,32768)
    assert 720*15*delta<F(1,2) and (1-delta)**3*(1-6*delta)**2>F(3,4)
    report={'pass':True,'profiles':len(seen),'cells':cells,'largest_upper_bound':str(worst),'seconds':time.monotonic()-start,'implementation':'power shift to Bernstein, 18 log terms, outward atom rounding'}
    Path('replay_results.json').write_text(json.dumps(report,indent=2)+'\n')
    assert cells==130172 and worst==F(-1,20000)
    delta=F(1,32768)
    ep=720*delta*15; en=(1-delta)**3*(1-6*delta)**2
    assert ep<F(1,2) and en>F(3,4)
    print('ENDPOINT PASS m=6: positive upper',ep,'negative lower',en)
    print(json.dumps(report),flush=True)

if __name__=='__main__':main()
