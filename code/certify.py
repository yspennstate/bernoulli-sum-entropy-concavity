"""Exact certificates using Bernstein convex-hull enclosures.

Enumeration comes from enumerate.cpp. Each
probability polynomial and its derivatives are represented in the Bernstein
basis, restricted to a rational cell by de Casteljau's algorithm. Rational
logarithm bounds use the integral/geometric expansion of atanh, with 14 terms.
All inequalities in the certificate are exact Fraction comparisons.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from pathlib import Path
import json, time

if not __debug__:
    raise RuntimeError('Run without Python optimization: assertions are proof checks.')

def canonical(cs):
    return tuple(sorted(tuple(c) for c in cs))

def read_enumeration(m):
    lines = Path(f'enumerated_{m}.txt').read_text().splitlines()
    mm, bound, total, n = map(int, lines[0].split())
    assert mm == m and total == comb(bound+m-1, m)
    result = []
    for line in lines[1:]:
        nums = list(map(int, line.split()))
        w, size, encoded = nums[:m], nums[m], nums[m+1:]
        assert size == len(encoded)
        counts = canonical([[(v >> (6*k)) & 63 for k in range(m+1)] for v in encoded])
        assert [sum(c[k] for c in counts) for k in range(m+1)] == [comb(m,k) for k in range(m+1)]
        # A third check of each representative: direct dictionary of bit strings.
        atoms = {}
        for mask in range(1<<m):
            s = sum(w[i] for i in range(m) if mask & (1<<i))
            atoms.setdefault(s,[0]*(m+1))[bin(mask).count('1')] += 1
        assert counts == canonical(atoms.values())
        result.append((w,counts))
    assert len(result) == n and len({cs for _,cs in result}) == n
    print(f'ENUMERATION PASS m={m}: {total} vectors, {n} profiles, representative counts checked', flush=True)
    return result

def split_coefficients(coeff,t):
    work = list(coeff)
    left, right = [work[0]], [work[-1]]
    while len(work)>1:
        work = [(1-t)*a+t*b for a,b in zip(work,work[1:])]
        left.append(work[0]); right.append(work[-1])
    return tuple(left),tuple(reversed(right))

def derivative(coeff):
    degree=len(coeff)-1
    return tuple(degree*(b-a) for a,b in zip(coeff,coeff[1:]))

def bounds(coeff,left,right):
    # Restrict first to [left,1], then [left,right]. Bernstein basis functions
    # are nonnegative and sum to one, hence the convex-hull enclosure.
    _, tail = split_coefficients(coeff,left)
    restricted, _ = split_coefficients(tail,(right-left)/(1-left))
    return min(restricted),max(restricted)

@lru_cache(maxsize=32768)
def ln_bounds(x):
    assert x>0
    # Normalize by repeated exact scaling, independently of bit-length code.
    y=x; k=0
    while y<1: y*=2; k-=1
    while y>=2: y/=2; k+=1
    def log_unit(z):
        t=(z-1)/(z+1); power=t; answer=F(0)
        for j in range(14):
            answer += 2*power/(2*j+1)
            power *= t*t
        error=2*power/(29*(1-t*t))
        return answer,answer+error
    lo,hi=log_unit(y); a,b=log_unit(F(2))
    return (lo+k*a,hi+k*b) if k>=0 else (lo+k*b,hi+k*a)

@lru_cache(maxsize=32768)
def atom_upper(counts,left,right):
    m=len(counts)-1
    coeff=tuple(F(counts[k],comb(m,k)) for k in range(m+1))
    first=derivative(coeff); second=derivative(first)
    qlo,qhi=bounds(coeff,left,right)
    dlo,dhi=bounds(first,left,right)
    elo,ehi=bounds(second,left,right)
    assert qlo>0 and qhi>0
    llo,lhi=ln_bounds(qlo)[0],ln_bounds(qhi)[1]
    product_lower=min(elo*llo,elo*lhi,ehi*llo,ehi*lhi)
    square_lower=F(0) if dlo<=0<=dhi else min(dlo*dlo,dhi*dhi)
    return -product_lower-square_lower/qhi

def upper(cs,left,right):
    return sum((atom_upper(c,left,right) for c in cs),F(0))

def certify(m):
    start=time.monotonic(); profiles=read_enumeration(m)
    delta=F(1,4096 if m==4 else 16384)
    endpoint_positive=4*m*m*(m-1)*delta*(12 if m==4 else 14)
    endpoint_negative=(1-delta)**(m-3)*(1-m*delta)**2
    assert endpoint_positive < (F(3,4)) and endpoint_negative > F(3,4)
    certificates=[]; worst=None; cells=0
    for idx,(w,cs) in enumerate(profiles,1):
        stack=[(delta,F(1,2),0)]; leaves=[]
        while stack:
            a,b,depth=stack.pop(); u=upper(cs,a,b)
            if u<0:
                # Outward rounding produces compact rational certificates.
                scaled=u*10**6
                rounded=F(-((-scaled.numerator)//scaled.denominator),10**6)
                # Small negative bounds could round to zero; refine instead.
                if rounded<0:
                    assert u<=rounded
                    leaves.append([str(a),str(b),str(rounded)]); cells+=1
                    worst=rounded if worst is None else max(worst,rounded)
                    continue
            assert depth<24,(m,w,a,b)
            c=(a+b)/2
            stack.extend([(c,b,depth+1),(a,c,depth+1)])
        assert F(leaves[0][0])==delta and F(leaves[-1][1])==F(1,2)
        assert all(x[1]==y[0] for x,y in zip(leaves,leaves[1:]))
        certificates.append({'weights':w,'counts':cs,'cells':leaves})
        if idx%25==0: print(f'm={m}: profiles={idx}, cells={cells}, seconds={time.monotonic()-start:.2f}',flush=True)
    report={'m':m,'profiles':len(profiles),'delta':str(delta),'cells':cells,'largest_upper_bound':str(worst),'endpoint_positive_upper':str(endpoint_positive),'endpoint_negative_lower':str(endpoint_negative),'seconds':time.monotonic()-start,'certificates':certificates,'pass':True}
    Path(f'certificate_{m}.json').write_text(json.dumps(report,indent=2)+'\n')
    print('CERTIFICATE PASS', {k:v for k,v in report.items() if k!='certificates'},flush=True)

if __name__=='__main__':
    certify(4)
    atom_upper.cache_clear(); ln_bounds.cache_clear()
    certify(5)
