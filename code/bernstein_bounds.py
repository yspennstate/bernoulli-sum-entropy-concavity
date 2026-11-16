"""Rational Bernstein and logarithm enclosures for certificate generation."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb

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
