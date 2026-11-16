"""Replay every stored rational bound with direct Bernstein restriction.

The restriction coefficients are evaluated by the polar-form formula rather
than the subdivision routine used by certify.py. Log bounds use 18 terms.
Only the Python standard library is required; run without -O.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from pathlib import Path
import json

if not __debug__:
    raise RuntimeError('Run without Python optimization.')


def profile(w):
    atoms = {}
    for mask in range(1 << len(w)):
        s = sum(w[i] for i in range(len(w)) if mask & (1 << i))
        atoms.setdefault(s, [0] * (len(w) + 1))[bin(mask).count('1')] += 1
    return tuple(sorted(tuple(c) for c in atoms.values()))


@lru_cache(maxsize=32768)
def logarithm(x):
    assert x > 0
    exponent = 0
    y = x
    while y < 1:
        y *= 2
        exponent -= 1
    while y >= 2:
        y /= 2
        exponent += 1
    def unit(z):
        t = (z - 1) / (z + 1)
        lo = sum((2 * t ** (2*j + 1) / (2*j + 1) for j in range(18)), F(0))
        return lo, lo + 2*t**37 / (37*(1-t*t))
    lo, hi = unit(y)
    a, b = unit(F(2))
    return (lo + exponent*a, hi + exponent*b) if exponent >= 0 else (lo + exponent*b, hi + exponent*a)


@lru_cache(maxsize=32768)
def restriction_matrix(degree, a, b):
    # The j-th restricted coefficient is the blossom at degree-j copies of
    # a and j copies of b. Group equal arguments by their number of successes.
    rows = []
    for j in range(degree+1):
        row = [F(0)]*(degree+1)
        for k in range(degree-j+1):
            for ell in range(j+1):
                row[k+ell] += (comb(degree-j,k)*a**k*(1-a)**(degree-j-k)
                                 *comb(j,ell)*b**ell*(1-b)**(j-ell))
        assert sum(row) == 1 and min(row) >= 0
        rows.append(tuple(row))
    return tuple(rows)


def enclosure(c, a, b):
    values = [sum((x*y for x,y in zip(row,c)), F(0))
              for row in restriction_matrix(len(c)-1,a,b)]
    return min(values), max(values)


@lru_cache(maxsize=32768)
def atom_bound(counts, a, b):
    m = len(counts)-1
    c = tuple(F(counts[k],comb(m,k)) for k in range(m+1))
    d = tuple(m*(c[k+1]-c[k]) for k in range(m))
    e = tuple((m-1)*(d[k+1]-d[k]) for k in range(m-1))
    qlo,qhi = enclosure(c,a,b)
    dlo,dhi = enclosure(d,a,b)
    elo,ehi = enclosure(e,a,b)
    assert 0 < qlo <= qhi
    llo,lhi = logarithm(qlo)[0],logarithm(qhi)[1]
    lower = min(elo*llo,elo*lhi,ehi*llo,ehi*lhi)
    square = F(0) if dlo <= 0 <= dhi else min(dlo*dlo,dhi*dhi)
    return -lower-square/qhi


def main(m):
    data = json.loads(Path('certificate_%d.json' % m).read_text())
    lines = Path('enumerated_%d.txt' % m).read_text().splitlines()
    mm,bound,total,n = map(int,lines[0].split())
    assert (mm,bound,total,n) == ((4,16,3876,23) if m==4 else (5,55,5006386,216))
    enumerated = set()
    for line in lines[1:]:
        values = list(map(int,line.split()))
        w,size,encoded = values[:m],values[m],values[m+1:]
        counts = tuple(sorted(tuple((x >> (6*k)) & 63 for k in range(m+1)) for x in encoded))
        assert size == len(encoded) and counts == profile(w)
        enumerated.add(counts)
    assert len(enumerated) == n == len(lines)-1
    delta = F(1,4096 if m==4 else 16384)
    assert F(data['delta']) == delta and data['pass'] and data['profiles'] == n
    cells = 0
    worst = None
    seen = set()
    for certificate in data['certificates']:
        counts = tuple(tuple(c) for c in certificate['counts'])
        assert counts == profile(certificate['weights']) and counts not in seen
        seen.add(counts)
        end = delta
        for left,right,stored in certificate['cells']:
            a,b,u = F(left),F(right),F(stored)
            assert a == end and a < b <= F(1,2) and u < 0
            actual = sum((atom_bound(c,a,b) for c in counts),F(0))
            assert actual <= u, (m,certificate['weights'],a,b)
            worst = u if worst is None else max(worst,u)
            cells += 1
            end = b
        assert end == F(1,2)
    assert seen == enumerated
    assert cells == data['cells'] and worst == F(data['largest_upper_bound'])
    positive = 4*m*m*(m-1)*delta*(12 if m==4 else 14)
    negative = (1-delta)**(m-3)*(1-m*delta)**2
    assert positive == F(data['endpoint_positive_upper'])
    assert negative == F(data['endpoint_negative_lower'])
    assert positive <= (F(9,16) if m==4 else F(1,2)) and negative > F(3,4)
    print('REPLAY PASS m=%d: %d vectors, %d profiles, %d cells, largest upper bound %s' % (m,total,n,cells,worst),flush=True)
    print('ENDPOINT PASS m=%d: positive upper %s, negative lower %s' % (m,positive,negative),flush=True)
    atom_bound.cache_clear()
    logarithm.cache_clear()
    restriction_matrix.cache_clear()

if __name__ == '__main__':
    main(4)
    main(5)
