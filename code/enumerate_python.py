"""Alternative exhaustive enumeration, using integer polynomial products."""
from itertools import combinations_with_replacement
from math import comb
from pathlib import Path
import sys

if not __debug__:
    raise RuntimeError('Run without Python optimization.')

def profile(w):
    radix = len(w)+1
    polynomial = {0:1}
    for a in w:
        shift = radix*a+1
        nxt = polynomial.copy()
        for e,c in polynomial.items():
            nxt[e+shift] = nxt.get(e+shift,0)+c
        polynomial = nxt
    fibres = {}
    for e,c in polynomial.items():
        s,k = divmod(e,radix)
        fibres.setdefault(s,[0]*radix)[k] = c
    return tuple(sorted(tuple(c) for c in fibres.values()))

def main(m):
    bound = 16 if m==4 else 55
    seen = set()
    total = 0
    for w in combinations_with_replacement(range(1,bound+1),m):
        seen.add(profile(w))
        total += 1
    lines = Path('enumerated_%d.txt' % m).read_text().splitlines()
    header = tuple(map(int,lines[0].split()))
    expected = set()
    for line in lines[1:]:
        values = list(map(int,line.split()))
        assert values[m] == len(values[m+1:])
        expected.add(tuple(sorted(tuple((v >> (6*k)) & 63 for k in range(m+1)) for v in values[m+1:])))
    assert header == (m,bound,total,len(seen))
    assert total == comb(bound+m-1,m) and seen == expected
    print('POLYNOMIAL ENUMERATION PASS m=%d: %d vectors, %d profiles' % (m,total,len(seen)),flush=True)

if __name__ == '__main__':
    if len(sys.argv)==2:
        m=int(sys.argv[1]); assert m in (4,5); main(m)
    else:
        main(4); main(5)
