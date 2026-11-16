"""Print the exact tables and constants quoted in the article."""
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
from math import comb, isqrt, factorial
from pathlib import Path
import json
from bernstein_bounds import ln_bounds

if not __debug__:
    raise RuntimeError('Run without Python optimization.')


def profile(w):
    atoms = {}
    for bits in product((0,1),repeat=len(w)):
        s = sum(a*b for a,b in zip(w,bits))
        atoms.setdefault(s,[0]*(len(w)+1))[sum(bits)] += 1
    return tuple(sorted(tuple(c) for c in atoms.values()))


def jets(weights, parameters, rates):
    atoms={}
    for bits in product((0,1),repeat=len(weights)):
        v,d,e=F(1),F(0),F(0)
        for bit,p,r in zip(bits,parameters,rates):
            b = p if bit else 1-p
            db = r if bit else -r
            v,d,e = v*b,d*b+v*db,e*b+2*d*db
        s=sum(w*b for w,b in zip(weights,bits))
        old=atoms.get(s,(F(0),)*3)
        atoms[s]=tuple(a+b for a,b in zip(old,(v,d,e)))
    return atoms


def main():
    witnesses={}
    for w in combinations_with_replacement(range(1,17),4):
        key=profile(w)
        if key not in witnesses or (sum(w),w)<(sum(witnesses[key]),witnesses[key]):
            witnesses[key]=w
    assert len(witnesses)==23
    for key,w in sorted(witnesses.items(),key=lambda item:(len(item[0]),item[1])):
        print('representative',w,'support',len(key))
    for m in (4,5):
        d=json.loads(Path('certificate_%d.json'%m).read_text())
        print('SUMMARY m=%d bound=%d vectors=%d profiles=%d cells=%d worst=%s cutoff=%s' %
              (m,isqrt(m**m),comb(isqrt(m**m)+m-1,m),d['profiles'],d['cells'],d['largest_upper_bound'],d['delta']))
    for m in (4,5,6):
        data=json.loads(Path(f'flats{m}_profiles.json').read_text())
        ps={tuple(sorted(map(tuple,p['counts']))) for p in data['profiles']}
        print('ORBITS m=%d by_rank=%s total=%d profiles=%d normals=%d' %
              (m,data['orbits_by_rank'],sum(data['orbits_by_rank']),len(ps),data['planes']))
        if m<6:
            base=json.loads(Path(f'certificate_{m}.json').read_text())
            assert ps=={tuple(sorted(map(tuple,p['counts']))) for p in base['certificates']}
            print('PROFILE SET AGREEMENT m=%d PASS'%m)
    d=json.loads(Path('certificate_summary.json').read_text())
    print('SUMMARY m=6 profiles=%d cells=%d worst=%s cutoff=%s' %
          (d['profiles_total'],d['cells'],d['largest_upper_bound'],d['delta']))
    maximum=[0]*6
    for line in Path('flats6_flats.txt').read_text().splitlines():
        r,*v=map(int,line.split());maximum[r]=max(maximum[r],max(map(abs,v)))
    print('MAXIMUM PRIMITIVE EXTERIOR COORDINATES',maximum)
    compact=json.loads(Path('profiles_six_compact.json').read_text())['profiles']
    print('COMPACT witnesses=%d maximum_coordinate=%d' %
          (len(compact),max(max(p['weights']) for p in compact)))
    geom=json.loads(Path('geometry_results.json').read_text())
    total=sum(sum(geom[k]) for k in ('feasible_extensions','infeasible_extensions','already_contained'))
    print('GEOMETRY candidate_extensions=%d feasible_extensions=%d' %
          (total,sum(geom['feasible_extensions'])))
    print('DETERMINANT m=6 Hadamard=216 divisibility=192 Ehlich_Wojtas=%d' % ((2*6-2)*(6-2)**(6//2-1)))
    for m in (6,7):
        raw=isqrt(m**m); divisor=1<<(m-1); refined=divisor*(raw//divisor)
        normals=(3**m-2*2**m+1)//2
        print('COST m=%d Hadamard_bound=%d sign_divisor=%d refined_bound=%d sorted_vectors=%d normals=%d permutations=%d' %
              (m,raw,divisor,refined,comb(refined+m-1,m),normals,factorial(m)))
    atoms=jets((1,1,2),(F(1,5),F(1,5),F(1,2)),(F(1),F(-1),F(0)))
    assert all(d==0 for v,d,e in atoms.values())
    argument=F(1)
    for v,d,e in atoms.values():
        assert e.denominator==1
        argument *= v**(-int(e))
    assert argument==F(289,256)>1
    print('HETEROGENEOUS probabilities',tuple(atoms[k][0] for k in sorted(atoms)))
    print('HETEROGENEOUS second derivatives',tuple(atoms[k][2] for k in sorted(atoms)))
    print('HETEROGENEOUS curvature log(%s)'%argument)
    p=F(1,2)
    q=p**4+(1-p)**4
    dq=4*p**3-4*(1-p)**3
    ddq=12*p*p+12*(1-p)**2
    assert (q,dq,ddq)==(F(1,8),F(0),F(6))
    print('ALL-EQUAL probability=%s derivative=%s second=%s curvature=6 log(7)'%(q,dq,ddq))
    joint={(0,0):(F(1,4),F(-1),F(8)),(0,1):(F(1,4),F(0),F(-8)),
           (1,0):(F(1,4),F(0),F(-8)),(1,1):(F(1,4),F(1),F(8))}
    assert tuple(sum(q[k] for q in joint.values()) for k in range(3))==(1,0,0)
    for i in (0,1):
        assert tuple(sum(q[k] for bits,q in joint.items() if bits[i]) for k in range(3))==(F(1,2),1,0)
    dep={}
    for bits,q in joint.items():
        old=dep.get(sum(bits),(F(0),)*3)
        dep[sum(bits)]=tuple(x+y for x,y in zip(old,q))
    print('DEPENDENT probability jets',tuple(dep[k] for k in sorted(dep)))
    print('DEPENDENT curvature=16 log(2)-8')
    assert ln_bounds(F(289,256))[0]>0
    assert 6*ln_bounds(F(7))[0]>0
    assert 16*ln_bounds(F(2))[0]-8>0
    print('COUNTEREXAMPLE SIGNS PASS (rational logarithm lower bounds)')

if __name__=='__main__':
    main()
