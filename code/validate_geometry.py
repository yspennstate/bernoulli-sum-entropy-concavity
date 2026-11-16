"""Independent checks of all listed flats using integer row reduction.

This validates every record, not a second exhaustive generator.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['VECLIB_MAXIMUM_THREADS']='1'
from itertools import product,combinations
from pathlib import Path
from math import gcd
from functools import reduce
from collections import Counter
import json,time


def rank(rows):
    a=[list(v) for v in rows];r=0
    for j in range(6):
        pivot=next((i for i in range(r,len(a)) if a[i][j]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]
        for i in range(r+1,len(a)):
            if not a[i][j]:continue
            u,v=a[r][j],a[i][j]
            a[i]=[u*x-v*y for x,y in zip(a[i],a[r])]
            g=reduce(gcd,a[i]);a[i]=[x//g for x in a[i]] if g else a[i]
        r+=1
        if r==len(a):break
    return r

def determinant(a):
    a=[list(v) for v in a];n=len(a)
    if not n:return 1
    sign=1;previous=1
    for j in range(n-1):
        pivot=next((i for i in range(j,n) if a[i][j]),None)
        if pivot is None:return 0
        if pivot!=j:a[pivot],a[j]=a[j],a[pivot];sign=-sign
        for i in range(j+1,n):
            for k in range(j+1,n):
                x=a[i][k]*a[j][j]-a[i][j]*a[j][k]
                assert x%previous==0;a[i][k]=x//previous
        previous=a[j][j]
        for i in range(j+1,n):a[i][j]=0
    return sign*a[-1][-1]

def main():
    if not __debug__:raise RuntimeError('Run with assertions enabled')
    start=time.monotonic();subsets=[list(combinations(range(6),r)) for r in range(7)]
    # Match the source's increasing bit-mask order without using its functions.
    for ss in subsets:ss.sort(key=lambda s:sum(1<<i for i in s))
    planes=[a for a in product((-1,0,1),repeat=6) if 1 in a and -1 in a and next(x for x in a if x)==1]
    expected={tuple(sorted(tuple(c) for c in p['counts'])) for p in json.loads(Path('flats6_profiles.json').read_text())['profiles']}
    found=set();ranks=Counter();maxcoord=Counter()
    for num,line in enumerate(Path('flats6_flats.txt').open(),1):
        r,*p=map(int,line.split());ranks[r]+=1;maxcoord[r]=max(maxcoord[r],max(map(abs,p)))
        assert len(p)==len(subsets[r])
        pivot=next(i for i,x in enumerate(p) if x);I=subsets[r][pivot];D=p[pivot]
        assert D>0
        lookup=dict(zip(subsets[r],p))
        A=[[D*int(k==j) for k in range(6)] for j in I]
        for h,i in enumerate(I):
            for j in range(6):
                if j in I:continue
                J=tuple(sorted(set(I)-{i}|{j}));t=J.index(j)
                A[h][j]=(-1)**(h+t)*lookup[J]
        # Check decomposability and every minor, independently of wedge growth.
        if r:
            assert [determinant([[row[j] for j in J] for row in A]) for J in subsets[r]]==[D**(r-1)*x for x in p]
        null=[]
        for j in range(6):
            if j in I:continue
            v=[0]*6;v[j]=D
            for h,i in enumerate(I):v[i]=-A[h][j]
            null.append(v)
        assert len(null)==6-r
        contained=[a for a in planes if all(sum(x*y for x,y in zip(a,v))==0 for v in null)]
        assert rank(contained)==r
        # Recover all nonnegative circuits directly from the minors.
        positive=[0]*6
        for J in subsets[r+1]:
            v=[0]*6
            for t,j in enumerate(J):v[j]=(-1)**t*lookup[tuple(k for k in J if k!=j)]
            if max(v)<=0:v=[-x for x in v]
            if min(v)>=0:
                positive=[x+y for x,y in zip(positive,v)]
        assert min(positive)>0
        assert all(sum(a*w for a,w in zip(row,positive))==0 for row in A)
        # Recover subset-sum fibres from a nullspace basis, rather than rays.
        atoms={}
        for bits in product((0,1),repeat=6):
            signature=tuple(sum(a*b for a,b in zip(bits,v)) for v in null)
            atoms.setdefault(signature,[0]*7)[sum(bits)]+=1
        found.add(tuple(sorted(tuple(c) for c in atoms.values())))
        if num%5000==0:print(num,'validated',time.monotonic()-start,flush=True)
    assert found==expected
    result={'pass':True,'flats':sum(ranks.values()),'orbits_by_rank':dict(ranks),'profiles':len(found),'maximum_exterior_coordinate_by_rank':dict(maxcoord),'seconds':time.monotonic()-start}
    Path('geometry_record_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
