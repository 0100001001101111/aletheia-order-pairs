"""Exact finite-order checks. Python standard library only; no network access."""
import json,itertools,math,time
from fractions import Fraction as F
from collections import Counter

def trim(p):
    p=list(map(F,p))
    while len(p)>1 and not p[-1]:p.pop()
    return p

def add(a,b):return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def scale(a,c):return trim([c*x for x in a])
def mul(a,b):
    p=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):p[i+j]+=x*y
    return trim(p)
def power(a,n):
    p=[F(1)]
    for _ in range(n):p=mul(p,a)
    return p
def compose(a,b):
    p=[F(0)]
    for x in reversed(a):p=add(mul(p,b),[x])
    return p

def charpoly(M):
    """Exact Newton trace recurrence, used only on matrices of size at most10."""
    n=len(M);P=[[int(i==j) for j in range(n)] for i in range(n)];co=[F(1)];tr=[]
    for k in range(1,n+1):
        P=[[sum(P[i][l]*M[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        tr.append(sum(P[i][i] for i in range(n)));co.append(-sum(co[k-i]*tr[i-1] for i in range(1,k+1))/k)
    return list(reversed(co))

def gram(v,cols):return [[sum(i in c and j in c for c in cols) for j in range(v)] for i in range(v)]
def eyej(v,a,b):return [[a*int(i==j)+b for j in range(v)] for i in range(v)]
def transform_spectrum(p,factor,offset=0):
    n=len(p)-1;return scale(compose(p,[F(-offset,factor),F(1,factor)]),factor**n)

def spectral_lists(v,cols,kind):
    m=len(cols);n=v+m;G=gram(v,cols);q=charpoly(G);z=[F(0),F(1)];t=[F(-1),F(1)];pad=power(z,m-v)
    # Singular-value block identities give all coefficients, including zero modes.
    sa=mul(pad,compose(q,[0,0,1]));asa=scale(mul(pad,compose(q,[0,0,-1])),(-1)**v)
    comm=scale(mul(mul(pad,q),compose(q,[0,-1])),(-1)**v)
    if kind==30:
        assert {len(c) for c in cols}=={2}
        assert all(sum(i in c for c in cols)==4 for i in range(v))
        L=mul(power([-2,1],m-v),compose(q,[8,-6,1]))
        normalized=scale(mul(power(t,m-v),compose(q,scale(power(t,2),8))),F(1,8**v))
        half=scale(mul(power(t,m-v),compose(q,scale(power(t,2),2))),F(1,2**v))
    else:
        assert v==7 and m==63
        classes=[]
        for k in [1,2,3]:
            c=[x for x in cols if len(x)==k];g=gram(v,c);b=g[0][1];a=g[0][0]-b
            assert g==eyej(v,a,b),'CLASS_GRAM_STRUCTURE';classes.append((k,len(c),a,b))
        assert all(sum(i in c for c in cols)==22 for i in range(v))
        factors=[[-k,1] for k,_,_,_ in classes];D=[F(1)]
        for f in factors:D=mul(D,f)
        pa=mul([-22,1],D);pc=list(pa);outside=[F(1)]
        for i,(k,count,a,b) in enumerate(classes):
            other=[F(1)]
            for j,f in enumerate(factors):
                if i!=j:other=mul(other,f)
            pa=add(pa,scale(other,-a));pc=add(pc,scale(other,-a-7*b));outside=mul(outside,power(factors[i],count-v))
        L=mul(outside,mul(power(pa,6),pc))
        a=sum(F(a,22*k) for k,count,a,b in classes);c=sum(F(a+7*b,22*k) for k,count,a,b in classes)
        normalized=mul(power(t,m-v),mul(power(add(power(t,2),[-a]),6),add(power(t,2),[-c])))
        half=mul(power(t,m-v),mul(power(add(power(t,2),[-4*a]),6),add(power(t,2),[-4*c])))
    polys={'link_laplacian':L,'normalized_link':normalized,'halfdegree_link':half,'normalized_active':normalized,'halfdegree_active':half,
      'C_sum':sa,'C_difference':asa,'C_commutator':comm,'B_sum':transform_spectrum(sa,2,-2),
      'B_difference':transform_spectrum(asa,2),'B_commutator':transform_spectrum(comm,4),'C_square_sum':power(z,n)}
    assert all(len(p)==n+1 and p[-1]==1 for p in polys.values())
    return {k:list(map(str,reversed(p))) for k,p in polys.items()}

def coordinate_check(points,covers):
    p=[[F(x) for x in row] for row in points];n=len(p);assert len(set(map(tuple,p)))==n,'DUPLICATE_POINT'
    expected=set(map(tuple,covers));margin=None
    for i in range(n):
        for j in range(n):
            if i==j:continue
            dt=p[j][0]-p[i][0];q=dt*dt-sum((p[j][k]-p[i][k])**2 for k in range(1,len(p[0])))
            assert ((dt>0 and q>0)==((i,j) in expected)),'COORDINATE_RELATION'
            assert q!=0,'NULL_PAIR';margin=abs(q) if margin is None else min(margin,abs(q))
    return {'ordered_pairs':n*(n-1),'minimum_absolute_quadratic_margin':str(margin)}

def invariants(v,cols,kind):
    m=len(cols);n=v+m;r=sum(map(len,cols));out=[sum(i in c for c in cols) for i in range(v)];inc=list(map(len,cols))
    assert all(out) and all(inc)
    return {'n':n,'chains':[1,n,r]+[0]*(n-2),'interval_counts':[r]+[0]*(n-2),
      'MM_input':str(F(r,n*n)),'ordering_fraction':str(F(r,n*(n-1)//2)),
      'interval_probabilities':['1']+['0']*(n-2),'scaled_interval_density':[n]+[0]*(n-2),
      'actions':[n-2*r,n-r],'link_degrees':sorted(out+inc),
      'half_degree_histogram':[[str(F(k,2)),str(F(c,n))] for k,c in sorted(Counter(out+inc).items())],
      'directed_link_degrees':sorted([[0,d] for d in out]+[[d,0] for d in inc]),
      'height_elements':2,'minimal_count':v,'minimal_heights':[0,v]+[0]*(n-2),'minimal_height_probabilities':['0','1']+['0']*(n-2),
      'polynomials':spectral_lists(v,cols,kind)}

def score(payload):
    start=time.monotonic();d=payload['data']
    assert charpoly([[1,-1],[-1,1]])==[0,-2,1]
    assert compose([-1,0,1],[1,1])==[0,2,1]
    assert transform_spectrum([-1,0,1],2,-2)==[0,4,1]
    assert coordinate_check([[0,0],[2,1]],[(0,1)])['minimum_absolute_quadratic_margin']=='3'
    caught=False
    try:coordinate_check([[0,0],[1,2]],[(0,1)])
    except AssertionError as e:caught=str(e)=='COORDINATE_RELATION'
    assert caught
    a=d['pair30'];graphs=[a['P_edges'],a['Q_edges']];q=[]
    for tag,edges in zip(['P','Q'],graphs):assert a[tag+'_covers']==[[v,10+k] for k,e in enumerate(edges) for v in e],'COVER_LIST_30'
    for edges in graphs:
        assert len(edges)==20 and len({tuple(sorted(e)) for e in edges})==20
        assert Counter(v for e in edges for v in e)=={i:4 for i in range(10)}
        adj=[[0]*10 for _ in range(10)]
        for u,v in edges:assert u!=v;adj[u][v]=adj[v][u]=1
        q.append(charpoly(adj))
    assert q[0]==q[1]==[80,200,9,-320,-180,136,110,-16,-20,0,1]
    paths=a['Q_subdivision_paths'];edges=set(map(tuple,graphs[1]));branch={1,2,7,3,8,9};inner=[]
    assert {(p[0],p[-1]) for p in paths}==set(itertools.product([1,2,7],[3,8,9]))
    for path in paths:
        assert all(tuple(sorted((u,v))) in edges for u,v in zip(path,path[1:]))
        assert not branch.intersection(path[1:-1]);inner+=path[1:-1]
    assert len(set(inner))==len(inner)
    covers=[(v,10+k) for k,e in enumerate(graphs[0]) for v in e];coord30=coordinate_check(a['P_coordinates'],covers)
    bad=[p[:] for p in a['P_coordinates']];bad[0][0]=1000000;caught=False
    try:coordinate_check(bad,covers)
    except AssertionError as e:caught=str(e)=='COORDINATE_RELATION'
    assert caught
    x=[invariants(10,e,30) for e in graphs];assert x[0]==x[1],'THIRTY_LISTS'
    A=[list(c) for k in [1,2,3] for c in itertools.combinations(range(7),k)]
    B=[list(c) for k in [1,2] for c in itertools.combinations(range(7),k)]+[sorted([i,(i+1)%7,(i+3)%7]) for _ in range(5) for i in range(7)]
    b=d['pair70'];assert b['B_subsets']==B,'COORDINATE_LABEL_ORDER'
    for tag,cols in [('A',A),('B',B)]:assert b[tag+'_covers']==[[v,7+k] for k,c in enumerate(cols) for v in c],'COVER_LIST_70'
    coord70=coordinate_check(b['B_coordinates'],[(v,7+k) for k,col in enumerate(B) for v in col])
    assert [len(set(map(tuple,z))) for z in [A,B]]==[63,35]
    y=[invariants(7,z,70) for z in [A,B]];assert y[0]==y[1],'SEVENTY_LISTS'
    changed=[c[:] for c in B];changed[0]=[1];assert gram(7,changed)!=gram(7,B)
    return {'status':'PASS','coordinate30':coord30,'coordinate70':coord70,'thirty_invariants':x[0],'seventy_invariants':y[0],
      'proof_limits':'Finite equalities, coordinates and obstruction data checked. Nonembedding also uses the written circle-order theorem and Radon argument. No formal proof certificate is claimed.',
      'controls':'two-point polynomial, polynomial substitution, timelike and spacelike controls, changed coordinate and incidence all passed','seconds':time.monotonic()-start}

if __name__=='__main__':
    import os,sys,signal,resource,hashlib
    from pathlib import Path
    for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
    try:os.nice(15)
    except OSError:pass
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('55-second limit')));signal.alarm(55)
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent.parent
    data={k:json.loads((root/'pairs'/(k+'.json')).read_text()) for k in ['pair30','pair70']}
    result=score({'data':data});peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024);assert peak<1000000000
    for line in ['known answers and changed-coordinate/changed-incidence controls','30-element pair: both fixed lists and strict2+1 coordinate certificate','30-element pair: K3,3 subdivision data','70-element pair: both fixed lists and strict3+1 coordinate certificate','70-element pair: all required subset tops for the written Radon argument']:print('PASS '+line)
    print('Nonembedding conclusions also depend on the written proofs and cited theorem. No sprinkled-set claim.')
    print('seconds: %.3f; peak memory bytes: %s'%(result['seconds'],peak))
    for p in [root/'pairs/pair30.json',root/'pairs/pair70.json']:print(p.name+' sha256:'+hashlib.sha256(p.read_bytes()).hexdigest())
