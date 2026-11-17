"""Deductions from the Part A rules (claim C44): stable states of a Region 1 block, computed exactly over all 2^25 states."""
import numpy as np, itertools
def fixed_points(p,q,c):
    """Exact stable states of one Region 1 block whose source is itself, for outside inputs p,q,c (0/1)."""
    out=[]
    N=1<<25; CH=1<<22
    for start in range(0,N,CH):
        x=np.arange(start,start+CH,dtype=np.uint32)
        S=lambda i,j: ((x>>((i-1)*5+(j-1)))&1).astype(bool)
        n=lambda a: ~a
        P,Q,C=bool(p),bool(q),bool(c)
        T=[n(S(3,1)), S(1,5)&n(S(1,3)), S(3,2), S(5,2)&n(S(4,5)), (S(1,4)|S(2,2)|S(2,3))&n(S(5,3)),
           S(2,5), np.full(CH,P&Q), S(4,4)&P, S(4,1), S(1,2)|S(3,1),
           S(4,2)&n(S(1,3)), S(2,4)&n(S(1,4)), n(S(1,1)), n(S(1,5))&(S(2,2)|S(4,4)), S(3,4),
           S(3,4)&S(5,4), n(S(3,2))&(S(3,1)|S(3,5)), S(2,5), np.full(CH,C), S(3,2)&n(S(4,3)),
           S(1,4), S(1,1)&S(2,1)&n(S(4,5)), S(3,2), n(S(4,3))&(S(3,2)|S(5,4)), S(3,2)]
        ok=np.ones(CH,bool)
        for k,t in enumerate(T): ok&= (((x>>k)&1).astype(bool)==t)
        out+= [int(v) for v in x[ok]]
    return out
