def get_irreducible(deg,corps):
    F=corps
    R.<z>=PolynomialRing(F)
    g=R.random_element(degree=deg)
    while not g.is_irreducible():
        g=R.random_element(degree=deg)
    return g
    
def get_distinct_elements(n,corps):
    L=[]
    for i in range(n):
        e=corps.random_element()
        while e in L:
            e=corps.random_element()
        L.append(e)
    return L
    
def gen_Goppa(corps,g,L):
    r=g.degree()
    n=len(L)
    M=MatrixSpace(corps,r,n)
    D=diagonal_matrix(corps,vector(corps,[1/g(s) for s in L]))
    V=M(0)
    for i in range(r):
        for j in range(n):
            V[i,j]=L[j]**i
    
    C=V*D
    G=C.right_kernel().basis_matrix()
    return G
    
def sum_inverse(inv,mot,polyn):
    s=0
    n=len(inv)
    for i in range(n):
        s+=mot[i]*inv[i]
    return s%polyn
    
    
    
def euclide_etendu_partiel(a,b,t,Ann):
    r0,u0,v0=a,Ann(1),Ann(0)
    r1,u1,v1 =b,Ann(0),Ann(1)
    while r1.degree()>t-1:
        r2=r0%r1
        q2=r0//r1
        u,v=u0,v0
        r0,u0,v0=r1,u1,v1
        r1,u1,v1=r2,u-q2*u1,v-q2*v1
    return (r1/v1(0),u1/v1(0), v1/v1(0))
    
    
def decode_Goppa(L,g,Inv,G,y):
    A=g.parent()
    x=A.gen()
    t=g.degree()//2
    n=len(L)
    S=A(0)
    for i in range(n):
        S+=y[i]*Inv[i]
    if S==0:
        return y
    else:
        w,nimp,sigma= euclide_etendu_partiel(g,S,t,A)
        P=y.list()
        for i in range(n):
            if sigma(L[i])==0:
                sigma=sigma//(x-L[i])
                P[i]=P[i]+GF(2)(1)
                
    message=(G.transpose()).solve_right(vector(P))
                
    return message
    

        
        

