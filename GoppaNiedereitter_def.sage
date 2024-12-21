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
    

def Goppa_control_Matrix_1(corps,g,L):
    V=corps.vector_space(map=False)
    dim=len(V.basis())
    r=g.degree()
    n=len(L)
    M=MatrixSpace(corps,r,n)
    M1=MatrixSpace(corps,r*dim,n)
    D=diagonal_matrix(corps,vector(corps,[1/g(s) for s in L]))
    V=M(0)
    for i in range(r):
        for j in range(n):
            V[i,j]=L[j]**i
    
    return V*D
    
    
    
def Goppa_control_Matrix_2(corps,g,L):
    V=corps.vector_space(map=False)
    dim=len(V.basis())
    r=g.degree()
    n=len(L)
    M=MatrixSpace(corps,r,n)
    M1=MatrixSpace(corps,r*dim,n)
    D=diagonal_matrix(corps,vector(corps,[1/g(s) for s in L]))
    V=M(0)
    for i in range(r):
        for j in range(n):
            V[i,j]=L[j]**i
    
    C=V*D
    C1=[]
    for h in C.columns():
        L=[]
        for i in h:
            L+=corps.vector_space(map=False)(i)
        C1.append(L)
    C1=Matrix(corps,C1)
    C1=C1.transpose()
    return C1
    
def gen_cle(m,n,t):
    field=GF(2**m)
    M=Matrix(field,m*t,n,0)
    l=[]
    g=0
    while M[:,0:m*t] != identity_matrix(field,m*t):
        g=get_irreducible(t,field)
        l=get_distinct_elements(n,field)
        H=Goppa_control_Matrix_2(field,g,l)
        M=H.echelon_form()
    return M[:,m*t:],l,g


def chiffre(e,T):
    l=T.nrows()
    Id=identity_matrix(l)
    H=block_matrix(1,2,[[Id,T]])
    S=H*e
    return S
    
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

def decode_Goppa_Niedereitter(L,g,y):
    A=g.parent()
    x=A.gen()
    t=g.degree()//2
    n=len(L)
    S=A(0)
    Inv=[inverse_mod(x-e,g) for e in L]
    for i in range(n):
        S+=y[i]*Inv[i]
    if S==0:
        print("no error")
        return (A^n)(0)
    else:
        w,nimp,sigma= euclide_etendu_partiel(g,S,t,A)
        P=y.list()
        for i in range(n):
            if sigma(L[i])==0:
                sigma=sigma//(x-L[i])
                P[i]=P[i]+GF(2)(1)
    
    
                
    return P


    
def dechiffre(s,l,g):
    n=len(l)
    k=len(list(s))
    corps=l[0].parent()
    
    c0=vector(corps, list(s)+[0 for i in range(n-k)])
    
    m= decode_Goppa_Niedereitter(l,g**2,c0)
    m=vector(corps,m)
    return c0-m
    
    

    

    
  
    

