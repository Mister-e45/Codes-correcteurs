def poids_H(w):
    W=w.list()
    c=0
    for i in W:
        if i!=0:
            c+=1
    return c
    
def ReedSolomon_control_Matrix(field,n,k,locs,mults):
    V=zero_matrix(field,n-k,n)
    D=diagonal_matrix(field,vector(field,mults))
    
    for i in range(n-k):
        for j in range(n):
            V[i,j]=locs[j]**i
    return V*D,V,D
    
def ReedSolomon_gen_Matrix(field,n,k,locs,mults):
    H,V,D=ReedSolomon_control_Matrix(field,n,k,locs,mults)
    return H.right_kernel().basis_matrix()
    
def decode_RS_exhaustif(field,n,k,alpha,r):
    
    R.<x>=PolynomialRing(field)
    t=ZZ(int((n-k)/2))
    space=field^k
    for l in space:
        f=R(l.list())
        v=vector(field,[ f(alpha**i) for i in range(n)])
        delta=v-r
        
        if poids_H(delta)<=t:
            return f.list()
        
    return "no closest word found"
    
    
def decode_RS_interpolation(field,n,k,alpha,r):
    t=ZZ(int((n-k)/2))
    R.<x>=PolynomialRing(field)
    lr=r.list()
    parts=Subsets([i for i in range(n)],n-t)
    for I in parts:
        L=[(alpha^i,lr[i]) for i in I]
        f=R.lagrange_polynomial(L)
        delta=vector(field,[f(alpha^i) for i in range(n)])-vector(field,lr)
        if f.degree()<k and poids_H(delta)<=t:
            return f.list()
    return "no closest word found"
    
    
    
def decode_RS_WB(field,n,k,alpha,r):
    R.<x>=PolynomialRing(field)
    t=ZZ(int((n-k)/2))

    M=MatrixSpace(field,n,2*t+k+1)
    lr=r.list()
    mat=M(0)
    for j in range(t+k):
        for i in range(n):
            mat[i,j]=alpha**(i*j)
            
    for j in range(t+1):
        for i in range(n):
            mat[i,j+t+k]=lr[i]*(alpha**(j*i))
            
    q=mat.right_kernel().basis_matrix()[0]
    q=q.list()
    q0=q[:t+k]
    q1=q[t+k:]
    Q0=R(q0)
    Q1=R(q1)
    
    f=-Q0//Q1
    
    
    return f.list()
