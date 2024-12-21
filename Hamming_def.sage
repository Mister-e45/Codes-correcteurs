import random



def hamming_controle(r):
    n=(2**r)-1
    MS=MatrixSpace(GF(2),n,r)
    L=[]
    for i in srange(1,n+1):
        L.append(i.digits(2,padto=r))
    H=MS.matrix(L)
    return H

def hamming(r):
    return hamming_controle(r).transpose()
    
def poids_H(w):
    W=w.list()
    c=0
    for i in W:
        if i!=GF(2)(0):
            c+=1
    return c

def distance(w1,w2):
    w=w1-w2
    return poids_H(w)
    
    
def perturbe(y, max_error):
    Y=y.list()
    n=len(Y)
    nb_err=random.randint(0,max_error)
    print("nombre d'erreurs: ",nb_err)
    L=[]
    for i in range(nb_err):
        e=random.randint(0,n-1)
        if not (i in L):
            L.append(e)
        else:
            e=random.randint(0,n-1)
            
            while e in L:
                e=random.randint(0,n-1)
                
            L.append(e)
            
    for i in L:
        Y[i]=Y[i]+GF(2)(1)
    
    return MatrixSpace(GF(2),1,n)(Y)
    
    
    

def dist_min(G):
    minimum=G.ncols() #la distance minimale est inferieure a
    C=G.row_space()
    for w in C:
        m=poids_H(w)
        if m!=0:
            if minimum>m:
                minimum=m
    return minimum
    
def decode_hamming(v,H):
    V=[]
    n=H.ncols()
    
    MS=MatrixSpace(GF(2),1,n)
    
    for i in v:
        V.append(GF(2)(i))
        
    V=MS.matrix(V)
    
    S=H*(V.transpose())
    r=H.nrows()
    c=0
    if S==0:
        c=V
    else:
        S=S.list()
        l=len(S)
        s=0
        
        for i in range(r):
            s+=(2**i)*ZZ(S[i])
        
        s=s-1
        
        e=[]
        for i in range(n):
            if i!=s :
                e.append(0)
            else:
                e.append(1)
        
        e=MS.matrix(e)
        
        
        c=V+e
    G=H.right_kernel().basis_matrix()
    L=(G.transpose()).solve_right(c.transpose())
    
    return L
    
    
    
def golay_gen(N):
    n=N/2
    M=matrix(GF(2),n,N,lambda i,j : 0)
    
    for i in range(n): #construction bloc identite
        M[i,i]=GF(2)(1)
        
    for i in range(1,n): # les bords de la matrice A dans M=(Id|A)
        M[i,n]=1
        M[0,i+n]=1
    
    carre=[] # les carres modulo 11
    for i in range(n-1):
        c=(i**2)%(n-1)
        if not (c in carre):
            carre.append(c)
    
    for j in range(0,n-1):
        for i in range(0,n-1):
            if (i+j)%(n-1) in carre:
                M[i+1,j+n+1]=1
                
    return M
    
    
    
        
    
def decode_golay24(r,G24):
    MS=MatrixSpace(GF(2),1,24)
    R=MS(r)
    T=G24*R.transpose()
    x=0
    if poids_H(T)<=3:
        x=T.transpose()
    else:
        A=G24[0:12,12:24]
        
        n=A.ncols()
        found=False
        for i in range(n):
            if poids_H(T+A[0:12,i:i+1])<=2:
                x=(T+A[0:12,i:i+1]).transpose()
                found=True
                break
                
        if not found:
            S=(A.transpose())*T
            if poids_H(S)<=3:
                x=0
            else:
                found=False
                for i in range(n):
                    if poids_H(S+A[0:12,i:i+1])<=2:
                        x=[]
                        for k in range(12):
                            if k==i:
                                x.append(GF(2)(1))
                            else:
                                x.append(GF(2)(0))
                        found=True
                        break
                x= MatrixSpace(GF(2),1,n)(x) #on transforme x en matrice ligne
                
                if not found:
                    return False
    return (R[0,0:12])+x
    
    
    
def decode_exhaustif(G,t,y):
    C=G.row_space()
    C=C.list()
    minimum=G.ncols()
    r=0
    for v in C:
        m=vector(y.list())-v
        if poids_H(m)<minimum:
            minimum=poids_H(m)
            r=v
    return (G.transpose()).solve_right(r)
    
    
    
def motsPoidsm(longueur,m):
    L=[]
    Ps=Subsets([i for i in range(longueur)],m)
    Ps=Ps.list()
    for s in Ps:
        c=s.list()
        a=0
        for i in c:
            a+=2**i
        w=a.digits(2,padto=longueur)
        W=[]
        for b in w:
            W.append(GF(2)(b))
        L.append(W)
        
    return L

def tableau_standard(H,t):
    n=H.ncols()
    s=H.nrows()
    L=[ [[GF(2)(0) for i in range(n)]] , [[GF(2)(0) for i in range(s)]] ]
    for w in range(1,t+1):
        P=motsPoidsm(n,w)
        for c in P:
            v=MatrixSpace(GF(2),n,1)(c)
            syndrome=H*v
            if not (syndrome.list() in L[1]):
                L[0].append(v.list())
                L[1].append(syndrome.list())
                
    return L
    
    
def decodage_sydrome(H,t,y):
    t_stand=tableau_standard(H,t) #deux listes: la premiere contient les mots de poids < à t
    n=len(t_stand)                #la deuxieme contient les syndromes de ces mots
    N=H.ncols()
    Y=MatrixSpace(GF(2),N,1)(y.list())
    syndrome=H*Y
    for i in range(n):
        if syndrome.list()==t_stand[1][i]:
            error=t_stand[0][i]
            error=MatrixSpace(GF(2),N,1)
            return (Y-error).list()
        
    return False #si le syndrome de y n'est pas dans t_stand c'est qu'il y a trop d'erreur dans le message recu y
    
    
def test_hamming(G,H,nb_test):
    VM=MatrixSpace(GF(2),1,11)
    for i in range(1,nb_test+1):
        print("test numero ",i)
        print()
        v=[]
        
        for i in range(11):
            b=random.randint(0,1)
            v.append(GF(2)(b))
        v=VM(v)
        m=v*G
        print("message initial:     ",v.list())
        print("encodage:            ",m.list())
        t=perturbe(m,1) #pour ce code de hamming on peut corriger au plus une erreur
        print("message recu:        ",t.list())
        d=decode_hamming(t.list(),H).list()
        print("correction/decodage: ",d)
        print("message initial egal au message decode?: ",d==v.list())
        print()


def test_Golay(G,nb_test):
    VM=MatrixSpace(GF(2),1,12)
    for i in range(1,nb_test+1):
        print("test numero ",i)
        print()
        v=[]
        
        for i in range(12):
            b=random.randint(0,1)
            v.append(GF(2)(b))
        v=VM(v)
        m=v*G
        print("message initial:     ",v.list())
        print("encodage:            ",m.list())
        t=perturbe(m,3) #pour ce code de hamming on peut corriger au plus 3 erreurs
        print("message recu:        ",t.list())
        d=decode_golay24(t.list(),G).list()
        print("correction/decodage: ",d)
        print("message initial egal au message decode?: ",d==v.list())
        print()
        
def test_decode_syndrome(H,t,nb_test): 
    G=G=H.right_kernel().basis_matrix()
    n=G.nrows()
    VM=MatrixSpace(GF(2),1,n)
    for i in range(1,nb_test+1):
        print("test numero ",i)
        print()
        v=[]
        
        for i in range(n):
            b=random.randint(0,1)
            v.append(GF(2)(b))
        v=VM(v)
        m=v*G
        print("message initial:     ",v.list())
        print("encodage:            ",m.list())
        r=perturbe(m,t) #pour ce code de hamming on peut corriger au plus 3 erreurs
        print("message recu:        ",r.list())
        d=decode_exhaustif(G,t,r).list()
        print("correction/decodage: ",d)
        print("message initial egal au message decode?: ",d==v.list())
        print()

    
    

    

