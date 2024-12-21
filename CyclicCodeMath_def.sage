def primitiveNroot(q,n):
    F=GF(q)
    o=q-1
    return F.multiplicative_generator()**(o//n)

print(primitiveNroot(3**4,1))


def classe_cyclotomique(q,n,s):
    D=IntegerModRing(n)
    S=D(s)
    Q=D(q)
    L=[S]
    i=1
    while S*(Q**i)!=S:
        L.append(S*(Q**i))
        i+=1
        
    return L
    

def qClasse_cyclotomiques(q,n):
    L=[[0]]
    for i in range(1,n):
        available=True
        for k in L:
            if i in k:
                available=False
        if available:
            L.append(classe_cyclotomique(q,n,i))
    return L
    
def irreducibleComponent(n,q,s):
    p=q
    i=1
    while ( (p-1)%n!=0 ):
        i+=1
        p=p*q
    
    F.<X>=PolynomialRing(GF(q**i))
    alpha=primitiveNroot(q**i,n)
    P=F(1)
    L=qClasse_cyclotomiques(q,n)
    for v in L:
        if s in v:
            for e in v:
                P=P*(X-alpha**e)
            break
        
        
    return P


def factorize(q,n): #gives the irreducible factors of X^n-1 in GF(q)[X]
    fact=set()
    for i in range(n):
        fact.add(irreducibleComponent(n,q,i))
    return list(fact)
