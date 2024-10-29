

def euclide_etendu_partiel(a,b,t):
    r0,u0,v0=a,1,0
    r1,u1,v1 =b,0,1
    while r1.degree()>t-1:
        r2=r0%r1
        q2=r0//r1
        u,v=u0,v0
        r0,u0,v0=r1,u1,v1
        r1,u1,v1=r2,u-q2*u1,v-q2*v1
    return (r1/v1(0),u1/v1(0), v1/v1(0)) 
    
def decode_BCH(message,g,alpha,t):
    
    S=0
    ordre=alpha.multiplicative_order()
    anneau=alpha.parent()
    A.<z>=PolynomialRing(anneau)
    
    Si=[message(alpha^k) for k in range(1,2*t+1)]
    N=len(Si)
    for i in range(N):
        S+=Si[i]*z^i
    
    
    w,uv,sigma=euclide_etendu_partiel(z^(2*t),S,t)
    racines=sigma.roots(multiplicities=False)
    N=len(racines)
    expo=[]
    for i in range(N):
        expo.append(ordre-racines[i].log(alpha))
    e=0
    for i in expo:
        e+=(message.parent().gen())^i #X^ij
    if (message+e)%g==0:
        return (message+e)//g
    else:
        print("probleme")
        return 0
