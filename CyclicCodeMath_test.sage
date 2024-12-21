#x^9-1 dans F2[x]
a=factorize(2,9)
print(a)

#x^15-1 dans F2[x]
a=factorize(2,15)
print(a)

#x^21-1 dans F2[x]
a=factorize(2,21)
print(a)

#x^23-1 dans F2[x]
a=factorize(2,23)
print(a)

#x^26-1 dans F3[x]
a=factorize(3,26)
print(a)


#partie 2

cl=qClasse_cyclotomiques(2,1023)
g=set()
for i in range(1,101):
    for p in cl:
        if i in p:
            g.add(irreducibleComponent(1023,2,i))
            break
P=1
L=list(g)
P=lcm(L)
print(P)

l=P.list()
M=MatrixSpace(GF(2),1023-450,1023)
G=M(0)
for i in range(1023-450):
    for j in range(450):
        G[i,i+j]=l[j]

print(G)
