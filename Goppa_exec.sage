n=25
s=5
t=3
F=GF(2**5)
L=get_distinct_elements(25,F)
g=get_irreducible(3,F)
G=gen_Goppa(F,g**2,L)

print(L)
print()
print(g)
print()
print(G.ncols())
print(G)


R.<z>=PolynomialRing(F)
Invs=[inverse_mod(z-e,g**2) for e in L]

Invs=[inverse_mod(z-e,g**2) for e in L]

v=vector(F,[GF(2).random_element() for i in range(19)])
err=vector(F,[GF(2)(0) for i in range(25)])
err[2]=GF(2)(1)
err[5]=GF(2)(1)
y=v*G+err
dec=decode_Goppa(L,g,Invs,G,y)


