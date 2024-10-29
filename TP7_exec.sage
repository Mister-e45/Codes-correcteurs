m=10
t=50
n=800
A=GF(2**m)
F.<x>=PolynomialRing(A)



g=get_irreducible(t,A)
gamma=get_distinct_elements(n,A)
print(g)
print()
print(gamma)

Htilde=Goppa_control_Matrix_1(A,g,gamma)
Hchapeau=Goppa_control_Matrix_2(A,g,gamma)

print(Htilde)
print()
print(Hchapeau)


public,L,G=gen_cle(m,n,t)
print(public)
print()
print(L)
print()
print(G)


#on essaye de chiffrer puis dechiffrer un message

indices=[randint(1,799) for i in range(t)]
E=GF(2)^n

message=E(0)

message[1]=GF(2)(1)
message[2]=GF(2)(1)
message[43]=GF(2)(1)
message[207]=GF(2)(1)
cryptogramme=chiffre(message,public)
messageBis=dechiffre(cryptogramme,gamma,g)

print(messageBis-message == E(0))
