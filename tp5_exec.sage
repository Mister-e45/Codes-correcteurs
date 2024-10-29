F.<a>=GF(2**4)
print("polynome minimal de a: ",a.minimal_polynomial())

H=ReedSolomon_control_Matrix(F,15,3,[a**k for k in range(15)],[a**k for k in range(15)])[0]
print(H)

G=ReedSolomon_gen_Matrix(F,15,3,[a**k for k in range(15)],[a**k for k in range(15)])
print(G)


#question 1.3.

G1,V1,D1=ReedSolomon_control_Matrix(F,15,12,[a**k for k in range(15)],[1 for i in range(15)])
print(G1)

#question 1.4.

G1Systematique=G1.echelon_form()

print(G1Systematique==G)

#on a donc que C=C'

m=vector(F,[1,a^3,1])*G1
e=vector(F,[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0])
M=m+e
print(decode_RS_exhaustif(F,15,3,a,M))

m=vector(F,[1,a^2,1])*G1
e=vector(F,[0,0,0,1,0,0,0,1,0,0,0,0,0,0,0])
M=m+e
print(decode_RS_interpolation(F,15,3,a,M))


m=vector(F,[1,a^2,1+a^3])*G1
e=vector(F,[0,0,0,1,0,0,0,1,0,0,0,0,0,0,0])
M=m+e
print(decode_RS_WB(F,15,3,a,M))



corps.<c>=GF(2**10)
print(c.minimal_polynomial())
Gen,V,D=ReedSolomon_control_Matrix(corps,1023,1023-512,[c**k for k in range(1023)],[1 for i in range(1023)])

m=[corps.random_element() for i in range(512)]
vectm=vector(corps,m)*Gen
e=[0 for i in range(1023)]
e[5]=1
e[100]=1
vecte=vector(corps,e)
M=vectm+vecte

res=decode_RS_WB(corps,1023,512,c,M)
print(res==m)
