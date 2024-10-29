H=hamming(4)
print("matrice de controle pour Hamming")
print(H)

G=H.right_kernel().basis_matrix()
print("matrice generatrice pour Hamming")
print(G)

print("distance minimale pour le code de hamming 4: ",dist_min(G))



G24=golay_gen(24)
print("matrice generatrice pour Golay24")
print(G24)
print("distance minimale pour code de Golay 24: ",dist_min(G24))


H_Golay24=G24.right_kernel().basis_matrix() 

test_hamming(G,H,5)

test_Golay(G24,5)

print("tests du decodage exhaustif pour hamming: ")
print()
test_decode_exhaustif(G,1,5)

print("tests du decodage exhaustif pour Golay: ")
print()
test_decode_exhaustif(G24,1,5)

print("tests du decodage par syndrome pour hamming: ")
print()
test_decode_syndrome(H,1,5)

print("tests du decodage pa syndrome pour Golay: ")
print()
test_decode_syndrome(H_Golay24,3,5)


