from QR_code import*

QR = QR_code('Q', 'optimal')

p = 5
m = 4
t = 7
k = 5
n = 5+2*t

generator = QR.makeGenerator(p, m, t)
# print(generator)

pp = galois.primitive_poly(p, m)
GF = galois.GF(p**m, irreducible_poly=pp)
message = GF(np.random.randint(0, 100, k))
GF.repr("int")
# print(message)
encoded = QR.encodeRS(message, p, m, n, k, generator)
# print(encoded)
error = np.random.randint(0, 100, n)
error[t:] = np.zeros(n-t, dtype=int)
# print(error)
received = encoded + GF(error)
# print(received)
decoded = QR.decodeRS(received, p, m, n, k, generator)
print(decoded)
print(np.all(decoded==message))

