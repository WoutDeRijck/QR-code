from QR_code import *

qr = QR_code('Q', 'optimal')

# Define the Galois field
pp = galois.primitive_poly(2,4)
field = galois.GF(2**4, irreducible_poly=pp)
# Define the Reed-Solomon code parameters
n = 15
k = 9

# Generate a random message
message = galois.Poly(np.array([11,3,6,5,12,0,0,4,2]),field=field)
message = field(message.coefficients())
print(message)
# Encode the message
generator = qr.makeGenerator(2,4,3)

rs = galois.ReedSolomon(n, k)
dummy = rs.generator_poly
test_encoded = rs.encode(message)
# print(generator)
encoded = qr.encodeRS(message, 2, 4, n, k,generator)
print(encoded)
error = field(np.array([1,0,3,0,0,0,0,0,0,0,0,0,0,0,9]))
received = encoded+error
print(received)


decoded = qr.decodeRS(received,2,4,n,k,generator)
print(decoded)

