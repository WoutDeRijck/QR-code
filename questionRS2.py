from QR_code import*


# Parameters
k = 50
N = 10
p = np.array([5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1])
qrgen = QR_code('Q', 'optimal')

# Simulated probability of decoding error
probs = np.zeros(p.shape)
for j, px in enumerate(p):
    print('j:'+str(j),end='\n')
    ber = 0
    for i in range(N):
        print('i:'+str(i),end='\r')
        # Encode 
        data = qrgen.generate_dataStream('KAAS')
        sent = qrgen.encodeData(data)

        # Introduce errors
        error = np.random.rand(sent.shape[0]) < px
        received = np.bitwise_xor(sent, error)

        # Decode
        received = qrgen.decodeData(received)

        # Remove padding
        received = received[:len(data)]

        # Check results
        ber += np.sum(np.bitwise_xor(received, data)) / k
    probs[j] = ber / N
plt.plot(p, probs)

# Plot layout
plt.xlabel("BSC $p$")
plt.ylabel("BER")
plt.show()
