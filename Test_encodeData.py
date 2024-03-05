from QR_code import QR_code
import numpy as np

qr = QR_code('Q', 'optimal')

datastream = qr.generate_dataStream('KAAS')

encodedData = qr.encodeData(datastream)

print(encodedData)
print(len(encodedData))

decodedData = qr.decodeData(encodedData)

print(decodedData)
print(len(decodedData))

output = qr.read_dataStream(decodedData)