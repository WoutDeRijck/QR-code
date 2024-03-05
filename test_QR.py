from QR_code import *

QR = QR_code('H', 'optimal')

QR_matrix = QR.generate("KAAS")

print(QR.read(QR_matrix))