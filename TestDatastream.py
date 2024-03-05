from QR_code import *

qr = QR_code('M', 'optimal')

qr.read_dataStream(qr.generate_dataStream('TESTDATASTREAM'))