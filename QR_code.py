import math

import galois
import numpy as np
import matplotlib.pyplot as plt

class QR_code:
    def __init__(self, level, mask):
        self.level = level # error correction level, can be 'L', 'M', 'Q', 'H'
        self.mask = mask # the mask pattern, can be 'optimal' or either the three bits representing the mask in a list
        self.version = 6 # the version number, range from 1 to 40, only version number 6 is implemented

        #the generator polynomial of the Reed-Solomon code
        if level=='L':
            self.generator=self.makeGenerator(2,8,9)
        elif level=='M':
            self.generator=self.makeGenerator(2,8,8)
        elif level=='Q':
            self.generator=self.makeGenerator(2,8,12)
        elif level=='H':
            self.generator=self.makeGenerator(2,8,14)
        else:
            Exception('Invalid error correction level!')

        self.NUM_MASKS = 8 #the number of masks


    def encodeData(self, bitstream):
        # first add padding to the bitstream obtained from generate_dataStream()
        # then split this datasequence in blocks and perform RS-coding
        # and apply interleaving on the bytes (see the specification
        # section 8.6 'Constructing the final message codeword sequence')
        # INPUT:
        #  -bitstream: bitstream to be encoded in QR code. 1D numpy array e.g. bitstream=np.array([1,0,0,1,0,1,0,0,...])
        # OUTPUT:
        #  -data_enc: encoded bits after interleaving. Length should be 172*8 (version 6). 1D numpy array e.g. data_enc=np.array([1,0,0,1,0,1,0,0,...])
        assert len(np.shape(bitstream))==1 and type(bitstream) is np.ndarray, 'bitstream must be a 1D numpy array'

        ################################################################################################################
        #insert your code here

        if self.level == 'L':
            num_codewords = 136
            capacity = 1088
            num_blocks = 2
            n = 86
        elif self.level == 'M':
            num_codewords = 108
            capacity = 864
            num_blocks = 4
            n = 43
        elif self.level == 'Q':
            num_codewords = 76
            capacity = 608
            num_blocks = 4
            n = 43
        else: #self.level == 'H'
            num_codewords = 60
            capacity = 480
            num_blocks = 4
            n = 43

        # Calculate the number of codewords in a block
        k = num_codewords // num_blocks

        # Check if bitstream valid
        assert capacity > len(bitstream), 'bitstream/input too long'

        # Add terminator '0000'
        terminator = min(4, capacity - len(bitstream))
        bitstream = np.concatenate((bitstream, np.zeros(terminator, dtype=int)))

        # Add zero's unit len(bitstream) multiple of 8
        padding = 8 - (len(bitstream) % 8)
        bitstream = np.concatenate((bitstream, np.zeros(padding, dtype=int)))

        # Add padding
        iterator = 0
        while len(bitstream) + 8 <= capacity:
            if (iterator % 2) == 0 :
                bitstream = np.concatenate((bitstream, np.array([1,1,1,0,1,1,0,0])))
            else:
                bitstream = np.concatenate((bitstream, np.array([0,0,0,1,0,0,0,1])))
            iterator += 1

        # Split in num_blocks
        blocks = np.split(bitstream, num_blocks)

        # Init Galois Field
        pp = galois.primitive_poly(2, 8)
        GF = galois.GF(2**8, irreducible_poly=pp)

        matrix = []
        for block in blocks:
            codewords = np.split(block, k)
            coef = []
            for codeword in codewords:
                value = int(''.join(map(str, codeword)), 2)
                coef.append(value)
            informationword = GF(coef)
            encoded_block = self.encodeRS(informationword, 2, 8, n, k, self.generator)
            matrix.append(encoded_block)
        matrix = np.array(matrix)
        
        ### Interleaving
        # Split encoded
        data_codewords, error_codewords = np.split(matrix, [k], axis=1)

        # Transpose and flatten to interleave
        data_codewords = data_codewords.transpose().flatten()
        error_codewords = error_codewords.transpose().flatten()

        # Append again
        total = np.concatenate((data_codewords, error_codewords))

        # Back to binary
        data_enc = ''
        for value in total:
            data_enc += format(value, '08b')

        data_enc = np.array(list(data_enc)).astype(int)
        ################################################################################################################

        assert len(np.shape(data_enc))==1 and type(data_enc) is np.ndarray, 'data_enc must be a 1D numpy array'
        return data_enc

    def decodeData(self, data_enc):
        # Function to decode data, this is the inverse function of encodeData
        # INPUT:
        #  -data_enc: encoded binary data with the bytes being interleaved. 1D numpy array e.g. data_enc=np.array([1,0,0,1,0,1,0,0,...])
        #   length is equal to 172*8
        # OUTPUT:
        #  -bitstream: a bitstream with the padding removed. 1D numpy array e.g. bitstream=np.array([1,0,0,1,0,1,0,0,...])
        assert len(np.shape(data_enc))==1 and type(data_enc) is np.ndarray, 'data_enc must be a 1D numpy array'

        ################################################################################################################
        #insert your code here

        # Make sure data_enc uses integers
        data_enc = np.array(data_enc, dtype=int)

        if self.level == 'L':
            num_codewords = 136
            capacity = 1088
            num_blocks = 2
            n = 86
        elif self.level == 'M':
            num_codewords = 108
            capacity = 864
            num_blocks = 4
            n = 43
        elif self.level == 'Q':
            num_codewords = 76
            capacity = 608
            num_blocks = 4
            n = 43
        else: #self.level == 'H'
            num_codewords = 60
            capacity = 480
            num_blocks = 4
            n = 43

        # Calculate the number of codewords in a block
        k = num_codewords // num_blocks

        # Reconstruct the total matrix, integer values
        total = []
        
        # Set data to string
        binaryData = ''.join(map(str, data_enc))

        # Set each 8 bit sequence to integer value
        for i in range(0, len(data_enc), 8):
            group = binaryData[i:i+8]
            value = int(group, 2)
            total.append(value)

        # Reverse interleaving
        data_codewords, error_codewords = np.split(total, [num_blocks*k])

        # Reverse transposition and reshape
        data_codewords = data_codewords.reshape(k, num_blocks).transpose()
        error_codewords = error_codewords.reshape(n-k, num_blocks).transpose()

        # Reconstruct the matrix
        matrix = np.concatenate((data_codewords, error_codewords), axis = 1)

        # Init Galois Field
        pp = galois.primitive_poly(2, 8)
        GF = galois.GF(2**8, irreducible_poly=pp)

        bitstream = ''
        for row in matrix:
            informationword = self.decodeRS(GF(row), 2, 8, n, k, self.generator)
            for coef in informationword:
                # form = '0' + str(k) + 'b'
                codeword = format(int(coef), '08b')
                bitstream += codeword
        
        # If padding was needed, then there is a terminator string to be removed
        padding = False

        # Remove padding
        if (bitstream[-8:] == '11101100'):
            bitstream = bitstream[:-8]
            padding = True
        
        while (bitstream[-16:] == '1110110000010001'):
            bitstream = bitstream[:-16]
            padding = True
            
        # Remove terminator
        if padding:
            bitstream = bitstream[:-4]

        # Remove additional zero's
        # -13 because of the Headers in the data_stream
        while not ((len(bitstream)-13) % 11 == 0 or (len(bitstream)-13) % 11 == 6):
            bitstream = bitstream[:-1]
        
        # Convert the bitstream to a np.array()
        bitstream = np.array(list(bitstream)).astype(int)

        ################################################################################################################

        assert len(np.shape(bitstream))==1 and type(bitstream) is np.ndarray, 'bitstream must be a 1D numpy array'
        return bitstream

    # QR-code generator/reader (do not change)
    def generate(self, data):
        # This function creates and displays a QR code matrix with either the optimal mask or a specific mask (depending on self.mask)
        # INPUT:
        #  -data: data to be encoded in the QR code. In this project a string with only characters from the alphanumeric mode
        #  e.g data="A1 %"
        # OUTPUT:
        #  -QRmatrix: a 41 by 41 numpy array with 0's and 1's
        assert type(data) is str, 'data must be a string'

        bitstream=self.generate_dataStream(data)
        data_bits=self.encodeData(bitstream)

        if self.mask=='optimal':
            # obtain optimal mask if mask=='optimal', otherwise use selected mask
            mask_code=[[int(x) for x in np.binary_repr(i,3)] for i in range(self.NUM_MASKS)]
            score=np.ones(self.NUM_MASKS)
            score[:]=float('inf')
            for m in range(self.NUM_MASKS):
                QRmatrix_m = self.construct(data_bits, mask_code[m], show=False)
                score[m] = self.evaluateMask(QRmatrix_m)
                if score[m]==np.min(score):
                    QRmatrix = QRmatrix_m.copy()
                    self.mask = mask_code[m]

        # create the QR-code using either the selected or the optimal mask
        QRmatrix = self.construct(data_bits, self.mask)

        return QRmatrix

    def construct(self, data, mask, show=True):
        # This function creates a QR code matrix with specified data and
        # mask (this might not be the optimal mask though)
        # INPUT:
        #  -data: the output from encodeData, i.e. encoded bits after interleaving. Length should be 172*8 (version 6).
        #  1D numpy array e.g. data=np.array([1,0,0,1,0,1,0,0,...])
        #  -mask: three bits that represent the mask. 1D list e.g. mask=[1,0,0]
        # OUTPUT:
        #  -QRmatrix: a 41 by 41 numpy array with 0's and 1's
        L = 17+4*self.version
        QRmatrix = np.zeros((L,L),dtype=int)

        PosPattern = np.ones((7,7),dtype=int)
        PosPattern[[1, 5],1:6] = 0
        PosPattern[1:6,[1, 5]] = 0

        QRmatrix[0:7,0:7] = PosPattern
        QRmatrix[-7:,0:7] = PosPattern
        QRmatrix[0:7, -7:] = PosPattern

        AlignPattern = np.ones((5,5),dtype=int)
        AlignPattern[[1,3],1:4] = 0
        AlignPattern[1:4,[1, 3]] = 0

        QRmatrix[32:37, L-9:L-4] = AlignPattern

        L_timing = L-2*8
        TimingPattern = np.zeros((1,L_timing),dtype=int)
        TimingPattern[0,0::2] = np.ones((1, (L_timing+1)//2),dtype=int)

        QRmatrix[6, 8:(L_timing+8)] = TimingPattern
        QRmatrix[8:(L_timing+8), 6] = TimingPattern

        FI = self.encodeFormat(self.level, mask)
        FI = np.flip(FI)

        QRmatrix[0:6, 8] = FI[0:6]
        QRmatrix[7:9, 8] = FI[6:8]
        QRmatrix[8, 7] = FI[8]
        QRmatrix[8, 5::-1]= FI[9:]
        QRmatrix[8, L-1:L-9:-1] = FI[0:8]
        QRmatrix[L-7:L, 8] = FI[8:]
        QRmatrix[L-8, 8] = 1

        nogo = np.zeros((L,L),dtype=int)
        nogo[0:9, 0:9] = np.ones((9,9),dtype=int)
        nogo[L-1:L-9:-1 , 0:9] = np.ones((8,9),dtype=int)
        nogo[0:9, L-1:L-9:-1] = np.ones((9,8),dtype=int)
        nogo[6, 8:(L_timing+8)] = np.ones(( L_timing),dtype=int)
        nogo[8:(L_timing+8), 6] = np.ones((1,L_timing),dtype=int)
        nogo[32:37, L-9:L-4] = np.ones((5,5),dtype=int)
        nogo =np.delete(nogo, 6, 1)
        nogo = nogo[ -1::-1, -1::-1];
        col1 = nogo[:, 0::2].copy()
        col2 = nogo[:, 1::2].copy()
        col1[:, 1::2] = col1[-1::-1, 1::2 ]
        col2[:, 1::2] = col2[-1::-1, 1::2 ]
        nogo_reshape = np.array([col1.flatten(order='F'), col2.flatten(order='F')])
        QR_reshape = np.zeros((2, np.shape(nogo_reshape)[1]),dtype=int)

        ind_col = 0
        ind_row = 0
        ind_data = 0

        for i in range (QR_reshape.size):
            if(nogo_reshape[ind_row, ind_col] == 0):
                QR_reshape[ind_row, ind_col] = data[ind_data]
                ind_data = ind_data + 1
                nogo_reshape[ind_row, ind_col] = 1

            ind_row = ind_row+1
            if ind_row > 1:
                ind_row = 0
                ind_col = ind_col + 1

            if ind_data >= len(data):
                break

        QR_data = np.zeros((L-1, L),dtype=int);
        colr = np.reshape(QR_reshape[0,:], (L, len(QR_reshape[0,:])//L),order='F')
        colr[:, 1::2] = colr[-1::-1, 1::2 ]
        QR_data[0::2, :] = np.transpose(colr)

        coll = np.reshape(QR_reshape[1,:], (L, len(QR_reshape[1,:])//L),order='F')
        coll[:, 1::2] = coll[-1::-1, 1::2]
        QR_data[1::2, :] = np.transpose(coll)

        QR_data = np.transpose(QR_data[-1::-1, -1::-1])
        QR_data = np.hstack((QR_data[:, 0:6], np.zeros((L,1),dtype=int), QR_data[:, 6:]))

        QRmatrix = QRmatrix + QR_data

        QRmatrix[30:33, 0:2] = np.ones((3,2),dtype=int)
        QRmatrix[29, 0] = 1

        nogo = nogo[ -1::-1, -1::-1]
        nogo = np.hstack((nogo[:, 0:6], np.ones((L,1),dtype=int), nogo[:, 6:]))

        QRmatrix = self.applyMask(mask, QRmatrix, nogo)

        if show == True:
            plt.matshow(QRmatrix,cmap='Greys')
            plt.show()

        return QRmatrix

    @staticmethod
    def read(QRmatrix):
        #function to read the encoded data from a QR code
        # INPUT:
        #  -QRmatrix: a 41 by 41 numpy array with 0's and 1's
        # OUTPUT:
        # -data_dec: data to be encoded in the QR code. In this project a string with only characters from the alphanumeric mode
        #  e.g data="A1 %"
        assert np.shape(QRmatrix)==(41,41) and type(QRmatrix) is np.ndarray, 'QRmatrix must be a 41 by numpy array'

        FI = np.zeros((15),dtype=int)
        FI[0:6] = QRmatrix[0:6, 8]
        FI[6:8] = QRmatrix[7:9, 8]
        FI[8] = QRmatrix[8, 7]
        FI[9:] = QRmatrix[8, 5::-1]
        FI = FI[-1::-1]

        L = np.shape(QRmatrix)[0]
        L_timing = L - 2*8

        [success, level, mask] = QR_code.decodeFormat(FI)

        if success:
            qr=QR_code(level, mask)
        else:
            FI = np.zeros((15),dtype=int)
            FI[0:8] = QRmatrix[8, L-1:L-9:-1]
            FI[8:] = QRmatrix[L-7:L, 8]

            [success, level, mask] = QR_code.decodeFormat(FI)
            if success:
                qr=QR_code(level, mask)
            else:
                # print('Format information was not decoded succesfully')
                exit(-1)

        nogo = np.zeros((L,L))
        nogo[0:9, 0:9] = np.ones((9,9),dtype=int)
        nogo[L-1:L-9:-1, 0:9] = np.ones((8,9),dtype=int)
        nogo[0:9, L-1:L-9:-1] = np.ones((9,8),dtype=int)

        nogo[6, 8:(L_timing+8)] = np.ones((1, L_timing),dtype=int)
        nogo[8:(L_timing+8), 6] = np.ones((L_timing),dtype=int)

        nogo[32:37, L-9:L-4] = np.ones((5,5),dtype=int)

        QRmatrix = QR_code.applyMask(mask, QRmatrix, nogo)

        nogo=np.delete(nogo,6,1)
        nogo = nogo[ -1::-1, -1::-1]
        col1 = nogo[:, 0::2]
        col2 = nogo[:, 1::2]
        col1[:, 1::2] = col1[-1::-1, 1::2 ]
        col2[:, 1::2] = col2[-1::-1, 1::2 ]

        nogo_reshape = np.vstack((np.transpose(col1.flatten(order='F')), np.transpose(col2.flatten(order='F'))))

        QRmatrix=np.delete(QRmatrix,6,1)
        QRmatrix = QRmatrix[ -1::-1, -1::-1]
        col1 = QRmatrix[:, 0::2]
        col2 = QRmatrix[:, 1::2]
        col1[:, 1::2] = col1[-1::-1, 1::2]
        col2[:, 1::2] = col2[-1::-1, 1::2]

        QR_reshape = np.vstack((np.transpose(col1.flatten(order='F')), np.transpose(col2.flatten(order='F'))))

        data = np.zeros((172*8, 1))
        ind_col = 0
        ind_row = 0
        ind_data =0
        for i in range(QR_reshape.size):
            if(nogo_reshape[ind_row, ind_col] == 0):
                data[ind_data] = QR_reshape[ind_row, ind_col]
                ind_data = ind_data + 1
                nogo_reshape[ind_row, ind_col] = 1

            ind_row = ind_row+1
            if ind_row > 1:
                ind_row = 0
                ind_col = ind_col + 1

            if ind_data >= len(data):
                break

        bitstream = qr.decodeData(data.flatten())
        data_dec = QR_code.read_dataStream(bitstream)

        assert type(data_dec) is str, 'data_dec must be a string'
        return  data_dec

    @staticmethod
    def generate_dataStream(data):
        # this function creates a bitstream from the user data.
        # ONLY IMPLEMENT ALPHANUMERIC MODE !!!!!!
        # INPUT:
        #  -data: the data string (for example 'ABC012')
        # OUTPUT:
        #  -bitstream: a 1D numpy array containing the bits that
        #  represent the input data, headers should be added, no padding must be added here.
        #  Add padding in EncodeData. e.g. bitstream=np.array([1,0,1,1,0,1,0,0,...])
        assert type(data) is str, 'data must be a string'

        ################################################################################################################
        #insert your code here

        # Alphanumeric mode
        encodingTable = {'0': 0, '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, 'A': 10, 'B': 11, 'C': 12, 'D': 13, 'E': 14, 'F': 15, 'G': 16, 'H': 17, 'I': 18, 'J': 19, 'K': 20, 'L': 21, 'M': 22, 'N': 23, 'O': 24, 'P': 25, 'Q': 26, 'R': 27, 'S': 28, 'T': 29, 'U': 30, 'V': 31, 'W': 32, 'X': 33, 'Y': 34, 'Z': 35, ' ': 36, '$': 37, '%': 38, '*': 39, '+': 40, '-': 41, '.': 42, '/': 43, ':': 44}
        
        # Character Count Indicator (9-bits): binary representation of the data length
        characterCountIndicator = np.array(list(format(len(data), '09b'))).astype(int)

        # Mode Indicator
        modeIndicator = np.array([0,0,1,0])

        # Encode data
        binaryData = ''
        for i in range(0, len(data), 2):
            # the character value of the final character is encoded to a 6-bit binary number
            if i == len(data) - 1:
                binaryGroup = format(encodingTable[data[i]], '06b')
            # convert each group to its 11-bit binary equivalent
            else:
                # The character value of the first character is multiplied by 45 and the character value of the second digit is added to the product. 
                sum = encodingTable[data[i]] * 45 + encodingTable[data[i+1]]
                # The sum is then converted to an 11 bit binary number.
                binaryGroup = format(sum, '011b')

            binaryData += binaryGroup

        dataStream = np.array(list(binaryData)).astype(int)

        # Add Mode Indicator and Character Count Indicator to binary data
        bitstream = np.concatenate((modeIndicator, characterCountIndicator, dataStream))
        print(bitstream)

        ################################################################################################################

        assert len(np.shape(bitstream))==1 and type(bitstream) is np.ndarray, 'bitstream must be a 1D numpy array'
        return bitstream

    @staticmethod
    def read_dataStream(bitstream):
        # inverse function of generate_dataStream: convert a bitstream to an alphanumeric string
        # INPUT:
        #  -bitstream: a 1D numpy array of bits (including the header bits) e.g. bitstream=np.array([1,0,1,1,0,1,0,0,...])
        # OUTPUT:
        #  -data: the encoded data string (for example 'ABC012')
        assert len(np.shape(bitstream))==1 and type(bitstream) is np.ndarray, 'bitstream must be a 1D numpy array'

        ################################################################################################################
        #insert your code here

        # Alphanumeric mode
        decodingTable = {0: '0', 1: '1', 2: '2', 3: '3', 4: '4', 5: '5', 6: '6', 7: '7', 8: '8', 9: '9', 10: 'A', 11: 'B', 12: 'C', 13: 'D', 14: 'E', 15: 'F', 16: 'G', 17: 'H', 18: 'I', 19: 'J', 20: 'K', 21: 'L', 22: 'M', 23: 'N', 24: 'O', 25: 'P', 26: 'Q', 27: 'R', 28: 'S', 29: 'T', 30: 'U', 31: 'V', 32: 'W', 33: 'X', 34: 'Y', 35: 'Z', 36: ' ', 37: '$', 38: '%', 39: '*', 40: '+', 41: '-', 42: '.', 43: '/', 44: ':'}

        # Mode Indicator
        modeIndicator = bitstream[:4]

        # Character Count Indicator (9-bits)
        characterCountIndicator = bitstream[4:13]
        dataLength = int(''.join(map(str, characterCountIndicator)), 2)

        # Extract encoded data
        encodedData = bitstream[13:(13+dataLength*6)]
        binaryData = ''.join(map(str, encodedData))

        # Decode data
        data = ''
        for i in range(0, len(binaryData), 11):
            # the final group, containing only 6 bits, is converted to its decimal value and represents the value of the final character
            if i + 6 >= len(binaryData):
                data += decodingTable[int(binaryData[i:i+6], 2)]
            # convert each group of 11 bits to its two corresponding characters.
            else:
                group = binaryData[i:i+11]
                value = int(group, 2)
                value1 = value // 45
                value2 = value % 45
                # value1 = int(group[:5], 2)
                # value2 = int(group[5:], 2)
                data += decodingTable[value1] + decodingTable[value2]
        print(data)

        ################################################################################################################

        assert type(data) is str, 'data must be a string'
        return data

    @staticmethod
    def encodeFormat(level, mask):
        # Encodes the 5 bit format to a 15 bit sequence using a BCH code
        # INPUT:
        #  -level: specified level 'L','M','Q' or'H'
        #  -mask: three bits that represent the mask. 1D list e.g. mask=[1,0,0]
        # OUTPUT:
        # format: 1D numpy array with the FI-codeword, with the special FI-mask applied (see specification)
        assert type(mask) is list and len(mask)==3, 'mask must be a list of length 3'

        ################################################################################################################
        #insert your code here

        # get the data bits by concatenating the level bits with the 3 mask bits
        level_bits = {"L": [0, 1], "M": [0, 0], "Q": [1, 1], "H": [1, 0]}
        data_bits = level_bits[level] + mask
        
        # create a polynomial that represents the data bits
        GF2 = galois.GF(2)
        data_poly = galois.Poly(data_bits + [0]*10, GF2) # already shifted 10 to the left
        
        # create the generator polynomial
        generator = galois.Poly([1,0,1,0,0,1,1,0,1,1,1], GF2)
        
        remainder = data_poly % generator
        
        # concatenate the original data bits with the BCH bits
        concat_bits = data_bits + list(remainder.coefficients(size=10))
        
        # XOR these with the XOR mask
        XOR_mask = [1,0,1,0,1,0,0,0,0,0,1,0,0,1,0]
        format = np.logical_xor(concat_bits, XOR_mask).astype(int)
        
        ################################################################################################################

        assert len(np.shape(format))==1 and type(format) is np.ndarray and format.size==15, 'format must be a 1D numpy array of length 15'
        return format

    @staticmethod
    def decodeFormat(Format):
        # Decode the format information
        # INPUT:
        # -format: 1D numpy array (15bits) with format information (with FI-mask applied)
        # OUTPUT:
        # -success: True if decodation succeeded, False if decodation failed
        # -level: being an element of {'L','M','Q','H'}
        # -mask: three bits that represent the mask. 1D list e.g. mask=[1,0,0]
        assert len(np.shape(Format))==1 and type(Format) is np.ndarray and Format.size==15, 'format must be a 1D numpy array of length 15'

        ################################################################################################################
        #insert your code here

        # define the galois field we will be working in (q=2, m=4)
        GF16 = galois.GF(2**4)
        # get primitive element, alpha
        alpha = GF16.primitive_element

        # apply mask again to get the FI data bits, get poly representation
        XOR_mask = [1,0,1,0,1,0,0,0,0,0,1,0,0,1,0]
        data = np.logical_xor(Format, XOR_mask).astype(int)
        data_poly = galois.Poly(data, GF16)
        data = list(data)
        
        # get syndrome numbers
        s1 = data_poly(alpha)
        s2 = data_poly(alpha ** 2)
        s3 = data_poly(alpha ** 3)
        s4 = data_poly(alpha ** 4)
        s5 = data_poly(alpha ** 5)
        s6 = data_poly(alpha ** 6)
        
        # starting matrix for PGZ algorithm
        A = GF16([[s1, s2, s3], [s2, s3, s4], [s3, s4, s5]])
        b = GF16([-s4, -s5, -s6])
        
        # initial rank estimate is three
        rank = 3
        while np.linalg.matrix_rank(A) != rank:
            # there is at least one less error in the data
            rank -= 1
            
            if rank == 0:
                break
            
            # make subselection for new rank calculation
            A = A[1:,1:]
            b = b[1:]
        
        if rank != 0:
            if rank != np.linalg.matrix_rank(A):
                # not decodable...
                return False, "M", [0,0,0]
            
            # solve set of equations
            sigmas = np.linalg.solve(A, b)
            
            lambdas = list(sigmas) + [1]  
            error_correcting = galois.Poly(lambdas, GF16)
        
            # now locate the errors
            roots = error_correcting.roots()

            # if not enough roots, we cannot correct all errors, return False
            if len(roots) < rank:
                return False, "M", [0,0,0]
            
            # calculate error locations based on roots
            errors = []
            for j in range(15):
                if alpha ** (-j) in roots:
                    errors += [14 - j]
            
            # correct the errors
            for error in errors:
                data[error] = (data[error] + 1) % 2
        
        
        # data_poly = galois.Poly(data, GF16)
        # remainder = data_poly % generator
        
        # if we are here, no errors are left (not present, or corrected already)
        success = True
        
        # bits 2 through 4 are the mask
        mask = data[2:5]
        
        # map the first two bits to the level character
        level_dict = {0: "M", 1: "L", 3 : "Q", 2 : "H"}
        level = level_dict[data[0]*2 + data[1]]
        
        ################################################################################################################

        assert type(mask) is list and len(mask)==3, 'mask must be a list of length 3'
        return success, level, mask


    @staticmethod
    def makeGenerator(p, m, t):
        # Generate the Reed-Solomon generator polynomial with error correcting capability t over GF(p^m)
        # INPUT:
        #  -p: field characteristic, prime number
        #  -m: positive integer
        #  -t: error correction capability of the Reed-Solomon code, positive integer > 1
        # OUTPUT:
        #  -generator: galois.Poly object representing the generator polynomial

        ################################################################################################################
        #insert your code here
        
        ### Creation of Galois field
        pp = galois.primitive_poly(p, m)
        GF = galois.GF(p**m, irreducible_poly=pp)
        GF.repr("power")
        
        ### extract primitive element
        prim = GF.primitive_element
        
        ### initialize generator with m0 = 0
        generator = galois.Poly([1], field=GF)
        
        ### generate generator
        for i in range(2*t):
            temp = galois.Poly([1, -prim**(i)], field=GF)
            generator = generator*temp
            
        ### Testing
        # print(f"generator for p = {p}, m = {m} and t = {t}: {generator}")

        ################################################################################################################

        assert type(generator) == type(galois.Poly([0], field=galois.GF(m))), 'generator must be a galois.Poly object'
        return generator

    @staticmethod
    def encodeRS(informationword, p, m, n, k, generator):
        # Encode the informationword
        # INPUT:
        #  -informationword: a 1D array of galois.GF elements that represents the information word coefficients in GF(p^m) (first element is the highest degree coefficient)
        #  -p: field characteristic, prime number
        #  -m: positive integer
        #  -n: codeword length, <= p^m-1
        #  -k: information word length
        #  -generator: galois.Poly object representing the generator polynomial
        # OUTPUT:
        #  -codeword: a 1D array of galois.GF elements that represents the codeword coefficients in GF(p^m) corresponding to systematic Reed-Solomon coding of the corresponding information word (first element is the highest degree coefficient)
        prim_poly = galois.primitive_poly(p,m)
        GF = galois.GF(p**m, irreducible_poly=prim_poly)
        assert type(informationword) is GF and len(np.shape(informationword))==1, 'each element of informationword(1D)  must be a galois.GF element'
        assert type(generator) == type(galois.Poly([0], field=galois.GF(m))), 'generator must be a galois.Poly object'

        ################################################################################################################
        #insert your code here
        
        ### turn code word into poly
        x = galois.Poly(GF([1, 0]))
        code_poly = galois.Poly(informationword, field=GF)*x**(n-k)
        
        ### remainder after division with generator
        rest_poly = code_poly % generator
        
        ### create codeword
        codeword = code_poly - rest_poly
        
        ### type casting
        codeword = GF(codeword.coefficients())
        
        ################################################################################################################

        assert type(codeword) is GF and len(np.shape(codeword)) == 1, 'each element of codeword(1D)  must be a galois.GF element'
        return codeword

    @staticmethod
    def decodeRS(codeword, p, m, n, k, generator):
        # Decode the codeword
        # INPUT:
        #  -codeword: a 1D array of galois.GF elements that represents the codeword coefficients in GF(p^m) corresponding to systematic Reed-Solomon coding of the corresponding information word (first element is the highest degree coefficient)
        #  -p: field characteristic, prime number
        #  -m: positive integer
        #  -n: codeword length, <= p^m-1
        #  -k: decoded word length
        #  -generator: galois.Poly object representing the generator polynomial
        # OUTPUT:
        #  -decoded: a 1D array of galois.GF elements that represents the decoded information word coefficients in GF(p^m) (first element is the highest degree coefficient)
        prim_poly = galois.primitive_poly(p,m)
        GF = galois.GF(p**m, irreducible_poly=prim_poly)
        assert type(codeword) is GF and len(np.shape(codeword))==1, 'each element of codeword(1D)  must be a galois.GF element'
        assert type(generator) == type(galois.Poly([0],field=galois.GF(m))), 'generator must be a galois.Poly object'

        ################################################################################################################
        #insert your code here
        t = int(np.floor(n-k)/2)
        alpha = GF.primitive_element
        codeword_poly = galois.Poly(codeword, field=GF)

        ##### PGZ
        ### calculating the syndrome
        S = GF(np.array([codeword_poly(alpha**j) for j in range(2*t)]))
        
        ### check if there are errors
        if np.all(S == GF.Zeros(2*t)):
            decoded = codeword[:k]
        else:
            ### set numbers of errors to max
            curr_nu = t
            
            ### construct S and M
            S_mat = np.negative(S[t:2*t])
            M = GF(np.array([S[i:i+t] for i in range(t)]))
            
            ### loop of the PGZ algorithm
            while(np.linalg.det(M) == 0):
                # remove first row and collumn
                M = M[1:,1:]
                # remove first row
                S_mat = S_mat[1:]
                # decrease number of errors
                curr_nu -= 1
            
            ### computation of lambda matrix
            L = np.dot(np.linalg.inv(M), S_mat)
            
            ### create error locator polynomial
            L_poly = galois.Poly(GF(1))+galois.Poly(L)*galois.Poly(GF([1,0]))
            
            ### extract location from roots with Shien method
            positions = []
            for j in range(n):
                if L_poly(alpha**(j+p**m-n)) == GF(0): # shortened code fix 
                    positions.append(j)      
            
            ### extract error evaluator polynomial
            div = np.zeros(2*t+1, int)
            div[0] = 1
            omega_poly = (galois.Poly(np.flip(S), field=GF)*L_poly)%(galois.Poly(div, field=GF))
            
            ### compute dirivative of error locator polynomial
            L_deriv = L_poly.derivative()
            
            ### Forney algorithm
            error_values = GF(np.array([GF((alpha**(n-1-pos)))*GF(omega_poly(alpha**(-(n-1-pos))))/GF(L_deriv(alpha**(-(n-1-pos)))) for pos in positions]).astype(np.int16)) # again fix for shortened codes
            ### error vector formatting
            e = GF.Zeros(n)
            for pos in enumerate(positions):
                e[pos[1]] = error_values[pos[0]]
            
            ### substract error values
            c = codeword + e
            
            ### systematic code
            decoded = c[:k]
    
        
        

        ################################################################################################################

        assert type(decoded) is GF and len(np.shape(decoded))==1, 'each element of decoded(1D)  must be a galois.GF element'
        return decoded

    # function to mask or unmask a QR_code matrix and to evaluate the masked QR symbol (do not change)
    @staticmethod
    def applyMask(mask, QRmatrix, nogo):
        #define all the masking functions
        maskfun1=lambda i, j : (i+j)%2==0
        maskfun2=lambda i, j : (i)%2==0
        maskfun3=lambda i, j : (j)%3==0
        maskfun4=lambda i, j : (i+j)%3==0
        maskfun5=lambda i, j : (math.floor(i/2)+math.floor(j/3))%2==0
        maskfun6=lambda i, j : (i*j)%2 + (i*j)%3==0
        maskfun7=lambda i, j : ((i*j)%2 + (i*j)%3)%2==0
        maskfun8=lambda i, j : ((i+j)%2 + (i*j)%3)%2==0

        maskfun=[maskfun1,maskfun2,maskfun3,maskfun4,maskfun5,maskfun6,maskfun7,maskfun8]

        L = len(QRmatrix)
        QRmatrix_masked = QRmatrix.copy()

        mask_number=int(''.join(str(el) for el in mask),2)

        maskfunction=maskfun[mask_number]

        for i in range(L):
            for j in range(L):
                if nogo[i,j]==0:
                    QRmatrix_masked[i,j] = (QRmatrix[i,j] + maskfunction(i,j))%2

        return QRmatrix_masked

    @staticmethod
    def evaluateMask(QRmatrix):
        Ni = [3, 3, 40, 10]
        L = len(QRmatrix)
        score = 0
        QRmatrix_temp=np.vstack((QRmatrix, 2*np.ones((1,L)),np.transpose(QRmatrix), 2*np.ones((1,L))  ))

        vector=QRmatrix_temp.flatten(order='F')
        splt=QR_code.SplitVec(vector)

        neighbours=np.array([len(x) for x in splt])
        temp=neighbours>5
        if (temp).any():
            score+= sum([x-5+Ni[0] for x in neighbours if x>5])

        QRmatrix_tmp = QRmatrix
        rec_sizes = np.array([[5, 2, 4, 4, 3, 4, 2, 3, 2, 3, 2], [2, 5, 4, 3, 4, 2, 4, 3, 3, 2, 2]])

        for i in range(np.shape(rec_sizes)[1]):

            QRmatrix_tmp, num=QR_code.find_rect(QRmatrix_tmp, rec_sizes[0, i], rec_sizes[1,i])
            score +=num*(rec_sizes[0, i]-1)*(rec_sizes[1,i]-1)*Ni[1]

        QRmatrix_tmp = np.vstack((QRmatrix, 2*np.ones((1, L)), np.transpose(QRmatrix), 2*np.ones((1, L))))
        temp=QRmatrix_tmp.flatten(order='F')
        temp2=[x for x in range(len(temp)-6) if (temp[x:x+7]==[1, 0, 1, 1, 1, 0, 1]).all()]
        score += Ni[2]*len(temp2)

        nDark = sum(sum(QRmatrix==1))/L**2
        k = math.floor(abs(nDark-0.5)/0.05)
        score += Ni[3]*k

        return score

    @staticmethod
    def SplitVec(vector):
        output=[]
        temp=np.where(np.diff(vector)!=0)[0]
        temp=temp+1
        temp=np.insert(temp,0,0)

        for i in range(len(temp)):
            if i==len(temp)-1:
                output.append(vector[temp[i]:])
            else:
                output.append(vector[temp[i]:temp[i+1]])

        return output

    @staticmethod
    def find_rect(A,nR,nC):

        Lx = np.shape(A)[0]
        Ly = np.shape(A)[1]
        num = 0
        A_new = A.copy()

        for x in range(Lx-nR+1):
            for y in range(Ly-nC+1):
                test = np.unique(A_new[x:x+nR, y:y+nC])

                if len(test) == 1:
                    num += 1
                    A_new[x:x+nR, y:y+nC] = np.reshape(np.arange(2+x*nR+y,2+nR*nC+x*nR+y),(nR, nC))

        return A_new, num
