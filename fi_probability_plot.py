import matplotlib.pyplot as plt
import numpy as np
from scipy.special import binom
from QR_code import QR_code
from random import random

def binomial(k, p):
    return binom(15, k) * p**k * (1-p)**(15 - k)

def cumul_binomial(p):
    sum = 0
    for i in range(4):
        sum += binomial(i, p)
    return 1 - sum

p = [5e-1, 2e-1, 1e-1, 5e-2, 2e-2, 1e-2]
dec_prob = [cumul_binomial(pr) for pr in p]

def mess_up_and_calc(orig_data : list, prob : float):
    data = orig_data.copy()
                
    errors = 0
    
    # simulate channel for all data bits
    for bit in range(15):
        if random() < prob:
            errors += 1
            data[bit] = (data[bit] + 1) % 2

    if errors < 4:
        pass
    
    result = QR_code.decodeFormat(np.array(data))
    return result

def simulate(inputs : list):
    results_one = []
    results_two = []
    
    for prob in p:
        count_one = 0
        count_two = 0
        
        for input in inputs:
            # count the # decoded words
            orig_data = list(QR_code.encodeFormat(input[0], input[1]))
            
            # do five times per input (thus 5 * len(input) for each prob)
            for _ in range(5):
                
                # simulate channel and decode
                result = mess_up_and_calc(orig_data, prob)

                # test if decoding worked
                actual_result = (result[1] == input[0]) and (result[2] == input[1])
                
                # count failed decodings
                if not actual_result:
                    count_one += 1
                    
                    # second format word
                    result = mess_up_and_calc(orig_data, prob)
                    actual_result = (result[1] == input[0]) and (result[2] == input[1])
                    if not actual_result:
                        count_two += 1

                    
                            
        results_one.append(count_one / (len(inputs) * 5))
        results_two.append(count_two / (len(inputs) * 5))
        
    return results_one, results_two

inputs = []
for level in ["M", "L", "H", "Q"]:
    for mask in [[i,j,k] for i in [0,1] for j in [0,1] for k in [0,1]]:
        inputs.append([level, mask])
        
results_one, results_two = simulate(inputs)


p_1 = p.copy()
p_2 = p.copy()
try: 
    while True:
        i = results_one.index(0)
        results_one = results_one[:i] + results_one[i+1:]
        p_1 = p_1[:i] + p_1[i+1:]
except:
    pass

try: 
    while True:
        i = results_two.index(0)
        results_two = results_two[:i] + results_two[i+1:]
        p_2 = p_2[:i] + p_2[i+1:]
except:
    pass


fig, ax = plt.subplots()

# plot results for one FI word
plt.plot(p, dec_prob, ":", linewidth=3, color="#479e6d", label="1 word calculated")
plt.plot(p_1, results_one, linewidth=3, color="#107d3f", label="1 word simulated")

name = "img/format_one_word.png"

# plot results for two FI words
if bool(str(input("Wanna show two plots?\n     "))[0] == "y"):
    name = "img/format_two_words.png"
    dec_prob_squared = [x**2 for x in dec_prob]
    plt.plot(p, dec_prob_squared, ":", linewidth=3, color="#f5da82", label="2 words calculated")
    plt.plot(p_2, results_two, linewidth=3, color="#edb705", label="2 words simulated")

plt.grid()
# plt.xlim(0, 0.55)
plt.xscale("log")
# plt.ylim(0, 1.10)
plt.yscale("log")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.legend()
plt.title("FI Codeword Decode Error Probability", fontsize=16)
plt.xlabel("bit error probability", fontsize=12)
# plt.show()

plt.savefig(name, dpi=1000)