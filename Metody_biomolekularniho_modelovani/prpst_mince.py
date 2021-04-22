#!/bin/python 3
#encoding='UTF8'

import math
#log() je prirozeny, log10() je dekadicky

N=10000
k= 5067

p = 0.5067
p_k = math.exp( N * math.log(N) + k * math.log(p) + (N-k) * math.log((1-p)) - k * math.log(k) - (N - k) * math.log(N - k) ) 
print(p_k)
print(1-p_k)

p = 0.5
p_k = math.exp( N * math.log(N) + k * math.log(p) + (N-k) * math.log((1-p)) - k * math.log(k) - (N - k) * math.log(N - k) ) 
print(p_k)