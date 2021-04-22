# #!/bin/python 3
# #encoding='UTF8'

import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns

import os
print("PYTHONPATH:", os.environ.get('PYTHONPATH'))
print("PATH:", os.environ.get('PATH'))

# PocetBodu = 70
# random.seed(43)
# theta = [random.uniform(-np.pi/2, np.pi/2) for i in range(PocetBodu)]
# d = 1 
# X = [round(d * np.tan(t),8) for t in theta]

# def cauchy_X(X, d, theta):
# 	return  1/np.pi * np.abs(d/(d**2 + X**2))  

# t = np.linspace(np.tan(-np.pi/2), np.tan(np.pi/2), 350) # 300 points between -30 and 30
# y = np.zeros(len(t))         # allocate y with float elements
# for i in range(len(t)):
#     y[i] = cauchy_X(t[i], d, 0)

# Density Plot and Histogram of all arrival delays
# sns.distplot(X, hist=True, kde=True, 
#              bins=int(180/5), color = 'darkblue', 
#              hist_kws={'edgecolor':'black'},
#              kde_kws={'linewidth': 4})

# # fig = plt.figure(figsize=(5,4))
# ax = fig.add_subplot(111)
# # plt.hist(X)
# plt.plot(t, y)
# # ax.margins(0.3,0.1)
# plt.title("Cauchyho rozdeleni mista dopadu castic")
# plt.show()