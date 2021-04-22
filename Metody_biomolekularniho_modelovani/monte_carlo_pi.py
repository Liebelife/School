#!/bin/python 3
#encoding='UTF8'

import matplotlib
import numpy as np
import matplotlib.pyplot as plt

import random


def NagenerujSouradnice(PocetBodu):
	kruh = 0
	celkem = 0
	x = [random.uniform(0, 1) for i in range(PocetBodu)]
	y = [random.uniform(0, 1) for i in range(PocetBodu)]

	for a, b in zip(x,y):
		r = a**2 + b**2 	#vzdalenost od stredu (r^2 == 1)
		if r <= 1:
			kruh += 1
		celkem += 1
	return(x,y, kruh, celkem)

# rovnice kruznice se stredem v pocatku: r^2 = x^2 + y^2
# vykresleni na efekt
body_kruznice = [random.uniform(0, 1) for i in range(10000)]
body_kruznice_y = [np.sqrt(1 - x**2) for x in body_kruznice]

chyby = {}
std = {}
for pocet_bodu in [100, 1000, 10000, 100000]:
	hodnoty_odchylek = []
	for opakovani in range(100):
		x, y, kruh, celkem = NagenerujSouradnice(pocet_bodu)
		hodnoty_odchylek.append( np.abs( np.pi - (4 * kruh / celkem) ) )
	chyby[pocet_bodu] = round( np.mean(hodnoty_odchylek), 4)
	std[pocet_bodu] = round( np.mean( np.std(hodnoty_odchylek) ), 4)

print(chyby.values())
print(std.values())

fig = plt.figure()
ax = fig.add_subplot(111)
plt.scatter(chyby.keys(), chyby.values(), color = "black")
plt.plot(chyby.keys(), chyby.values(), color = "black")
plt.errorbar(chyby.keys(), chyby.values(), std.values(), marker='', color = "black")
plt.hlines(y= 10e-3, xmin = 0, xmax = 10e4, color="red")
plt.xscale(value="log")
plt.title("Monte Carlo: prumerna chyba predikce Pi")
plt.xlabel("Pocet bodu")
plt.ylabel("Prumerna chyba")
plt.show()

# plot vizualizace ulohy

# fig = plt.figure()
# ax = fig.add_subplot(111)
# plt.xlim(0,1)
# plt.ylim(0,1)
# ax.set_aspect('equal', adjustable='box')
# plt.scatter(x, y, color = "blue")
# plt.scatter(body_kruznice, body_kruznice_y, marker = ".", color ="black")
# plt.show()
