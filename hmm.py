#!/bin/python 3
#encoding: utf-8


#### http://www.vyuka.ookami.cz/materialy/bioinfo/motifs/hmm.xml
import math
from pprint import pprint

zarovnani = \
"""\
A C A – – – A T G
T C A A C T A T C
A C A C – – A G C
A G A – – – A T C
A C C G – – A T C\
"""
zarovnani = zarovnani.split('\n')

# prpst vyskytu baze na pozici 1 az 7 (poradi v radku)
EmisniMatice = {
					'A':[0.8, 0.0, 0.8, 0.2, 1.0, 0.0, 0.0],
					'T':[0.2, 0.0, 0.0, 0.2, 0.0, 0.8, 0.0],
					'G':[0.0, 0.2, 0.0, 0.2, 0.0, 0.2, 0.2],
					'C':[0.0, 0.8, 0.2, 0.4, 0.0, 0.0, 0.8]
				}

# prpst prechodu mezi stavy
Transmise = {	
				'1>2': 1,
				'2>3': 1,
				'3>4': 0.6,
				'3>5': 0.4,
				'4>4': 0.4,
				'4>5': 0.6,
				'5>6': 1,
				'6>7': 1,
			}
# navstivene stavy modelu pri urcite delce sekvence a provedene transmise
MoznePrechody = {
					"6": [[1, 2, 3, 5, 6, 7], ['1>2', '2>3', '3>5', '5>6', '6>7']],
					"7": [[1, 2, 3, 4, 5, 6, 7], ['1>2', '2>3', '3>4', '4>5', '5>6', '6>7']],
					"8": [[1, 2, 3, 4, 4, 5, 6, 7], ['1>2', '2>3', '3>4', '4>4', '4>5', '5>6', '6>7']],
					"9": [[1, 2, 3, 4, 4, 4, 5, 6, 7], ['1>2', '2>3', '3>4', '4>4', '4>4', '4>5', '5>6', '6>7']]
				}
					
"""
UKOL 1

Spocitat prpsti a log-odds vsech zadanych sekvenci (zde "zarovnani")
"""
def GenerujPrechodyProSekvenci(sekvence):
	Sekvence = "".join(sekvence.split())
	DelkaSekvence = 0
	for znak in Sekvence:
		if znak != '–':
			DelkaSekvence += 1
	NavstiveneStavy = MoznePrechody[str(DelkaSekvence)][0]
	Prechody = MoznePrechody[str(DelkaSekvence)][1]
	return [Sekvence, NavstiveneStavy, Prechody]

def SpocitejPravdepodobnost(sekvence, EmisniMatice, Transmise):
	sekvence, stavy, prechody = GenerujPrechodyProSekvenci(sekvence)
	znaky = []
	for znak in sekvence:
		if znak != '–':
			znaky.append(znak)
	ZnakStav = [dvojice for dvojice in zip(znaky, stavy)]

	prpsti = [EmisniMatice[i[0]][i[1]-1] for i in ZnakStav]
	prchd = [Transmise[i] for i in prechody]
	prpst_prchd = prpsti + prchd

	# slouceni pravdepodobnosti znaku ve stavech a pravdepodobnosti prechodu do jednoho vektoru násobků
	Prpst = 1
	for i in prpst_prchd:
		Prpst *= i
	return(round(Prpst,4))

def SpocitejLogOdd(pravdepodobnost, sekvence):
	PocetZnaku = 0
	for i in sekvence:
		if i != '–':
			PocetZnaku += 1
	Pozadi = 0.25**PocetZnaku
	logodd = math.log(pravdepodobnost/Pozadi)
	return round(logodd, 2)

## Řešení úkolu číslo 1
Vysledky = []
for sekvence in zarovnani:
	Seq = "".join(sekvence.split())
	Pst = SpocitejPravdepodobnost(sekvence, EmisniMatice, Transmise)
	Logodd = SpocitejLogOdd(Pst, Seq)
	Vysledky.append([Seq, Pst, Logodd])
pprint(Vysledky)

"""
UKOL 2

Se stejnym modelem prohledej zadanou sekvenci. Kde v textu se nejpravdepodobneji (a s jakou prpsti) hledany vzor nachazi?
"""

ProhledavanaSekvence = 'AGATCCATTGACCGTTACACATCAGATTGATAGATTGATTTTGATCGACAAAGTG'

# převod Emisní matice na pseudocounts přičtením jedničky ke všem hodnotám matice
EmisniMaticePseudoCounts = {}
for key, value in EmisniMatice.items():
	EmisniMaticePseudoCounts[key] = [i + 1 for i in value]
# pprint(EmisniMaticePseudoCounts)

# vsechny moznosti klouzaveho okenka
## v delsich nez 6 se muze vyskytnout insert => od 4. do 6. pozice
"""z toho plyne, ze kdyz je sekvence dlouha 7, ma insert na 4. pozici, kdyz 8, ma insert na 4 a 5. pozici (tak byl vytvoren model)"""
moznosti = {}
for delka in range(6,10): # 6, 7, 8, 9 
	moznost = []
	for i in range(len(ProhledavanaSekvence)-delka+1):
		usek = ProhledavanaSekvence[i:i+delka]
		moznost.append(usek)
	moznosti[delka] = moznost
# pprint(moznosti)

## prpsti modelu
PrpstiModelu = {}
for delka in moznosti:
	hodnoty = []
	for sekvence in moznosti[delka]:
		Prpst = SpocitejPravdepodobnost(sekvence, EmisniMaticePseudoCounts, Transmise)
		hodnoty.append(Prpst)
	PrpstiModelu[delka] = hodnoty

# nejpravdepodobnejsi mezi delkami patternu
MaxPrpsti = []
for delka in PrpstiModelu:
	MaxPrpsti.append([delka, max(PrpstiModelu[delka])])

# nejpravdepodobnejsi vychazi pattern s delkou 7 a prpsti 19.0468:
# sort sestupne podle druheho elementu listu (pravdepodobnosti)
MaxPrpsti.sort(key = lambda x: x[1], reverse = True) 
NejvyssiPrpstPatternu = MaxPrpsti[0]
# print(NejvyssiPrpstPatternu)

# nalezneme na ktere pozici vzor zacina:
DelkaPatternu, HledanaPrpst = NejvyssiPrpstPatternu
PocatekPatternu = 0
KontrolaDuplikatu = 0
PoziceNalezena = False
for pozice in PrpstiModelu[DelkaPatternu]:
	if PoziceNalezena == False:
		PocatekPatternu += 1
	if pozice == HledanaPrpst:
		KontrolaDuplikatu += 1
		PoziceNalezena = True
	
## vypis vysledku
if KontrolaDuplikatu == 1:
	print("Pattern zacina na", PocatekPatternu, "pismenku, je dlouhy", DelkaPatternu,\
		 	"znaku. Pravdepodobnost jeho vyskytu je", HledanaPrpst)
else:	
	print("Bylo nalezeno vice pozic se stejnou pravdepodobnosti, pattern se vyskytuje asi na vice mistech.")

## prpsti pozadi
PrpstPozadi = {}
for delka in range(6,10):
	pozadi = round(0.25**delka, 6)
	PrpstPozadi[delka] = pozadi
PozadiPatternu = PrpstPozadi[DelkaPatternu]

print("Pravdepodobnost pozadi pro pattern delky", DelkaPatternu, "je", PozadiPatternu)

## vizualizace vysledku
print(ProhledavanaSekvence)
DelkaSekvence = len(ProhledavanaSekvence)

print( (PocatekPatternu - 1) * "b" + DelkaPatternu * "m" + (DelkaSekvence - (PocatekPatternu - 1 + DelkaPatternu) ) * "b")

print("Pro zitrek bych zvolila cheesecake")
