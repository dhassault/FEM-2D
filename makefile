#
# Makefile du TP1 de MEF: le mailleur
#
mailleur.exe:	maillage.o	fonctions.o
	gcc maillage.o fonctions.o -o mailleur.exe

maillage.o:	maillage.c
	gcc -o maillage.o -c maillage.c -W -Wall

fonctions.o:	fonctions.c
	gcc -o fonctions.o -c fonctions.c -W -Wall
