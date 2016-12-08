
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


float a00(float* x);
float a11(float* x);
float a12(float* x);
float a21(float* x);
float a22(float* x);
void adwdw(int nbneel, float** dfcbas, float** ijacob, float eltdif, float** cofvar, float** matelm);
float **alloctab(int dim1, int dim2);
int **alloctabint(int dim1, int dim2);
void assemb(int nutyge, int nbtng, int nbtel, int nbneel, int nbaret, float** coord, int** ngnel, int** nRefAr, int nrefdom, int nbrefD0, int* nurefD0, int nbrefD1, int* nurefD1, int nbrefF1, int* nurefF1, FILE* ficsmd);

float bn(float* x);

void calelb(int nutygear, int nbnear, float** cooreb, float** matela, float* scmela);
void calelm(int nutyge, int nbneel, float** coorel, float** matelm, float *scmelm);
void coel2b(int nbnear, int* nunear, float** coorel, float** cooreb);
void coma2el(int k, int nbneel, float** coord, int** ngnel, float** coorel);

void dmdamo(FILE* ficsmd, FILE* ficsmo);

void elemkl(int nrefdom, int nbrefD0, int* nurefD0, int nbrefD1, int* nurefD1, int nbrefF1, int* nurefF1, int nutyge, int nbneel, float** coorel, int nbaret, int* nrefak,float** matelm, float* scmelm, int* nudlel, float* udel);

float fomega(float* x);
float fn(float* x);
void freetab(void *ptr);

void lecfima(FILE* ficmai, int* nutyge, int* nbtng, float*** pcoord, int* nbtel, int* nbneel, int*** pngnel, int* nbaret, int*** pnRefAr);
void lecnuref(FILE* firef, int* nrefdom, int* nbrefD0, int* nurefD0, int* nbrefD1, int* nurefD1, int* nbrefF1, int* nurefF1);
void lecsmd(FILE* ficsmd);
void lecsmo(FILE* ficsmo);

int max(int a, int b);
int min(int a, int b);

float ud(float* x);

void w(int nbneel, float *fctbas, float eltdif, float cofvar, float *scmelm);
void ww(int nbneel, float *fctbas, float eltdif, float cofvar, float **matelm);
