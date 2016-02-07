#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void etiqar(int typel, int n1, int n2, int nrefdom, const int *nrefcot, int nbtel, int nbaret, int **nRefAr);
int **alloctabINT(int dim1, int dim2);
float **alloctabFLOAT(int dim1, int dim2);
int lecfima(char *ficmai, int *nutyge, int *nbtng, float ***pcoord, int *nbtel, int ***pngnel, int *nbneel, int *nbaret, int ***pnRefAr);
void createMeshFile();
