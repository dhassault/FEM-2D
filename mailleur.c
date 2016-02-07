#include <stdio.h>
#include <stdlib.h>
#include "fonctions.h"

int main() {
  char a;
  printf("Do you want to create a custom mesh file ? [y]es or [n]o\n");
  scanf("%c",&a);
  if (a=="y") {
    void createMeshFile();
  }
  /****************************
  Reading of a mesh file
  *****************************/
  float **coord;
  //char *ficmai = filename;
  int nbtng;
  float ***pcoord=&coord;
  int ***pnRefAr=&nRefAr;

  //pnRefAr=lecfima(ficmai,nbtng,pcoord,nbtel,pngnel,nbneel,nbaret,pnRefAr); //remember to add nutyge on university's computers

  return 0;
}
