#include <stdio.h>
#include <stdlib.h>
#include "fonctions.h"

int main() {
  char a;
  char filename[20]="mesh.dat";
  int i,m,n,p,q,t;
  int j; //pour tester sur printf avec lecfima, peut être supprimé après
  int **nRefAr,**ngnel;
  float **coord;
  printf("Do you want to create a custom mesh file ? [y]es or [n]o\n");
  scanf("%c",&a);
  if (a=='y') {
  createMeshFile();
  }

  /****************************
  Reading of a mesh file
  *****************************/
  char *ficmai=filename;
  int *nbtel=&m;
  int *nbneel=&p;
  int *nbtng=&n;
  int *nbaret=&q;
  int *nutyge=&t;
  int ***pnRefAr=&nRefAr;
  int ***pngnel=&ngnel;
  float ***pcoord=&coord;

  lecfima(ficmai,nutyge,nbtng,pcoord,nbtel,pngnel,nbneel,nbaret,pnRefAr); //remember to add nutyge on university's computers
  printf("There is %d geometrical nodes.\n",n);
  printf("m=%d t=%d p=%d q=%d\n",m,t,p,q);
  ///************************************************
  printf("The coordiantes are:\n");
  for (i=0;i<n;i++){
    printf("%f %f\n",coord[i][0],coord[i][1]);
  }
  printf("ngnel:-----------\n");
  for (i=0;i<m;i++){
    for (j=0;j<p;j++){
      printf("%d ",ngnel[i][j]);
    }
    printf("\n");
  }
  printf("-----------\n");
  printf("nRefAr:------------\n");
  for (i=0;i<m;i++){
    for (j=0;j<q;j++){
      printf("%d ",nRefAr[i][j]);
    }
    printf("\n");
  }
  printf("-----------\n");
  //*******************************************/
  return 0;
}
