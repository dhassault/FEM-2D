
#include "fonctions.h"

int **alloctabINT(int dim1, int dim2) {
  /********************************************
   create a table with adjacent memory spaces
  **********************************************/
  int **ptr;

  ptr = malloc(dim1*sizeof(int *));
  if (ptr != NULL) {
    int i, taille_ligne = dim2*sizeof(int);
    int *tmp = malloc(dim1*taille_ligne);
    if (tmp != NULL) {
      for (i=0; i<dim1; i++) {
  	ptr[i] = tmp;
  	tmp += dim2;
  	}
      }
    else
      ptr = NULL;
    }
  return(ptr);
}

float **alloctabFLOAT(int dim1, int dim2) {
  /********************************************
   create a table with adjacent memory spaces
  **********************************************/
  float **ptr;

  ptr = malloc(dim1*sizeof(float *));
  if (ptr != NULL) {
    int i, taille_ligne = dim2*sizeof(float);
    float *tmp = malloc(dim1*taille_ligne);
    if (tmp != NULL) {
      for (i=0; i<dim1; i++) {
  	ptr[i] = tmp;
  	tmp += dim2;
  	}
      }
    else
      ptr = NULL;
    }
  return(ptr);
}

void createMeshFile(){
  /******************************
  Creation of a custom mesh file
  ******************************/
  int i, j;
  int n1,n2,typel,nbaret,nbtel;
  double a,b,c,d,x,y,dx,dy;
  int nrefcot[4] = {1, 2, 3, 4}; //number of sides
  int nrefdom=0; //the central domain is 0
  printf("Choose the number of nodes n1 (on [a,b]) and n2 (on [c,d])...\n");
  printf("n1: ");
  scanf("%d",&n1);
  printf("n2: ");
  scanf("%d",&n2);
  printf("Choose the coordinate of [a,b]x[c,d] in (O;x,y)...\n");
  printf("a: ");
  scanf("%lf",&a);
  printf("b: ");
  scanf("%lf",&b);
  printf("c: ");
  scanf("%lf",&c);
  printf("d: ");
  scanf("%lf",&d);

  dx=(b-a)/n2; //step on x
  dy=(d-c)/n1; //step on y

  FILE *fp = NULL;
  char filename[20]="mesh.dat";
  fp = fopen (filename,"w");
  if (fp == NULL) {
      fprintf (stderr,"Cannot open the file\n");
  }

  fprintf(fp,"%d\n",n1*n2); //number of nodes
  //we write coordinates of each points in the file mesh.dat
  for (y=c; y<d; y+=dy){
    for (x=a; x<b; x+=dx){
      fprintf (fp,"%3.1lf %3.1lf\n",x,y);
    }
  }

  printf("Do you want squares (1) or triangles (2)?\n");
  scanf("%d",&typel);
  if (typel==1) {  //for squares
    nbaret=4;
    nbtel=(n1-1)*(n2-1); //number of squares
    fprintf(fp,"%d %d %d %d\n",nbtel,typel,nbaret,nbaret);
  }
  else if (typel==2) { //for triangles
    nbaret=3;
    nbtel=(n1-1)*(n2-1)*2; //number of triangles
    fprintf(fp,"%d %d %d %d\n",nbtel,typel,nbaret,nbaret);
  }
  else {
    printf("Please restart the program and choose (1) for squares or (2) for triangles...\n");
  }

  int **nRefAr=alloctabINT(nbtel, nbaret);
  etiqar(typel,n1,n2,nrefdom,nrefcot,nbtel,nbaret,nRefAr);

  for (i=0;i<nbtel;i++){
    for (j=1;j<=nbaret;j++){
      fprintf(fp,"%d ",j);
    }
    for (j=0;j<nbaret;j++){
      fprintf(fp,"%d ",nRefAr[i][j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  printf("The file created is named mesh.dat...\n");
}
void etiqar(int typel, int n1, int n2, int nrefdom, const int *nrefcot, int nbtel, int nbaret, int **nRefAr) {
  /******************************
        Fill nRefAr table
  ******************************/
  int i,j;
  for (i=0;i<nbtel;i++){ //the whole nRefAr is firstly initilized to domain 0 for convinience
    for (j=0;j<nbaret;j++){
      nRefAr[i][j]=nrefdom;
    }
  }
  if (typel==1) {  //for squares
    for (i=0;i<n1-1;i++){
      nRefAr[i][0]=1;
    }
    for (i=n1-2;i<nbtel;i+=n1-1){
      nRefAr[i][1]=2;
    }
    for (i=nbtel-(n2-1);i<nbtel;i++){
      nRefAr[i][2]=3;
    }
    for (i=0;i<nbtel-(n1-3);i+=n1-1){
      nRefAr[i][3]=4;
    }
}
  else if (typel==2) { //for triangles
    for (i=0;i<((n1-1)*2)-1;i+=2){
      nRefAr[i][0]=1;
    }
    for (i=((n1-1)*2)-1;i<nbtel;i+=(n1*2)-2){
      nRefAr[i][0]=2;
    }
    for (i=(nbtel-1)-(((n1-1)*2)-2);i<nbtel;i+=2){
      nRefAr[i][1]=3;
    }
    for (i=0;i<nbtel-(n1-3);i+=(n1-1)*2){
      nRefAr[i][2]=4;
    }
  }
}

void lecfima(char *ficmai, int *nutyge, int *nbtng, float ***pcoord, int *nbtel, int ***pngnel, int *nbneel, int *nbaret, int ***pnRefAr){
  /******************************************
   Read an existing mesh file named mesh.dat
  ******************************************/
  int i,j;
  FILE *fp = NULL;
  fp = fopen (ficmai,"r");
  if (fp == NULL) {
      fprintf (stderr,"Cannot open the file\n");
  }
  fscanf(fp,"%d",nbtng);
  *pcoord=alloctabFLOAT(*nbtng,2);
  for (i=0;i<*nbtng;i++){
      fscanf(fp,"%f %f",&(*pcoord)[i][0],&(*pcoord)[i][1]);
  }
  fscanf(fp,"%d %d %d %d",nbtel,nutyge,nbneel,nbaret);
  *pngnel=alloctabINT(*nbtel,*nbneel);
  *pnRefAr=alloctabINT(*nbtel,*nbaret);
  for (i=0;i<*nbtel;i++){
    for (j=0;j<*nbneel;j++){
      fscanf(fp,"%d",&(*pngnel)[i][j]);
    }
    for (j=0;j<*nbaret;j++){
      fscanf(fp,"%d",&(*pnRefAr)[i][j]);
    }
  }
}
