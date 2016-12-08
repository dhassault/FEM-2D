#include <stdio.h>
#include <stdlib.h>
#include "fonctions.h"
#include "melina.h"

int main()
{
  int i,j,k,kp,l;

  printf("\n");
  printf("--TP2,3 and 4: 2D finite element method--\n\n");

  char filename1[20];
  char filename2[20];
  printf("Please, enter the name of your mesh file:");
  scanf("%s",&filename1);

  int nbtng, nbtel, nbneel, nbaret;
  int nutyge;

  float** coord;
  int** ngnel;
  int** nRefAr;

  FILE* ficmai = NULL;
  ficmai = fopen(filename1,"r");
  if (ficmai != NULL)
    {
      lecfima(ficmai,&nutyge,&nbtng,&coord,&nbtel,&nbneel,&ngnel,&nbaret,&nRefAr);
      fclose(ficmai);
    }
  else
    {
      printf("Cannot read the file %s.\n", filename1);
      exit(0);
    }

//__________________________________________________________________________

  printf("\n\n");
  printf("-- Reference numbers --\n\n");

  int nrefdom;
  int nbrefD0, nbrefD1, nbrefF1;
  int nurefD0[4], nurefD1[4], nurefF1[4];

  printf("Please, enter the name of the file that contain the reference numbers: ");
  scanf("%s",&filename1);
  printf("\n\n");

  FILE* firef = NULL;
  firef = fopen(filename1,"r");
  if (firef != NULL)
    {
      lecnuref(firef,&nrefdom,&nbrefD0,nurefD0,&nbrefD1,nurefD1,&nbrefF1,nurefF1);
      fclose(firef);
    }
  else
    {
      printf("Cannot read the file %s.\n", filename1);
      exit(1);
    }

//___________________________________________________________________

  int ndim = 2;
  int nblmx = nbneel;

  float** coorel = alloctab(nbneel,ndim);
  int* nrefak = malloc(nbaret*sizeof(int));
  float** matelm = alloctab(nbneel,nbneel);
  float* scmelm = malloc(nbneel*sizeof(float));
  int* nudlel = malloc(nbneel*sizeof(int));
  float* udel = malloc(nbneel*sizeof(float));

  for (k=0;k<nbtel;k++)
    {

      coma2el(k,nbneel,coord,ngnel,coorel);
      elemkl(nrefdom,nbrefD0,nurefD0,nbrefD1,nurefD1,nbrefF1,nurefF1,nutyge,nbneel,coorel,nbaret,nRefAr[k],matelm,scmelm,nudlel,udel);
      kp = k+1;
      impmat_(&kp,&nutyge,&nbneel,&nblmx,matelm[0],scmelm,nudlel,udel);
    }
  printf("\n\n");

//____________________________________________________________________________


  printf("Please, enter the name of the SMD file (fsmd):");
  scanf("%s",&filename1);
  printf("\n");


  FILE* ficsmd = NULL;
  ficsmd = fopen(filename1,"wb");
  if (ficsmd != NULL)
    {
      assemb(nutyge,nbtng,nbtel,nbneel,nbaret,coord,ngnel,nRefAr,nrefdom,nbrefD0,nurefD0,nbrefD1,nurefD1,nbrefF1,nurefF1,ficsmd);
      fclose(ficsmd);
    }
  else
    {
      printf("Cannot read the file %s.\n", filename1);
      exit(2);
    }


  ficsmd = NULL;
  ficsmd = fopen(filename1,"rb");
  if (ficsmd != NULL)
    {
      lecsmd(ficsmd);
      fclose(ficsmd);
    }
  else
    {
      printf("Cannot read the file %s.\n", filename1);
      exit(3);
    }


//___________________________________________________________________________

  ficsmd = NULL;
  FILE* ficsmo = NULL;

  printf("Please, enter the name of the SMD file (fsmd):");
  scanf("%s",&filename1);
  ficsmd = fopen(filename1,"rb");

  printf("Please, enter the name of the SMO file (fsmo):");
  scanf("%s",&filename2);
  printf("\n");
  ficsmo = fopen(filename2,"wb");


  if (ficsmd != NULL)
    {
      if (ficsmo != NULL)
	{
	  dmdamo(ficsmd,ficsmo);
	  fclose(ficsmo);
	  fclose(ficsmd);
	}
      else
	{
	  fclose(ficsmd);
	  printf("Cannot read the file %s.\n", filename2);
	  exit(5);
	}
    }
  else
    {
      printf("Cannot read the file %s.\n", filename1);
      exit(4);
    }

  ficsmo = NULL;
  ficsmo = fopen(filename2,"rb");
  if (ficsmo != NULL)
    {
      lecsmo(ficsmo);
      fclose(ficsmo);
    }
  else
    {
      printf("Cannot read the file %s.\n", filename2);
      exit(6);
    }


  freetab(coorel);
  free(nrefak);
  freetab(matelm);
  free(scmelm);
  free(nudlel);
  free(udel);

  freetab(coord);
  freetab(ngnel);
  freetab(nRefAr);

  return 0;
}
