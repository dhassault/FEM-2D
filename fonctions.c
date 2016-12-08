
#include "fonctions.h"
#include "melina.h"


void adwdw(int nbneel, float** dfcbas, float** ijacob, float eltdif, float** cofvar, float** matelm) {

  int i;
  for (i=0;i<nbneel;i++) {
    int j;
    for (j=0;j<nbneel;j++) {
      int a;
      for (a=0;a<2;a++) {
        float ci;
        ci = dfcbas[i][0]*ijacob[a][0] + dfcbas[i][1]*ijacob[a][1];
        int b;
        for (b=0;b<2;b++) {
          float cj;
          cj = dfcbas[j][0]*ijacob[b][0] + dfcbas[j][1]*ijacob[b][1];
          matelm[i][j] += eltdif*cofvar[a][b]*ci*cj;
        }
      }
    }
  }

}




float **alloctab(int dim1, int dim2) {
  /*--------------------------------------------------------------------------------
    Cette fonction alloue de la memoire pour stocker une matrice de dimensions
    dim1 x dim2 de type float. La fonction alloue un tableau de pointeurs (ptr)
    dont chacun des dim1 elements pointe vers une zone de dim2 elements.

    La fonction renvoie NULL en cas d'erreur lors de l'allocation.

    La liberation de la memoire allouee peut etre faite par appel a freetab.
    Remarque :
      Les elements utiles de la matrice sont ranges consecutivement en memoire,
      en dim1 blocs de dim2 elements chacun. Il s'ensuit que cet ensemble peut
      etre transmis a une procedure Fortran via la valeur de ptr[0], comme un seul
      tableau de dimensions (dim2, dim1). On a donc la correspondance d'adressage
      suivante (ptr[i][j] et tab(j,i) designent le meme element) :
        - langage C : ptr[i][j], i=0,dim1-1, j=0,dim2-1,
        - Fortran   : tab(j,i),  i=1,dim1,   j=1,dim2.
      En pratique, etant donnees le sous-programme Fortran :
        subroutine trucF (tab,nbli,nbco)
        integer nbli, nbco
        real tab(nbli, nbco)
      et la fonction C :
        void trucC(float **ptr)
      l'allocation de memoire et les appels se feront de la maniere suivante :
        ptr = alloctab(nbco, nbli);
        FORTRANNAME(trucF) (ptr[0], nbli, nbco);
        trucC (ptr);
      La reference aux elements de la matrice se font alors par la notation a
      deux indices indiquee ci-dessus.
  --------------------------------------------------------------------------------*/
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


int **alloctabint(int dim1, int dim2) {
  /*--------------------------------------------------------------------------------
    Cette fonction alloue de la memoire pour stocker une matrice de dimensions
    dim1 x dim2 de type int. La fonction alloue un tableau de pointeurs (ptr)
    dont chacun des dim1 elements pointe vers une zone de dim2 elements.

    La fonction renvoie NULL en cas d'erreur lors de l'allocation.

    La liberation de la memoire allouee peut etre faite par appel a freetab.
    Remarque :
      Les elements utiles de la matrice sont ranges consecutivement en memoire,
      en dim1 blocs de dim2 elements chacun. Il s'ensuit que cet ensemble peut
      etre transmis a une procedure Fortran via la valeur de ptr[0], comme un seul
      tableau de dimensions (dim2, dim1). On a donc la correspondance d'adressage
      suivante (ptr[i][j] et tab(j,i) designent le meme element) :
        - langage C : ptr[i][j], i=0,dim1-1, j=0,dim2-1,
        - Fortran   : tab(j,i),  i=1,dim1,   j=1,dim2.
      En pratique, etant donnees le sous-programme Fortran :
        subroutine trucF (tab,nbli,nbco)
        integer nbli, nbco
        integer tab(nbli, nbco)
      et la fonction C :
        void trucC(int **ptr)
      l'allocation de memoire et les appels se feront de la maniere suivante :
        ptr = alloctab(nbco, nbli);
        FORTRANNAME(trucF) (ptr[0], nbli, nbco);
        trucC (ptr);
      La reference aux elements de la matrice se font alors par la notation a
      deux indices indiquee ci-dessus.
  --------------------------------------------------------------------------------*/

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

void assemb(int nutyge, int nbtng, int nbtel, int nbneel, int nbaret, float** coord,
            int** ngnel, int** nRefAr, int nrefdom, int nbrefD0, int* nurefD0, int nbrefD1,
            int* nurefD1, int nbrefF1, int* nurefF1, FILE* ficsmd) {

  int ig,j,jg,ngmax,ngmin;

  int ndim = 2;
  int nblign = nbtng;
  int nbcoef, nbcoefmax = 3*nblign;
  int nextad;


  float** coorel = alloctab(nbneel,ndim);
  int* nrefak = malloc(nbaret*sizeof(int));
  float** matelm = alloctab(nbneel,nbneel);
  float* scmelm = malloc(nbneel*sizeof(float));
  int* nudlel = malloc(nbneel*sizeof(int));
  float* udel = malloc(nbneel*sizeof(float));

  float* secmbr = malloc(nblign*sizeof(float));
  int* nuddir = malloc(nblign*sizeof(int));
  float* valdir = malloc(nblign*sizeof(float));
  int* adprcl = malloc(nblign*sizeof(int));
  float* matris = malloc((nblign+nbcoefmax)*sizeof(float));
  int* numcol = malloc(nbcoefmax*sizeof(int));
  int* adsucl = malloc(nbcoefmax*sizeof(int));

  int i;
  for (i=0;i<nblign;i++) {
      secmbr[i] = 0.;
      nuddir[i] = i+1;
      adprcl[i] = 0;
      valdir[i] = 0.;
    }
  for (i=0;i<(nblign+nbcoefmax);i++) matris[i] = 0.;
  for (i=0;i<nbcoefmax;i++) {
      numcol[i] = 0;
      adsucl[i] = 0;
    }


  nextad = 1;
  int k;
  for (k=0;k<nbtel;k++) {
    coma2el(k,nbneel,coord,ngnel,coorel);
    elemkl(nrefdom,nbrefD0,nurefD0,nbrefD1,nurefD1,nbrefF1,nurefF1,nutyge,nbneel,coorel,nbaret,nRefAr[k],matelm,scmelm,nudlel,udel);

    for (i=0;i<nbneel;i++) {
      ig = ngnel[k][i];

      matris[ig-1] += matelm[i][i];

      secmbr[ig-1] += scmelm[i];

      // Note: conditions de Dirichlet non homogene (-1) et homogene (0)
      if (nudlel[i] == -1) {
        nuddir[ig-1] = -ig;
        valdir[ig-1] = udel[i];
      }
      else if (nudlel[i] == 0) {
        nuddir[ig-1] = 0;
      }
    }

      // Lower triangle part of the matrix
      for (i=1;i<nbneel;i++) {
        ig = ngnel[k][i];

        for (j=0;j<i;j++) {
          jg = ngnel[k][j];
          ngmax = max(ig,jg);
          ngmin = min(ig,jg);
          assmat_(&ngmax,&ngmin,&matelm[i][j],adprcl,numcol,adsucl,&matris[nblign],&nextad);
        }
      }
    }

  adprcl[nblign-1] = nextad;
  nbcoef = adprcl[nblign-1]-1;

  // Binary file ficsmd
  fwrite(&nblign,sizeof(int),1,ficsmd);
  fwrite(secmbr,sizeof(float),nblign,ficsmd);
  fwrite(nuddir,sizeof(int),nblign,ficsmd);
  fwrite(valdir,sizeof(float),nblign,ficsmd);
  fwrite(adprcl,sizeof(int),nblign,ficsmd);
  fwrite(matris,sizeof(float),(nblign+nbcoef),ficsmd);
  fwrite(numcol,sizeof(int),nbcoef,ficsmd);
  fwrite(adsucl,sizeof(int),nbcoef,ficsmd);

  free(secmbr);
  free(nuddir);
  free(valdir);
  free(numcol);
  free(adprcl);
  free(matris);
  free(adsucl);

  freetab(coorel);
  free(nrefak);
  freetab(matelm);
  free(scmelm);
  free(nudlel);
  free(udel);
}

int max(int a, int b){
  if (a < b) {
    a = b;
  }
  return a;
}

int min(int a, int b){
  if (b < a) {
    a = b;
  }
  return a;
}

void calelb(int nutygear, int nbnear, float** cooreb, float** matela, float* scmela) {
  int imp = 1;
  int ndim = 2;
  int ndmdom = ndim-1;
  int nuschq = -2;
  int nbptq, nbptqmax;
  nbptqmax = ptxqua_(&nuschq,&ndmdom); // max points


  float* poidsq = malloc(nbptqmax*sizeof(float));
  float** coorq = alloctab(nbptqmax,ndmdom);
  elemqd_(&nutygear,&ndmdom,&nuschq,&nbptq,poidsq,coorq[0],&imp);


  float* P = malloc(nbnear*sizeof(float));
  float** DP = alloctab(nbnear,ndim);
  float** DDP = alloctab(nbnear,ndim);
  char kelfct[3] = "YYN";


  int isgdet = 1;
  float** jacob = alloctab(ndmdom,ndim);
  float* normal = malloc(ndim*sizeof(float));


  float* coord_point = malloc(ndim*sizeof(float));

  float eltdif;
  float cofvar;

  int i;
  for (i=0;i<nbptq;i++) {

      elemfb_(&nutygear,kelfct,coorq[i],P,DP[0],DDP[0],&imp);


      delafo_(&ndim,&ndmdom,&nbnear,cooreb[0],DP[0],&imp,jacob[0]);


      eldiff_(&ndim,&ndmdom,jacob[0],&isgdet,&eltdif,normal,&imp);
      eltdif *= poidsq[i];


      transf_(&ndim,&nbnear,cooreb[0],P,coord_point);


      cofvar = bn(coord_point);
      ww(nbnear,P,eltdif,cofvar,matela);


      cofvar = fn(coord_point);
      w(nbnear,P,eltdif,cofvar,scmela);
    }

  free(coord_point);
  freetab(jacob);
  free(P);
  freetab(DP);
  freetab(DDP);
  free(poidsq);
  freetab(coorq);
}

void calelm(int nutyge, int nbneel, float** coorel, float** matelm, float *scmelm)
{
  int i,j;
  int imp = 1;
  int ndim = 2;
  int ndmdom = ndim;
  int nuschq = -2;
  int nbptq, nbptqmax;
  nbptqmax = ptxqua_(&nuschq,&ndmdom);


  float* poidsq;
  poidsq = malloc(nbptqmax*sizeof(float));
  float** coorq;
  coorq = alloctab(nbptqmax,ndmdom);
  elemqd_(&nutyge,&ndmdom,&nuschq,&nbptq,poidsq,coorq[0],&imp);


  float *P;
  float **DP, **DDP;
  P = malloc(nbneel*sizeof(float));
  DP = alloctab(nbneel,ndim);
  DDP = alloctab(nbneel,ndim);
  char kelfct[3] = "YYN";

  float **jacob, **ijacob;
  jacob = alloctab(ndmdom,ndim);
  ijacob = alloctab(ndmdom,ndim);
  float detjacob;
  float epsilon = 1e-5;


  float* coord_point = malloc(ndim*sizeof(float));

  float eltdif;
  float** cofvar_adwdw;
  cofvar_adwdw = alloctab(ndim,ndim);
  float cofvar;

  for (i=0;i<nbptq;i++)
    {

      elemfb_(&nutyge,kelfct,coorq[i],P,DP[0],DDP[0],&imp);


      delafo_(&ndim,&ndmdom,&nbneel,coorel[0],DP[0],&imp,jacob[0]);


      ofaled_(&ndim,jacob[0],&detjacob,ijacob[0],&epsilon,&imp);
      if (detjacob < 0) detjacob = -detjacob;


      transf_(&ndim,&nbneel,coorel[0],P,coord_point);


      eltdif = poidsq[i]*detjacob;
      cofvar_adwdw[0][0] = a11(coord_point);
      cofvar_adwdw[0][1] = a12(coord_point);
      cofvar_adwdw[1][0] = a21(coord_point);
      cofvar_adwdw[1][1] = a22(coord_point);
      adwdw(nbneel,DP,ijacob,eltdif,cofvar_adwdw,matelm);


      cofvar = a00(coord_point);
      ww(nbneel,P,eltdif,cofvar,matelm);


      cofvar = fomega(coord_point);
      w(nbneel,P,eltdif,cofvar,scmelm);
    }

  freetab(cofvar_adwdw);
  free(coord_point);
  freetab(jacob);
  freetab(ijacob);
  free(P);
  freetab(DP);
  freetab(DDP);
  free(poidsq);
  freetab(coorq);
}

void coel2b(int nbnear, int* nunear, float** coorel, float** cooreb)
{
  int i,j,ni;
  int ndim = 2;

  for (i=0;i<nbnear;i++)
    {
      ni = nunear[i]-1;
      for (j=0;j<ndim;j++) cooreb[i][j] = coorel[ni][j];
    }
}

void coma2el(int k, int nbneel, float** coord, int** ngnel, float** coorel) {
  int i, numglb;

  for (i=0; i<nbneel; i++)
    {
      numglb = ngnel[k][i] - 1;
      coorel[i][0] = coord[numglb][0];
      coorel[i][1] = coord[numglb][1];
    }
}

float a11(float* x) {
  return 1.;
}

float a12(float* x) {
  return 0.;
}

float a21(float* x) {
  return 0.;
}

float a22(float* x) {
  return 1.;
}

float a00(float* x) {
  return 1.;
}

float bn(float* x) {
  return 1.;
}

float fomega(float* x) {
  return 1.;
}

float fn(float* x) {
  return 1.;
}

float ud(float* x) {
  return (100*x[0]+x[1]);
}

void dmdamo(FILE* ficsmd, FILE* ficsmo)
{

  int nblign, nbcoef;
  fread(&nblign,sizeof(int),1,ficsmd);
  float* secmbr = malloc(nblign*sizeof(float));
  int* nuddir = malloc(nblign*sizeof(int));
  float* valdir = malloc(nblign*sizeof(float));
  int* adprcl = malloc(nblign*sizeof(int));
  fread(secmbr,sizeof(float),nblign,ficsmd);
  fread(nuddir,sizeof(int),nblign,ficsmd);
  fread(valdir,sizeof(float),nblign,ficsmd);
  fread(adprcl,sizeof(int),nblign,ficsmd);
  nbcoef = adprcl[nblign-1] -1;
  float* matris = malloc((nblign+nbcoef)*sizeof(float));
  int* numcol = malloc(nbcoef*sizeof(int));
  int* adsucl = malloc(nbcoef*sizeof(int));
  fread(matris,sizeof(float),(nblign+nbcoef),ficsmd);
  fread(numcol,sizeof(int),nbcoef,ficsmd);
  fread(adsucl,sizeof(int),nbcoef,ficsmd);


  int nbcoeo, nbcoeomax = nbcoef;
  float* secmbo = malloc(nblign*sizeof(float));
  int* adprco = malloc(nblign*sizeof(int));
  float* matrio = malloc((nblign+nbcoeomax)*sizeof(float));
  int* numcoo = malloc(nbcoeomax*sizeof(int));

  int i;
  for (i=0;i<nblign;i++)
    {
      secmbo[i] = 0.;
      adprco[i] = 0;
    }
  for (i=0;i<(nblign+nbcoeomax);i++) matrio[i] = 0.;
  for (i=0;i<nbcoeomax;i++) numcoo[i] = 0;

  cdesse_(&nblign,adprcl,numcol,adsucl,matris,secmbr,nuddir,valdir,adprco,numcoo,matrio,secmbo);
  nbcoeo = adprco[nblign-1] - 1;


  fwrite(&nblign,sizeof(int),1,ficsmo);
  fwrite(secmbo,sizeof(float),nblign,ficsmo);
  fwrite(adprco,sizeof(int),nblign,ficsmo);
  fwrite(matrio,sizeof(float),(nblign+nbcoeo),ficsmo);
  fwrite(numcoo,sizeof(int),nbcoeo,ficsmo);

  free(secmbr);
  free(nuddir);
  free(valdir);
  free(numcol);
  free(adprcl);
  free(matris);
  free(adsucl);

  free(secmbo);
  free(adprco);
  free(matrio);
  free(numcoo);
}

void elemkl(int nrefdom, int nbrefD0, int* nurefD0, int nbrefD1, int* nurefD1,
  int nbrefF1, int* nurefF1, int nutyge, int nbneel, float** coorel, int nbaret,
  int* nrefak,float** matelm, float* scmelm, int* nudlel, float* udel) {

  int imp = 1;
  int ip,l,nj,nl;
  int ndim = 2;
  int ndmdom = ndim -1;


  int nutygear;
  int nbnear, nbnearmax = 2;
  int nbdlar,nbdlarmax = nbneel;
  char kelfct[3] = "YNN";
  float udx;


  int* nunear = malloc(nbnearmax*sizeof(int));
  int* nudlar = malloc(nbdlarmax*sizeof(int));
  float** cooreb = alloctab(nbnearmax,ndim);
  float* corloc = malloc(nbnearmax*ndmdom*sizeof(float));
  float* P = malloc(nbnearmax*sizeof(float));
  float** DP = alloctab(nbnearmax,ndim);
  float** DDP = alloctab(nbnearmax,ndim);
  float* coord_point = malloc(ndim*sizeof(float));
  float** matela = alloctab(nbnearmax,nbnearmax);
  float* scmela = malloc(nbnearmax*sizeof(float));

  int i;
  for (i=0;i<nbneel;i++) {
    int j;
    for (j=0;j<nbneel;j++) matelm[i][j] = 0.;
    scmelm[i] = 0.;
    nudlel[i] = 1;
    udel[i] = 0.;
  }

  calelm(nutyge,nbneel,coorel,matelm,scmelm);

  int cas = 0;
  for (i=0;i<nbaret;i++) {

    if (nrefak[i] != nrefdom) {

  	ip = i+1;
  	elemar_(&nutyge,&ip,&nutygear,&nbnear,nunear,&nbdlar,nudlar);
  	coel2b(nbnear,nunear,coorel,cooreb);

    int j;
    for (j=0;j<nbrefD0;j++) {
  	  if (nrefak[i] == nurefD0[j]) {
	      cas = 1;
	      break;
	    }
    }


	for (j=0;j<nbrefD1;j++) {
	  if (nrefak[i] == nurefD1[j]) {
	      cas = 2;
	      break;
      }
    }


	for (j=0;j<nbrefF1;j++) {
	    if (nrefak[i] == nurefF1[j]) {
        cas = 3;
        break;
      }
    }


	if (cas == 1) {
    for (l=0;l<nbnear;l++) {
      nl = nunear[l]-1;
      nudlel[nl] = 0;
    }
  }


	else if (cas == 2) {
    elemco_(&nutygear,&ndmdom,corloc);
    for (l=0;l<nbnear;l++) {
      nl = nunear[l] - 1;
      nudlel[nl] = -1;
      elemfb_(&nutygear,kelfct,&corloc[l],P,DP[0],DDP[0],&imp);
      transf_(&ndim,&nbnear,cooreb[0],P,coord_point);
      udel[nl] = ud(coord_point);
    }
  }


  else if (cas == 3) {
    for (j=0;j<nbnearmax;j++) {
      for (l=0;l<nbnearmax;l++) matela[j][l] = 0;
      scmela[j] = 0;
    }

    calelb(nutygear,nbnear,cooreb,matela,scmela);

    for (j=0;j<nbnear;j++) {
      nj = nunear[j]-1;
      for (l=0;l<nbnear;l++) {
        nl = nunear[l]-1;
        matelm[nj][nl] += matela[j][l];
      }
      scmelm[nj] += scmela[j];
    }
  }
  }
  }

  free(nunear);
  free(nudlar);
  freetab(cooreb);
  free(corloc);
  free(P);
  freetab(DP);
  freetab(DDP);
  free(coord_point);
  freetab(matela);
  free(scmela);
}

void freetab(void *ptr) {
  void **ptrT=ptr;
  free(ptrT[0]);
  free(ptr);
}

void lecfima(FILE* ficmai, int* nutyge, int* nbtng, float*** pcoord, int* nbtel,
  int* nbneel, int*** pngnel, int* nbaret, int*** pnRefAr) {

  fscanf(ficmai,"%d", nbtng);

  *pcoord = alloctab(*nbtng,2);
  int i;
  for (i=0;i<(*nbtng);i++) fscanf(ficmai,"%f %f", &((*pcoord)[i][0]), &((*pcoord)[i][1]));

  int t;
  fscanf(ficmai,"%d %d",nbtel,&t);

  int inter = 1;
  *nutyge = codeel_(&t,&inter);
  *nbneel = 5-t;
  *nbaret = 5-t;

  *pngnel = alloctabint(*nbtel,*nbneel);
  *pnRefAr = alloctabint(*nbtel,*nbaret);
  int k;
  for (k=0;k<(*nbtel);k++)
    {
      int l;
      for (l=0;l<(*nbneel);l++) fscanf(ficmai,"%d", &((*pngnel)[k][l]));
      for (l=0;l<(*nbaret);l++) fscanf(ficmai,"%d", &((*pnRefAr)[k][l]));
    }
}

void lecnuref(FILE* firef, int* nrefdom, int* nbrefD0, int* nurefD0, int* nbrefD1,
  int* nurefD1, int* nbrefF1, int* nurefF1) {

  int i;
  fscanf(firef, "%d", nrefdom);
  fscanf(firef, "%d", nbrefD0);
  for (i=0;i<(*nbrefD0);i++) fscanf(firef, "%d", &(nurefD0[i]));
  fscanf(firef, "%d", nbrefD1);
  for (i=0;i<(*nbrefD1);i++) fscanf(firef, "%d", &(nurefD1[i]));
  fscanf(firef, "%d", nbrefF1);
  for (i=0;i<(*nbrefF1);i++) fscanf(firef, "%d", &(nurefF1[i]));
}

void lecsmd(FILE* ficsmd) {

	int nblign, nbcoef;
	fread(&nblign,sizeof(int),1,ficsmd);
	float* secmbr = malloc(nblign*sizeof(float));
	int* nuddir = malloc(nblign*sizeof(int));
	float* valdir = malloc(nblign*sizeof(float));
	int* adprcl = malloc(nblign*sizeof(int));
	fread(secmbr,sizeof(float),nblign,ficsmd);
	fread(nuddir,sizeof(int),nblign,ficsmd);
	fread(valdir,sizeof(float),nblign,ficsmd);
	fread(adprcl,sizeof(int),nblign,ficsmd);
	nbcoef = adprcl[nblign-1]-1;
	float* matris = malloc((nblign+nbcoef)*sizeof(float));
	int* numcol = malloc(nbcoef*sizeof(int));
	int* adsucl = malloc(nbcoef*sizeof(int));
	fread(matris,sizeof(float),(nblign+nbcoef),ficsmd);
	fread(numcol,sizeof(int),nbcoef,ficsmd);
	fread(adsucl,sizeof(int),nbcoef,ficsmd);


	affsmd_(&nblign,adprcl,numcol,adsucl,matris,secmbr,nuddir,valdir);
	printf("\n");

	free(secmbr);
	free(nuddir);
	free(valdir);
	free(numcol);
	free(adprcl);
	free(matris);
	free(adsucl);
}

void lecsmo(FILE* ficsmo) {

	int nblign, nbcoeo;
	fread(&nblign,sizeof(int),1,ficsmo);
	float* secmbo = malloc(nblign*sizeof(float));
	int* adprco = malloc(nblign*sizeof(int));
	fread(secmbo,sizeof(float),nblign,ficsmo);
	fread(adprco,sizeof(int),nblign,ficsmo);
	nbcoeo = adprco[nblign-1] - 1;
	float* matrio = malloc((nblign+nbcoeo)*sizeof(float));
	int* numcoo = malloc(nbcoeo*sizeof(int));
	fread(matrio,sizeof(float),(nblign+nbcoeo),ficsmo);
	fread(numcoo,sizeof(int),nbcoeo,ficsmo);


	affsmo_(&nblign,adprco,numcoo,matrio,secmbo);
	printf("\n");

	free(secmbo);
	free(adprco);
	free(matrio);
	free(numcoo);
}

void w(int nbneel, float *fctbas, float eltdif, float cofvar, float *scmelm) {

  float coeff;

  coeff = eltdif*cofvar;
  int j;
  for (j=0;j<nbneel;j++) scmelm[j] += coeff*fctbas[j];
}

void ww(int nbneel, float *fctbas, float eltdif, float cofvar, float **matelm) {

  float coeff;

  int i;
  for (i=0; i<nbneel; i++) {
    coeff = eltdif*cofvar*fctbas[i];
    int j;
    for (j=0; j<nbneel; j++)
      matelm[i][j] += coeff*fctbas[j];
    }
}
