#ifndef MELINA_H
#define MELINA_H

/* ---------- Macro pour gerer le nom des sous-programmes Fortran ----------- */
#ifdef C_NO_TRAILING_UNDERSCORE
     /* Cas du compilateur xlf (IBM, Apple) */
     #define FORTRANNAME(x) x
#else 
     /* Cas de la plupart des autres compilateurs */
     #define FORTRANNAME(x) x##_
#endif
/* -------------------------------------------------------------------------- */

/*
           Interface avec la bibliotheque Fortran Melina
           =============================================
*/

void FORTRANNAME(delafo)
     (const int *NDIM, const int *NDMDOM, const int *NBDLG, const float *CORDLG,
      const float *DFCTBA, const int *IMPFCH, float *JACOB);

void FORTRANNAME(eldiff)
     (const int *NDIM, const int *NDMDOM, const float *TANGTS, const int *ISGDET,
      float *ELTDIF, float *NORMAL, const int *IMPFCH);

void FORTRANNAME(elemar)
     (const int *NUTYPE, const int *NUARET, int *NUTYAR, int *NBNEAR,
      int *NUNEAR, int *NBDLAR, int *NUDLAR);

void FORTRANNAME(elemco)
     (const int *NUTYPE, const int *NDIM, float *CORLOC);

void FORTRANNAME(elemfb)
     (const int *NUTYPE, const char *KELFCT, const float *COR, float *P,
      float *DP, float *DDP, const int *IMPFCH);

void FORTRANNAME(elemnb)
     (const int *NUTYPE, int *LAGHER, int *CLASSE, int *NBSOMT, int *NBARET,
      int *NBFACE, int *NBNEEL, int *NBDLEL, int *NBNEAX, int *NBDLAX,
      int *NBNEFX, int *NBDLFX);

void FORTRANNAME(elemqd)
     (const int *NUTYPE, const int *NDMDOM, const int *NUSCHQ, int *NBPTQ,
      float *POIDQ, float *COORQ, const int *IMPFCH);

int FORTRANNAME(nutyel)
    (const char *ENTREE);

void FORTRANNAME(ofaled)
     (const int *NBLIGN, const float *JACF, float *DETDF, float *JACM,
      const float *EPSMAC, const int *IMPFCH);

int FORTRANNAME(ptxqua)
    (const int *K, const int *N);

void FORTRANNAME(transf)
     (const int *NDIM, const int *NBPOIN, const float *CORPTL,
      const float *FCTBAS, float *POINT);


/*
           Interface avec les procedures Fortran fournies
           ==============================================
*/

void FORTRANNAME(affsmd)
     (const int *nblign, const int *adprcl, const int *numcol,
      const int *adsucl, const float *matris, const float *secmbr,
      const int *nuddir, const float *valdir);

void FORTRANNAME(affsmo)
     (const int *NBLIGN, const int *ADPRCO, const int *NUMCO0,
      const float *MATRI0, const float *SECMB0);

void FORTRANNAME(affsol)
     (const int *nblign, const float*coor, const float*u, const float*uex,
      const int *impfch);

void FORTRANNAME(assmat)
     (const int *I, const int *J, const float *X, int *ADPRCL, int *NUMCOL,
      int *ADSUCL,  float *LMATRI, int *NEXTAD);

void FORTRANNAME(cdesse)
     (const int *NBLIGN, const int *ADPRCL, const int *NUMCOL,
      const int *ADSUCL, const float *MATRIS, const float *SECMBR,
      const int *NUDDIR, const float *VALDIR,
      int *ADPRC0, int *NUMCO0, float *MATRI0, float *SECMB0);

int FORTRANNAME(codeel)
    (const int *itypma, const int *inter);

void FORTRANNAME(impmat)
     (const int *K, const int *NUTYEL, const int *NBNEEL, const int *NBLMX,
      const float *MATELM, const float *SCMELM,
      const int *NUDLEL, const float *UDEL);

void FORTRANNAME(impmpr)
     (const int*IMPFCH, const int *RANG, const int *PROFIL, const float *D,
      const float *L);

void FORTRANNAME(ltlpr)
     (const int*rang, const int *profil, const float *ad, const float *al,
      const float *eps, float *ld, float *ll);

void FORTRANNAME(rsprl)
     (const int*rang, const int *profil, const float *d, const float *l,
      const float *b, float *x);

void FORTRANNAME(rspru)
     (const int*rang, const int *profil, const float *d, const float *l,
      const float *b, float *x);

void FORTRANNAME(tri)
     (const int *N, int *NTAB, float *RTAB);

#endif
