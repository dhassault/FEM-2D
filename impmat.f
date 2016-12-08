************************************************************************
      SUBROUTINE IMPMAT(K,NUTYEL,NBNEEL,NBLMX,MATELM,SCMELM,NUDLEL,UDEL)
************************************************************************
* Imprime les resultats de la matrice et du second membre elementaires
* ainsi que les conditions Dirichlet en  chaque noeud 
* et les valeurs des conditions Dirichlet non homogene 
* 
*** Arguments *** 
*  K      : Numero de l'element  
*  NUTYEL : Numero de type de l'element
*  NBNEEL : Nombre de noeuds de l'element
*  NBLMX  : Nombre de lignes de MATELM (leading dimension >= NBNEEL)
*  MATELM : Matrice elementaire de masse
*  SCMELM : Second membre elementaire
*  NUDLEL : Tableau des conditions de Dirichlet pour les noeuds
*  UDEL   : Tableau des valeurs de blocage 
*           pour les noeuds Dirichlet non homogene
************************************************************************
*
      INTEGER K,NUTYEL,NBNEEL,NUDLEL(*)
      REAL MATELM(NBLMX,*),SCMELM(*),UDEL(*)
*
      INTEGER I,J
*-----------------------------------------------------------------------
      WRITE(*,*)
      WRITE(*,10000)K,NUTYEL,NBNEEL
      WRITE(*,10001)
      DO 10 I=1,NBNEEL
        WRITE(*,10002)NUDLEL(I),UDEL(I),SCMELM(I),(MATELM(I,J),J=1,I)
10    CONTINUE
*
*  FORMATS
10000 FORMAT(T2,'ELEMENT=',I3,'    DE TYPE=',I5,'    NB NOEUDS=',I2)
10001 FORMAT(T2,'NUDLEL',TR4,'UDEL',TR6,'SCMELM',TR6,'MATELM')
10002 FORMAT(T2,I6,11(' ',E10.4))
      END
