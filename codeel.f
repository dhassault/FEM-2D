      integer function CODEEL(itypma,inter)
      integer itypma, inter
*  Point d'entree 
      integer NUTYEL
*  Variables locales
      character*5 chcode
*     ------------------------------------------------------------------------
*     Renvoie le code interne "Melina" correspondant au type de l'element
*     souhaite caracterise par les donnees d'entree.
*     Arguments d'entree  : itypma, inter
*      itypma : type des elements du maillage (1 = quadrangles, 2 = triangles)
*      inter  : degre de l'interpolation de Lagrange (1 ou 2)
*     ------------------------------------------------------------------------
      if (itypma.eq.2) then
        chcode = 'TR00x'
      else
        chcode = 'QU00x'
      endif
      write(chcode(5:5),'(i1)') inter
      CODEEL = NUTYEL(chcode)
      end
