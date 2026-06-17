      subroutine xislfi (nsnind, lindxg, invsup)
 
c
c     update global indices with inverse permutation generated
c     during the factorization
c
c---------------------------------------------------------------------
 
      integer    nsnind
      integer    lindxg (*), invsup (*)
 
      integer     i
 
c---------------------------------------------------------------------
 
      do i = 1, nsnind
         lindxg (i) = invsup ( lindxg (i) )
      enddo
 
c---------------------------------------------------------------------
 
      return
      end
