      subroutine xdslf7( lson  , sonind, sonstk, lpar  , parstk )
 
c
c  purpose -- to add a frontal matrix of a son into that of the father
c
c  created               10-feb-87, cca
c  last modifications -- 10-feb-87, cca
c
c  input variables --
c
c      lson   -- size of the son's frontal matrix
c      sonind -- local index vector for the son
c      sonstk -- son's frontal matrix on the stack
c      lpar   -- size of parent's frontal matrix
c
c  output variable
c
c      parstk -- parent's frontal matrix
c
c  subprograms called --
c
c      none
c
c  =====================================================================
 
      integer           lson, sonind(*), lpar
 
      double precision  sonstk(*), parstk(*)
 
      integer           ic    , j     , jc    , kpar  , kson  , parbeg
 
c  =====================================================================
 
      kson = 1
      do jc = 1,lson
          j = sonind(jc)
          parbeg = lpar*(j-1) - (j*(j-1))/2
cdir$ ivdep
          do ic = jc,lson
              kpar         = parbeg + sonind(ic)
              parstk(kpar) = parstk(kpar) + sonstk(kson)
              kson         = kson + 1
          enddo
      enddo
 
c  =====================================================================
 
      return
      end
