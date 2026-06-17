      subroutine xdslf5( lson  , sonind, lpar  , stack, sonoff,
     1                   paroff )
 
c
c  purpose -- to slide a frontal matrix of a son into that of the father
c
c  created               10-feb-87, cca
c  last modifications -- 10-feb-87, cca
c                        12-sep-91, rgg, modified to remove hidden
c                                        equivalencing in stack array
c
c  input variables --
c
c      lson   -- size of the son's frontal matrix
c      sonind -- local index vector for the son
c      lpar   -- size of parent's frontal matrix
c      stack  -- array containing both the son update matrix and
c                the parent frontal matrix which overlap.
c      sonoff -- offset in stack of son update matrix
c      paroff -- offset in stack of parent frontal matrix
c
c  output variable
c
c      stack  -- array containing both the son update matrix and
c                the parent frontal matrix which overlap.
c
c  subprograms called --
c
c      none
c
c  =====================================================================
 
      integer           lson, sonind(*), lpar, sonoff, paroff
 
      double precision  stack(*)
 
      integer           ic    , j     , jc    , kpar  , kson  , parbeg
 
      double precision  temp
 
c  =====================================================================
 
      kson = (lson*(lson+1))/2
 
      do jc = lson,1,-1
 
          j = sonind(jc)
          parbeg = lpar*(j-1) - (j*(j-1))/2
 
          do ic = lson,jc,-1
              temp               = stack(sonoff+kson)
              stack(sonoff+kson) = 0.0
              kson               = kson - 1
              kpar         =       parbeg + sonind(ic)
              stack(paroff+kpar) = temp
          enddo
 
      enddo
 
c  =====================================================================
 
      return
      end
