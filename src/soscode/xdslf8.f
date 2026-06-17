      subroutine xdslf8( lson  , sonind, wafil2, iopos2, ocbufr,
     1                   bfrpos, bfrlen, locbfr, lpar  , parstk,
     2                   ierr )
 
c
c  purpose -- to add a frontal matrix of a son into that of the father
c             where son's matrix is out-of-core.  derived from xdslf7.
c
c  created               13-sep-89, rgg
c  last modifications    04-nov-92, rgg -- reduce number of i/o accesses
c                                          to wafil2 by use of ocbufr
c                                          array.
c
c  input variables --
c
c      lson   -- size of the son's frontal matrix
c      sonind -- local index vector for the son
c      wafil2 -- i/o file where son's frontal matrix is located.
c      iopos2 -- position on wafil2 where son'f frontal matrix starts.
c      bfrpos -- current position into out-of-core buffer array, ocbufr.
c      bfrlen -- current length of out-of-core buffer array, ocbufr.
c      locbfr -- dimensioned length of out-of-core buffer array, ocbufr.
c      lpar   -- size of parent's frontal matrix
c
c  working storage --
c
c      ocbufr -- buffer for reading in son's frontal matrix
c
c  output variable
c
c      parstk -- parent's frontal matrix
c
c  subprograms called --
c
c      xdslw1
c
c  =====================================================================
 
      integer           lson, sonind(*), wafil2, iopos2, lpar,
     1                  bfrpos, bfrlen, locbfr
 
      double precision  ocbufr(*), parstk(*)
 
      integer           ic    , ierr  , j     , jc    , kpar  ,
     1                  len   , ncol  , nrow  , parbeg

      integer           idummy(1)
 
      double precision  temp
 
c  =====================================================================
 
      do jc = 1,lson
 
          if ( bfrpos .gt. bfrlen ) then
 
              nrow   = lson - jc + 1
              len    = locbfr
              temp   = max ( ( .5 + nrow ) ** 2 - 2. * len, 0. )
              ncol   = .5 + nrow - sqrt ( temp )
              bfrlen = nrow * ncol - ncol * ( ncol - 1 ) / 2
 
              call xdslw1 ( wafil2, 2, idummy, idummy, ocbufr, 
     1                      iopos2, bfrlen, ierr)
              if ( ierr .ne. 0 ) return
 
              iopos2 = iopos2 + bfrlen
              bfrpos = 1
 
          end if
 
          j = sonind(jc)
          parbeg = lpar*(j-1) - (j*(j-1))/2
cdir$ ivdep
          do ic = jc,lson
              kpar         = parbeg + sonind(ic)
              parstk(kpar) = parstk(kpar) + ocbufr(bfrpos)
              bfrpos       = bfrpos + 1
          enddo
      enddo
 
c  =====================================================================
 
      return
      end
