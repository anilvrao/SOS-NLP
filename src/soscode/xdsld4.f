      subroutine xdsld4( zpcntl, npcntl, pnlcol, pnlrow, panell,
     1                   panelu, pvtblk, fctops, slvops, error  )
 
c
c  purpose -- perform block elimination step for the current panel of the
c             front, without pivoting.
c             unsymmetric version
c 
c             note that the panel has the columns of the front stored
c             as columns.
c
c  created            -- 21-may-98, rgg, derived from xdsle4 and hdsluf
c  last modifications    31-jan-01, jgl  error handling modified
c
c  input variables --
c
c      zpcntl -- array for zero pivot control
c      npcntl -- array for negative pivot control
c      pnlcol -- number of columns in panel
c      pnlrow -- number of rows    in panel
c      panell -- panel of the lower tri. of the frontal matrix
c      panelu -- panel of the lower tri. of the frontal matrix
c      fctops -- factor operation count
c      slvops -- solve  operation count
c
c  output variable --
c
c      a      -- the columns of the factor and the update in
c                packed storage
c      pvtblk -- type of pivot vector
c      fctops -- factor operation count
c      slvops -- solve  operation count
c      error  -- error return,
c                if error =  0, success,
c                         =  -1 exact zero diagonal encountered as pivot
c                         =  -2 negative diagonal encountered while
c                               prohibited by user
c                         = -99 undocumented error in lower level routine
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           pnlcol, pnlrow, error
 
      integer           pvtblk(*)
 
      double precision  fctops, slvops
 
      double precision  zpcntl(*), npcntl(*), 
     1                  panell(pnlrow,pnlcol),
     2                  panelu(pnlrow,pnlcol)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer           i     , j     , m     , ierr
 
      integer           inrtia(3)

      double precision  diag1 , 
     1                  fcol  , frow  ,
     2                  ops   , rdiag1
 
c  =====================================================================
 
      error  = 0

      fcol   = pnlcol
      frow   = pnlrow - pnlcol

      ops    = frow*fcol*fcol
     1       + fcol*(fcol+1)*(2*fcol+1)/6
      fctops = fctops + 2 * ops

      ops    = fcol + 4 * ( (fcol+1)*fcol/2 + fcol*frow )
      slvops = slvops + ops
 
c     ------------------------
c     loop through the columns
c     ------------------------
 
      do 50 m = 1, pnlcol
c.debug
c     write(6,'("in xdsld4 do 50 loop - m = ", i8)') m
c     call xdslp5('panell', pnlcol*pnlrow, panell, 6 )
c     call xdslp5('panelu', pnlcol*pnlrow, panelu, 6 )
c.debug
 
c         ----------------------------------------
c         scale the m-th row of the upper triangle
c         ----------------------------------------
 
          call xdslep ( m, panell(m,m), pnlrow-m, panell(m+1,m),
     1                  1, zpcntl, npcntl, inrtia, ierr )
          if ( ierr .ne. 0 )  then
             if  ( ( ierr .eq. -1 ) .or. ( ierr .eq. -2 ) )  then
                error = ierr
             else
                error = -99
             end if
             return
          end if

          diag1  = panell(m,m)
          rdiag1 = 1. / diag1
 
          pvtblk(m) = 1
 
cdir$ ivdep
          do i = m+1, pnlrow
              panelu(i,m) = rdiag1*panelu(i,m)
          enddo
 
c         ----------------------------
c         update the remaining columns
c         ----------------------------
 
          do i = m+1, pnlcol

              panell(i,i) = panell(i,i)
     1                    - panell(i,m) * panelu(i,m)

cdir$ ivdep
              do j = i+1, pnlrow 
                  panell(j,i) = panell(j,i)  
     1                        - panell(j,m) * panelu(i,m)
                  panelu(j,i) = panelu(j,i)  
     1                        - panell(i,m) * panelu(j,m)
              enddo
          enddo
 
c         --------------------------
c         ... scale m-th column of l
c         --------------------------
 
cdir$ ivdep
          do i = m+1, pnlrow
              panell(i,m) = rdiag1*panell(i,m)
          enddo
 
   50 continue
 
      return
 
c  =====================================================================
 
      end
