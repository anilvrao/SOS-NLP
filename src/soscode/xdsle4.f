      subroutine xdsle4 ( zpcntl, npcntl, pnlcol, pnlrow, panel ,
     1                    pvtblk, inrtia, fctops, slvops, error  )
 
c
c  purpose -- perform block elimination step for the current panel of the
c             front, without pivoting.
c 
c             note that the panel has the columns of the front stored
c             as columns.
c
c             this routine performs rank three update on the
c             frontal matrix using horizontal unrolling
c             of degree two.
c
c  created            -- 07-feb-97, rgg
c  last modifications -- 28-oct-97, rgg, added zero and negative pivot
c                                        controls
c
c                        31-jan-01, jgl  error handling modified
c  input variables --
c
c      zpcntl -- array for zero pivot control
c      npcntl -- array for negative pivot control
c      pnlcol -- number of columns in panel
c      pnlrow -- number of rows    in panel
c      panel  -- rectangular the symmetric frontal matrix in packed storage
c      inrtia -- matrix inertia
c      fctops -- factor operation count
c      slvops -- solve  operation count
c
c  output variable --
c
c      a      -- the columns of the factor and the update in
c                packed storage
c      pvtblk -- type of pivot vector
c      inrtia -- matrix inertia
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
 
      integer           pvtblk(*), inrtia(3)
 
      double precision  fctops, slvops
 
      double precision  zpcntl(*), npcntl(*), panel(pnlrow,pnlcol)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer           i     , j     , m     , ierr

      double precision  diag1 , diag2 , diag3 , fac11 , fac12 , 
     1                  fac21 , fac22 , fac31 , fac32 , fcol  ,
     2                  frow  , ops   , rdiag1, rdiag2, rdiag3
 
c  =====================================================================
 
      error  = 0
 
      fcol   = pnlcol
      frow   = pnlrow - pnlcol
 
      ops    = frow*fcol*fcol
     1       + fcol*(fcol+1)*(2*fcol+1)/6
      fctops = fctops + ops
 
      ops    = fcol + 4 * ( (fcol+1)*fcol/2 + fcol*frow )
      slvops = slvops + ops

c  ----------------------------------------------------------
c  loop through the columns belonging to the supernode by 2's
c  ----------------------------------------------------------
 
      do 50 m = 3, pnlcol, 3
 
c         -------------------------
c         compute the three columns
c         -------------------------
 
c  -------------------------------------------
c.debug
c     write(6,'("in xdsle4-m, pnlcol, pnlrow, panel(m-2,m-2) = "
c    1      3i8,1pd15.5)') m, pnlcol, pnlrow, panel(m-2,m-2)
c.debug
          call xdslep ( m-2, panel(m-2,m-2), 
     1                  pnlrow-m+2, panel(m-1,m-2), 1,
     2                  zpcntl, npcntl, inrtia, ierr )
          if ( ierr .ne. 0 )  then
             if  ( ( ierr .eq. -1 ) .or. ( ierr .eq. -2 ) )  then
                error = ierr
             else
                error = -99
             end if
             return
          end if

          diag1  = panel(m-2,m-2)
          rdiag1 = 1. / diag1
c  -------------------------------------------
          fac21          = panel(m-1,m-2)
          panel(m-1,m-2) = rdiag1*panel(m-1,m-2)
          panel(m-1,m-1) = panel(m-1,m-1) - fac21*panel(m-1,m-2)

c.debug
c     write(6,'("in xdsle4-m, pnlcol, pnlrow, panel(m-1,m-1) = "
c    1      3i8,1pd15.5)') m, pnlcol, pnlrow, panel(m-1,m-1)
c.debug
          call xdslep ( m-1, panel(m-1,m-1), 
     1                  pnlrow-m+1, panel(m,m-1), 1,
     2                  zpcntl, npcntl, inrtia, ierr )
          if ( ierr .ne. 0 )  then
             if  ( ( ierr .eq. -1 ) .or. ( ierr .eq. -2 ) )  then
                error = ierr
             else
                error = -99
             end if
             return
          end if

          diag2  = panel(m-1,m-1)
          rdiag2 = 1. / diag2
c  ------------------------------------------
          fac31        = panel(m,m-2)
          panel(m,m-2) = rdiag1*panel(m,m-2)
          fac32        = panel(m,m-1) - fac21*panel(m,m-2)
          panel(m,m-1) = rdiag2*fac32
          panel(m  ,m) = panel(m  ,m) - fac31*panel(m,m-2) 
     1                                - fac32*panel(m,m-1)

c.debug
c     write(6,'("in xdsle4-m, pnlcol, pnlrow, panel(m,m) = "
c    1      3i8,1pd15.5)') m, pnlcol, pnlrow, panel(m,m)
c.debug
          call xdslep ( m, panel(m,m), 
     1                  pnlrow-m, panel(m+1,m), 1,
     2                  zpcntl, npcntl, inrtia, ierr )
          if ( ierr .ne. 0 )  then
             if  ( ( ierr .eq. -1 ) .or. ( ierr .eq. -2 ) )  then
                error = ierr
             else
                error = -99
             end if
             return
          end if

          diag3  = panel(m,m)
          rdiag3 = 1. / diag3
c  -------------------------------------------
          pvtblk(m-2) = 1
          pvtblk(m-1) = 1
          pvtblk(m  ) = 1
c  -------------------------------------------
cdir$ ivdep
          do i = m+1, pnlrow
              panel(i,m-2) = rdiag1* panel(i,m-2)
              panel(i,m-1) = rdiag2*(panel(i,m-1) - fac21*panel(i,m-2))
              panel(i,m  ) = rdiag3*(panel(i,m  ) - fac31*panel(i,m-2) 
     1                                            - fac32*panel(i,m-1))
          enddo
 
c         -------------------------------------
c         update the remaining columns of panel
c         -------------------------------------
 
          do j = m+1, pnlcol-1, 2

              fac11 = diag1*panel(j  ,m-2)
              fac12 = diag1*panel(j+1,m-2)
              fac21 = diag2*panel(j  ,m-1)
              fac22 = diag2*panel(j+1,m-1)
              fac31 = diag3*panel(j  ,m  )
              fac32 = diag3*panel(j+1,m  )

              panel(j,j) = panel(j,j) - fac11*panel(j,m-2) 
     1                                - fac21*panel(j,m-1) 
     2                                - fac31*panel(j,m  )

cdir$ ivdep
              do i = j+1, pnlrow
                  panel(i,j  ) = panel(i,j  ) - fac11*panel(i,m-2)
     1                                        - fac21*panel(i,m-1) 
     2                                        - fac31*panel(i,m  )
                  panel(i,j+1) = panel(i,j+1) - fac12*panel(i,m-2)
     1                                        - fac22*panel(i,m-1) 
     2                                        - fac32*panel(i,m  )
              enddo

          enddo
c.debug
c     write(6,'("after 30 - j, pnlcol = ", 2i8)') j, pnlcol
c.debug

          if ( j .eq. pnlcol ) then

              fac11 = diag1*panel(j,m-2)
              fac21 = diag2*panel(j,m-1)
              fac31 = diag3*panel(j,m  )

cdir$ ivdep
              do i = j, pnlrow
                  panel(i,j  ) = panel(i,j  ) - fac11*panel(i,m-2)
     1                                        - fac21*panel(i,m-1) 
     2                                        - fac31*panel(i,m  )
              enddo

          endif

c.debug
c     call xdslp5('panel after 40', pnlrow*pnlcol, panel, 6 )
c.debug
 
   50 continue
 
      if ( m - 1 .eq. pnlcol ) then
 
c         -----------------------
c         compute the two columns
c         -----------------------
 
c  -------------------------------------------
c.debug
c     write(6,'("in xdsle4-m, pnlcol, pnlrow, panel(m-2,m-2) = "
c    1      3i8,1pd15.5)') m, pnlcol, pnlrow, panel(m-2,m-2)
c.debug
          call xdslep ( m-2, panel(m-2,m-2), 
     1                  pnlrow-m+2, panel(m-1,m-2), 1,
     2                  zpcntl, npcntl, inrtia, ierr )
          if ( ierr .ne. 0 )  then
             if  ( ( ierr .eq. -1 ) .or. ( ierr .eq. -2 ) )  then
                error = ierr
             else
                error = -99
             end if
             return
          end if

          diag1  = panel(m-2,m-2)
          rdiag1 = 1. / diag1
c  -------------------------------------------
          fac21          = panel(m-1,m-2)
          panel(m-1,m-2) = rdiag1*panel(m-1,m-2)
          panel(m-1,m-1) = panel(m-1,m-1) - fac21*panel(m-1,m-2)

c.debug
c     write(6,'("in xdsle4-m, pnlcol, pnlrow, panel(m-1,m-1) = "
c    1      3i8,1pd15.5)') m, pnlcol, pnlrow, panel(m-1,m-1)
c.debug
cdbg   looks like reference to panel(m,m-1) violates bounds, jtb.
cdbg      replace this call
cdbg          call xdslep ( m-1, panel(m-1,m-1), 
cdbg     1                  pnlrow-m+1, panel(m,m-1), 1,
cdbg     2                  zpcntl, npcntl, inrtia, ierr )
cdbg      with this call
          if(zpcntl(1).eq.0.) then
            if(panel(m-1,m-1).eq.0.) ierr = -1
          else
c
            call xdslep ( m-1, panel(m-1,m-1), 
     1                  pnlrow-m+1, panel(m,m-1), 1,
     2                  zpcntl, npcntl, inrtia, ierr )
c
          endif
c
          if ( ierr .ne. 0 )  then
             if  ( ( ierr .eq. -1 ) .or. ( ierr .eq. -2 ) )  then
                error = ierr
             else
                error = -99
             end if
             return
          end if

          diag2  = panel(m-1,m-1)
          rdiag2 = 1. / diag2
c  -------------------------------------------
          pvtblk(m-2) = 1
          pvtblk(m-1) = 1
c  -------------------------------------------
cdir$ ivdep
          do i = m, pnlrow
              panel(i,m-2) = rdiag1* panel(i,m-2)
              panel(i,m-1) = rdiag2*(panel(i,m-1) - fac21*panel(i,m-2))
          enddo
 
      elseif ( m - 2 .eq. pnlcol ) then
 
c         -----------------------
c         compute the last column
c         -----------------------
 
c  -------------------------------------------
c.debug
c     write(6,'("in xdsle4-m, pnlcol, pnlrow, panel(m-2,m-2) = "
c    1      3i8,1pd15.5)') m, pnlcol, pnlrow, panel(m-2,m-2)
c.debug
          call xdslep ( m-2, panel(m-2,m-2), 
     1                  pnlrow-m+2, panel(m-1,m-2), 1,
     2                  zpcntl, npcntl, inrtia, ierr )
          if ( ierr .ne. 0 )  then
             if  ( ( ierr .eq. -1 ) .or. ( ierr .eq. -2 ) )  then
                error = ierr
             else
                error = -99
             end if
             return
          end if

          diag1  = panel(m-2,m-2)
          rdiag1 = 1. / diag1
c  -------------------------------------------
          pvtblk(m-2) = 1
c  -------------------------------------------
cdir$ ivdep
          do i = m-1, pnlrow
              panel(i,m-2) = rdiag1* panel(i,m-2)
          enddo
 
      endif

c.debug
c     call xdslp5('panel before return', pnlrow*pnlcol, panel, 6 )
c.debug

c  =====================================================================
 
      return
      end
