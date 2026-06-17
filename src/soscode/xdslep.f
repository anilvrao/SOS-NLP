      subroutine xdslep ( jcol  , diag  , lofdia, offdia, incrmt,
     1                    zpcntl, npcntl, inrtia, error  )
 
c
c  purpose -- perform pivot control tests and inertia computations.
c 
c  created            -- 28-oct-97, rgg
c  last modifications -- 31-jan-01, jgl -- error codes changed
c                        18-mar-02, jgl -- protected against special
c                                          case -- last column of root
c
c  input variables --
c
c      jcol   -- column in the current front that is being tested
c      lofdia -- length of the offdiagonal array
c      offdia -- array of offdiagonal entries
c      incrmt -- increment used for the offdia array
c
c  input/output variables --
c
c      diag   -- diagonal entry being tested
c      zpcntl -- array for zero pivot control
c                zpcntl(1) = 0 no control
c                       .ne. 0 fudge small pivots by this amount
c                zpcntl(2) = number of near zero pivots fudged
c                zpcntl(3) = first near zero pivot fudged
c                zpcntl(4) = minimum fudging of near zero pivots
c                zpcntl(5) = maximum fudging of near zero pivots 
c      npcntl -- array for negative pivot control
c                npcntl(1) = 0 no control
c                          = 1 monitor negative pivots
c                          = 2 abort on a negative pivots
c                npcntl(2) = number of negative pivots found 
c                npcntl(3) = first negative pivot found 
c                npcntl(4) = minimum negative pivots
c                npcntl(5) = maximum negative pivots 
c      inrtia -- matrix inertia
c
c  output variable --
c
c      error   -- error return,
c                if error =  0, success,
c                         = -1, zero pivot detected
c                         = -2, negative pivot encountered with 
c                               npcntl(1) .eq. 2
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           jcol  , lofdia, incrmt, error
 
      integer           inrtia(3)
 
      double precision  diag  
 
      double precision  offdia(*), zpcntl(*), npcntl(*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer           k

      double precision  big   , fudge , olddia, one, damax

      external  damax
 
c  =====================================================================

      error = 0
      fudge = zpcntl(1)
      one   = 1.0
c.debug
c     write(6,'("xdslep-zpcntl(1), npcntl(1), diag = ", 1p3d15.5)')
c    1                   zpcntl(1), npcntl(1), diag 
c.debug
 
c  =====================================================================

c     ----------------------------
c     ... control near zero pivots
c     ----------------------------

      if ( zpcntl(1) .eq. 0. ) then

c         ----------------------------------
c         ... no control, just test for zero
c         ----------------------------------

          if ( diag .eq. 0. ) then
              error = -1
              return
          end if

      else

          if ( abs ( diag ) .lt. fudge ) then

c             -----------------------------------------
c             ... near zero pivot that should be fudged
c             -----------------------------------------
c.debug
c     write(6,'("in xdslep adjusting zero pivot")')
c     write(6,'("diag  , fudge , lofdia, incrmt, jcol = ", 
c    1    1p2d15.5, 3i8)') diag, fudge, lofdia, incrmt, jcol
c     call xdslp5 ( 'offdia', lofdia*incrmt, offdia, 6 )
c.debug
              olddia = diag

              if ( lofdia .ge. 1 ) then

                 big   = damax ( lofdia, offdia, incrmt )
                 big   = max ( big, one )
                 diag  = fudge * big
c.debug
c                write(6,'("diag  , big   , k              = ", 
c    1                     1p2d15.5, 2i8)') diag, big   , k             
c.debug
              else

                 diag  = fudge
c.debug
c                write(6,'("diag  = ", 1p2d15.5, 2i8)') diag
c.debug

              end if

              zpcntl(2) = zpcntl(2) + 1.

              if ( zpcntl(2) .eq. 1. ) then

                  zpcntl(3) = jcol
                  zpcntl(4) = abs ( diag - olddia )
                  zpcntl(5) = abs ( diag - olddia )

              else

                  zpcntl(4) = min ( zpcntl(4), abs ( diag - olddia ) )
                  zpcntl(5) = max ( zpcntl(5), abs ( diag - olddia ) )

              end if

          end if

      end if
 
c  =====================================================================

c     ---------------------------
c     ... control negative pivots
c     ---------------------------

      if ( npcntl(1) .gt. 0. ) then

          if ( diag .lt. 0. ) then

              npcntl(2) = npcntl(2) + 1.
c.debug
c     write(6,'("in xdslep - negative - jcol = ", i8)') jcol
c.debug

              if ( npcntl(2) .eq. 1. ) then

                  npcntl(3) = jcol
                  npcntl(4) = diag
                  npcntl(5) = diag 

                  if ( npcntl(1) .eq. 2 ) then
                      error = -2
                      return
                  end if

              else

                  npcntl(4) = min ( npcntl(4), diag ) 
                  npcntl(5) = max ( npcntl(5), diag  ) 

              end if

          end if

      end if
 
c  =====================================================================

c     ---------------------------
c     ... track the inertia count
c     ---------------------------

      if ( diag .gt. 0. ) then
          inrtia(1) = inrtia(1) + 1
      else
          inrtia(2) = inrtia(2) + 1
      end if
 
c  =====================================================================
 
      return
      end
