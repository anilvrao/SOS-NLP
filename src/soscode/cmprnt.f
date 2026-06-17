      subroutine   CMPRNT   ( msglvl, iPrint, n, nclin, nctotl, bigbnd,
     1                        named, names, istate,
     2                        bl, bu, clamda, featol, r )
c     ==================================================================
c     ==================================================================
c     ====  cmprnt / CMPRNT -- common iteration output routine     =====
c     ==================================================================
c     ==================================================================

      integer            msglvl, iPrint, n, nclin, nctotl

      double precision   bigbnd

      character(len=16)  names(*)
      
      logical            named
      
      integer            istate (nctotl)
      
      double precision   bl (nctotl), bu (nctotl),
     1                   clamda (nctotl), featol (nctotl), r (nctotl)

c     ==================================================================
c     derived from qpopt version 1.0
c     last modification -- 25-March-1996
c
c         Original Fortran 77 version written  October 1984.
c        This version of  cmprnt dated  11-May-95.
c
c     cmprnt/CMPRNT  prints r(x) (x,  A*x and c(x)), the bounds, the
c     multipliers, and the slacks (distance to the nearer bound).
c
c     ==================================================================

      integer            is, j, nplin, number

      double precision   b1, b2, rj, slk, slk1, slk2, tol, wlam
      
      character(len=1)   key
      
      character(len=8)   name
      
      character(len=102) line

      double precision   zero
      
      parameter        ( zero  = 0.0d0 )
      
      character(len=2)   lstate(-2:4), state

      data               lstate(-2) / '--' /, lstate(-1) / '++' /
      data               lstate( 0) / 'FR' /, lstate( 1) / 'LL' /
      data               lstate( 2) / 'UL' /, lstate( 3) / 'EQ' /
      data               lstate( 4) / 'TF' /

c     ==================================================================
      
      if  ( iPrint .eq. 0 .or.  (msglvl .lt. 10  .and.  msglvl .ne. 1) )
     1then
         return
      endif

      write (iPrint, 1000) 'Variable       '
      name   = 'variable'
      nplin  = n + nclin

      do j = 1, nctotl
         b1     = bl(j)
         b2     = bu(j)
         wlam   = clamda(j)
         rj     = r(j)

         if  ( j .le. n )  then
            number = j
         else
     1   if  ( j .le. nplin )  then
            number = j - n
            if  ( number .eq. 1 )  then
               write (iPrint, 1000) 'Linear constrnt'
               name = 'lincon  '
            end if
         else
c           (else clause cannot occur in qp setting)
            number = j - nplin
            if  ( number .eq. 1 )  then
               write (iPrint, 1000) 'Nonlin constrnt'
               name = 'nlncon  '
            end if
         end if

c        Print a line for the jth variable or constraint.
c        ------------------------------------------------
         is     = istate(j)
         state  = lstate(is)
         tol    = featol(j)
   
         slk1   = rj - b1
         slk2   = b2 - rj
         if  ( abs(slk1) .lt. abs(slk2) )  then
            slk = slk1
            if  ( b1 .le. - bigbnd) slk = slk2 
         else
            slk = slk2
            if  ( b2 .ge.   bigbnd) slk = slk1 
         end if

c        Flag infeasibilities, primal and dual degeneracies, 
c        and active QP constraints that are loose in NP.
c      
         key    = ' ' 
         if  ( slk1 .lt. -tol  .or.       slk2  .lt. -tol)  key = 'I'
         if  ( is   .eq.  0    .and.  abs(slk ) .le.  tol)  key = 'D'
         if  ( is   .ge.  1    .and.  abs(wlam) .le.  tol)  key = 'A'

         write (line, 2000) name, number, key, state, 
     $                     rj, b1, b2, wlam, slk

c        Reset special cases:
c           Infinite bounds
c           Zero bounds
c           Lagrange multipliers for inactive constraints
c           Lagrange multipliers for infinite bounds
c           Infinite slacks
c           Zero slacks

         if  (       named      )  line( 2: 17) = names(j)
         if  ( b1  .le. - bigbnd)  line(39: 54) = '      None      '
         if  ( b2  .ge.   bigbnd)  line(55: 70) = '      None      '
         if  ( b1  .eq.   zero  )  line(39: 54) = '        .       '
         if  ( b2  .eq.   zero  )  line(55: 70) = '        .       '
         if  ( is  .eq.   0       .or.    
     1         wlam.eq.   zero  )  then
                                   line(71: 86) = '        .       '
         endif
         if  ( b1  .le. - bigbnd  .and. 
     1         b2  .ge.   bigbnd )  then
                                   line(71: 86) = '                '
                                   line(87:102) = '                '
         end if
         if  ( slk .eq.   zero  )  line(87:102) = '        .       '

         write (iPrint, '(a)') line
      enddo

      return

 1000 format(//  1x,  a15, 2x, 'State', 6x, 'Value',
     $           7x, 'Lower bound', 5x, 'Upper bound',
     $           3x, 'Lagr multiplier', 4x, '   Slack' / )
 2000 format( 1x, a8, i6, 3x, a1, 1x, a2, 4g16.7, g16.4 )

c     end of CMPRNT (cmprnt)
      end
