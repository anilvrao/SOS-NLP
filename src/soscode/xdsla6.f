      subroutine xdsla6 ( alpha , n     , x     , index , 
     1                    poffst, lpan  , y     )
 
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------

      integer             n     , poffst, lpan

      integer             index(*)

      double precision    alpha

      double precision    x(*)  , y(*)
 
c     -------------------
c     ... local variables
c     -------------------

      integer             i     , ibgn  , iend  , j
 
c  =====================================================================

c     -------------------------------------------------------------
c     ... find active section of y (fron poffst+1 to poffst + lpan)
c     -------------------------------------------------------------
c.debug
c     write(6,'("in xdsla6 before 100 - poffst, lpan = ", 2i8)')
c    1                                   poffst, lpan
c     call xislp3 ( 'index', n, index, 6 )
c.debug

      do i = 1, n
          if ( index(i) .gt. poffst ) then
              ibgn = i
              go to 200
          end if
      enddo

      return

  200 continue
c.debug
c     write(6,'("in xdsla6 before 300 - ibgn         = ", 2i8)')
c    1                                   ibgn        
c.debug

      do i = n, ibgn, -1
          if ( index(i) .le. poffst + lpan ) then
              iend = i
              go to 400
          end if
      enddo

      return
 
c  =====================================================================

c     ---------------------------------------
c     ... add contributions to active section
c         branch if alpha = 1.0 or not
c     ---------------------------------------

  400 continue
      if ( alpha .eq. 1.0d0 ) then

c.debug
c     write(6,'("in xdsla6 before 500 - ibgn, iend = ", 2i8)')
c    1                                   ibgn, iend
c.debug

          do i = ibgn, iend
              j    = index(i) - poffst
              y(j) = y(j) + x(i)
          enddo

      else

c.debug
c     write(6,'("in xdsla6 before 600 - ibgn, iend = ", 2i8)')
c    1                                   ibgn, iend
c.debug

          do i = ibgn, iend
              j    = index(i) - poffst
              y(j) = y(j) + alpha * x(i)
          enddo

      end if
 
c  =====================================================================
 
      return
      end
