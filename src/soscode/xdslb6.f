      subroutine xdslb6 ( alpha , n     , x     , index , 
     1                    poffst, lpan  , matsiz, panell, panelu )
 
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------

      integer             n     , poffst, lpan, matsiz

      integer             index(*)

      double precision    alpha

      double precision    x(*)  , panell(*), panelu(*)
 
c     -------------------
c     ... local variables
c     -------------------

      integer             i     , j, uoffst
 
c  =====================================================================

      uoffst = poffst + matsiz

c     ---------------------------------------
c         branch if alpha = 1.0 or not
c     ---------------------------------------

      if ( alpha .eq. 1.0d0 ) then

c.debug
c     write(6,'("                       poffst, uoffst = ", 2i8)')
c    1                                   poffst, uoffst
c     write(6,'("                       matsiz, lpan   = ", 2i8)')
c    1                                   matsiz, lpan  
c.debug

          do i = 1, n
              if ( index(i) .lt. poffst ) cycle
              if ( index(i) - poffst .gt. lpan ) cycle
              j    = index(i) - poffst
c.debug
c             write(6,'("in 500 - i, j = ", 2i8)') i, j
c.debug
              panell(j) = panell(j) + x(i)
          enddo

          do i = 1, n
              if ( index(i) .lt. uoffst ) cycle
              if ( index(i) - uoffst .gt. lpan ) cycle
              j    = index(i) - uoffst
c.debug
c             write(6,'("in 510 - i, j = ", 2i8)') i, j
c.debug
              panelu(j) = panelu(j) + x(i)
          enddo

      else

          do i = 1, n
              if ( index(i) .lt. poffst ) cycle
              if ( index(i) - poffst .gt. lpan ) cycle
              j    = index(i) - poffst
              panell(j) = panell(j) + alpha * x(i)
          enddo

          do i = 1, n
              if ( index(i) .lt. uoffst ) cycle
              if ( index(i) - uoffst .gt. lpan ) cycle
              j    = index(i) - uoffst
              panelu(j) = panelu(j) + alpha * x(i)
          enddo

      end if
 
c  =====================================================================
 
      return
      end
