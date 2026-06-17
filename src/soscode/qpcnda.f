      subroutine  QPCNDA   ( n, ldT, T, sminT, smaxT, sigma, work,
     1                       lwork )

c     ==================================================================
c     ==================================================================
c     ====  QPCNDA /                                                ====
c     ====  qpcnda -- compute condition number of T                 ====
c     ==================================================================
c     ==================================================================

      integer           n, ldT, lwork

      double precision  sminT, smaxT

      double precision  T (ldT, n), sigma (n), work (lwork)

c     ... local variables

      integer           i, ierror, j, job, ldV, nerror

      parameter       ( job = 00 )

      double precision  dummy, zero

      parameter        ( zero = 0.0d0 )

c     ==================================================================

c     ... last modification  26-March-1996
      
c     ==================================================================

c     ... make matrix explicitly upper triangle
       
      do j = 1, n-1
         do i = j+1, n
            T (i,j) = zero
         enddo
      enddo

c     ... brute force compute SVD of T

      ldV = n

      call hdsvdd ( T     , ldT   , n     , n     , job   , ldV   ,
     1              work  , lwork , nerror, sigma , dummy , dummy ,
     2              ierror )

      if  ( ierror .ne. 0 )   then
c
c        condition number of active constraint matrix cannot be 
c        computed --- set it to (-one)
c
         sigma(1) = 1.d0
         sigma(n) = -1.d0

      endif

      sminT = sigma (n)
      smaxT = sigma (1)

      return

c     end of  QPCNDA (qpcnda)
      
      end
