      subroutine  xislex  ( n, length, clsfit, lnused )

c
c     purpose -- determine how many leading columns of a triangular
c                matrix of order  n  will fit into a space of 
c                predetermined length
c
c     created            -- 15-apr-03, jgl
c
c    input variables --
c
c        n      -- order of the triangular matrix (> 0)
c        length -- amount of space available (> 0)
c
c    output variables --
c
c        clsfit -- the number of columns that fit, in the range (0,n)
c        lnused -- space actually occupied by these columns, zero if
c                  no columns fit       
c
c  =====================================================================
 
c     -------------
c     ... arguments
c     -------------
 
      integer           n, length, clsfit,lnused
 
c  =====================================================================
 
      if ( n .le. length ) then

         if  ( ( n*(n+1) ) / 2 .le. length )  then

c           -------------------------
c           ... the whole matrix fits
c           -------------------------

            clsfit = n
            lnused = ( n*(n+1) ) / 2

         else

c           -------------------------------------------------------------
c           ... Otherwise, a block of  k  columns will occupy  
c                   nk - ((k)(k-1)/2) space
c               So, solve the quadratic equation
c                    k**2 - (2n+1) k + 2 length = 0
c               for its smaller root and then verify that floating 
c               point arithmetic didn't slightly miss the discrete answer
c           -------------------------------------------------------------

            clsfit = ( (2.0*n + 1.0) - 
     1                  sqrt ( (2.0*n + 1.0)**2 - 8.0*length ) )
     2               * 0.5

            lnused = n*clsfit - ( (clsfit)*(clsfit-1) ) / 2

            if  ( lnused .gt. length )  then

               clsfit = clsfit - 1
               lnused = lnused - (n - clsfit)

            elseif  ( lnused + n - clsfit  .le. length )  then

               lnused = lnused + (n - clsfit)
               clsfit = clsfit + 1

            end if

         endif

      else

c        --------------------------------
c        ... n is so large no columns fit
c        --------------------------------

         clsfit = 0
         lnused = 0

      end if

      return

      end
