      subroutine  xdslmv(n,x,loc1,loc2)
 
c
c     moves a subset of a vector to another location in the same
c     vector
c
c     ------------------------------------------------------------------

      integer              n, loc1, loc2

      double precision     x(*)

      integer              i, i1, i2, nm1
 
c     ------------------------------------------------------------------
 
      if ( n .le. 0 ) return

      nm1 = n - 1

      if ( ( loc2 .gt. loc1 + nm1 ) .or. ( loc1 .gt. loc2 + nm1 ) ) then

c         -----------------------------------------------
c         ... the two subsets of the array do not overlap.
c             it is safe to allowed an compiler to optimize
c             this loop.  
c         -----------------------------------------------

          i1 = loc1
          i2 = loc2

cdir$ ivdep 
          do i = 1, n
              x(i2) = x(i1)
              i1    = i1 + 1
              i2    = i2 + 1
          enddo

      else if ( loc2 .gt. loc1 ) then

c         -------------------------------------------------
c         ... the two subsets of the array overlap and
c             the data is being moved to a greater location.
c             data must be moved in reverse manner and the
c             loop can not be optimized.
c         -------------------------------------------------

          i1 = loc1 + nm1
          i2 = loc2 + nm1

          do i = n, 1, -1
              x(i2) = x(i1)
              i1    = i1 - 1
              i2    = i2 - 1
          enddo

      else if ( loc1 .gt. loc2 ) then

c         -------------------------------------------------
c         ... the two subsets of the array overlap and
c             the data is being moved to a lesser location.
c             note that the loop is the same as the 
c             non-overlapping case but it is not safe to
c             optimize this loop.
c         -------------------------------------------------

          i1 = loc1
          i2 = loc2

          do i = 1, n
              x(i2) = x(i1)
              i1    = i1 + 1
              i2    = i2 + 1
          enddo

c     else 

c         ----------------------------------------------
c         ... the two subsets of the array are identical
c             so nothing is moved
c         ----------------------------------------------

      end if
 
c     ------------------------------------------------------------------
 
      return
      end
