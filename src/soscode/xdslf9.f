      subroutine xdslf9 ( n    , x    , index , y      )
c
c  purpose -- to scatter/add in a vector into another
c
c  created            -- 22-jun-87, cca
c  last modifications -- 22-jun-87, cca
c
c  input variables --
c
c      n      -- length of the operation
c      x      -- vector to be scatter/added into y
c      index  -- index vector showing locations to receive entries
c
c  output variables --
c
c      y      -- vector to receive entries from x
c
c  subprograms called --
c
c      none
c
c  =====================================================================
 
      integer           n, index(*)
      double precision  x(*), y(*)
 
      integer           i     , j
 
c  =====================================================================
 
cdir$ ivdep
      do i = 1,n
          j    = index(i)
          y(j) = y(j) + x(i)
      enddo
 
c  =====================================================================
 
      return
      end
