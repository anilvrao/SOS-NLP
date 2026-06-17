      integer function   xdslil   (n)
c
c
c     ==================================================================
c     ====  xdslil  --  number of integers in a floating point vector ==
c     ==================================================================
c     ==================================================================
c
c     created       13-mar-91  -- rgg --
c
c     purpose - compute length the number of long integers in  n
c               floating point words
c
c     ==================================================================
 
c     -------------
c     ... parameter
c     -------------
 
      integer           n
 
c     ==================================================================
 
      xdslil = 2 * n
 
      return
 
      end
