      integer function   xdslni   (n)
 
c
c     ==================================================================
c     ====  xdslni  --  length of long integer vector               ====
c     ==================================================================
c     ==================================================================
c
c     created       07-jun-89  -- jgl --
c
c     purpose - compute length, in floating point words, of a
c               vector of  n  long integers
c
c     ==================================================================
 
c     -------------
c     ... parameter
c     -------------
 
      integer           n
 
c     ==================================================================
 
      xdslni = ( n + 1 ) / 2
 
      return
 
      end
