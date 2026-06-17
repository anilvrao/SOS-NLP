      double precision function   DDIVFN   ( a, b, fail )

c     ==================================================================
c     ==================================================================
c     ====  DDIVFN /  f06blf /                                      ====
c     ====  ddiv   -- cautious division with overflow check         ====
c     ==================================================================
c     ==================================================================


      double precision                  a, b
      
      logical                           fail

c     ==================================================================

c
c     ddiv / DDIVFN  returns the value div given by
c
c     div = ( a/b                 if a/b does not overflow,
c           (
c           ( 0.0                 if a .eq. 0.0,
c           (
c           ( sign( a/b )*flmax   if a .ne. 0.0  and a/b would overflow,
c
c  where  flmax  is a large value, via the function name. In addition if
c  a/b would overflow then  fail is returned as true, otherwise  fail is
c  returned as false.
c
c  Note that when  a and b  are both zero, fail is returned as true, but
c  div  is returned as  0.0. In all other cases of overflow  div is such
c  that  abs( div ) = flmax.
c
c  When  b = 0  then  sign( a/b )  is taken as  sign( a ).
c
c  Nag Fortran 77 O( 1 ) basic linear algebra routine.
c
c  -- Written on 26-October-1982.
c     Sven Hammarling, Nag Central Office.
c
c     ==================================================================

      double precision      absb, div, flmax, flmin
      
      logical               first
      
      intrinsic             abs, sign

      save                  first, flmin, flmax

      double precision      one , zero
      
      parameter           ( one = 1.0d0, zero = 0.0d0 )

      double precision      hdmcon

      external              hdmcon
      
      data                  first/ .true. /

c     ==================================================================

      if  ( a .eq. zero )  then
         div = zero
         if  ( b .eq. zero )  then
            fail = .true.
         else
            fail = .false.
         end if
      else

         if  ( first )  then
            first  = .false.
            flmin  = hdmcon (4)
            flmax  = hdmcon (2)
         end if

         if  ( b .eq. zero )  then
            div  =  sign( flmax, a )
            fail = .true.
         else
            absb = abs( b )
            if  ( absb .ge. one )  then
               fail = .false.
               if  ( abs( a ) .ge. absb*flmin )  then
                  div = a/b
               else
                  div = zero
               end if
            else
               if  ( abs( a ) .le. absb*flmax )  then
                  fail = .false.
                  div  =  a/b
               else
                  fail = .true.
                  div  = flmax
                  if  ( ( ( a .lt .zero )  .and.  ( b .gt .zero ) )
     1                  .or.
     2                  ( ( a .gt .zero )  .and.  ( b .lt .zero ) ) )
     3            then
                     div = -div
                  endif
               end if
            end if
         end if
      end if

      DDIVFN   = div

c     end of DDIVFN /  f06blf / ddiv
      
      end
