      subroutine  xislog ( n, perm, invp )
 
c
c     forms the inverse permutation vector, invp, given the
c     permutation vector, perm.
c
c     roger grimes, january, 1989.
c
c---------------------------------------------------------------------
 
      integer     n
      integer     perm(*), invp(*)
 
      integer     i, j
 
c---------------------------------------------------------------------
 
      if ( n .le. 0 ) return
 
      do i = 1,n
          j         = perm ( i )
          invp (j ) = i
      enddo
 
c---------------------------------------------------------------------
 
      return
      end
