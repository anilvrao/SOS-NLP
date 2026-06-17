      subroutine  xdslwb ( mtxunt, qascii, mtxtyp, mtxstr, neqns , 
     1                     perm  , diag  , colstr, rowind, mtxval )
c
c     ==================================================================
c     ====  xdslwb -- writes matrices in either ascii or binary     ====
c     ====            format to i/o unit, in memory version         ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xdslwb writes a matrix stored in memory to i/o unit mtxunt
c
c     created         20-apr-99   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     mtxunt      i   i/o unit to use.
c     qascii      l   ascii or binary logical flag
c     mtxtyp      i   1 for symmetric, 2 for unsymmetric
c     mtxstr      i   1 for general, 3 for diagonal
c     neqns       i   order of matrix
c     perm        i   new-to-old permutation
c     diag        d   diagonal array
c     colstr      i   column start array
c     rowind      i   row indices array
c     mtxval      d   matrix off-diagonal value array
c
c     output arguments
c     ----------------
c
c     error       i   error flag
c                     = -1    i/o error
c
c     ==================================================================
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------

      logical              qascii

      integer              mtxunt, mtxtyp, mtxstr, neqns 

      integer              perm  (*),      colstr(*),      rowind(*)

      double precision     diag  (*),      mtxval(*)
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------

      integer              inew  ,  iold  , jchv  , jnew  , jold  ,
     1                     k     , kbgn  , kend  

      double precision     val
 
c     ==================================================================

      do inew = 1, neqns

          iold = perm(inew)

          if ( qascii ) then
              write(mtxunt,'(2i10,d25.17)') iold, iold, diag(inew)
          else
              write(mtxunt                ) iold, iold, diag(inew)
          end if

      enddo

      if ( mtxstr .ne. 1 ) return

      do 300 jchv = 1, neqns

          kbgn  = colstr(jchv)
          kend  = colstr(jchv+1) - 1

          do k = kbgn, kend

              if ( mtxtyp .eq. 1 ) then

                  jnew = jchv
                  inew = rowind(k)

                  jold = min ( perm(jnew), perm(inew) )
                  iold = max ( perm(jnew), perm(inew) ) 

                  val  = mtxval(k)

              else

                  if ( rowind(k) .le. 0 ) then
                      jnew = jchv
                      inew = jchv - rowind(k)
                  else
                      inew = jchv
                      jnew = jchv + rowind(k)
                  end if

                  jold = perm(jnew)
                  iold = perm(inew)

                  val  = mtxval(k)

              end if

              if ( qascii ) then
                  write(mtxunt,'(2i10,d25.17)') iold, jold, val
              else
                  write(mtxunt                ) iold, jold, val
              end if

          enddo

  300 continue
 
c     ==================================================================

      return
      end
