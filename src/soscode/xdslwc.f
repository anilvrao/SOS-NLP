      subroutine  xdslwc ( mtxunt, qascii, mtxtyp, mtxstr, neqns , 
     1                     nsuper, sqfile, perm  , xsup  ,
     1                     diag  , colstr, mxvlib, mxvlrb, error )
c
c     ==================================================================
c     ====  xdslwc -- writes matrices in either ascii or binary     ====
c     ====            format to i/o unit, out of memory version     ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xdslwc writes a matrix stored on disk to i/o unit mtxunt
c
c     created         20-apr-99   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     mtxunt      i   i/o unit to use.
c     qascii      l   ascii or binary logical flag
c     mtxtyp      i   = 1 if symmetric, = 2 if unsymmetric
c     mtxstr      i   = 1 if general, = 3 if diagonal
c     neqns       i   order of matrix
c     nsuper      i   number of super nodes
c     sqfile      i   i/o file holding matrix data
c     perm        i   new-to-old permutation
c     xsup        i   super node partition array
c     diag        d   diagonal array
c     colstr      i   column start array
c
c     working storage
c     ---------------
c
c     mxvlib      i   buffer to hold integer matrix data 
c     mxvlrb      d   buffer to hold floating point matrix values
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

      integer              mtxunt, mtxtyp, mtxstr, neqns , nsuper, 
     1                     sqfile, error

      integer              perm(*),        xsup(*),        colstr(*),      
     1                     mxvlib(*)

      double precision     diag(*),        mxvlrb(*)
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------

      integer              inew  , iold  , isuper, 
     1                     jbgn  , jchv  , jend  , jnew  , jold  ,
     1                     k     , kbgn  , kend  , klong , koff

      double precision     val
 
c     ==================================================================

      call xislrw ( sqfile, error ) 
      if ( error .ne. 0 ) go to 8000

      do inew = 1, neqns

          iold = perm(inew)

          if ( qascii ) then
              write(mtxunt,'(2i10,d25.17)') iold, iold, diag(inew)
          else
              write(mtxunt                ) iold, iold, diag(inew)
          end if

      enddo

      if ( mtxstr .ne. 1 ) return

      do 400 isuper = 1, nsuper

          jbgn = xsup(isuper)
          jend = xsup(isuper+1) - 1

          kbgn  = colstr(jbgn)
          klong = colstr(jend+1) - kbgn
          koff  = kbgn - 1

c         ---------------------------
c         ... read in the matrix data
c         ---------------------------

          if ( klong .gt. 0 ) then 

              call xislvr ( sqfile, klong, mxvlib, error )
              if ( error .ne. 0 ) go to 8000

              call xislvs ( sqfile, error )
              if ( error .ne. 0 ) go to 8000

              call xdslvr ( sqfile, klong, mxvlrb, error )
              if ( error .ne. 0 ) go to 8000

          end if

c         --------------------------------------------
c         ... write out the data one chevron at a time
c         --------------------------------------------

          do jchv = jbgn, jend

              kbgn  = colstr(jchv  ) - koff
              kend  = colstr(jchv+1) - koff - 1

              do k = kbgn, kend

                  if ( mtxtyp .eq. 1 ) then

                      jnew = jchv
                      inew = mxvlib(k)

                      jold = min ( perm(jnew), perm(inew) )
                      iold = max ( perm(jnew), perm(inew) )

                      val  = mxvlrb(k)

                  else

                      if ( mxvlib(k) .le. 0 ) then
                          jnew = jchv
                          inew = jchv - mxvlib(k)
                      else
                          inew = jchv
                          jnew = jchv + mxvlib(k)
                      end if

                      jold = perm(jnew)
                      iold = perm(inew)

                      val  = mxvlrb(k)

                  end if

                  if ( qascii ) then
                      write(mtxunt,'(2i10,d25.17)') iold, jold, val
                  else
                      write(mtxunt                ) iold, jold, val
                  end if

              enddo

          enddo

  400 continue
 
c     ==================================================================

      return
 
c     ==================================================================
 
c     -------------
c     ... i/o error
c     -------------

 8000 continue
      error = -1

      return
      end
