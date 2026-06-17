      subroutine xdslf3 ( nstack, stknod, rstbeg, xlindx, istack,
     1                    mstack, iopos2, wafil2, walen2, watrn2,
     2                    rstack, ierr )
 
c
c     created            -- 24-feb-89, rgg
c     last modifications -- 12-feb-91, rgg - added i/o stats for wafil2
c                           31-jan-01, jgl - error handling documented
c
c     purpose
c     -------
c
c     spill the stack to i/o file wafil2
c
c     input variables
c     ---------------
c
c     nstack i   number of entries in the stack, that if not already
c                out of core, to put on wafil2.
c     stknod i   array containing the super node number of the matrix
c                on the stack
c     rstbeg i   array containing the starting location in rstack of
c                matrix on the stack
c     xlindx i   supernodal index pointer array -- used here to
c                determine size of update matrix.
c     wafil2 i   word addressable i/o unit
c     rstack d   real stack.
c
c     input/output variables
c     ----------------------
c
c     iopos2 i   next free position on i/o file wafil2
c     walen2 d   length of i/o file wafil2
c     watrn2 d   amount of i/o transfer to and from wafil2
c
c     output variables
c     ----------------
c
c     ierr   i   error flag
c                 0 -- success
c                -1 -- error writing to wafil2
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           nstack, iopos2, wafil2, mstack, ierr
 
      integer           stknod(*),     rstbeg(*),
     1                  xlindx(*),          istack(4,*)
 
      double precision  walen2, watrn2
 
      double precision  rstack(*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer           isuper, jstack, l,      mtxbgn, ncolpp, nupdat,
     1                  iobgn , iolast

      integer           idummy(1)
 
      double precision  temp
 
c  =====================================================================

c.debug
c     write(6,'(/,"in xdslf3 - nstack = ",i8)') nstack
c     call xislp2 ( 'stknod',   nstack, stknod, 6 )
c     call xislp2 ( 'rstbeg',   nstack, rstbeg, 6 )
c     call xislp2 ( 'istack', 4*nstack, istack, 6 )
c.debug

c     ------------------------
c     ... find end of i/o file
c     ------------------------

      iolast = 1

      do jstack = 1, nstack

          if ( rstbeg(jstack) .eq. 0 ) then

              iobgn  = istack(4,jstack)
              ncolpp = istack(3,jstack)
              nupdat = istack(1,jstack) + ncolpp
              l      = ncolpp + ( nupdat * ( nupdat + 1 ) ) / 2
              iolast = max ( iolast, iobgn + l )

          end if

      enddo
c.debug
c     write(6,'("resetting iopos2 - old and new    =  ", 2i8)')
c    1            iopos2, iolast
c.debug
      
      iopos2 = iolast
 
      do jstack = 1, nstack
 
          mtxbgn = rstbeg(jstack)
 
          if ( mtxbgn .gt. 0 ) then

              isuper = stknod(jstack)
              ncolpp = istack(3,jstack)
              nupdat = istack(1,jstack) + ncolpp
              l      = ncolpp + ( nupdat * ( nupdat + 1 ) ) / 2
 
              call xdslw2 ( wafil2, 2, idummy, idummy, rstack(mtxbgn), 
     1                      iopos2, l, ierr )
c.debug
c     write(6,'("jstack, rstack(mtxbgn), rstack(mtxbgn+l-1) = ", 
c    1        i8,1p,2d25.15)') 
c    2            jstack, rstack(mtxbgn), rstack(mtxbgn+l-1)
c.debug
 
              istack(4,jstack) = iopos2
              rstbeg(jstack)   = 0
              iopos2           = iopos2 + l
              temp             = iopos2 - 1
              walen2           = max ( walen2, temp )
              watrn2           = watrn2 + l
c.debug 
c     write(6,'("jstack, l, iopos2, walen2, watrn2 =  ", 3i8,2f15.3)')
c    1            jstack, l, iopos2, walen2, watrn2 
c.debug 
 
          end if
 
      enddo
 
c  =====================================================================
 
      return
      end
