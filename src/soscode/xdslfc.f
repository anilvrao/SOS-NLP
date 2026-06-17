      subroutine xdslfc ( unsym , nsuper, npanel, xsup  , xlindx,
     1                    xpanel, wafil1, wafil4, lnz   , error )
c
c     purpose
c     -------
c
c     xdslfc reads the out-of-core representation of lnz and unz from 
c     i/o files wafil1 and wafil4 and stores into array lnz
c
c     created         30-jul-98   -- rgg --
c     modified        
c
c     input arguments
c     ---------------
c
c     unsym       l   logical indicating unsym or sym problem
c     nsuper      i   number of super nodes
c     npanel      i   number of panels in the factorization
c     xsup        i   supernode pointer array
c     xlindx      i   supernodal indicies pointer array
c     xpanel      i   panel point array
c     wafil1      i   word addressable i/o file holding lnz
c     wafil4      i   word addressable i/o file holding unz
c
c     output arguments
c     ----------------
c
c     lnz         d   array holding lnz and unz
c     error       i   error flag
c                     =    0  normal return
c                     =    1  i/o error on wafil1
c                     =    2  i/o error on wafil4
c
c---------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------

      logical             unsym
 
      integer             nsuper, npanel, wafil1, wafil4, error
 
      integer             xsup(*),        xlindx(*),
     1                    xpanel(*)
 
      double precision    lnz(*)
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer             fctpnt, ipanel, ipbgn , ipend , 
     1                    istart, isuper,
     2                    l     , n1    , n2    , nodext

      integer             ipos(2)
 
c     ================================================================
 
c     -----------------------------------------------------
c     ... perform the forward elimination step and diagonal
c         scaling for each supernode
c     -----------------------------------------------------
 
      error  = 0
      fctpnt = 1

      ipos(1) = 1
      ipos(2) = 1
 
      do 200 isuper = 1, nsuper

          ipbgn  = xsup(isuper)
          ipend  = xsup(isuper+1) - 1
          n2     = xpanel(ipend+1) - xpanel(ipbgn)
          
          nodext = xlindx ( isuper+1 ) - xlindx ( isuper )

          do ipanel = ipbgn, ipend

c             -----------------------------------------
c             ... read this portion of lnz from wafil1.
c             -----------------------------------------
 
              istart = xpanel ( ipanel )
              n1     = xpanel ( ipanel+1 ) - istart

              n2     = n2 - n1

              l      = n1 * ( n1 - 1 ) / 2
     1               + n1 * ( n2 + nodext )

              call xdslw5 ( wafil1, lnz(fctpnt),
     1                      ipos, l, error )

              if ( error .ne.0 ) then
                  error = 1
                  return
              end if

              fctpnt = fctpnt + l

c             -----------------------------------------
c             ... read this portion of unz from wafil4.
c             -----------------------------------------

              if ( unsym ) then
 
                  call xdslw5 ( wafil4, lnz(fctpnt),
     1                          ipos, l, error )

                  if ( error .ne.0 ) then
                      error = 2
                      return
                  end if

                  fctpnt = fctpnt + l

              end if

c             ------------------------------------------
c             ... set i/o position for wafil1 and wafil4
c             ------------------------------------------

              call xdslw9 ( ipos, l )

          enddo

  200 continue

c     ================================================================
 
c     ------------------------
c     ... end of module xdslfc
c     ------------------------
 
      return
      end
