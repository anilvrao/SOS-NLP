      subroutine xislit ( maxint, neqns , nnzero, offset, jcolo , 
     1                    rowlst, collst, colptr, itemp1, itemp2,
     2                    itemp3, itemp4, itemp5, itemp6  )
c
c
c     ==================================================================
c     ====  xislit -- gathers statistics for graph compression      ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xislit gathers the statistics from the adjacency structure
c     to be used in computing the graph compression>
c
c     created         28-apr-97   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     maxint      i   largest integer allowed.
c     neqns       i   order of uncompressed problem.
c     nnzero      i   the number of row and columns entries currently
c                     in memory
c     jcolo       i   last column number encountered.
c     rowlst,
c       collst    i   arrays to hold sections of row and column entries 
c                     of the strict lower triangle
c     itemp1      i   array of length neqns holding first random
c                     permutation for use in computing second checksum
c     itemp2      i   array of length neqns holding second random
c                     permutation for use in computing first checksum
c
c     output arguments
c     ----------------
c
c     colptr      i   pointer to the start of each column in 
c                     rowlst/collst on i/o file wafile
c     itemp3      i   holds first  set of checksums
c     itemp4      i   holds second set of checksums
c     itemp5      i   holds third  set of checksums
c     itemp6      i   holds first  row index for each column
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------

      integer             maxint, neqns , nnzero, offset, jcolo
 
      integer             rowlst (*), collst (*), colptr (*),
     1                    itemp1 (*), itemp2 (*), itemp3 (*),
     2                    itemp4 (*), itemp5 (*), itemp6 (*)

      real                f
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             i     , irow  , irow1 , irow2 , 
     1                    jcol  , jcol1 , jcol2 , k

c     ==================================================================

      f = maxint

      if ( sqrt(f) .gt. neqns ) then 

c         -----------------------------------------------------
c         ... maxint is big enough relative to neqns to assure 
c             no integer overflow in checksum computation.  use 
c             version of loops without overflow protection.
c         -----------------------------------------------------

          do i = 1, nnzero

              irow = rowlst(i)
              jcol = collst(i)

              if ( jcol .ne. jcolo ) then

                  do k = jcolo+1, jcol
                      colptr(k) = i + offset
                  enddo

                  jcolo     = jcol

              end if

              irow1 = itemp1(irow)
              irow2 = itemp2(irow)
              jcol1 = itemp1(jcol)
              jcol2 = itemp2(jcol)

              itemp3(jcol) = itemp3(jcol) + irow 
              itemp4(jcol) = itemp4(jcol) + irow1
              itemp5(jcol) = itemp5(jcol) + irow2

              itemp6(jcol) = min ( itemp6(jcol), irow )
    
              if ( irow .eq. jcol ) cycle
    
              itemp3(irow) = itemp3(irow) + jcol 
              itemp4(irow) = itemp4(irow) + jcol1
              itemp5(irow) = itemp5(irow) + jcol2
    
              itemp6(irow) = min ( itemp6(irow), jcol )

          enddo

      else

c         --------------------------------------------------------
c         ... maxint is not big enough relative to neqns to assure 
c             no integer overflow in checksum computation.  use 
c             version of loops with overflow protection.  since
c             the checksums are additive and maxint is set to be
c             neqns less than largest integer, a simple check
c             is sufficient.
c         --------------------------------------------------------

          do 220 i = 1, nnzero

              irow = rowlst(i)
              jcol = collst(i)

              if ( jcol .ne. jcolo ) then

                  do k = jcolo+1, jcol
                      colptr(k) = i + offset
                  enddo

                  jcolo     = jcol

              end if

              irow1 = itemp1(irow)
              irow2 = itemp2(irow)
              jcol1 = itemp1(jcol)
              jcol2 = itemp2(jcol)

              itemp3(jcol) = itemp3(jcol) + irow

              if ( itemp3(jcol) .gt. maxint ) then
                  itemp3(jcol) = itemp3(jcol) - maxint
              end if

              itemp4(jcol) = itemp4(jcol) + irow1

              if ( itemp4(jcol) .gt. maxint ) then
                  itemp4(jcol) = itemp4(jcol) - maxint
              end if

              itemp5(jcol) = itemp5(jcol) + irow2

              if ( itemp5(jcol) .gt. maxint ) then
                  itemp5(jcol) = itemp5(jcol) - maxint
              end if

              itemp6(jcol) = min ( itemp6(jcol), irow )
    
              if ( irow .eq. jcol ) go to 220
    
              itemp3(irow) = mod ( itemp3(irow) + jcol , maxint )

              if ( itemp3(irow) .gt. maxint ) then
                  itemp3(irow) = itemp3(irow) - maxint
              end if

              itemp4(irow) = mod ( itemp4(irow) + jcol1, maxint )

              if ( itemp4(irow) .gt. maxint ) then
                  itemp4(irow) = itemp4(irow) - maxint
              end if

              itemp5(irow) = mod ( itemp5(irow) + jcol2, maxint )

              if ( itemp5(irow) .gt. maxint ) then
                  itemp5(irow) = itemp5(irow) - maxint
              end if

              itemp6(irow) = min ( itemp6(irow), jcol )

  220     continue

      end if

c     ==================================================================
          
      return
      end
