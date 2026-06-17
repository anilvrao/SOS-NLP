      subroutine  xisli7  ( neqns , neqp1 , xladj , adjncy,
     1                      xadj  , last  , maxrow )
 
c
c     creation date: 15-apr-92 (rgg)
c     last updated:
c
c     purpose -- converts the lower triangular adjacency structure
c                of a sym. matrix, such as that recorded for any matrix
c                in the boeing-harwell format, to the corresponding full
c                adjacency structure required by sparse ordering
c                algorithms.  this routine is based on a routine of the
c                same name written by cleve ashcraft.  it differs in
c                requiring one less work vector.  another difference is
c                that it handles lower adjacency structures with or
c                without entries for the diagonal.
c
c                this version starts with the lower adjacency structure
c                in the array adjncy.  the lower adjacency structure
c                is overwritten with the full adjacency structure and
c                returned.  diagonals must not be present.
c
c     input parameters:
c
c        neqns  - number of equations
c        neqp1  - neqns+1, dimension of xladj, xadj, last
c       (xladj
c        adjncy)- lower adjacency structure
c
c     output parameters:
c
c       (xadj
c        adjncy) - full adjacency structure
c        maxrow  - max. no. of entries in a row a.
c
c     work parameters:
c
c        last   - pointers for full adjacency structure
c
c     external subprograms:
c
c        none
c
c     ==================================================================
 
c     ----------
c     parameters
c     ----------
 
      integer             neqns , neqp1 , maxrow
 
      integer             xladj(*)  ,     xadj(*)   ,
     1                    adjncy(*) ,     last(*)
 
c     ---------------
c     local variables
c     ---------------
 
      integer             j     , jp1   , k     , locatn, rowind
 
c     ==================================================================
 
c     --------------
c     initialization
c     --------------
 
      do j = 1, neqns
          last(j) = 0
      enddo
 
c     --------------------------------------
c     compute individual adjacency set sizes
c     for full adjacency structure
c     --------------------------------------
 
      do j = 1, neqns
          jp1 = j + 1
          do  k = xladj(j), xladj(jp1)-1
              rowind = adjncy(k)
              if  ( rowind .ne. j )  then
                  last(j)      = last(j)    + 1
                  last(rowind) = last(rowind) + 1
              endif
          enddo
      enddo
 
c     ---------------------------------------------------
c     compute pointers to the full adjacency structure
c     xadj(j) will point to the beginning of the j-th row
c     last(j) will point to the end of the j-th row
c     ---------------------------------------------------
 
      xadj(1) = 1
      maxrow   = 1
 
      do j = 1, neqns
          maxrow = max ( maxrow, last(j) )
          xadj(j+1) = xadj(j) + last(j)
          last(j)   = xadj(j+1) - 1
      enddo
 
c     ----------------------------------------------------------
c     copy the lower adjacency structure into the full adjacency
c     structure.  process the rows in reverse order so that the
c     full structure can overwrite the lower only structure.
c     ----------------------------------------------------------
 
      do j = neqns, 1, -1
 
          do k = xladj(j+1)-1, xladj(j), -1
 
              rowind = adjncy(k)
 
              if  ( rowind .ne. j )  then
 
                  locatn         = last(j)
                  adjncy(locatn) = rowind
                  last(j)        = last(j) - 1
 
                  locatn         = last(rowind)
                  adjncy(locatn) = j
                  last(rowind)   = last(rowind) - 1
 
              endif
 
          enddo
 
      enddo
 
c     ------------------------------------------
c     normal return - no error detection for now
c     ------------------------------------------
 
      return
      end
