      subroutine   xisli3   ( jcol  , nlist,  rowind, neqns , nnzero, 
     1                        maxnz , nrecrd, maxrec, lstcol, sort  ,
     2                        wafile, rowlst, collst, srtlst, recpos,
     3                        reclen, colstr, error  )
c
c
c     ==================================================================
c     ====  xisli3 -- add all of the entries from a column          ====
c     ====            to the row and column lists                   ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xisli3 adds all of the entries from a column to the row 
c     and columns lists representing the matrix structure.  
c     when the list is full it is sorted and compressed.  If 
c     necessary, some of list is written to file wafile.
c
c     created         12-dec-96   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     jcol        i   column number
c     nlist       i   number of row indicies
c     rowind      i   array holding row indicies
c     neqns       i   order of the problem
c     maxnz       i   size of rowlst and collst arrays.
c     maxrec      i   size of recpos and reclen arrays.
c     wafile      i   unit number for i/o file
c
c     input/output arguments
c     ----------------------
c
c     nnzero      i   the number of entries in rowlst and collst
c     nrecrd      i   number of records written to wafile
c     lstcol      i   last column written to wafile
c     sort        i   flag on whether to perform sort/merge.
c     rowlst,
c       collst    i   lists of row and column entries
c     srtlst      i   array used in sorting rowlst and collst
c     recpos      i   array of length nrecrd holding starting i/o
c                     positions
c     reclen      i   array of length nrecrd holding length of
c                     each i/o record
c     colstr      i   array used in sorting rowlst and collst
c
c     output arguments
c     ----------------
c
c     error       i   error flag, = -1 if i/o error occurred.
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      integer             jcol  , nlist , neqns , nnzero, maxnz ,
     1                    nrecrd, maxrec, lstcol, sort  , wafile, 
     2                    error
 
      integer             rowind (*), rowlst (*), collst (*),
     1                    srtlst (*), recpos (*), reclen (*),
     2                    colstr (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             i
 
c     ==================================================================

      error = 0

      if ( nlist .eq. 0 ) return

c.debug
c     write(6,'("in xisli3 - nnzero, maxnz, jcol = ", 3i8)')
c    1                        nnzero, maxnz, jcol
c     call xislp3 ( 'rowind', nlist, rowind, 6 )
c.debug
 

c     ==================================================================

c     -------------------------------------------------------
c     ... start entering the entries into rowlst and 
c         collst until something happens.
c     -------------------------------------------------------
 
      do i = 1, nlist

          if ( rowind(i) .eq. 0 .or. rowind(i) .eq. jcol ) cycle

          if ( nnzero+1 .gt. maxnz ) then

              call xisliz ( .false., neqns , nnzero, maxnz , nrecrd, 
     1                      maxrec, lstcol , sort  , wafile, 0,
     2                      rowlst, collst , srtlst, recpos, reclen, 
     3                      colstr, error )

              if ( error .ne. 0 ) return

          end if

          nnzero = nnzero + 1

          collst(nnzero) = min ( jcol, rowind(i) )
          rowlst(nnzero) = max ( jcol, rowind(i) )

      enddo
 
c     ==================================================================

      return
      end
