      subroutine   xisli2   ( nlist , lrowls, lcolls, neqns , nnzero,
     1                        maxnz , nrecrd, maxrec, lstcol, sort  , 
     2                        wafile, rowlst, collst, srtlst, recpos,
     3                        reclen, colstr, error  )
c
c
c     ==================================================================
c     ====  xisli2 -- add a list of triples to the row and column   ====
c     ====            lists                                         ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xisli2 adds a list of triples to the row and columns lists 
c     representing the matrix structure.  when the list is full
c     it is sorted and compressed.  If necessary, some of list
c     is written to file wafile.
c
c     created         03-dec-96   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     nlist       i   number of entries in the local list
c     lrowls      i   array holding row entries to insert
c     lcolls      i   array holding column entries to insert
c     maxnz       i   size of rowlst and collst arrays.
c     maxrec      i   size of recpos and reclen arrays.
c     wafile      i   unit number for i/o file
c
c     input/output arguments
c     ----------------------
c
c     neqns       i   number of equations
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
 
      integer             nlist , neqns , nnzero, maxnz , nrecrd, 
     1                    maxrec, lstcol, sort  , wafile, error
 
      integer             lrowls (*), lcolls (*), 
     1                    rowlst (*), collst (*), srtlst(*),
     2                    recpos (*), reclen (*), colstr(*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             i
 
c     ==================================================================

      error = 0

      if ( nlist .eq. 0 ) return

c.debug
c     write(6,'("in xisli2 - nnzero, maxnz = ", 2i8)')
c    1                        nnzero, maxnz
c     call xislp3 ( 'local row    list', nlist, lrowls, 6 )
c     call xislp3 ( 'local column list', nlist, lcolls, 6 )
c.debug
 
c     ==================================================================

c     -------------------------------------------------------
c     ... start entering the entries into rowlst and 
c         collst until something happens.
c     -------------------------------------------------------
 
      do i = 1, nlist

          if ( lcolls(i) .eq. lrowls(i) ) cycle

          if ( nnzero+1 .gt. maxnz ) then

              call xisliz ( .false., neqns, nnzero, maxnz, nrecrd,
     1                      maxrec, lstcol, sort  , wafile, 0, 
     2                      rowlst, collst, srtlst, recpos, reclen, 
     3                      colstr, error )

              if ( error .ne. 0 ) return

          end if

          nnzero = nnzero + 1

          collst(nnzero) = min ( lcolls(i), lrowls(i) )
          rowlst(nnzero) = max ( lcolls(i), lrowls(i) )

      enddo
 
c     ==================================================================

      return
      end
