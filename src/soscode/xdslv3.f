      subroutine   xdslv3   ( jcol  , nzcol , jrowin, lvalue, unsym , 
     1                        neqns , nnzero, maxnz , nrecrd, maxrec,
     2                        wafile, invp  , coord1, coord2, diag  , 
     3                        values, recpos, reclen, srtlst, srtval,
     4                        colstr, error  )
c
c
c     ==================================================================
c     ====  xdslv3 -- add a list of row entries for a given column  ====
c     ====            to the row and column and value lists         ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xdslv3 adds a list of row entries to the coordinate and value 
c     lists representing the matrix structure.  when the lists are full
c     they are sorted and compressed.  If necessary, some of lists
c     are written to file wafile.
c
c     created         20-dec-96   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c     jcol        i   column number.
c     nzcol       i   no. of nonzeroes in the column.
c     jrowin      i   row indicies for the column.
c     lvalue      d   values for the column.
c     unsym       l   symmetric/unsymmetric logical flag
c     neqns       i   number of equations
c     maxnz       i   size of coord2 and coord1 arrays.
c     maxrec      i   length of recpos and reclen arrays.
c     wafile      i   unit number for i/o file
c     invp        i   old to new permutation
c
c     input/output arguments
c     ----------------------
c
c     nnzero      i   the number of entries in coord2 and coord1
c     nrecrd      i   number of records written to wafile
c     diag        d   values of the permuted diagonal
c     coord1,
c       coord2    i   lists of coordinates
c     values      d   corresponding value entries
c     recpos      i   array of length nrecrd holding starting i/o
c                     positions
c     reclen      i   array of length nrecrd holding length of
c                     each i/o record
c
c     working storage
c     ---------------
c
c     srtlst      i   array used for integer sorting
c     srtval      i   array used for floating point sorting
c     colstr      i   array used for sorting
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
 
      integer             jcol  , nzcol , neqns , nnzero,
     1                    maxnz , nrecrd, maxrec, wafile, error

      logical             unsym
 
      integer             jrowin (*), 
     1                    invp   (*), coord1 (*), coord2 (*),
     2                    recpos (*), reclen (*),
     3                    srtlst (*), colstr (*)

      double precision    lvalue (*), diag   (*), values (*),
     1                    srtval (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             i     , newcol, newrow, oldrow 

      double precision    value
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer             xdslni
 
      external            xdslni
 
c     ==================================================================

      error = 0

c     --------------------------------
c     ... for each entry of the column 
c     --------------------------------

      newcol = invp(jcol)

      do i = 1, nzcol                        

          oldrow = jrowin(i)

c         --------------------           
c         ... get new indicies
c         --------------------

          newrow = invp(oldrow)

          value  = lvalue(i)

c         ----------------------------
c         ... process a diagonal entry
c         ----------------------------

          if ( newcol .eq. newrow ) then
              diag(newcol) = diag(newcol) + value
              cycle
          end if

c         ------------------------------------------------
c         ... check to see if entry can be added to lists.
c             if not, sort and compress.
c         ------------------------------------------------

          if ( nnzero+1 .gt. maxnz ) then 

c.debug
c     write(6,'("before xdslvz")')
c.debug
              call xdslvz ( .false., neqns, nnzero, maxnz, nrecrd,
     1                      maxrec, wafile, coord1, coord2, values,
     2                      recpos, reclen, srtlst, srtval, colstr,
     3                      error )
 
              if ( error .ne. 0 ) return
 
          end if

c         ----------------------------------
c         ... convert to chevron coordinates
c         ----------------------------------
 
          nnzero = nnzero + 1
 
          coord1(nnzero) = min ( newcol, newrow )

          if ( unsym ) then
              coord2(nnzero) = newcol - newrow
          else
              coord2(nnzero) = abs ( newcol - newrow )
          end if

          values(nnzero) = value

      enddo    

c-----------------------------------------------------------------------

c     -----------------                                                  
c     ... normal return
c     -----------------

      return            
      end
