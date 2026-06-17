      subroutine   xisliz   ( qlast , neqns , nnzero, maxnz , nrecrd, 
     1                        maxrec, lstcol, sort  , wafile, length,
     2                        rowlst, collst, srtlst, recpos, reclen,
     3                        colstr, error  )
c
c
c     ==================================================================
c     ====  xisliz -- sorts and compresses rowlst and collst        ====
c     ====            if necessary parts are written to wafile      ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xisliz sorts and compress rowlst and collst when the list is full.
c     If necessary, some of list is written to file wafile.
c
c     created         03-dec-96   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     qlast       l   logical flag.  if .true. then this is the
c                     last call so do not spill to disk.
c     neqns       i   number of equations
c     maxnz       i   size of rowlst and collst arrays.
c     maxrec      i   size of recpos and reclen arrays.
c     wafile      i   unit number for i/o file
c     length      i   length of working storage in integers.  used
c                     when qlast = .true. to determine if data should
c                     remain in memory.
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

      logical             qlast
 
      integer             neqns , nnzero, maxnz , nrecrd, maxrec,
     1                    lstcol, sort  , wafile, length, error
 
      integer             rowlst (*), collst (*), srtlst (*),
     1                    recpos (*), reclen (*), colstr (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             count , i,      ibgn,   iend,   ipos, 
     1                    irow  , j,      jbgn,   jcol,   k,
     2                    len   , locerr, orgnnz
c.timer
c     double precision    t1, t2, t3, w1, w2, w3
c.timer

c     ==================================================================

c.timer
c     call xdslt1 ( t1, w1 )
c.timer
c.debug
c     write(6,'("at start of xisliz")')
c     write(6,'("qlast  = ", l8)') qlast
c     write(6,'("nnzero = ", i8)') nnzero
c     write(6,'("maxnz  = ", i8)') maxnz 
c     write(6,'("nrecrd = ", i8)') nrecrd
c     write(6,'("wafile = ", i8)') wafile
c     write(6,'("error  = ", i8)') error 
c     call xislp3 ( 'row    list', nnzero, rowlst, 6 )
c     call xislp3 ( 'column list', nnzero, collst, 6 )
c     call xislp3 ( 'row    list', nnzero, rowlst, 6 )
c     call xislp3 ( 'column list', nnzero, collst, 6 )
c     if ( nrecrd .gt. 0 ) then
c         call xislp3 ( 'record pos.', nrecrd, recpos, 6 )
c         call xislp3 ( 'record len.', nrecrd, reclen, 6 )
c     end if
c.debug

c     ---------------------------------------------------
c     ... rowlst and collst are full.  sort and compress.
c         first accumulate counts of nonzeroes in each
c         column.  build colstr array.
c     ---------------------------------------------------

      orgnnz = nnzero

      do i = 1, neqns
          colstr(i) = 0
      enddo

      do k = 1, nnzero
          jcol = collst(k)
          colstr(jcol) = colstr(jcol) + 1
      enddo

      j     = 1
      count = 0
 
      do i = 1, neqns
          len       = colstr(i)
          count     = count + len
          colstr(i) = j
          j         = j + len
      enddo

c     --------------------------------------
c     ... move row indices into column order
c     --------------------------------------

      do k = 1, nnzero
          jcol         = collst(k)
          j            = colstr(jcol)
          srtlst(j)    = rowlst(k) 
          colstr(jcol) = j + 1
      enddo

c     ---------------------------------------------
c     ... sort indices for each column and compress
c     ---------------------------------------------

      iend   = 0 
      nnzero = 0

      do jcol = 1, neqns
 
          ibgn = iend + 1
          iend = colstr(jcol) - 1
          len  = iend - ibgn + 1

          if ( len .eq. 0 ) cycle

          call xislq1 ( len, srtlst(ibgn), locerr )

          irow = 0
          jbgn = nnzero + 1

          do i = ibgn, iend

              if ( srtlst(i) .ne. irow ) then

                irow           = srtlst(i)
                nnzero         = nnzero + 1
                collst(nnzero) = jcol
                rowlst(nnzero) = irow
              endif

          enddo
c.debug
c     write(6,'("after 50")')
c     call xislp3 ( 'column list', nnzero-jbgn+1, collst(jbgn), 6 )
c     call xislp3 ( 'row    list', nnzero-jbgn+1, rowlst(jbgn), 6 )
c.debug

      enddo

c.debug
c     write(6,'("in xislvz - nonzeroes - org. & final = ", 
c    1          2i15)') orgnnz, nnzero
c     call xislp3 ( 'row    list', nnzero, rowlst, 6 )
c     call xislp3 ( 'column list', nnzero, collst, 6 )
c.debug
c.timer
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'("in xislvz - time to sort = ", f15.3)') t2
c.timer

c     ==================================================================

c     ---------------------------------------------------
c     ... rowlst and collst are now sorted and compressed.
c         if rowlst and collst are still nearly full then
c         write part of them to wafile.
c     ---------------------------------------------------

      k = 8 * ( neqns + 1 ) + 2 * nnzero
c.debug
c     k = length + 1
c     write(6,'("in xisliz - qlast, nrecrd, length = ", l12, 2i12)')
c    1                        qlast, nrecrd, length
c     write(6,'("            neqns, nnzero, k      = ", 3i12)')
c    1                        neqns, nnzero, k     
c.debug

      if ( ( .not. qlast  .and.  nnzero .gt. .8 * maxnz ) 
     1                 .or.
     2     ( qlast .and. nrecrd .gt. 0 )
     3                 .or.
     4     ( qlast .and. k .gt. length ) ) then

c         --------------------------------
c         ... if not opened, open i/o file
c         --------------------------------

          if ( nrecrd .ge. maxrec+1 ) then 
              error = -1
              return
          end if

          if ( nrecrd .eq. 0 ) then 

              call xislw4 ( wafile, error )
              if ( error .ne. 0 ) return

          end if

c         ---------------------------------------------
c         ... test if sort/merge needs to be performed.
c         ---------------------------------------------

          if ( collst(1) .le. lstcol ) sort = 1

c         ------------------------------------------------------
c         ... select the amount of rowlst and collst to write to
c             wafile.  Start at 70 percent and find a nearby 
c             column boundary.
c         ------------------------------------------------------
 
          ipos    = 1

          do i = 1, nrecrd
              ipos = ipos + reclen (i)
          enddo

          if ( qlast ) then

              len = nnzero

          else

              len    = .7 * maxnz
              j      = collst(len)
  
              do i = len-1, 1, -1
                  if ( collst(i) .ne. j ) then
                      len    = i
                      lstcol = collst(i)
                      go to 340
                  end if
              enddo
  
              do i = len+1, maxnz
                  if ( collst(i) .ne. j ) then
                      len    = i-1
                      lstcol = j
                      go to 340
                  end if
              enddo

c.debug
              write(6,'("oops after 321 in xisliz")')
              stop
c.debug

          end if

c         ----------------------------
c         ... perform the actual write
c         ----------------------------

  340     continue
          call xislw2 ( wafile, rowlst, collst, ipos, len, error )
          if ( error .ne. 0 ) return

          nrecrd         = nrecrd + 1
          reclen(nrecrd) = len
          recpos(nrecrd) = ipos
c.debug
c     write(6,'("after write - nrecrd, len, ipos = ", 3i8)')
c    1                          nrecrd, len, ipos 
c     i = min ( 50, len )
c     j = len - i + 1
c     call xislp3 ( 'column list - first 50', i, collst, 6 )
c     call xislp3 ( 'row    list - first 50', i, rowlst, 6 )
c     call xislp3 ( 'column list - last  50', i, collst(j), 6 )
c     call xislp3 ( 'row    list - last  50', i, rowlst(j), 6 )
c.debug

c        ----------------------------------------------------
c         ... collapse in-memory storage for rowlst and collst
c         ----------------------------------------------------

          j      = len + 1 
          len    = nnzero - len
 
          call xislmv ( len, rowlst, j, 1 )
          call xislmv ( len, collst, j, 1 )
 
          nnzero = len

c.debug
c         write(6,'("after dumping to disk in xisliz")')
c         write(6,'("nnzero = ", i8)') nnzero
c     if ( nnzero .gt. 0 ) then
c         call xislp3 ( 'row    list', nnzero, rowlst, 6 )
c         call xislp3 ( 'column list', nnzero, collst, 6 )
c     end if
c     if ( nrecrd .gt. 0 ) then
c         call xislp3 ( 'record pos.', nrecrd, recpos, 6 )
c         call xislp3 ( 'record len.', nrecrd, reclen, 6 )
c     end if
c.debug
c.timer
c     call xdslt2 ( t1, w1, t3, w3 )
c     write(6,'("in xislvz - time to dump = ", f15.3)') t3-t2
c.timer
 
      end if
 
c     ==================================================================

      return

      end
