      subroutine   xdslvz   ( qlast , neqns , nnzero, maxnz , nrecrd, 
     1                        maxrec, wafile, coord1, coord2, values,
     2                        recpos, reclen, srtlst, srtval, colstr,
     3                        error  )
c
c
c     ==================================================================
c     ====  xdslvz -- sorts and compresses coord2, coord1, and      ====
c     ====            values.  if necessary parts are written to    ====
c     ====            wafile.                                       ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xdslvz sorts and compresses the  coordinate and values arrays when 
c     the list is full.  If necessary, some of list is written to file 
c     wafile.
c
c     created         19-dec-96   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     qlast       l   logical flag.  if .true. then this is the
c                     last call so do not spill to disk.
c     neqns       i   number of equations
c     maxnz       i   size of coord2 and coord1 arrays.
c     maxrec      i   size of recpos and reclen arrays.
c     wafile      i   unit number for i/o file
c
c     input/output arguments
c     ----------------------
c
c     nnzero      i   the number of entries in coord2 and coord1
c     nrecrd      i   number of records written to wafile
c     coord1,
c       coord2    i   lists of coordinates
c     values      d   lists of value entries
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

      logical             qlast
 
      integer             neqns , nnzero, maxnz , nrecrd, maxrec, 
     1                    wafile, error
 
      integer             coord1 (*), coord2 (*),
     1                    recpos (*), reclen (*),
     2                    srtlst (*), colstr (*)

      double precision    values (*), srtval (*)
c.timer
c     double precision    t1, t2, t3, w1, w2, w3
c.timer
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             count , i     , ibgn  , iend  , ipos  ,
     1                    irow  , j     , jbgn  , jchv  , k     ,
     2                    len   , lstchv, orgnnz

c     ==================================================================

c.timer
c     call xdslt1 ( t1, w1 )
c.timer
c.debug
c     write(6,'("at start of xdslvz")')
c     write(6,'("qlast  = ", l8)') qlast
c     write(6,'("nnzero = ", i8)') nnzero
c     write(6,'("maxnz  = ", i8)') maxnz 
c     write(6,'("nrecrd = ", i8)') nrecrd
c     write(6,'("wafile = ", i8)') wafile
c     write(6,'("error  = ", i8)') error 
c     call xislp3 ( 'coord1 list', nnzero, coord1, 6 )
c     call xislp3 ( 'coord2 list', nnzero, coord2, 6 )
c     call xdslp5 ( 'values list', nnzero, values, 6 )
c     if ( nrecrd .gt. 0 ) then
c         call xislp3 ( 'record pos.', nrecrd, recpos, 6 )
c         call xislp3 ( 'record len.', nrecrd, reclen, 6 )
c     end if
c.debug

c     -----------------------------------------------------------
c     ... coord1, coord2, and values are full.  sort and compress.
c         first accumulate counts of nonzeroes in each
c         column.  build colstr array.
c     -----------------------------------------------------------

      orgnnz = nnzero

      do i = 1, neqns
          colstr(i) = 0
      enddo

      do k = 1, nnzero
          jchv = coord1(k)
          colstr(jchv) = colstr(jchv) + 1
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
          jchv         = coord1(k)
          j            = colstr(jchv)
          srtlst(j)    = coord2(k) 
          srtval(j)    = values(k) 
          colstr(jchv) = j + 1
      enddo
c.debug
c     call xislp3 ( 'column start after 40', neqns, colstr, 6 )
c.debug

c     ---------------------------------------------
c     ... sort indices for each column and compress
c     ---------------------------------------------

      iend   = 0 
      nnzero = 0

      do jchv = 1, neqns
 
          ibgn = iend + 1
          iend = colstr(jchv) - 1
          len  = iend - ibgn + 1

          if ( len .eq. 0 ) cycle

          call xdslq2 ( len, srtlst(ibgn), srtval(ibgn) )

          irow = 0
          jbgn = nnzero + 1

          do i = ibgn, iend

              if ( srtlst(i) .eq. irow ) then

                  values(nnzero) = values(nnzero) + srtval(i)

              else

                  irow           = srtlst(i)
                  nnzero         = nnzero + 1
                  coord1(nnzero) = jchv
                  coord2(nnzero) = irow
                  values(nnzero) = srtval(i)

              end if

          enddo
c.debug
c     write(6,'("after 50 - jchv, jbgn, nnzero = ",3i8)')
c    1                       jchv, jbgn, nnzero 
c     call xislp3 ( 'coord1 list', nnzero-jbgn+1, coord1(jbgn), 6 )
c     call xislp3 ( 'coord2 list', nnzero-jbgn+1, coord2(jbgn), 6 )
c     call xdslp5 ( 'values list', nnzero-jbgn+1, values(jbgn), 6 )
c.debug

      enddo

c.debug
c     write(6,'("in xdslvz - nonzeroes - org. & final = ",
c    1          2i15)') orgnnz, nnzero
c     call xislp3 ( 'coord1 list', nnzero, coord1, 6 )
c     call xislp3 ( 'coord2 list', nnzero, coord2, 6 )
c     call xdslp5 ( 'values list', nnzero, values, 6 )
c.debug
c.timer
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'("in xdslvz - time to sort = ", f15.3)') t2
c.timer

c     ==================================================================

c     -------------------------------------------------------------
c     ... coord1, coord2, and values are now sorted and compressed.
c         if they are still nearly full then write part of them 
c         to wafile.
c     -------------------------------------------------------------

      if ( ( .not. qlast  .and.  nnzero .gt. .8 * maxnz ) 
     1                 .or.
     2     ( qlast .and. nrecrd .gt. 0 ) ) then

c         --------------------------------
c         ... if not opened, open i/o file
c         --------------------------------

          if ( nrecrd .gt. maxrec ) then 
              error = -1
              return
          end if

          if ( nrecrd .eq. 0 ) then 

              call xdslw4 ( wafile, 3, error )
              if ( error .ne. 0 ) return

          end if

c         ------------------------------------------------------------
c         ... select the amount of coord1, coord2, and values to write 
c             to wafile.  Start at 70 percent and find a nearby 
c             column boundary.
c         ------------------------------------------------------------
 
          ipos    = 1

          do i = 1, nrecrd
              ipos = ipos + reclen (i)
          enddo

          if ( qlast ) then

              len = nnzero

          else

              len    = .7 * maxnz
              j      = coord1(len)
  
              do i = len-1, 1, -1
                  if ( coord1(i) .ne. j ) then
                      len    = i
                      lstchv = coord1(i)
                      go to 340
                  end if
              enddo
  
              do i = len+1, maxnz
                  if ( coord1(i) .ne. j ) then
                      len    = i-1
                      lstchv = j
                      go to 340
                  end if
              enddo

c.debug
              write(6,'("oops after 321 in xdslvz")')
              RETURN
c.debug

          end if

c         ----------------------------
c         ... perform the actual write
c         ----------------------------

  340     continue
c.debug
c     write(6,'("before write - nrecrd, len, ipos = ", 3i8)')
c    1                           nrecrd, len, ipos 
c.debug

          call xdslw2 ( wafile, 3, coord1, coord2, values,
     1                  ipos, len, error )
          if ( error .ne. 0 ) return

          nrecrd         = nrecrd + 1
          reclen(nrecrd) = len
          recpos(nrecrd) = ipos
c.debug
c     write(6,'("after  write - nrecrd, len, ipos = ", 3i8)')
c    1                           nrecrd, len, ipos 
c.debug

c        ----------------------------------------------------
c         ... collapse in-memory storage for rowlst and collst
c         ----------------------------------------------------

          j      = len + 1 
          len    = nnzero - len
 
          call xislmv ( len, coord1, j, 1 )
          call xislmv ( len, coord2, j, 1 )
          call xdslmv ( len, values, j, 1 )
 
          nnzero = len

c.debug
c     write(6,'("after dumping to disk in xdslvz")')
c     write(6,'("nnzero = ", i8)') nnzero
c     if ( nnzero .gt. 0 ) then
c         call xislp3 ( 'coord1 list', nnzero, coord1, 6 )
c         call xislp3 ( 'coord2 list', nnzero, coord2, 6 )
c         call xdslp5 ( 'values list', nnzero, values, 6 )
c     end if
c     if ( nrecrd .gt. 0 ) then
c         call xislp3 ( 'record pos.', nrecrd, recpos, 6 )
c         call xislp3 ( 'record len.', nrecrd, reclen, 6 )
c     end if
c.debug
c.timer
c     call xdslt2 ( t1, w1, t3, w3 )
c     write(6,'("in xdslvz - time to dump = ", f15.3)') t3-t2
c.timer
 
      end if
 
c     ==================================================================
c.debug
c     write(6,'("at end   of xdslvz")')
c     write(6,'("qlast  = ", l8)') qlast
c     write(6,'("nnzero = ", i8)') nnzero
c     write(6,'("maxnz  = ", i8)') maxnz 
c     write(6,'("nrecrd = ", i8)') nrecrd
c     write(6,'("wafile = ", i8)') wafile
c     write(6,'("error  = ", i8)') error 
c     call xislp3 ( 'coord1 list', nnzero, coord1, 6 )
c     call xislp3 ( 'coord2 list', nnzero, coord2, 6 )
c     call xdslp5 ( 'values list', nnzero, values, 6 )
c     if ( nrecrd .gt. 0 ) then
c         call xislp3 ( 'record pos.', nrecrd, recpos, 6 )
c         call xislp3 ( 'record len.', nrecrd, reclen, 6 )
c     end if
c.debug

      return

      end
