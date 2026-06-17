      subroutine xdslel ( rowoff, pnlcol, pnlrow, nswap , 
     1                    swap  , rectan )
      
 
c
c  purpose -- to apply permutations from later panels to earlier panels
c
c  created            -- 19-nov-97, rgg
c                        08-dec-97, rgg
c
c  input variables --
c
c      rowoff -- row offset between swap and local storage
c      pnlcol -- number of columns in the panel
c      pnlrow -- number of rows in the rectangular portion of the 
c                factorization for this panel 
c      nswap  -- number of swaps to apply
c      swap   -- swapping record array.
c
c  output variable --
c
c      rectan -- array holding the rectangular portion of the 
c                factorization for this panel 
c     
c     ------------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------

      integer            rowoff, pnlcol, pnlrow, nswap

      integer            swap(*)

      double precision   rectan(pnlrow,pnlcol)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer            irow  , jrow
      
c     ------------------------------------------------------------------
c.debug
c     write(6,'("rowoff on entry to xdslel = ", 3i8)') rowoff
c     write(6,'("nswap  on entry to xdslel = ", 3i8)') nswap 
c     call xislp3 ( 'swap in xdslel', nswap, swap, 6 )
c.debug

      do irow = 1, nswap

          jrow = swap(irow) - rowoff
c.debug
c         write(6,'("irow, jrow = ", 3i8)') irow, jrow              
c         if ( jrow .lt. 1 .or. jrow .gt. pnlrow ) stop
c.debug

          if ( irow .ne. jrow ) then

              call dswap ( pnlcol, rectan(irow,1), pnlrow,
     1                             rectan(jrow,1), pnlrow )

          end if

      enddo
      
c     ------------------------------------------------------------------
 
      return
      end 
