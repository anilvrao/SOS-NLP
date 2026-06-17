      subroutine xdsle2 ( cmajor, colbgn, colend, loclfr, front,
     1                    pnlcol, pnlrow, panel )
 
c
c  purpose -- extract the current panel from the front.
c
c  created            -- 07-feb-97, rgg
c  last modifications -- 
c
c  input variables --
c
c      cmajor -- column major flag
c                .eq. 0 row    major form for panel is used
c                .ne. 0 column major form for panel is used
c      colbgn -- first column to extract
c      colend -- last  column to extract
c      loclfr -- size of front
c      front  -- frontal matrix
c      pnlcol -- number of columns in panel
c      pnlrow -- number of rows    in panel
c
c  output variables --
c      panel  -- rectangular the symmetric frontal matrix in packed storage
c
c     ------------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer            cmajor, colbgn, colend, loclfr, pnlcol, 
     1                   pnlrow
 
      double precision   front (*),      panel (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer            jcol  , k     , kk    , kzbgn ,
     1                   l     , length, stride
 
c     --------------------
c     ... subprograms used
c     --------------------
 
c     ------------------------------------------------------------------

      k      = 1
      kk     = 1
      length = loclfr

      if ( cmajor .eq. 0 ) then
          stride = pnlcol
      else
          stride = 1
      end if

      do jcol = colbgn, colend

c         ------------------------------------------------
c         ... zero out the leading portion of the column 
c             down to but not including the diagonal.
c             this required due to the use of dgemm during
c             the factorization.
c         ------------------------------------------------

          l = jcol - colbgn
          kzbgn = kk - stride * l
          do i=0,l-1
            panel(kzbgn+i*stride) = 0.d0
          enddo

c         ------------------------------------------------
c         ... move column jcol from front into panel
c             note that k points to the diagonal entry in
c             the front and kk points to the corresponding
c             the panel.
c         ------------------------------------------------

c.debug
c     write(6,'("in xdsle2 - jcol, colbgn, colend  = ", 4i8)')
c    1                        jcol, colbgn, colend  
c     write(6,'("            length, k, kk, stride = ", 4i8)')
c    1                        length, k, kk, stride
c.debug
          do j=0,length-1
            panel(kk+j*stride) = front(k+j)
          enddo

          k      = k      + length
          length = length - 1

          if ( cmajor .eq. 0 ) then
              kk = kk + pnlcol + 1
          else
              kk = kk + pnlrow + 1
          end if
          
      enddo

c     ------------------------------------------------------------------
 
      return
      end 
