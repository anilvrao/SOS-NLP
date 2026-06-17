      subroutine xdslek ( cmajor, ncolf , colbgn, colend,
     1                    loclfr, front , pnlcol, pnlrow, panel )
c
c  purpose -- stuff postponed columns back into the front.
c
c  created            -- 09-apr-97, rgg
c  last modifications -- 
c
c  input variables --
c
c      cmajor -- column major flag
c                .eq. 0 row    major form for panel is used
c                .ne. 0 column major form for panel is used
c      ncolf  -- number of columns factored
c      colbgn -- first column of panel  
c      colend -- last  column of panel  
c      loclfr -- size of front
c      pnlcol -- number of columns in panel
c      pnlrow -- number of rows    in panel
c      panel  -- rectangular the symmetric frontal matrix in packed storage
c
c  output variables --
c      front  -- frontal matrix
c
c     ------------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer            cmajor, ncolf , colbgn, colend, loclfr, pnlcol, 
     1                   pnlrow
 
      double precision   front (*),      panel (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer            jcol  , k     , kk    , length, stride
 
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
c         ... move column jcol from panel back into front
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
          if ( jcol .gt. ncolf ) then
              do j=0,length-1
                front(k+j) = panel(kk+j*stride)
              enddo
          end if

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
