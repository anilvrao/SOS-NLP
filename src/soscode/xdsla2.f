      subroutine xdsla2 ( nnode, pospon, loclfr, front ) 
c
c  purpose -- expand frontal matrix to make room for postponed columns.
c             in-core version.
c
c             the postponed columns are to be inserted after the first
c             nnode columns.
c
c  created   -- 05-feb-97, rgg
c  revisions -- 
c
c  input variables
c  ---------------
c
c      nnode  -- number of original columns that can be eliminated for
c                this front
c      pospon -- number of postponed columns 
c      loclfr -- size of the original front
c
c  input/output variables
c  ----------------------
c 
c      front  -- frontal matrix
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           nnode,  pospon, loclfr
 
      double precision  front(*)
 
c     -------------------
c     ... local variables
c     -------------------

      integer           jcol  , kcur  , kfnl  , l     , lfront, lupdat,
     1                  sizcur, sizfnl, sizupd, space
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      external          xdslmv
 
c  =====================================================================

c     ----------------------------------------------------------------
c     ... note:  that since this subroutine expands the size of the
c                frontal matrix but the initial position in memory
c                stays fix, the data movement starts at the end of the
c                current front and finishes at the beginning.
c                kcur is the begining position of the org. data to be
c                moved.
c                kfnl is the begining position where the data is to be
c                moved to.
c     ----------------------------------------------------------------

      lupdat = loclfr - nnode
      lfront = loclfr + pospon

c     ---------------------------
c     ... move update matrix down
c     ---------------------------

      sizupd = ( lupdat * ( lupdat + 1 ) ) / 2
      sizcur = ( loclfr * ( loclfr + 1 ) ) / 2
      sizfnl = ( lfront * ( lfront + 1 ) ) / 2

      kcur   = sizcur - sizupd + 1
      kfnl   = sizfnl - sizupd + 1

      call xdslmv ( sizupd, front, kcur, kfnl )

c     -----------------------------------------
c     ... adjust kfnl for postponed columns and
c         zero out that portion of the front
c     -----------------------------------------

      space = ( lupdat + pospon ) * pospon - ( pospon*(pospon-1) ) / 2 
      kfnl  = kfnl - space
    
      front(kfnl:kfnl+space-1) = 0.d0

c     ------------------------------------------------------
c     ... expand each of the original nnode columns to allow
c         for the extra pospon rows in each
c     ------------------------------------------------------

      do jcol = nnode, 1, -1

c         ---------------------------------------------------
c         ... move part of column associated with update down
c         ---------------------------------------------------

          kfnl = kfnl - lupdat
          kcur = kcur - lupdat

          call xdslmv ( lupdat, front, kcur, kfnl )

c         -----------------------------------------
c         ... adjust kfnl for postponed columns and
c             zero out that portion of the front
c         -----------------------------------------

          kfnl = kfnl - pospon

          front(kfnl:kfnl+pospon-1) = 0.d0

c         ---------------------------------------------------
c         ... move part of column associated with first nnode
c             columns
c         ---------------------------------------------------

          l    = nnode - jcol + 1
          kfnl = kfnl - l
          kcur = kcur - l

          call xdslmv ( l, front, kcur, kfnl )

      enddo

c.debug
      if ( kcur .ne. 1 .or. kfnl .ne. 1 ) then
          write(6,'("oops at end of xdsla2 - kcur, kfnl = ", 2i8)')
     1                                        kcur, kfnl
          stop
      end if
c.debug
 
c  =====================================================================

      return
      end
