      subroutine xislo8( nsuper, parsup, fstsup, brosup, snroot, snsize,
     1                   uleng , fleng , maxzer, zeroes, fellow, lstfel,
     2                   lstson                                        )
 
c
c  purpose -- to obtain the relaxed supernode partition when a supernode
c             may contain logical zeroes in its data structure. the
c             maximum number for any given supernode is an input
c             parameter.
c             the tree of strict supernodes is input, the tree of
c             relaxed supernodes is returned.
c
c  created            -- 20-jul-87, cca
c  last modifications -- 20-jul-87, cca
c
c  input variables --
c
c      nsuper -- number of strict supernodes
c      parsup -- parent supernode vector
c      fstsup -- first son supernode vector
c      brosup -- brother supernode vector
c      snroot -- root of the supernode elimination forest
c      snsize -- number of nodes in each supernode
c      uleng  -- size of the update matrices
c      fleng  -- size of the frontal matrices
c      maxzer -- maximum number of zeroes in each supernode
c
c  updated variables
c
c      nsuper -- number of relaxed supernodes
c      parsup -- parent supernode vector
c      fstsup -- first son supernode vector
c      brosup -- brother supernode vector
c      snroot -- root of the supernode elimination forest
c      snsize -- number of nodes in each supernode
c      uleng  -- size of the update matrices
c      fleng  -- size of the frontal matrices
c
c  working variables --
c
c      lstfel -- last fellow link vector
c      lstson -- last son link vector
c
c  output variables --
c
c      zeroes -- number of zeroes in each supernode's factor entries
c      fellow -- link vector for strict supernodes merged into
c                a relaxed supernode.
c
c  subprograms called
c
c  =====================================================================
 
      integer   nsuper, parsup(*), fstsup(*), brosup(*),
     1          snroot, snsize(*), uleng(*) , fleng(*),
     2          maxzer
 
      integer   lstfel(*), lstson(*)
 
      integer   zeroes(*), fellow(*)
 
      integer   extra , fstgsn, gsn   , isuper, least , lstgsn, maxsiz,
     1          mrgson, son   , jhmcon
 
      external  jhmcon
 
c  =====================================================================
 
c  --------------------------------------------
c  initialize zeroes, fellow, lstfel and lstson
c  --------------------------------------------
 
      do isuper = 1,nsuper
          zeroes(isuper) = 0
          fellow(isuper) = 0
          if ( ( brosup(isuper) .eq. 0 ) .and.
     1         ( parsup(isuper) .ne. 0 ) ) then
              lstson(parsup(isuper)) = isuper
          endif
          lstfel(isuper) = isuper
      enddo
 
c  -----------------------------------------------------------
c  perform a post-order traversal of the strict supernode tree
c  -----------------------------------------------------------
 
      isuper = snroot
      fstgsn = 0
 
   10 continue
 
      if ( fstsup(isuper) .gt. 0 ) then
          isuper = fstsup(isuper)
          go to 10
      endif
 
   20 continue
 
c  -----------------------------------------
c  visiting a strict supernode, check to see
c  if it can be merged with any of its sons.
c  -----------------------------------------
 
      mrgson = 0
      least  = jhmcon(3)
      maxsiz = 0
 
      son = fstsup(isuper)
   30 continue
      if ( son .gt. 0 ) then
 
          extra = snsize(son)*(fleng(isuper) - uleng(son))
          if ( extra+zeroes(son)+zeroes(isuper) .le. maxzer ) then
 
c             ---------------------------------------
c             son is a candidate to merge with isuper
c             ---------------------------------------
 
              if ( extra .lt. least ) then
                  mrgson = son
                  least  = extra
              elseif ( extra .eq. least ) then
                  if ( snsize(son) .gt. maxsiz ) then
                      mrgson = son
                      maxsiz = snsize(son)
                  endif
              endif
 
          endif
 
          son = brosup(son)
          go to 30
 
      endif
 
      if ( mrgson .gt. 0 ) then
 
c         ------------------------
c         merge isuper with mrgson
c         ------------------------
 
c         -----------------------
c         update the  merge links
c         -----------------------
 
          if ( lstfel(isuper) .eq. isuper) then
              fellow(isuper) = mrgson
          else
              fellow(lstfel(isuper)) = mrgson
          endif
          lstfel(isuper) = lstfel(mrgson)
 
c         -------------------------------------
c         remove mrgson from the supernode tree
c         -------------------------------------
 
          if ( mrgson .eq. fstsup(isuper) ) then
 
              fstsup(isuper) = brosup(mrgson)
              if ( mrgson .eq. lstson(isuper) ) lstson(isuper) = 0
 
          else
 
              son = fstsup(isuper)
   40         continue
              if ( brosup(son) .ne. mrgson ) then
                  son = brosup(son)
                  go to 40
              endif
              brosup(son) = brosup(mrgson)
              if ( mrgson .eq. lstson(isuper) ) lstson(isuper) = son
 
          endif
 
c         ---------------------------------------------
c         add the sons of mrgson to the grandchild list
c         ---------------------------------------------
 
          if ( fstsup(mrgson) .ne. 0 ) then
 
              if ( fstgsn .eq. 0 ) then
                   fstgsn = fstsup(mrgson)
              else
                   brosup(lstgsn) = fstsup(mrgson)
              endif
              lstgsn = lstson(mrgson)
 
           endif
 
c         -----------------------------------
c         adjust the information about isuper
c         -----------------------------------
 
          fleng(isuper)  = fleng(isuper)  + snsize(mrgson)
          zeroes(isuper) = zeroes(mrgson) + zeroes(isuper) + least
          snsize(isuper) = snsize(isuper) + snsize(mrgson)
 
c         ----------------------------------------
c         decrease the number of supernodes by one
c         ----------------------------------------
 
          nsuper = nsuper - 1
 
c         -----------------------------
c         try to merge with another son
c         -----------------------------
 
          go to 20
 
      else
 
c         --------------------------------------------
c         no more sons to merge with, attach grandsons
c         to the list of sons for isuper
c         --------------------------------------------
 
          if ( fstgsn .ne. 0 ) then
              if ( lstson(isuper) .eq. 0 ) then
                  fstsup(isuper) = fstgsn
              else
                  brosup(lstson(isuper)) = fstgsn
              endif
              lstson(isuper) = lstgsn
              gsn = fstgsn
   50         continue
              parsup(gsn) = isuper
              if ( brosup(gsn) .ne. 0 ) then
                  gsn = brosup(gsn)
                  go to 50
              endif
              fstgsn = 0
          endif
 
      endif
 
c  ------------------------------------------------
c  proceed to the next strict supernode in the tree
c  ------------------------------------------------
 
      if ( brosup(isuper) .ne. 0 ) then
          isuper = brosup(isuper)
          go to 10
      endif
 
      if ( parsup(isuper) .ne. 0 ) then
          isuper = parsup(isuper)
          go to 20
      endif
 
c  =====================================================================
 
      return
      end
