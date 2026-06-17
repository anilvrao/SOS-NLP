      subroutine  xisls2  ( nnodes, nladj , xladj , ladjcy, nsuper,
     .                      xsup  , nsubs , xlindx, perm  , 
     .                      cmpmap, xrowls, rowlst, lgindx, lrindx,
     .                      qarndx, arindx, mrglnk, collnk, marker,
     .                      mxacol, alist , supmap, relloc, error   )
 
 
c    ==================================================================
c    ==================================================================
c    ======                                                      ======
c    ====== xisls2 - multifrontal symbolic factorization         ======
c    ======                                                      ======
c    ==================================================================
c    ==================================================================
c
c     creation date: 10-07-87 (bwp)
c     last updated: 10-22-87 (bwp)
c     last updated: 4-88 (bwp)   jess and kees modifications,
c                                 subsequently removed except for
c                                 'right-side up' subscripts
c     last updated: 05-10-88 (bwp) noetic multifrontal modifications
c     last updated: 06-17-88 (bwp) modified for general purpose
c                                  multifrontal code
c     last updated: 01-30-89 (rgg) kept xladj in nodal form instead
c                                  of supernodal form.
c     last updated: 11-02-90 (mlc) changed nsubs usage in dimensioning
c                                  to use of *
c     last updated: 12-17-96 (rgg) converted to use compressed form of
c                                  (xladj, ladjcy) 
c
c     purpose -- to perform a symbolic factorization in preparation for
c                a multifrontal numerical factorization and solve.
c                it generates global and local supernodal row indices
c                stored column-wise in ascending order for the factor
c                matrix. it generates local indices for the lower
c                triangle of pap'.
c                for some finite element matrices each node may
c                actually represent more than one node in the
c                original graph.  in this case, the nodes are
c                a priori supernodes of the same size. this
c                has no effect on this routine at all, but it does
c                mean that the output from this routine must be expanded
c                by the companion routine symexp to be use by the
c                the multifrontal factorization and solve modules.
c
c     ***** note that the indices are stored in ascending
c           order--the usual arrangement
c
c     ***** also note that collnk(*) requires nnodes+1 words of storage
c           and is dimensioned as collnk(0:nnodes)
c
c     ***** note that each set of supernodal indices includes
c           only the 'external' indices, i.e. those not actually in
c           the supernode itself.
c
c     input parameters:
c
c        nnodes  - number of nodes (a priori supernodes)
c        nladj   - dimension of ladjcy
c       (xladj,
c        ladjcy) - adjacency structure of the lower triangle of
c                  the coefficient matrix pap', modulo some uniform
c                  a priori supernode size
c        xsup    - supernode partition of the columns of l modulo some
c                  uniform a priori supernode size
c        nsubs   - number of supernodal row subscripts for l, modulo
c                  some uniform a priori supernode size.
c                  dimension of lindx.
c        nsuper  - number of supernodes
c        xlindx  - pointers to the indices of l
c        perm    - new to old permutation
c        cmpmap  - compression map - rows to compressed nodes
c       (xrowls,
c        rowlst) - row lists for each compressed node
c 
c        mxacol  - maximum number of nonzeros in the strict lower
c                  triangle of pap'
c
c     output parameters:
c
c        lgindx  - supernodal global indices for factor l, modulo
c                  some uniform a priori supernode size
c        lrindx  - supernodal relative indices of factor l, modulo
c                  some uniform a priori supernode size
c        error   - error code
c                     0 = no error
c                     1 = called sort routine xislq2 with a negative
c                         list length
c
c     work parameters:
c
c        mrglnk - maintains list of columns to merge
c        collnk - linked list in which current column is built
c        marker - marks nodes as merged into current column
c        alist  - collects indices from pap' for sorting
c        supmap - maps each node to its supernode
c        relloc - relative locations of each index in a frontal matrix
c
c     external subprograms:
c
c        xislq2 - sorts integer list in descending order
c
c     ==================================================================
 
c     ----------
c     parameters
c     ----------
 
      integer             nnodes, nladj , nsubs , nsuper,
     .                    error , mxacol
 
      integer             xladj(nnodes+1) , ladjcy(*) ,
     1                    xsup(*)         , arindx(*) ,
     2                    xlindx(*)       , lgindx(*) ,
     3                    mrglnk(*)       , collnk(0:nnodes),
     4                    marker(nnodes)  , alist(*) ,
     5                    supmap(nnodes)  , lrindx(*) ,
     6                    relloc(nnodes)

      integer             perm(*)         , cmpmap(*)       ,
     2                    xrowls(*)       , rowlst(*)
 
      logical             qarndx
 
c     ---------------
c     local variables
c     ---------------
 
      integer             colsiz, frtlen, fstloc, fstnod, inode ,
     1                    isuper, jsuper, ksuper, locatn, lstlen,
     2                    lstloc, lstnod, newind, nindex, nnodp1,
     3                    nxtind, nxtloc, nxtsup, offset, prvind,
     4                    rloc  , rowind, supsiz

      integer             icomp , k     , kbgn  , kend  , newnod
 
c     --------------------
c     external subprograms
c     --------------------
 
      external            xislq2
 
c     ==================================================================
 
c     --------------
c     initialization
c     --------------
 
      nnodp1    = nnodes + 1
      error     = 0
 
      do inode = 1, nnodes
          mrglnk(inode) = 0
          marker(inode) = 0
      enddo
 
c     --------------------------------
c     initialize node to supernode map
c     --------------------------------
 
      do isuper = 1, nsuper
          do inode = xsup(isuper), xsup(isuper+1)-1
              supmap(inode) = isuper
          enddo
      enddo
 
c     ------------------------
c     for supernode jsuper ...
c     ------------------------
 
      do  2000  jsuper = 1, nsuper
 
c         -----------------------------
c         get supernode info for jsuper
c         -----------------------------
 
          fstnod = xsup(jsuper)
          lstnod = xsup(jsuper+1) - 1
          supsiz = lstnod - fstnod + 1
          colsiz = xlindx(jsuper+1) - xlindx(jsuper)
 
c         ------------------------------------------
c         create null list for row indices of jsuper
c         and initialize index counter
c         ------------------------------------------
 
          collnk(0) = nnodp1
          nindex    = 0
 
c         --------------------------------------------
c         skip merges if there are no external indices
c         for supernode jsuper
c         --------------------------------------------
 
          if  ( colsiz .eq. 0 )  go to 1150
 
c         -----------------------------------------------------
c         first merge: copy the row indices of supernode ksuper
c                      into the null list for jsuper
c         -----------------------------------------------------
 
          ksuper = mrglnk(jsuper)
 
          if  ( ksuper .eq. 0 )  go to 800
 
          fstloc = xlindx(ksuper)
          lstloc = xlindx(ksuper+1) - 1
          nindex = 0
 
          do locatn = lstloc, fstloc, -1
              rowind = lgindx(locatn)
              if  ( rowind .gt. lstnod )  then
                  nindex         = nindex + 1
                  collnk(rowind) = collnk(0)
                  collnk(0)      = rowind
                  marker(rowind) = jsuper
              else
                  go to 200
              endif
          enddo
 
c         ----------------------------------
c         quit if no more indices are needed
c         ----------------------------------
 
  200     continue
          if  ( nindex .eq. colsiz )  go to 1150
 
c         -----------------------------------------------------
c         subsequent merges: merge the row indices of supernode
c                            ksuper into the current list for
c                            jsuper
c         -----------------------------------------------------
 
          ksuper = mrglnk(ksuper)
          if  ( ksuper .eq. 0 )  go to 800
 
          fstloc = xlindx(ksuper)
          lstloc = xlindx(ksuper+1) - 1
 
          prvind = 0
          nxtind = collnk(0)
 
          do  700  locatn = fstloc, lstloc
 
              newind = lgindx(locatn)
 
              if  ( newind .gt. lstnod )  then
 
c.obs
c 300             if  ( newind - nxtind ) 400, 500, 600
c.obs

  300             continue
                  if  ( newind .eq. nxtind ) go to 500
                  if  ( newind .gt. nxtind ) go to 600
 
c                 --------------------------------------------------
c                 new index lt next index: link in new index so that
c                 the sequence -- prvind, newind, nxtind -- appears
c                 in the list, then proceed to the next new index
c                 if the list is not yet complete
c                 --------------------------------------------------
 
  400             continue
                  collnk(newind) = nxtind
                  collnk(prvind) = newind
                  prvind         = newind
                  marker(newind) = jsuper
 
                  nindex         = nindex + 1
                  if  ( nindex .eq. colsiz )  go to 1150
                  go to 700
 
c                 ---------------------------------------------------
c                 new index eq next index: bump the linked list, then
c                 proceed to the next new index
c                 ---------------------------------------------------
 
  500             continue
                  prvind = nxtind
                  nxtind = collnk(nxtind)
                  go to 700
 
c                 ---------------------------------------------------
c                 new index gt next index: bump the linked list, then
c                 try the current new index again
c                 ---------------------------------------------------
 
  600             continue
                  prvind = nxtind
                  nxtind = collnk(nxtind)
                  go to 300
 
              endif
 
  700     continue
 
          go to 200
 
c         ---------------------------------------------------------
c         merge into the list any row indices from the columns of
c         the lower triangle of pap' (in supernode jsuper) that are
c         not already in the list
c         ---------------------------------------------------------
 
c         -----------------
c         1. build the list
c         -----------------
 
  800     continue
 
          nxtloc = 1
 
          fstloc = xladj(fstnod)
          lstloc = xladj(lstnod+1) - 1
c.debug
c     write(6,'("before 900 - fstloc, lstloc = ", 2i8)')
c    1                         fstloc, lstloc 
c     call xislp3 ( 'ladjcy', lstloc-fstloc+1, ladjcy(fstloc), 6 )
c.debug
 
          do  locatn = fstloc, lstloc
              newnod = ladjcy(locatn)
c.debug
c     write(6,'("in 900 - jsuper, nxtloc, locatn, newnod = ", 4i8)')
c    1                     jsuper, nxtloc, locatn, newnod 
c.debug
              icomp  = cmpmap(newnod)
              kbgn   = xrowls(icomp)
              kend   = xrowls(icomp+1) - 1
c.debug
c     write(6,'("pre 850- oldnod, icomp, kbgn, kend      = ", 4i8)')
c    1                     oldnod, icomp, kbgn, kend      
c.debug

              do k = kbgn, kend

                  newind = rowlst(k)
c.debug
c     write(6,'("in 850 - k, newind = ", 3i8)')
c    1                     k, newind 
c.debug

                  if ( marker(newind) .eq. jsuper .or.
     1                 newind         .le. lstnod     )  cycle

                      marker(newind) = jsuper
                      alist(nxtloc)  = newind
                      nxtloc         = nxtloc + 1
                      nindex         = nindex + 1
                      if  ( nindex .eq. colsiz )  go to 950

              enddo

          enddo
 
  950     continue
          lstlen = nxtloc - 1
c.debug
c     write(6,'("after 950 - lstlen = ", i8)')
c     call xislp3 ( 'alist', lstlen, alist, 6 )
c.debug
 
c         ----------------
c         2. sort the list
c         ----------------
 
          if  ( lstlen .ge. 2 )  then
 
c.debug
c     write(6,'("before sort - lstlen = ", 2i8)')
c     call xislp3 ( 'alist', lstlen, alist, 6 )
c.debug
              call xislq2 ( lstlen, alist, error )
c.debug
c     write(6,'("after  sort - lstlen = ", 2i8)')
c     call xislp3 ( 'alist', lstlen, alist, 6 )
c.debug
 
              if  ( error .ne. 0 )  then
                  error = 1
                  return
              endif
 
          endif
 
c         -----------------
c         3. merge the list
c         -----------------
c.debug
c     write(6,'("before merge - lstlen = ", 2i8)')
c     call xislp3 ( 'alist', lstlen, alist, 6 )
c.debug
 
          if  ( lstlen .ge. 1 )  then
 
              prvind = 0
              nxtind = collnk(0)
 
              do locatn = lstlen, 1, -1
 
                  newind = alist(locatn)
 
 1000             continue
                  if  ( newind .lt. nxtind )  then
 
c                     ----------------------------------
c                     newind lt nxtind: insert new index
c                     in front of next index
c                     ----------------------------------
 
                      collnk(newind) = nxtind
                      collnk(prvind) = newind
                      prvind         = newind
                      marker(newind) = jsuper
c.debug
c     write(6,'("insert  - prvind, newind, nxtind = ", 3i8)') 
c    1                      prvind, newind, nxtind
c.debug
 
                  else
 
c                     ----------------------------------
c                     newind gt nxtind: bump linked list
c                     and try to insert new index again
c                     ----------------------------------
 
                      prvind = nxtind
                      nxtind = collnk(nxtind)
c.debug
c     write(6,'("bumping - prvind, newind, nxtind = ", 3i8)') 
c    1                      prvind, newind, nxtind
c.debug
                      go to 1000
 
                  endif
 
              enddo
 
          endif
 
c         -------------------------------------
c         scatter the relative locations of the
c         internal indices
c         -------------------------------------
 
 1150     continue
          rloc   = 1
          do nxtind = fstnod, lstnod
              relloc(nxtind) = rloc
              rloc           = rloc + 1
          enddo
 
c         ------------------------------------
c         copy column's index list from linked
c         list to permanent data structure
c         and scatter relative locations
c         ------------------------------------
 
          nxtloc = xlindx(jsuper)
          nxtind = collnk(0)
 
 1200     continue
          if  ( nxtind .lt. nnodp1 )  then
c.debug
c     write(6,'("after 1200 - nxtloc, nxtind = ", 2i8)')
c    1                         nxtloc, nxtind 
c.debug
              lgindx(nxtloc) = nxtind
              relloc(nxtind) = rloc
              nxtloc         = nxtloc + 1
              rloc           = rloc   + 1
              nxtind         = collnk(nxtind)
              go to 1200
          endif
 
c         ------------------------------------------------
c         compute the relative row indices for each column
c         of l merged to form the column just constructed
c         ------------------------------------------------
 
          ksuper = mrglnk(jsuper)
 1400     continue
          if  ( ksuper .gt. 0 )  then
              fstloc = xlindx(ksuper)
              lstloc = xlindx(ksuper+1) - 1
cdir$         ivdep
              do locatn = fstloc, lstloc
                  newind = lgindx(locatn)
                  lrindx(locatn) = relloc(newind)
              enddo
              ksuper = mrglnk(ksuper)
              go to 1400
          endif
 
c         ----------------------------------------------------
c         compute the relative indices for the columns of pap'
c         ----------------------------------------------------
 
          if ( qarndx ) then
 
              offset = 0
              frtlen = colsiz + supsiz - 1
              do inode = fstnod, lstnod
                  fstloc = xladj(inode)
                  lstloc = xladj(inode+1) - 1
cdir$             ivdep
                  do locatn = fstloc, lstloc
                      newind         = ladjcy(locatn)
                      arindx(locatn) = relloc(newind) + offset
c.debug 
c     write(6,'("in 1550 - locatn, newind, relloc(newind), offset = "
c    1  4i8)')              locatn, newind, relloc(newind), offset 
c.debug
                  enddo
                  offset = offset + frtlen
                  frtlen = frtlen - 1
              enddo
 
          end if
 
c         --------------------------------------------
c         place supernode jsuper in the merge list for
c         the supernode into whose index list its list
c         be merged
c         --------------------------------------------
 
          if  ( xlindx(jsuper+1) .gt. xlindx(jsuper) ) then
 
              nxtsup         = supmap( lgindx(xlindx(jsuper)) )
              mrglnk(jsuper) = mrglnk(nxtsup)
              mrglnk(nxtsup) = jsuper
 
          endif
 
 2000 continue
 
c     ------------------
c     normal termination
c     ------------------
 
      return
      end
