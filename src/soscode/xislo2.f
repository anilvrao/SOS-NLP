      subroutine   xislo2   ( mdnode, xadj,   adjncy, dhead,  dforw,
     1                        dbakw,  xrowls, rowlst, cmpmap, qsize, 
     2                        llist,  marker, maxint, tag  ,  neqns, 
     3                        nadj ,  neqp1 , xcliqu, clqsiz, sparnt  )
 
c
c     ==================================================================
c     ==================================================================
c     ====  xislo2 -- multiple min deg node elimination             ====
c     ==================================================================
c     ==================================================================
c
c     created by joseph w. h. liu, york university
c     last modification   mar. 05, 1986
c     last modification   mar. 25, 1993, cca
c
c     purpose: this routine eliminates the node 'mdnode' of
c              minimum degree from the adjacency structure, which
c              is stored in the quotient graph format.  it also
c              transforms the quotient graph representation of the
c              elimination graph.
c
c     input parameters:
c         mdnode - node of minimum degree.
c         maxint - estimate of maximum representable integer.
c         tag    - tag value.
c
c     updated parameters:
c         (xadj, adjncy) - updated adjacency structure.
c         (dhead, dforw, dbakw) - degree doubly linked structure.
c         qsize - size of supernode.
c         marker - marker vector.
c         llist - temporary linked list of eliminated nbrs.
c         xcliqu - pointers to the first clique stored in each
c                  uneliminated node's adjacency set
c         clqsiz - elimination clique sizes
c         sparnt - eventually transformed into the parent vector
c                  for the supernodal elimination tree
c
c     reason for 93-mar-25 modification :
c     we added the capability to compress the graph prior to
c     the minimum degree ordering. this means that indistinguishable
c     nodes with respect to the original data structure are
c     recognized and the adjancency structures for each representative
c     node was compressed. therefore the original indices are found
c     in adjncy(xadj(node):xadj(node+1)-1), but they may not be
c     consume the entire array fragment. their end is denoted
c     by a zero adjncy(i) value, as usual. prior to this mod
c     the code assumed that all indices in the range were valid
c
c     the rest of the code was slightly modified for clarity, imho,
c     in case i ever have to deal with it again.
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      integer             mdnode, maxint, tag, neqns, nadj, neqp1
 
      integer             adjncy (*), dhead  (*), dforw  (*),
     1                    dbakw  (*), qsize  (*), llist  (*),
     2                    marker (*), xcliqu (*), clqsiz (*),
     3                    sparnt (*), xrowls (*), rowlst (*),
     4                    cmpmap (*)
 
      integer             xadj (*)
 
c     --------------------
c     ... global variables
c     --------------------
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             clqnbr, elmnt , i     , istop ,
     1                    istrt , j     , jstop , jstrt , k     ,
     2                    kstrt , kstop , link  , lstloc, node  ,
     3                    nodnbr, nqnbrs, nxnode, oldpnt, pvnode,
     4                    rloc  , rlmt  , rnode , xqnbr , k1    ,
     5                    k3

c.debug
c     integer             k2
c.debug
 
c     ==================================================================
 
c     ------------------------------------------------------------------
c     ... data structures
 
c     adjncy -- quotient graph representation of the reduced matrix
c               (see george and liu for details)
c         adjncy -- > 0 -- node number
c                   = 0 -- end of list
c                   < 0 -- pointer to another list
 
c     dhead, dforw and dbakw  -- doubly linked list of nodes by degree
c                                with overloaded meanings as below
 
c         dforw -- > 0 -- next node in list
c                  = 0 --
c                  < 0 -- -inverse permutation position of node
c                         (means node has been eliminated)
 
c         dbakw -- > 0 -- previous node
c                  = 0 -- node needs degree updating
c                  < 0 -- -degree of node (pointer to head of list)
c                  = -maxint -- node does not require degree update
c                               and is outmatched by another node
 
c         tag   -- used to note stage at which node has been marked
c                  (multiple values used to reduce amount of resetting
c                   required)
 
c         marker -- = maxint -- node has been combined into a clique
c                               either an eliminated superelement or
c                               a clique of indistinguishable nodes
c                               these cause further overloadings of
c                               dforw and dbakw
c     ------------------------------------------------------------------
c.debug
c     write(6,'("in xislo2")')
c.debug
 
c     =====================================
c     tag mdnode's quotient graph neighbors
c     =====================================
c
c     ------------------
c     set up loop limits
c     ------------------
c
      istrt = xadj(mdnode)
      istop = xadj(mdnode+1) - 1
      kstrt = xcliqu(mdnode)
      kstop = xadj(mdnode+1) - 1
      if  ( kstrt .ne. 0 ) then
         istop = kstrt - 1
      endif
c
c     ---------------------------------------------------------
c     tag node neighbors that do not belong to adjacent cliques
c     ---------------------------------------------------------
c
      marker(mdnode) = tag
      do i = istrt, istop
         nodnbr = adjncy(i)
         marker(nodnbr) = tag
      enddo
c
c     ---------------------------------------------------
c     tag clique neighbors and add each to linked list so
c     they can be merged into mdnode's elimination clique
c     ---------------------------------------------------
c
      elmnt = 0
      if  ( kstrt .ne. 0 )  then
         k = kstrt
  200    continue
         if ( k .le. kstop ) then
            clqnbr = adjncy(k)
            if ( clqnbr .gt. 0 ) then
               marker(clqnbr) = tag
               llist(clqnbr)  = elmnt
               elmnt          = clqnbr
               sparnt(clqnbr) = mdnode
               k = k + 1
               go to 200
            endif
         endif
c
c        ====================================
c        generate mdnode's elimination clique
c        (not necessary if mdnode has not yet
c        appeared in an elimination clique)
c        ====================================
c
c        ---------------------------------------------------------
c        ... merge into the list all nodes reachable from already
c            eliminated generalized elements adjacent to  'mdnode'
c            if we run out of space for the list, we will continue
c            the list in the space occupied by these elements.
c        ---------------------------------------------------------
c
         rlmt = kstop
         rloc = kstrt
  400    continue
         if  ( elmnt .gt. 0 )  then
            adjncy(rlmt) = - elmnt
            link         = elmnt
  500       continue
            jstrt = xadj(link)
            jstop = xadj(link+1) - 1
            j = jstrt
  800       continue
            if ( j .le. jstop ) then
               node = adjncy(j)
               if ( node .gt. 0 ) k1 = cmpmap(node)
c.debug
c             if ( node .gt. 0 ) then
c             k2 = rowlst(xrowls(k1))
c             if ( k2 .ne. node ) then
c                 write(6,'("oops no. 04 - node, k1, k2 = ", 3i8)')
c    1                                      node, k1, k2
c             end if
c             end if
c.debug
               if ( node .lt. 0 ) then
c
c                 ----------------------------
c                 move to other array fragment
c                 ----------------------------
c
                  link = - node
                  go to 500
               else if ( node .eq. 0 ) then
c
c                 ------------
c                 end of array
c                 ------------
c
                  go to 600
               else if ( marker(node) .eq. tag
     1              .or. qsize(k1) .eq. 0 ) then
c
c                 -------------------------------------
c                 node already present or else inactive
c                 -------------------------------------
c
                  j = j + 1
                  go to 800
               else
c
c                 ---------------------------------------------
c                 ... new node -- neither marked nor eliminated
c                     add it to list, using storage from
c                     eliminated nodes if necessary
c                 ---------------------------------------------
 
                  marker(node) = tag
  700             continue
                  if  ( rloc .ge. rlmt ) then
c
c                    ------------------------------------
c                    end of array fragment to hold
c                    boundary nodes, move to new fragment
c                    ------------------------------------
c
                     link = - adjncy(rlmt)
                     rloc = xadj(link)
                     rlmt = xadj(link+1) - 1
                     go to 700
                 endif
                 adjncy(rloc) = node
                 rloc = rloc + 1
                 j = j + 1
                 go to 800
               endif
            endif
c
c           ----------------------------------------------------
c           ... adjacency list for generalized element exhausted
c               move to next element
c           ----------------------------------------------------
c
  600       continue
            elmnt = llist(elmnt)
            go to 400
         endif
c
c        ------------------------------------
c        tie up the tail of the boundary list
c        ------------------------------------
c
         if ( rloc .le. rlmt ) then
             adjncy(rloc) = 0
         endif
      endif
c
c=======================================================================
c
c     -----------------------------------------------------
c     ... reachable set has been built.  for each node in
c         reachable set, update its adjacency list to point
c         to the new superelement and remove links to nodes
c         within the superelement.  flag all nodes in
c         reachable set for degree update.
c     -----------------------------------------------------
c
      link = mdnode
 1100 continue
      istrt = xadj(link)
      istop = xadj(link+1) - 1
      i = istrt
 1700 continue
      if ( i .le. istop ) then
         rnode = adjncy(i)
         if ( rnode .lt. 0 ) then
c
c           ---------------------------
c           move to next array fragment
c           ---------------------------
c
            link = - rnode
            go to 1100
         else if ( rnode .eq. 0 ) then
c
c           -----------------
c           end of reach list
c           -----------------
c
            continue
         else
c
c           --------------------
c           ... an ordinary node
c           --------------------
c
            pvnode = dbakw(rnode)
            if  ( pvnode .ne. 0  .and.  pvnode .ne. - maxint ) then
c
c              ------------------------------------------------
c              ... rnode is still in the degree list, remove it
c                  because its degree is no longer correct.
c              ------------------------------------------------
c
               nxnode = dforw(rnode)
               if ( nxnode .gt. 0 ) then
                  dbakw(nxnode) = pvnode
               endif
               if ( pvnode .gt. 0 ) then
                  dforw(pvnode) = nxnode
               else
                  dhead(-pvnode) = nxnode
               endif
            endif
 
c           -------------------------------------------
c           remove uneliminated neighbors of rnode that
c           also lie in mdnode's elimination clique
c           -------------------------------------------
 
 1300       continue
            jstrt  = xadj(rnode)
            jstop  = xadj(rnode+1) - 1
            lstloc = jstop
            oldpnt = xcliqu(rnode)
            if ( oldpnt .ne. 0 ) then
               jstop = oldpnt - 1
            endif
            xqnbr = jstrt
            if ( jstrt .le. jstop ) then
               do j = jstrt, jstop
                  nodnbr = adjncy(j)
                  if ( marker(nodnbr) .lt. tag ) then
c
c                    --------------------------------------
c                    nodnbr not part of new clique, keep it
c                    --------------------------------------
c
                     adjncy(xqnbr) = nodnbr
                     xqnbr = xqnbr + 1
                  endif
               enddo
            endif
c
c           ------------------------------------------
c           reset pointer to clique neighbors of rnode
c           ------------------------------------------
c
            xcliqu(rnode) = xqnbr
c
c           -------------------------------------------------
c           remove from rnode's adjacency list any cliques
c           that were merged into mdnode's elimination clique
c           -------------------------------------------------
 
            if  ( oldpnt .ne. 0 )  then
               kstrt  = oldpnt
               kstop  = lstloc
               k = kstrt
 1400          continue
               if ( k .le. kstop ) then
                  clqnbr = adjncy(k)
                  if ( clqnbr .gt. 0 ) then
                     if ( marker(clqnbr) .lt. tag ) then
c
c                       ------------------------------------------
c                       clqnbr not absorbed by new clique, keep it
c                       ------------------------------------------
c
                        adjncy(xqnbr) = clqnbr
                        xqnbr = xqnbr + 1
                     endif
                     k = k + 1
                     go to 1400
                  endif
               endif
            endif
c
c           ------------------------------------------
c           add mdnode to the list of adjacent cliques
c           ------------------------------------------
c
            adjncy(xqnbr) = mdnode
            xqnbr = xqnbr + 1
            if ( xqnbr .le. lstloc ) then
               adjncy(xqnbr) = 0
            endif
            nqnbrs = xqnbr - xadj(rnode)
            if  ( nqnbrs .eq. 1 )  then
c
c              ----------------------------------------------
c              rnode is connected to only one clique and no other nodes
c              therefore it is indistinguishable to mdnode,
c              ----------------------------------------------

              k1 = cmpmap(mdnode)
              k3 = cmpmap(rnode)

c.debug
c             k2 = rowlst(xrowls(k1))
c             if ( k2 .ne. mdnode ) then
c                 write(6,'("oops no. 05 - mdnode, k1, k2 = ", 3i8)')
c    1                                      mdnode, k1, k2
c             end if
c.debug
c.debug
c             k2 = rowlst(xrowls(k1))
c             if ( k2 .ne. rnode ) then
c                 write(6,'("oops no. 06 - rnode, k3, k2 = ", 3i8)')
c    1                                      rnode, k3, k2
c             end if
c.debug
               qsize(k1)      = qsize(k1)      + qsize(k3)
               clqsiz(mdnode) = clqsiz(mdnode) - qsize(k3)
               qsize(k3)      = 0
               marker(rnode)  = maxint
               dforw(rnode)   = - mdnode
               dbakw(rnode)   = - maxint
            else
c
c           ------------------------------
c           flag 'rnode' for degree update
c           ------------------------------
c
               dforw(rnode)  = nqnbrs
               dbakw(rnode)  = 0
            endif
c
c           -------------------------------------
c           move onto next node in the reach list
c           -------------------------------------
c
            i = i + 1
            go to 1700
         endif
      endif
c
c     ------
c     return
c     ------
c
      return
      end
