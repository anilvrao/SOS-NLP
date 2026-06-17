      subroutine   xislo6   ( ehead,  neqns,  xadj,   adjncy, delta,
     1                        mdeg,   dhead,  dforw,  dbakw,  xrowls, 
     2                        rowlst, cmpmap, qsize,  llist,  marker,
     3                        maxint, tag,    nadj,   neqp1,  xcliqu, 
     4                        clqsiz )
 
c
c     ==================================================================
c     ==================================================================
c     ====  xislo6 -- multiple min degree degree update             ====
c     ==================================================================
c     ==================================================================
c
c     created by joseph w. h. liu, york university
c     last modification   mar. 05, 1986
c     last modification   mar. 25, 1993, cca
c
c     purpose: this routine updates the degrees of nodes
c              after a multiple elimination step.
c
c     input parameters:
c         ehead - the beginning of the list of eliminated
c                 nodes (i.e. newly formed elements).
c         neqns - number of equations.
c         (xadj, adjncy) - adjacency structure.
c         xcliqu - pointers to the first clique in each uneliminated
c                  node's adjacency set
c         clqsiz - elimination clique size
c
c     updated parameters:
c         mdeg   - new minimum degree after degree update.
c         (dhead, dforw, dbakw) - degree doubly linked structure.
c         qsize  - size of supernode.
c         llist  - working linked list.
c         marker - marker vector for degree update.
c         tag    - tag value.
c
c     modification history :
c
c       april 1993 : some debug statements inserted,
c          but most important edges not yet covered by an
c          eliminated element are weighed individually
c          instead of assuming they have unit weight,
c          which is true if the graph has not been compressed.
c
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      integer             ehead, neqns, delta, mdeg, maxint, tag , nadj,
     1                    neqp1
 
      integer             adjncy (*), dhead  (*), dforw  (*),
     1                    dbakw  (*), qsize  (*), llist  (*),
     2                    marker (*), xcliqu (*), clqsiz (*),
     3                    xrowls (*), rowlst (*), cmpmap (*)
 
      integer             xadj (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             deg0  , deg   , elmnt , enode , fnode ,
     1                    fstclq, fstnbr, i     , istrt , istop ,
     2                    j     , jstrt , jstop , link  , mdeg0 ,
     3                    mtag  , nabor , node  , q2head, qxhead,
     4                    k1    , k3

c.debug
c     integer             k2
c.debug
 
      logical             iq2
 
c.d   integer             msglvl
c.d   parameter           (msglvl = 0)
 
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
c     write(6,'("in xislo6")')
c.debug
 
c     -----------------------------------------------------------
c     ... one or more nodes have been eliminated, creating one or
c         more new superelements, each containing the adjacency
c         lists for the eliminated node.  for each superelement,
c         update the degrees of every uneliminated neighbor.
c         also look for certain cases of indistinguishable nodes
c         and of outmatching nodes
c     -----------------------------------------------------------
c
      i = xadj(neqns+1) - 1
cdbg      print *,'Entry into xislo6:'
cdbg      print *,' ehead, neqns, delta, mdeg, maxint, tag, nadj, neqp1',
cdbg     1          ehead, neqns, delta, mdeg, maxint, tag, nadj, neqp1
cdbg      call xislp3 ( 'adjncy', i      , adjncy, 6 )
cdbg      call xislp3 ( 'dhead ', neqns, dhead , 6 )
cdbg      call xislp3 ( 'dforw ', neqns, dforw  , 6 )
cdbg      call xislp3 ( 'dbakw ', neqns, dbakw  , 6 )
cdbg      call xislp3 ( 'qsize ', 5, qsize , 6 )
cdbg      call xislp3 ( 'llist ', neqns, llist , 6 )
cdbg      call xislp3 ( 'marker', neqns, marker, 6 )
cdbg      call xislp3 ( 'xcliqu', neqns, xcliqu, 6 )
cdbg      call xislp3 ( 'clqsiz', neqns, clqsiz, 6 )
cdbg      call xislp3 ( 'xrowls', 5+1, xrowls, 6 )
cdbg      call xislp3 ( 'rowlst', neqns  , rowlst, 6 )
cdbg      call xislp3 ( 'cmpmap', neqns  , cmpmap, 6 )
cdbg      call xislp3 ( 'xadj  ', neqns+1, xadj  , 6 )

      
      mdeg0 = mdeg + delta
      elmnt = ehead
 
  100 continue
      if  ( elmnt .le. 0 )  go to 2400
 
      mtag = tag + mdeg0
      if  ( tag .lt. maxint - mdeg0 )  go to 200
 
c         ---------------------------------------------------
c         ... each member of a new superelement will have its
c             own tag value.  reset tag value if necessary to
c             avoid integer overflow.  reset markers of all
c             uneliminated nodes in this case.
c         ---------------------------------------------------
 
          tag = 1
cdbg      print *,' xislo6 before 150 neqns, maxint=',neqns,maxint
cdbg      call xislp3 ( 'marker', neqns, marker, 6 )
          do i = 1, neqns
             if ( marker(i) .lt. maxint ) then
                marker(i) = 0
             endif
          enddo
          mtag = tag + mdeg0
cdbg      print *,' xislo6 after 150 mtag,tag,mdeg0=',mtag,tag,mdeg0
cdbg      call xislp3 ( 'marker', neqns, marker, 6 )
 
c         --------------------------------------------------
c         create 2 linked lists from nodes associated with
c         elmnt': one with 2 nbrs (q2head) in adjacency
c         structure, and the other with more than 2 (qxhead)
c         neighbors.  we can afford to search the 2 nbr list
c         for special cases. while we create the lists, also
c         compute 'deg0', no of nodes in this element.
c         --------------------------------------------------
 
  200     continue
          q2head = 0
          qxhead = 0
          deg0   = 1 + clqsiz(elmnt)
c         write(*,*) 'elmnt ', elmnt, ', mtag', mtag, ', deg0 ', deg0
 
          link = elmnt
 
c         -----------------------------------
c         ... search through one superelement
c         -----------------------------------
 
  300     continue
          istrt = xadj(link)
          istop = xadj(link+1) - 1
          do 700 i = istrt, istop
cdbg      original code
cdbg              enode = adjncy(i)
cdbg              k1    = cmpmap(enode)
cdbg              link  = - enode
cdbgc.obs         if  ( enode )  300, 800, 400
cdbg              if  ( enode .lt. 0 ) go to 300
cdbg              if  ( enode .eq. 0 ) go to 800
cdbg      jtb version
              enode = adjncy(i)
              if  ( enode .lt. 0 ) then
                link = - enode
                go to 300
              elseif( enode .eq. 0 ) then
                link = 0
                go to 800
              else
                k1    = cmpmap(enode)
              endif
cdbg      jtb version
 
  400         continue
c.d           if ( msglvl .ge. 2 ) then
c.d              write(*,*) '   enode ', enode, ' qsize ', qsize(k1),
c.d  1                      ', dbakw ', dbakw(enode),
c.d  2                      ', dforw ', dforw(enode)
c.d           endif
c.debug
c             k2 = rowlst(xrowls(k1))
c             if ( k2 .ne. enode ) then
c                 write(6,'("oops no. 08 - enode, k1, k2 = ", 3i8)')
c    1                                      enode, k1, k2
c             end if
c.debug
cdbg      print *,' xislo6 after 400 k1=',k1
cdbg      call xislp3 ( 'qsize ', 5, qsize , 6 )
              if  ( qsize(k1) .eq. 0 )  go to 700
cdbg      print *,' xislo6 after 1f go to 700 enode,mtag=',enode,mtag
cdbg      call xislp3 ( 'marker', neqns, marker, 6 )
                  marker(enode) = mtag
 
c                 ------------------------------------------------
c                 ... when can the next test ever be taken????????
c                     when a node is outmatched?
c                 (see comment below, cca)
c                 ------------------------------------------------
cdbg      print *,' xislo6 before go to 700 enode=',enode
cdbg      call xislp3 ( 'dbakw ', neqns, dbakw  , 6 )
                  if  ( dbakw(enode) .ne. 0 )  go to 700
c
c                     --------------------------------------
c                     enode has not yet been outmatched
c                     nor has its degree been computed yet
c                     place either in qxhead or q2head lists
c                     --------------------------------------
c
                      if  ( dforw(enode) .eq. 2 )  go to 600
c.d                       if ( msglvl .ge. 2 ) then
c.d                          write(*,*) 'placing in qxhead'
c.d                       endif
                          llist(enode) = qxhead
                          qxhead       = enode
                          go to 700
 
  600                     continue
c.d                       if ( msglvl .ge. 2 ) then
c.d                          write(*,*) 'placing in q2head'
c.d                       endif
                          llist(enode) = q2head
                          q2head       = enode
 
  700     continue
 
c         --------------------------------------------
c         ... nodes in reachable set that need degree
c             updates have been divided into two lists
c         --------------------------------------------
 
  800     continue
          enode = q2head
          iq2   = .true.
 
c         ---------------------------------------------------
c         ... process all nodes with only two neighbors.
c             since one of the neighbors is the new
c             superelement, the node has only one non-trivial
c             neighbor.  look for superelements becoming
c             indistinguishable or nodes being outmatched by
c             another for these cheap cases
c         ---------------------------------------------------
 
  900     continue
          if  ( enode .le. 0 )  go to 1500
c.d       if ( msglvl .ge. 2 ) then
c.d          write(*,*) ' 2-list, enode ', enode,
c.d  1                  ', dbakw ', dbakw(enode)
c.d       endif
cdbg      print *,' xislo6 before go to 2200 enode=',enode
cdbg      call xislp3 ( 'dbakw ', neqns, dbakw  , 6 )
          if  ( dbakw(enode) .ne. 0 )  go to 2200
              tag = tag + 1
              deg = deg0
 
c             -----------------------------------------------
c             ... identify the other adjacent element 'nabor'
c             -----------------------------------------------
 
              istrt = xadj(enode)
              nabor = adjncy(istrt)
              if  ( nabor .eq. elmnt )  nabor = adjncy(istrt+1)
c.d           if ( msglvl .ge. 2 ) then
c.d              write(*,*) ' nabor ', nabor, ', dforw ', dforw(nabor)
c.d           endif
 
              link = nabor
              k3   = cmpmap(nabor)
              if  ( dforw(nabor) .lt. 0 )  go to 1000
 
c                 --------------------------------------------------
c                 ... if 'nabor' is uneliminated, increase deg count
c                 --------------------------------------------------
 
c                 write(*,*) ' adding ', qsize(k3), ' to degree'
c.debug
c             k2 = rowlst(xrowls(k3))
c             if ( k2 .ne. nabor ) then
c                 write(6,'("oops no. 09 - nabor, k3, k2 = ", 3i8)')
c    1                                      nabor, k3, k2
c             end if
c.debug
cdbg      print *,' xislo6 before go to 2100 k3=',k3
cdbg      call xislp3 ( 'qsize ', 5, qsize , 6 )
                  deg = deg + qsize(k3)
 
                  go to 2100
 
c                 -----------------------------------------------
c                 ... otherwise, neighbor is another superelement
c                     for each node in the 2nd element see if it
c                     is indistinguishable from  'enode'  or
c                     outmatched by  'enode'
c                 -----------------------------------------------
 
 1000             continue
                  istrt = xadj(link)
                  istop = xadj(link+1) - 1
                  do 1400 i = istrt, istop
cdbg      original code
cdbg                      node = adjncy(i)
cdbg                      k3   = cmpmap(node)
cdbg                      link = - node
cdbg                      if  ( node .eq. enode )  go to 1400
cdbgc.obs                 if  ( node )  1000, 2100, 1100
cdbg                      if  ( node .lt. 0 )  go to 1000 
cdbg                      if  ( node .eq. 0 )  go to 2100 
cdbg      jtb version
                       node = adjncy(i)
                       if ( node .eq. enode ) go to 1400
                       if  ( node .lt. 0 ) then
                         link = - node
                         go to 1000
                       elseif( node .eq. 0 ) then
                         link = 0
                         go to 2100
                       else
                         k3    = cmpmap(node)
                       endif
cdbg      jtb version
 
 1100                     continue
c.d                   if ( msglvl .ge. 2 ) then
c.d                       write(*,*) ' node ', node, ', qsize ',
c.d  1                    qsize(node), ', dbakw ', dbakw(node),
c.d  2                    ', dforw ', dforw(node)
c.d                   endif
c.debug
c             k2 = rowlst(xrowls(k3))
c             if ( k2 .ne. node ) then
c                 write(6,'("oops no. 10 - node, k3, k2 = ", 3i8)')
c    1                                      node, k3, k2
c             end if
c.debug
cdbg      print *,' xislo6 before if go to 1400 k3=',k3
cdbg      call xislp3 ( 'qsize ', 5, qsize , 6 )
cdbg                          if  ( qsize(k3) .eq. 0 )  go to 1400
cdbg      print *,' xislo6 before if goto 1200 node,tag=',node,tag
cdbg      call xislp3 ( 'marker', neqns, marker, 6 )
                          if  ( marker(node) .ge. tag)  go to 1200
 
c                             ------------------------------------
c                             ... a new node not previously marked
c                                 adjacent to  enode's  clique
c                             ------------------------------------
 
c.d                           if ( msglvl .ge. 2 ) then
c.d                              write(*,*) ' adding ', qsize(k3),
c.d  1                                      ' to degree'
c.d                           endif
cdbg      print *,' xislo6 b4 goto 1400 k3,node,tag,deg=',k3,node,tag,deg
cdbg      call xislp3 ( 'qsize ', 5, qsize , 6 )
cdbg      call xislp3 ( 'marker', neqns, marker, 6 )
                              marker(node) = tag
                              deg = deg + qsize(k3)
 
                              go to 1400
 
 1200                 continue
cdbg      print *,' xislo6 before go to 1400 node=',node
cdbg      call xislp3 ( 'dbakw ', neqns, dbakw  , 6 )
cdbg      call xislp3 ( 'dforw ', neqns, dforw  , 6 )
                      if  ( dbakw(node) .ne. 0 )  go to 1400
                      if  ( dforw(node) .ne. 2 )  go to 1300
 
c                         -------------------------------------------
c                         ... new node has only two neighbors, its
c                             own superelement and  enode's
c                             superelement.  since these are the
c                             same as  enode's  neighbors, the two
c                             nodes are indistinguishable. merge them
c                             into a new supernode.
c                         -------------------------------------------
c
c.d                           if ( msglvl .ge. 2 ) then
c.d                              write(*,*) ' node ', node,
c.d  1                               ' is indistinguishable to ', enode
c.d                           endif

                              k1 = cmpmap(enode)
                              k3 = cmpmap(node)

c.debug
c             k2 = rowlst(xrowls(k1))
c             if ( k2 .ne. enode ) then
c                 write(6,'("oops no. 11 - enode, k1, k2 = ", 3i8)')
c    1                                      enode, k1, k2
c             end if
c.debug
c.debug
c             k2 = rowlst(xrowls(k3))
c             if ( k2 .ne. node ) then
c                 write(6,'("oops no. 12 - node, k3, k2 = ", 3i8)')
c    1                                      node, k3, k2
c             end if
c.debug
cdbg      print *,' xislo6 b4 sets goto 1400 k1,k3,node,enode,-maxint=',
cdbg     1           k1,k3,node,enode,-maxint
cdbg      call xislp3 ( 'qsize ', 5, qsize , 6 )
cdbg      call xislp3 ( 'dbakw ', neqns, dbakw  , 6 )
cdbg      call xislp3 ( 'marker', neqns, marker, 6 )
                              qsize(k1)    = qsize(k1) + qsize(k3)
                              qsize(k3)    = 0
                              marker(node) = maxint
                              dforw(node)  = -enode
                              dbakw(node)  = -maxint
cdbg      print *,' xislo6 after sets go to 1400 k1, k3, node=',
cdbg     1           k1,k3,node
cdbg      call xislp3 ( 'qsize ', 5, qsize , 6 )
cdbg      call xislp3 ( 'dbakw ', neqns, dbakw  , 6 )
cdbg      call xislp3 ( 'marker', neqns, marker, 6 )
                              go to 1400
 
c                         ------------------------------------------
c                         ... new node's neighbors are a superset of
c                             enode's  neighbors.
c                             it is outmatched by enode.
c                         ------------------------------------------
 
 1300                     continue
                          if  ( dbakw(node).eq.0 )  then
c.d                           if ( msglvl .ge. 2 ) then
c.d                              write(*,*) ' node ', enode,
c.d  1                               ' is outmatched by ', node
c.d                           endif
cdbg      print *,' xislo6 after 1300 node, maxint, -maxint=',
cdbg     1                            node,maxint,-maxint
cdbg      call xislp3 ( 'dbakw ', neqns, dbakw  , 6 )
cdbg                              dbakw(node) = -maxint
                          endif
 
 1400             continue
                  go to 2100
 
c         -----------------------------------------------------
c         ... process each  'enode'  in the qx list,
c             the nodes with more than one non-trivial neighbor
c         -----------------------------------------------------
 
 1500     continue
          enode = qxhead
          iq2   = .false.
 
 1600     continue
          if  ( enode .le. 0 )  go to 2300
c.d       if ( msglvl .ge. 2 ) then
c.d          write(*,*) ' x-list, enode ', enode,
c.d  1                  ', dforw ', dforw(enode),
c.d  2                  ', dbakw ', dbakw(enode),
c.d  3                  ', llist ', llist(enode)
c.d       endif
cdbg      print *,' xislo6 after 1600 enode=',enode
cdbg      call xislp3 ( 'dbakw ', neqns, dbakw  , 6 )
          if  ( dbakw(enode) .ne. 0 )  go to 2200
              tag    = tag + 1
              deg    = deg0
              fstnbr = xadj(enode)
              fstclq = xcliqu(enode)
c
c             -------------------------------------------------
c             here is a 1993/cca modification.
c             without graph compression the leading entries
c             in an adjacency set are singletons, so their
c             size is just their number. with graph compression
c             the weight of each must be added to the degree
c             -------------------------------------------------
c
              do i = fstnbr, fstclq-1
                 k1 = cmpmap(adjncy(i))
c.d              if ( msglvl .ge. 2 ) then
c.d                 write(*,*) ' nbr ', adjncy(i), ' adds ',
c.d  1                         qsize(k1), ' to degree'
c.d              endif
c.d              if ( adjncy(i) .le. 0 ) then
c.d                 write(*,*) ' hey, fatal error'
c.d                 stop
c.d              endif
c.debug
c             k2 = rowlst(xrowls(k1))
c             if ( k2 .ne. adjncy(i) ) then
c                 write(6,'("oops no. 13 - adjncy(i), k1, k2 = ", 
c    1                3i8)')                adjncy(i), k1, k2
c             end if
c.debug
cdbg      print *,' xislo6 before 1600 k1,deg=',k1,deg
cdbg      call xislp3 ( 'qsize ', 5, qsize , 6 )
                 deg = deg + qsize(k1)
              enddo
 
c             ------------------------------------------------------
c             ... search through list of neighbors of  enode, adding
c                 unmarked neighbors into degree count for  enode
c             ------------------------------------------------------
 
              istrt = xcliqu(enode)
              istop = xadj(enode+1) - 1
              do 2000 i = istrt, istop
                  nabor = adjncy(i)
                  if  ( nabor .eq. elmnt )  go to 2100
cdbg      print *,' xislo6 after if goto 2100 nabor,tag=',nabor,tag
cdbg      call xislp3 ( 'marker', neqns, marker, 6 )
                      marker(nabor) = tag
                      link          = nabor
 
 1700                 continue
                      jstrt = xadj(link)
                      jstop = xadj(link+1) - 1
                      do j = jstrt, jstop
                          node = adjncy(j)
cdbg    original code
cdbg                          node = adjncy(j)
cdbg                          k1   = cmpmap(node)
cdbg                          link = - node
cdbgc.obs                     if  ( node )  1700, 2000, 1800
cdbg                          if  ( node .lt. 0 )  go to 1700
cdbg                          if  ( node .eq. 0 )  go to 2000
cdbg    new version, jtb 
                          node = adjncy(j)
                          link = - node
c.obs                     if  ( node )  1700, 2000, 1800
                          if  ( node .lt. 0 ) then
                            go to 1700
                          elseif  ( node .eq. 0 ) then
                            go to 2000
                          else
                            k1   = cmpmap(node)
                          endif
cdbg    new version, jtb 

 1800                         continue
c.d                       if ( msglvl .ge. 2 ) then
c.d                          write(*,*) ' node ', node, ' marker ',
c.d  1                                  marker(node)
c.d                       endif
cdbg      print *,' xislo6 after 1800 node,tag=',node,tag
cdbg      call xislp3 ( 'marker', neqns, marker, 6 )
                              if  ( marker(node) .ge. tag ) cycle
                                  marker(node) = tag
c.d                               if ( msglvl .ge. 2 ) then
c.d                                  write(*,*) ' adding ',
c.d  1                                     qsize(k1), ' to degree'
c.d                               endif
c.debug
c             k2 = rowlst(xrowls(k1))
c             if ( k2 .ne. node ) then
c                 write(6,'("oops no. 14 - node, k1, k2 = ", 3i8)')
c    1                                      node, k1, k2
c             end if
c.debug
cdbg      print *,' xislo6 before 1600 k1,deg=',k1,deg
cdbg      call xislp3 ( 'qsize ', 5, qsize , 6 )
                                  deg = deg + qsize(k1)
 
                      enddo
 
 2000         continue
 
c         ----------------------------------------------------
c         ... common processing for all nodes in reachable set
c             update external degree of 'enode' in degree
c             structure, and 'mdeg' (min deg) if necessary.
c         ----------------------------------------------------
 
 2100     continue
          k1   = cmpmap(enode)
cdbg      print *,' xislo6 after 2100 enode,k1,deg=',enode,k1,deg
cdbg      call xislp3 ( 'qsize ', 5, qsize , 6 )
          deg  = deg - qsize(k1)
c.debug
c             k2 = rowlst(xrowls(k1))
c             if ( k2 .ne. enode ) then
c                 write(6,'("oops no. 15 - enode, k1, k2 = ", 3i8)')
c    1                                      enode, k1, k2
c             end if
c.debug
c.d       if ( deg .le. 0 ) then
c.d          write(*,*) ' fatal error, element ', elmnt,
c.d  1                  ' enode ', enode, ', deg = ', deg
c.d          stop
c.d       endif
          fnode        = dhead(deg)
          dforw(enode) = fnode
cdbg      print *,' xislo6 after 2100 enode, -deg=',enode, -deg
cdbg      call xislp3 ( 'dbakw ', neqns, dbakw  , 6 )
          dbakw(enode) = -deg
          if  ( fnode .gt. 0 )  dbakw(fnode) = enode
cdbg      print *,' xislo6 after if and set enode, fnode=',enode,fnode
cdbg      call xislp3 ( 'dbakw ', neqns, dbakw  , 6 )
          dhead(deg)   = enode
          if  ( deg .lt. mdeg )  mdeg = deg
 
c         ---------------------------------------
c         ... return to either  q2  or  qx  list,
c             get next 'enode' in current element
c         ---------------------------------------
 
 2200     continue
          enode = llist(enode)
          if  ( iq2 )  then
              go to 900
          else
              go to 1600
          endif
 
c     --------------------------------------------------------------
c     ... finished processing all nodes in the reachable set for one
c         eliminated node.  now process the next element in the list
c     --------------------------------------------------------------
 
 2300 continue
      tag   = mtag
      elmnt = llist(elmnt)
      go to 100
 
 2400 continue
      return
 
      end
