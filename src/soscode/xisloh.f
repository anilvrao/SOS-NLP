      subroutine xisloh ( ncomp , xrowls, rowlst, cmpmap, neqns,
     1                    xadj,   adjncy, dhead,  dforw,  dbakw,
     2                    maxint, qsize,  llist,  marker, xcliqu, 
     3                    clqsiz )
 
c
c     ==================================================================
c     ==================================================================
c     ====  xisloh -- multiple mininum degree initialization        ====
c     ==================================================================
c     ==================================================================
c
c
c     purpose : this routine performs initialization for the multiple
c               elimination version of the minimum degree algorithm.
c
c     created : march 93, cca
c     modified: dec.  96, rgg, to use input compression 
c
c     input parameters:
c         ncomp  - number of compressed nodes
c         (xrowls, rowlst) - list of rows for each compressed node
c         cmpmap - map from row numbers to compressed nodes
c         neqns - number of equations.
c         (xadj, adjncy) - adjacency structure.
c
c     output parameters:
c         (dhead, dforw, dbakw) -- degree doubly linked structure.
c         maxint -- maximum integer
c         qsize  -- size of supernode (initialized to one).
c         llist  -- linked list.
c         marker -- marker vector.
c         xcliqu -- pointers to cliques in adjacency lists
c         clqsiz -- elimination clique sizes
c
c     ==================================================================
 
      integer             ncomp , neqns , maxint 
 
      integer             xrowls (*), rowlst (*), cmpmap (*), 
     1                    xadj   (*), adjncy (*),
     1                    dhead  (*), dforw  (*), dbakw  (*),
     2                    qsize  (*), llist  (*),
     3                    marker (*), xcliqu (*), clqsiz (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             fnode , i     , icomp , istart, 
     1                    istop,  jrep  , k1    , ndeg  , node
 
c     ==================================================================
 
c   ----------------------------
c   initialize some array fields
c   ----------------------------
c
c.debug
c     write(6,'("in xisloh")')
c.debug
      do node = 1, neqns

         dhead(node)  = 0
         marker(node) = 0
         llist(node)  = 0
         xcliqu(node) = 0
         clqsiz(node) = 0

         icomp  = cmpmap(node)
         jrep   = rowlst(xrowls(icomp))
c.debug
c        write(6,'("node, icomp, jrep              = ", 3i8)')
c    1               node, icomp, jrep 
c        write(6,'("xrowls(icomp), xrowls(icomp+1) = ", 2i8)')
c    1               xrowls(icomp), xrowls(icomp+1)
c.debug

         if ( jrep .eq. node ) then

c            -------------------------------------
c            ... this is the representative column
c            -------------------------------------

             qsize (icomp) = xrowls(icomp+1) - xrowls(icomp)
c.debug
c     write(6,'("setting qsize - icomp, qsize(icomp) =", 2i8)')
c    1                            icomp, qsize(icomp)
c     call xislp3 ( 'qsize ', icomp, qsize , 6 )
c.debug
             dforw (node)  = 0
             dbakw (node)  = 0

         else

c            ---------------------------------------
c            ... this is an indistinguishable column
c            ---------------------------------------

             dforw (node) = -jrep
             dbakw (node) = -maxint
             marker(node) = maxint

         end if

      enddo
c.debug
c     write(6,'("in xisloh after 100")')
c     call xislp3 ( 'qsize ', ncomp, qsize , 6 )
c     call xislp3 ( 'dforw ', neqns, dforw , 6 )
c     call xislp3 ( 'dbakw ', neqns, dbakw , 6 )
c.debug
 
c     ------------------------------------------
c     initialize the degree doubly linked lists.
c     ------------------------------------------
 
      do node = 1, neqns

          icomp  = cmpmap(node)
          jrep   = rowlst(xrowls(icomp))

          if ( jrep .ne. node ) cycle

          if ( qsize(icomp) .eq. 0 ) cycle
c
c         -----------------------------------------------
c         node is representative, compute external degree
c         plus one and insert into the degree list
c         -----------------------------------------------
c
          istart = xadj(node)
          istop  = xadj(node+1) - 1
c         ndeg   = qsize(node)
          ndeg   = 1 
          do i = istart, istop
             k1   = cmpmap(adjncy(i))
             ndeg = ndeg + qsize(k1)
          enddo
c
c         -----------------------
c         add node to degree list
c         -----------------------

          fnode       = dhead(ndeg)
          dforw(node) = fnode
          dhead(ndeg) = node
          if ( fnode .gt. 0 ) then
             dbakw(fnode) = node
          endif
          dbakw(node) = - ndeg
c.debug
c     write(6,'("before 200 - fnode, node, ndeg        = ", 3i8)')
c    1                         fnode, node, ndeg 
c     write(6,'("dforw(node), dhead(ndeg), dbakw(node) = ", 3i8)')
c    1            dforw(node), dhead(ndeg), dbakw(node)
c.debug

      enddo
c.debug
c     write(6,'("in xisloh after 200")')
c     call xislp3 ( 'dforw ', neqns, dforw , 6 )
c     call xislp3 ( 'dbakw ', neqns, dbakw , 6 )
c.debug
 
      return
      end
