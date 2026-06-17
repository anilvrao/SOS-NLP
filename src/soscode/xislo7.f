      subroutine xislo7 ( nsuper, nnodes, parsup, xsup  , lnzcol,
     .                    fstsup, brosup, parent, nodlst, snsize,
     .                    uleng , fleng , mxsons, nsons   )
 
c
c     ==================================================================
c     ==================================================================
c     ====  xislo7 -- initialize for partition routines             ====
c     ==================================================================
c     ==================================================================
c
c     originally this routine was betree, created by joseph w.h. liu
c     modification   apr. 24, 1987
c     modification for noetic  may 16, 1988 (bwp)
c     modification for the out-of-core multifrontal code occured
c                        jun 10, 1988
c
c     purpose - determine a binary tree representation of a
c               supernodal elimination tree represented by the
c               parsup vector.
c               the returned representation will be given by
c               the  'first-son'  and  'brother'  vectors.
c               the root of the binary tree is always 'nsuper'.
c               code has been modified for noetic to return
c               maximum number of sons occuring in the tree.
c               code has been further modified for the out-of-core
c               multifrontal solver to return the parent vector
c               of the nodal elimination tree, node lists for
c               each supernode---also, supernode sizes, update matrix
c               dimensions, and frontal matrix dimensions.
c
c     input parameters -
c         nsuper - number of supernodes
c         nnodes = number of nodes
c         parsup - the supernodal elimination tree parent vector
c                  it is assumed that parsup(i) > i except at the roots.
c         xsup   - the supernode partition
c         lnzcol - the number of nonzero entries in each column of l
c                  excluding the diagonal entry
c
c     output parameters -
c         fstsup - the first son vector for the supernodal elimination
c                  tree
c         brosup - the next brother vector for the supenodal elimination
c                  tree
c         parent - the parent vector for the nodal elimination tree
c         nodlst - the nodes in each supernode
c         snsize - supernode sizes
c         uleng  - dimensions of the update matrices
c         fleng  - dimensions of the frontal matrices
c         mxsons - maximum number of sons
c
c     work parameters -
c         nsons  - number of sons
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      integer         nsuper, nnodes, mxsons
 
      integer         parsup(*), xsup(*), lnzcol(*),
     .                fstsup(*), brosup(*), parent(*),
     .                nodlst(*), snsize(*), uleng (*),
     .                fleng (*), nsons (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer         fstnod, inode , isuper, lroot , lstnod,
     .                par
 
c     ==================================================================
 
c     --------------
c     initialization
c     --------------
 
      mxsons = 0
 
      if  ( nsuper .le. 0 )  return
 
      do isuper = 1, nsuper
          fstsup(isuper) = 0
          brosup(isuper) = 0
          nsons (isuper) = 0
      enddo
 
      do inode = 1, nnodes
          parent(inode) = 0
      enddo
 
      lroot = nsuper
 
c     --------------------------------------------
c     for each isuper := nsuper-1 step -1 downto 1
c     --------------------------------------------
 
      if  ( nsuper .gt. 1 )  then
 
          do isuper = nsuper-1, 1, -1
 
              par = parsup(isuper)
 
              if  ( par .gt. 0 .and. par .ne. isuper )  then
 
c                 ----------------------------------------
c                 ... isuper becomes first son of its parent
c                 ----------------------------------------
 
                  brosup(isuper) = fstsup(par)
                  fstsup(par)    = isuper
                  nsons (par)    = nsons(par) + 1
                  mxsons         = max ( mxsons, nsons(par) )
 
              else
 
c                 -----------------------------------------------------
c                 ... isuper has no parent. given structure is a forest
c                     set isuper to be one of the roots of the trees
c                 -----------------------------------------------------
 
                  brosup(lroot) = isuper
                  lroot         = isuper
 
              endif
 
          enddo
 
      endif
 
      brosup(lroot) = 0
 
c     ------------------------------------
c     generate the rest of the information
c     ------------------------------------
 
      do isuper = 1, nsuper
 
          fstnod = xsup(isuper)
          lstnod = xsup(isuper+1) - 1
          snsize(isuper) = lstnod - fstnod + 1
          fleng (isuper) = lnzcol(fstnod) + 1
          uleng (isuper) = fleng(isuper) - snsize(isuper)
          par = parsup(isuper)
 
          do inode = fstnod, lstnod-1
              nodlst(inode) = inode
              parent(inode) = inode+1
          enddo
 
          if ( par .gt. 0 )  parent(lstnod) = xsup(par)
          nodlst(lstnod) = lstnod
 
      enddo
 
      return
 
      end
