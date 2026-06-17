      subroutine  xislo5  ( neqns , nsuper, nsupp1, xsup  , invp  ,
     1                      sparnt, supmap, temp   )
 
c
c     ==================================================================
c     ==================================================================
c     =====                                                        =====
c     ===== xislo5 - generates supernodal parent vector            =====
c     =====                                                        =====
c     ==================================================================
c     ==================================================================
c
c     creation date: 09-28-87 (bwp)
c     last updated: 09-28-87 (bwp)
c
c     purpose -- generates the parent vector for the supernodal
c                elimination tree from data gathered while computing
c                the minimum degree ordering.  this routine is for use
c                by the minimum degree driver after the call to the
c                numbering routine, mmdnum has been completed.
c
c     input parameters:
c
c        neqns  - number of equations
c        nsuper - number of supernodes
c        nsupp1 - nsuper+1, dimension of xsup
c        xsup   - the supernode partition
c        invp   - the inverse permutation vector, maps old labels to
c                 new labels
c        sparnt - maps a single node (using old labels) from each
c                 supernode to a single node in its parent supernode.
c                 overwritten on output.
c
c     output parameters:
c
c        sparnt - on output, sparnt contains the parent vector for the
c                 supernodal elimination tree
c
c     work parameters:
c
c        supmap - maps each node to the supernode containing it
c        temp   - temporary n-vector in which the parent vector for the
c                 supernodal elimination tree is computed.
c
c     external subprograms:
c
c        none
c
c     ==================================================================
 
c     ----------
c     parameters
c     ----------
 
      integer             neqns , nsuper, nsupp1
 
      integer             xsup  (*), invp  (* ),
     1                    sparnt(* ), supmap(* ),
     2                    temp  (*)
 
c     ---------------
c     local variables
c     ---------------
 
      integer             child , i     , j     , parent, schild
 
c     ==================================================================
 
c     ------------------
c     initialize temp(*)
c     ------------------
 
      do i = 1, nsuper
          temp(i) = 0
      enddo
 
c     ---------------------------
c     generate node-supernode map
c     ---------------------------
 
      do i = 1, nsuper
          do j = xsup(i), xsup(i+1)-1
              supmap(j) = i
          enddo
      enddo
 
c     ----------------------------------------------------
c     transform the
c           old node--to--old node
c     child-parent edges currently found in sparnt(*) into
c           supernode--to--supernode
c     child-parent edges
c     ----------------------------------------------------
 
      do i = 1, neqns
 
          if  ( sparnt(i) .ne. 0 )  then
 
              child        = invp(i)
              parent       = invp(sparnt(i))
 
              schild       = supmap(child)
              temp(schild) = supmap(parent)
 
          endif
 
      enddo
 
c     ----------------------------------------------
c     copy parent vector from temp(*) into sparnt(*)
c     ----------------------------------------------
 
      do i = 1, nsuper
          sparnt(i) = temp(i)
      enddo
 
      return
      end
