
      subroutine dualin(veclam,istatc,maxcon,mcon,vecnu,istatv,ndim)
c
      implicit double precision (a-h,o-z)
c
c         purpose: given estimates for the dual variables (lagrange 
c                  multipliers veclam, and vecnu), set the constraint
c                  status flags to insure active set estimate
c                  is reasonable.   Equality and ignored constraints
c                  are unchanged.   Inequality constraints are selected
c                  based on the absolute value of the multiplier---
c                  largest first.   Inequalities are not added to 
c                  the active set if the active set is full.
c
c         arguments:
c
c           veclam   Constraint multipliers (mcon)
c           istatc   constraint status (mcon)
c           maxcon   Max dimension of veclam, istatc.
c           vecnu    variable multipliers (ndim)
c           istatv   variable status (ndim)
c
c
      double precision, allocatable :: vmlt(:)
      integer, allocatable :: imlt(:)
      integer, allocatable :: iprm(:)
c
      dimension veclam(maxcon),istatc(maxcon),vecnu(ndim),istatv(ndim)
c
      parameter (zero=0.d0)
c
      lnmlt = ndim + maxcon
      allocate(vmlt(1:lnmlt))
      allocate(imlt(1:lnmlt))
      allocate(iprm(1:lnmlt))
c
      mactiv = 0
      ineq = 0
      varloop: do ii = 1,ndim
        if(istatv(ii).eq.3) then
          mactiv = mactiv + 1
          cycle varloop
        endif
        ineq = ineq + 1
        imlt(ineq) = ii
        vmlt(ineq) = abs(vecnu(ii))
      enddo varloop
c
      conloop: do ii = 1,mcon
        if(istatc(ii).ge.3) then
          mactiv = mactiv + 1
          cycle conloop
        endif
        ineq = ineq + 1
        imlt(ineq) = -ii
        vmlt(ineq) = abs(veclam(ii))
      enddo conloop
c
c         sort absolute value of inequality multipliers
c
      call hdsrtn(vmlt,ineq,0,0,iprm,ierp)
      if(ierp.ne.0) then
        RETURN
      endif
c
c         reorder row indices to correspond
c
      call hjprmx(imlt,ineq,iprm,ierp)
      if(ierp.ne.0) then
        RETURN
      endif
c
c         loop through the inequalities to reset multipliers
c         and status arrays
c
      ineqloop: do ii = 1,ineq
c
        jj = imlt(ii)
        if(jj.gt.0) then
c
c         variable multiplier
c
          if(vmlt(ii).gt.zero) then
            mactiv = mactiv + 1
            if(mactiv.gt.ndim) then
              istatv(jj) = 0
              vecnu(jj) = zero
            elseif(vecnu(jj).gt.zero) then
              istatv(jj) = 1
            else
              istatv(jj) = 2
            endif
          else
            istatv(jj) = 0
          endif
c
        else
c
          jj = -jj
c
c         constraint multiplier
c
          if(vmlt(ii).gt.zero) then
            mactiv = mactiv + 1
            if(mactiv.gt.ndim) then
              istatc(jj) = 0
              veclam(jj) = zero
            elseif(veclam(jj).gt.zero) then
              istatc(jj) = 1
            else
              istatc(jj) = 2
            endif
          else
            istatc(jj) = 0
          endif
c
        endif
c
      enddo ineqloop
c
      deallocate(vmlt,imlt,iprm)
c
      return
      end
