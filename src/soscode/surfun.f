
      double precision function surfun(xcoef,ndim,mcon)
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension xcoef(10)
c
      fun = xcoef(1) 
      fun = fun + xcoef(2)*dble(ndim)
      fun = fun + xcoef(3)*dble(ndim**2)
      fun = fun + xcoef(4)*dble(ndim**3)
      fun = fun + xcoef(5)*dble(mcon)
      fun = fun + xcoef(6)*dble(mcon**2)
      fun = fun + xcoef(7)*dble(mcon**3)
      fun = fun + xcoef(8)*dble(ndim)*dble(mcon)
      fun = fun + xcoef(9)*dble(ndim**2)*dble(mcon)
      fun = fun + xcoef(10)*dble(ndim)*dble(mcon**2)
c
      surfun = fun
c
      return
      end
