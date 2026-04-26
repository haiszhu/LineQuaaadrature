module lq_adaptive_mod
  use linequaaadrature_mod, only: r64
  implicit none
  private

  public :: line_quad_root_initial_guess_r64
  public :: line_quad_root_refine_r64

contains

  ! Copy/adapt the body of lqa_root_initial_guess_r64 here,
  ! but rename it:
  subroutine line_quad_root_initial_guess_r64(tj, xj, yj, zj, n, x0, y0, z0, tinit)
    ! Find complex initial guess for the closest-point root-finder.
    ! Returns tinit = (t_proj + i*perp_dist) in parameter space near the
    ! closest point on the curve to (x0,y0,z0).
    ! Replaces legacy rootfinderinitialguess() from linequad.f.
    integer(8), intent(in) :: n
    real(r64), intent(in) :: tj(n), xj(n), yj(n), zj(n), x0, y0, z0
    complex(8), intent(out) :: tinit
    real(r64) :: rsqvec(n), rsqmin1, rsqmin2, p(3), r(3)
    real(r64) :: pnorm2, pnorm2i, pnormi, rnormsq, rdotp, a, b
    integer(8) :: i, imin1, imin2
    ! Find two nearest GL nodes
    rsqvec  = (xj-x0)**2 + (yj-y0)**2 + (zj-z0)**2
    rsqmin1 = 1.0e16_r64;  rsqmin2 = 1.0e16_r64
    imin1 = 1;  imin2 = 1
    do i = 1, n
      if (rsqvec(i) < rsqmin1) then
        rsqmin2 = rsqmin1;  imin2 = imin1
        rsqmin1 = rsqvec(i);  imin1 = i
      else if (rsqvec(i) < rsqmin2) then
        rsqmin2 = rsqvec(i);  imin2 = i
      end if
    end do
    ! Complex initial guess via projection of r0 onto the chord i1→i2
    p(1) = xj(imin1) - xj(imin2)
    p(2) = yj(imin1) - yj(imin2)
    p(3) = zj(imin1) - zj(imin2)
    pnorm2  = p(1)**2 + p(2)**2 + p(3)**2
    pnorm2i = 1.0_r64 / pnorm2
    pnormi  = sqrt(pnorm2i)
    r(1) = x0 - xj(imin1)
    r(2) = y0 - yj(imin1)
    r(3) = z0 - zj(imin1)
    rnormsq = r(1)**2 + r(2)**2 + r(3)**2
    rdotp   = r(1)*p(1) + r(2)*p(2) + r(3)*p(3)
    a = (tj(imin1) - tj(imin2)) * rdotp * pnorm2i
    b = sqrt(max(0.0_r64, rnormsq - rdotp*rdotp*pnorm2i)) * (tj(imin1) - tj(imin2)) * pnormi
    tinit = cmplx(tj(imin1) + a, b, kind=8)
  end subroutine line_quad_root_initial_guess_r64


  ! Copy/adapt the body of lqa_root_refine_r64 here,
  ! but rename it:
  subroutine line_quad_root_refine_r64(xhat, yhat, zhat, n, x0, y0, z0, tinit, troot, ifconv)
    ! Complex root refinement: find t* such that |r(t*) - r0|^2 = 0 in complex plane.
    ! Uses Newton (20 iter) then Muller (20 iter) with complex Legendre evaluations.
    ! Replaces legacy rootfinder() from linequad.f.
    integer(8), intent(in) :: n
    real(r64), intent(in) :: xhat(n), yhat(n), zhat(n), x0, y0, z0
    complex(8), intent(in) :: tinit
    complex(8), intent(out) :: troot
    integer(8), intent(out) :: ifconv
    complex(8) :: P(n), D(n)
    complex(8) :: t, tp, tpp, dt
    complex(8) :: cx, cy, cz, cxp, cyp, czp, cdx, cdy, cdz
    complex(8) :: F, Fprime, Fp, Fpp
    complex(8) :: q, A, B, C, d1, d2
    real(r64), parameter :: tol = 1.0e-15_r64
    integer(8) :: iter
    t      = tinit
    troot  = tinit
    ifconv = 0_8
    Fp     = cmplx(0.0_r64, 0.0_r64, kind=8)
    tp     = t
    tpp    = t
    Fpp    = cmplx(0.0_r64, 0.0_r64, kind=8)
    ! Newton iterations
    do iter = 1, 20
      call line_quad_legendrederiv_cx8(n-1_8, t, P, D)
      cx  = sum(xhat*P);  cy  = sum(yhat*P);  cz  = sum(zhat*P)
      cxp = sum(xhat*D);  cyp = sum(yhat*D);  czp = sum(zhat*D)
      cdx = cx - x0;  cdy = cy - y0;  cdz = cz - z0
      F      = cdx*cdx + cdy*cdy + cdz*cdz
      Fprime = 2.0_r64*(cdx*cxp + cdy*cyp + cdz*czp)
      dt  = -F / Fprime
      tpp = tp;  Fpp = Fp;  Fp = F;  tp = t
      t   = t + dt
      if (abs(dt) < tol) then;  ifconv = 1_8;  troot = t;  return;  end if
    end do
    ! Muller iterations: seed F from the last Newton P (already in P) to save one eval
    ifconv = 0_8
    cx = sum(xhat*P);  cy = sum(yhat*P);  cz = sum(zhat*P)
    cdx = cx - x0;  cdy = cy - y0;  cdz = cz - z0
    F = cdx*cdx + cdy*cdy + cdz*cdz
    ! tp, Fp, tpp, Fpp already hold the last two Newton history points
    do iter = 1, 20
      q  = (t-tp) / (tp-tpp)
      A  = q*F - q*(q+1.0_r64)*Fp + q*q*Fpp
      B  = (2.0_r64*q+1.0_r64)*F - (1.0_r64+q)**2*Fp + q*q*Fpp
      C  = (1.0_r64+q)*F
      d1 = B + sqrt(B*B - 4.0_r64*A*C)
      d2 = B - sqrt(B*B - 4.0_r64*A*C)
      if (abs(d1) > abs(d2)) then
        dt = -(t-tp)*2.0_r64*C/d1
      else
        dt = -(t-tp)*2.0_r64*C/d2
      end if
      tpp = tp;  Fpp = Fp;  Fp = F;  tp = t
      t   = t + dt
      if (abs(dt) < tol) then;  ifconv = 1_8;  troot = t;  return;  end if
      call line_quad_legendre_cx8(n-1_8, t, P)
      cx = sum(xhat*P);  cy = sum(yhat*P);  cz = sum(zhat*P)
      cdx = cx - x0;  cdy = cy - y0;  cdz = cz - z0
      F = cdx*cdx + cdy*cdy + cdz*cdz
    end do
  end subroutine line_quad_root_refine_r64

  subroutine line_quad_legendre_cx8(n, x, P)
    ! Complex Legendre polynomial values P(0)..P(n) at complex x. (n+1 values)
    integer(8), intent(in) :: n
    complex(8), intent(in) :: x
    complex(8), intent(out) :: P(n+1)
    integer(8) :: l
    P(1) = 1.0_r64
    if (n < 1) return
    P(2) = x
    do l = 1, n-1
      P(l+2) = ((2*l+1)*x*P(l+1) - l*P(l)) / real(l+1, r64)
    end do
  end subroutine line_quad_legendre_cx8

  subroutine line_quad_legendrederiv_cx8(n, x, P, D)
    ! Complex Legendre values P(0)..P(n) and derivatives D(0)..D(n) at complex x.
    integer(8), intent(in) :: n
    complex(8), intent(in) :: x
    complex(8), intent(out) :: P(n+1), D(n+1)
    integer(8) :: l
    P(1) = 1.0_r64;  D(1) = 0.0_r64
    if (n < 1) return
    P(2) = x;  D(2) = 1.0_r64
    do l = 1, n-1
      P(l+2) = ((2*l+1)*x*P(l+1) - l*P(l)) / real(l+1, r64)
      D(l+2) = ((2*l+1)*(P(l+1) + x*D(l+1)) - l*D(l)) / real(l+1, r64)
    end do
  end subroutine line_quad_legendrederiv_cx8
end module lq_adaptive_mod