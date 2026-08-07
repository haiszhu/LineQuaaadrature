module lq_adaptive_mod
  use linequaaadrature_mod, only: r64, r128
  implicit none
  private

  public :: line_quad_root_initial_guess_r64
  public :: line_quad_root_refine_r64
#ifndef BIESOLVER_R64_ONLY
  public :: line_quad_root_initial_guess_r128
  public :: line_quad_root_refine_r128
#endif

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
    complex(8), intent(inout) :: tinit
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


#ifndef BIESOLVER_R64_ONLY
  subroutine line_quad_root_initial_guess_r128(tj, xj, yj, zj, n, x0, y0, z0, tinit)
    integer(8),  intent(in)  :: n
    real(r128),  intent(in)  :: tj(n), xj(n), yj(n), zj(n), x0, y0, z0
    complex(16), intent(out) :: tinit
    real(r128) :: rsqvec(n), rsqmin1, rsqmin2, p(3), r(3)
    real(r128) :: pnorm2, pnorm2i, pnormi, rnormsq, rdotp, a, b
    integer(8) :: i, imin1, imin2
    rsqvec  = (xj-x0)**2 + (yj-y0)**2 + (zj-z0)**2
    rsqmin1 = 1.0e32_r128;  rsqmin2 = 1.0e32_r128
    imin1 = 1;  imin2 = 1
    do i = 1, n
      if (rsqvec(i) < rsqmin1) then
        rsqmin2 = rsqmin1;  imin2 = imin1
        rsqmin1 = rsqvec(i);  imin1 = i
      else if (rsqvec(i) < rsqmin2) then
        rsqmin2 = rsqvec(i);  imin2 = i
      end if
    end do
    p(1) = xj(imin1) - xj(imin2)
    p(2) = yj(imin1) - yj(imin2)
    p(3) = zj(imin1) - zj(imin2)
    pnorm2  = p(1)**2 + p(2)**2 + p(3)**2
    pnorm2i = 1.0_r128 / pnorm2
    pnormi  = sqrt(pnorm2i)
    r(1) = x0 - xj(imin1)
    r(2) = y0 - yj(imin1)
    r(3) = z0 - zj(imin1)
    rnormsq = r(1)**2 + r(2)**2 + r(3)**2
    rdotp   = r(1)*p(1) + r(2)*p(2) + r(3)*p(3)
    a = (tj(imin1) - tj(imin2)) * rdotp * pnorm2i
    b = sqrt(max(0.0_r128, rnormsq - rdotp*rdotp*pnorm2i)) * (tj(imin1) - tj(imin2)) * pnormi
    tinit = cmplx(tj(imin1) + a, b, kind=16)
  end subroutine line_quad_root_initial_guess_r128
#endif


  ! Copy/adapt the body of lqa_root_refine_r64 here,
  ! but rename it:
  subroutine line_quad_root_refine_r64(xhat, yhat, zhat, n, x0, y0, z0, tinit, troot, ifconv)
    ! Complex root refinement: find t* such that |r(t*) - r0|^2 = 0 in complex plane.
    ! Uses Newton (20 iter) then Muller (20 iter) with complex Legendre evaluations.
    !
    ! Implementation note (perf): the Newton/Muller hot loops use the FUSED
    ! scalar 3-term Legendre recurrence pattern from the legacy rootfinder()
    ! in rrq-legacy/src/linequad.f. The Legendre values P_k(t) and
    ! derivatives D_k(t) are streamed into the six scalar accumulators
    ! (cx, cy, cz, cxp, cyp, czp) on the fly via the pkm2/pkm1/dkm2/dkm1
    ! scalar history, so we never materialize a length-n P/D array and
    ! never pay six dot-product passes per iteration. This matches the
    ! legacy timing on n=16 problems; an earlier "fill array, then
    ! sum(xhat*P)" form was ~30-40% slower because of the extra array
    ! stores plus six xhat*P/D temporaries created by the array intrinsic.
    integer(8), intent(in) :: n
    real(r64), intent(in) :: xhat(n), yhat(n), zhat(n), x0, y0, z0
    complex(8), intent(in) :: tinit
    complex(8), intent(out) :: troot
    integer(8), intent(out) :: ifconv
    complex(8) :: t, tp, tpp, dt
    complex(8) :: cx, cy, cz, cxp, cyp, czp, cdx, cdy, cdz
    complex(8) :: F, Fprime, Fp, Fpp
    complex(8) :: q, A, B, C, d1, d2
    complex(8) :: pkm2, pkm1, pk, dkm2, dkm1, dk
    ! logical, parameter :: use_r128_tail = .true.
    logical, parameter :: use_r128_tail = .false.
    real(r64), parameter :: tol = 1.0e-15_r64
    integer(8) :: iter, ell, i

    t      = tinit
    troot  = tinit
    ifconv = 0_8
    Fp     = cmplx(0.0_r64, 0.0_r64, kind=8)
    tp     = t
    tpp    = t
    Fpp    = cmplx(0.0_r64, 0.0_r64, kind=8)

    ! ---------------- Newton ----------------
    do iter = 1, 20
      ! Fused 3-term recurrence with on-the-fly accumulation (no P/D arrays).
      pkm2 = cmplx(1.0_r64, 0.0_r64, kind=8)
      dkm2 = cmplx(0.0_r64, 0.0_r64, kind=8)
      cx   = xhat(1)*pkm2;  cy   = yhat(1)*pkm2;  cz   = zhat(1)*pkm2
      cxp  = xhat(1)*dkm2;  cyp  = yhat(1)*dkm2;  czp  = zhat(1)*dkm2
      if (n > 1_8) then
        pkm1 = t
        dkm1 = cmplx(1.0_r64, 0.0_r64, kind=8)
        cx  = cx  + xhat(2)*pkm1;  cy  = cy  + yhat(2)*pkm1;  cz  = cz  + zhat(2)*pkm1
        cxp = cxp + xhat(2)*dkm1;  cyp = cyp + yhat(2)*dkm1;  czp = czp + zhat(2)*dkm1
        do ell = 1, n - 2
          i  = ell + 2
          pk = ((2*ell+1)*t*pkm1 - ell*pkm2) / real(ell+1, r64)
          dk = ((2*ell+1)*(pkm1 + t*dkm1) - ell*dkm2) / real(ell+1, r64)
          cx  = cx  + xhat(i)*pk;  cy  = cy  + yhat(i)*pk;  cz  = cz  + zhat(i)*pk
          cxp = cxp + xhat(i)*dk;  cyp = cyp + yhat(i)*dk;  czp = czp + zhat(i)*dk
          pkm2 = pkm1;  pkm1 = pk
          dkm2 = dkm1;  dkm1 = dk
        end do
      end if
      cdx = cx - x0;  cdy = cy - y0;  cdz = cz - z0
      F      = cdx*cdx + cdy*cdy + cdz*cdz
      Fprime = 2.0_r64*(cdx*cxp + cdy*cyp + cdz*czp)
      dt     = -F / Fprime
      tpp = tp;  Fpp = Fp;  Fp = F;  tp = t
      t   = t + dt
      if (abs(dt) < tol) then
        ifconv = 1_8
        troot  = t
#ifndef BIESOLVER_R64_ONLY
        if (use_r128_tail) then
          block
            real(r128) :: xhat128(n), yhat128(n), zhat128(n)
            real(r128) :: x0128, y0128, z0128
            complex(16) :: t128, dt128
            complex(16) :: P128(n), D128(n)
            complex(16) :: cx128, cy128, cz128, cxp128, cyp128, czp128
            complex(16) :: cdx128, cdy128, cdz128, F128, Fprime128
            real(r128), parameter :: tol_tail = 1.0e-16_r128
            integer(8) :: iter128

            xhat128 = real(xhat, r128)
            yhat128 = real(yhat, r128)
            zhat128 = real(zhat, r128)
            x0128 = real(x0, r128)
            y0128 = real(y0, r128)
            z0128 = real(z0, r128)
            t128 = cmplx(real(troot, r128), real(aimag(troot), r128), kind=16)

            do iter128 = 1, 20
              call line_quad_legendrederiv_cx16(n-1_8, t128, P128, D128)
              cx128  = sum(xhat128*P128);  cy128  = sum(yhat128*P128);  cz128  = sum(zhat128*P128)
              cxp128 = sum(xhat128*D128);  cyp128 = sum(yhat128*D128);  czp128 = sum(zhat128*D128)
              cdx128 = cx128 - x0128;  cdy128 = cy128 - y0128;  cdz128 = cz128 - z0128
              F128      = cdx128*cdx128 + cdy128*cdy128 + cdz128*cdz128
              Fprime128 = 2.0_r128*(cdx128*cxp128 + cdy128*cyp128 + cdz128*czp128)
              dt128 = -F128 / Fprime128
              t128 = t128 + dt128
              if (abs(dt128) < tol_tail) exit
            end do

            troot = cmplx(real(t128, r64), real(aimag(t128), r64), kind=8)
          end block
        end if
#endif
        return
      end if
    end do

    ! ---------------- Muller ----------------
    ! Each iteration computes F(t) at the *current* t using the fused
    ! recurrence (derivatives not needed). This matches legacy ordering;
    ! tp, tpp, Fp, Fpp already hold the last two Newton history points.
    ifconv = 0_8
    do iter = 1, 20
      pkm2 = cmplx(1.0_r64, 0.0_r64, kind=8)
      cx   = xhat(1)*pkm2;  cy   = yhat(1)*pkm2;  cz   = zhat(1)*pkm2
      if (n > 1_8) then
        pkm1 = t
        cx = cx + xhat(2)*pkm1;  cy = cy + yhat(2)*pkm1;  cz = cz + zhat(2)*pkm1
        do ell = 1, n - 2
          i  = ell + 2
          pk = ((2*ell+1)*t*pkm1 - ell*pkm2) / real(ell+1, r64)
          cx = cx + xhat(i)*pk;  cy = cy + yhat(i)*pk;  cz = cz + zhat(i)*pk
          pkm2 = pkm1;  pkm1 = pk
        end do
      end if
      cdx = cx - x0;  cdy = cy - y0;  cdz = cz - z0
      F   = cdx*cdx + cdy*cdy + cdz*cdz
      ! Muller step on the parabola through (tpp,Fpp), (tp,Fp), (t,F).
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
      if (abs(dt) < tol) then
        ifconv = 1_8
        troot  = t
#ifndef BIESOLVER_R64_ONLY
        if (use_r128_tail) then
          block
            real(r128) :: xhat128(n), yhat128(n), zhat128(n)
            real(r128) :: x0128, y0128, z0128
            complex(16) :: t128, dt128
            complex(16) :: P128(n), D128(n)
            complex(16) :: cx128, cy128, cz128, cxp128, cyp128, czp128
            complex(16) :: cdx128, cdy128, cdz128, F128, Fprime128
            real(r128), parameter :: tol_tail = 1.0e-16_r128
            integer(8) :: iter128

            xhat128 = real(xhat, r128)
            yhat128 = real(yhat, r128)
            zhat128 = real(zhat, r128)
            x0128 = real(x0, r128)
            y0128 = real(y0, r128)
            z0128 = real(z0, r128)
            t128 = cmplx(real(troot, r128), real(aimag(troot), r128), kind=16)

            do iter128 = 1, 20
              call line_quad_legendrederiv_cx16(n-1_8, t128, P128, D128)
              cx128  = sum(xhat128*P128);  cy128  = sum(yhat128*P128);  cz128  = sum(zhat128*P128)
              cxp128 = sum(xhat128*D128);  cyp128 = sum(yhat128*D128);  czp128 = sum(zhat128*D128)
              cdx128 = cx128 - x0128;  cdy128 = cy128 - y0128;  cdz128 = cz128 - z0128
              F128      = cdx128*cdx128 + cdy128*cdy128 + cdz128*cdz128
              Fprime128 = 2.0_r128*(cdx128*cxp128 + cdy128*cyp128 + cdz128*czp128)
              dt128 = -F128 / Fprime128
              t128 = t128 + dt128
              if (abs(dt128) < tol_tail) exit
            end do

            troot = cmplx(real(t128, r64), real(aimag(t128), r64), kind=8)
          end block
        end if
#endif
        return
      end if
    end do
  end subroutine line_quad_root_refine_r64

#ifndef BIESOLVER_R64_ONLY
  subroutine line_quad_root_refine_r128(xhat, yhat, zhat, n, x0, y0, z0, tinit, troot, ifconv)
    integer(8),  intent(in)  :: n
    real(r128),  intent(in)  :: xhat(n), yhat(n), zhat(n), x0, y0, z0
    complex(16), intent(in)  :: tinit
    complex(16), intent(out) :: troot
    integer(8),  intent(out) :: ifconv
    complex(16) :: P(n), D(n)
    complex(16) :: t, tp, tpp, dt
    complex(16) :: cx, cy, cz, cxp, cyp, czp, cdx, cdy, cdz
    complex(16) :: F, Fprime, Fp, Fpp
    complex(16) :: q, A, B, C, d1, d2
    real(r128), parameter :: tol = 1.0e-30_r128
    integer(8) :: iter
    t      = tinit
    troot  = tinit
    ifconv = 0_8
    Fp     = cmplx(0.0_r128, 0.0_r128, kind=16)
    tp     = t
    tpp    = t
    Fpp    = cmplx(0.0_r128, 0.0_r128, kind=16)
    do iter = 1, 20
      call line_quad_legendrederiv_cx16(n-1_8, t, P, D)
      cx  = sum(xhat*P);  cy  = sum(yhat*P);  cz  = sum(zhat*P)
      cxp = sum(xhat*D);  cyp = sum(yhat*D);  czp = sum(zhat*D)
      cdx = cx - x0;  cdy = cy - y0;  cdz = cz - z0
      F      = cdx*cdx + cdy*cdy + cdz*cdz
      Fprime = 2.0_r128*(cdx*cxp + cdy*cyp + cdz*czp)
      dt  = -F / Fprime
      tpp = tp;  Fpp = Fp;  Fp = F;  tp = t
      t   = t + dt
      if (abs(dt) < tol) then;  ifconv = 1_8;  troot = t;  return;  end if
    end do
    ifconv = 0_8
    cx = sum(xhat*P);  cy = sum(yhat*P);  cz = sum(zhat*P)
    cdx = cx - x0;  cdy = cy - y0;  cdz = cz - z0
    F = cdx*cdx + cdy*cdy + cdz*cdz
    do iter = 1, 20
      q  = (t-tp) / (tp-tpp)
      A  = q*F - q*(q+1.0_r128)*Fp + q*q*Fpp
      B  = (2.0_r128*q+1.0_r128)*F - (1.0_r128+q)**2*Fp + q*q*Fpp
      C  = (1.0_r128+q)*F
      d1 = B + sqrt(B*B - 4.0_r128*A*C)
      d2 = B - sqrt(B*B - 4.0_r128*A*C)
      if (abs(d1) > abs(d2)) then
        dt = -(t-tp)*2.0_r128*C/d1
      else
        dt = -(t-tp)*2.0_r128*C/d2
      end if
      tpp = tp;  Fpp = Fp;  Fp = F;  tp = t
      t   = t + dt
      if (abs(dt) < tol) then;  ifconv = 1_8;  troot = t;  return;  end if
      call line_quad_legendre_cx16(n-1_8, t, P)
      cx = sum(xhat*P);  cy = sum(yhat*P);  cz = sum(zhat*P)
      cdx = cx - x0;  cdy = cy - y0;  cdz = cz - z0
      F = cdx*cdx + cdy*cdy + cdz*cdz
    end do
  end subroutine line_quad_root_refine_r128
#endif

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

#ifndef BIESOLVER_R64_ONLY
  subroutine line_quad_legendre_cx16(n, x, P)
    ! Complex Legendre polynomial values P(0)..P(n) at complex x. (n+1 values)
    integer(8),  intent(in)  :: n
    complex(16), intent(in)  :: x
    complex(16), intent(out) :: P(n+1)
    integer(8) :: l
    P(1) = 1.0_r128
    if (n < 1) return
    P(2) = x
    do l = 1, n-1
      P(l+2) = ((2*l+1)*x*P(l+1) - l*P(l)) / real(l+1, r128)
    end do
  end subroutine line_quad_legendre_cx16

  subroutine line_quad_legendrederiv_cx16(n, x, P, D)
    ! Complex Legendre values P(0)..P(n) and derivatives D(0)..D(n) at complex x.
    integer(8),  intent(in)  :: n
    complex(16), intent(in)  :: x
    complex(16), intent(out) :: P(n+1), D(n+1)
    integer(8) :: l
    P(1) = 1.0_r128;  D(1) = 0.0_r128
    if (n < 1) return
    P(2) = x;  D(2) = 1.0_r128
    do l = 1, n-1
      P(l+2) = ((2*l+1)*x*P(l+1) - l*P(l)) / real(l+1, r128)
      D(l+2) = ((2*l+1)*(P(l+1) + x*D(l+1)) - l*D(l)) / real(l+1, r128)
    end do
  end subroutine line_quad_legendrederiv_cx16
#endif
end module lq_adaptive_mod
