! test_interpolation.f90
!
! Pure-Fortran port of test/interpolation/test_interpolation9.m.
!
! Mirrors, end-to-end, what the MATLAB script does:
!
!   1. Build a Gauss-Legendre panel of order n on [-1,1] (legacy `gauss` +
!      `legeexps`), define the toy curve r(t) = (t, k_par*t^2, t), expand
!      its components in Legendre coefficients, and lay down a grid of
!      targets X(t,d) = r(t) + d*(0, 1, -2*k_par*t).
!
!   2. For every target, compute two estimates of the near-singular
!      integral I = sum_k funvals0(k) * weights(k):
!
!        (a) "local" path, in r64 — a single call to the module routine
!            lq_kernel_mod::build_target_nearroot_weights_local_r64,
!            which is the module-level equivalent of the MATLAB MEX
!            wrapper lqk_build_target_nearroot_weights_local_mex. It
!            returns funvals0, weights, troot, accepted, and the
!            uncompressed local-panel reference I_local.
!
!        (b) "r128" path — orchestrated inline (mirroring the MATLAB
!            utility utils/lqk_build_target_nearroot_weights_r128.m) on
!            top of lq_adaptive_mod (root finder) and lq_kernel_mod
!            (kernel eval + compress_nearroot). Edge structures
!            (sxbd128, sxpbd128, ..., w_bclag128, bclagmatlr128,
!            lq_legmat128) are built once outside the per-target loop.
!
!   3. Reports three error scalars over the full grid:
!
!        max abs diff           = max_p | I_compress(p) - I_local(p) |
!        max rel diff           = same / max(|I_local|, eps)
!        max rel diff vs r128   = max_p | I_compress(p) - I_compress128(p) |
!                                 / max(|I_compress128|, eps)
!
! 05/04/26 Hai

module test_interpolation_callbacks
  use, intrinsic :: iso_c_binding, only: c_double
  use linequaaadrature_mod, only: r128
  implicit none
contains
  subroutine invr_kernel_r64(r_s, tau_s, r0j, kdata, val) bind(C, name="test_interpolation_invr_kernel_r64")
    real(c_double), intent(in) :: r_s(3), tau_s(3), r0j(3), kdata(3)
    real(c_double), intent(inout) :: val
    real(c_double) :: dx, dy, dz, r2, rinv
    integer(8) :: power

    dx = r_s(1) - r0j(1)
    dy = r_s(2) - r0j(2)
    dz = r_s(3) - r0j(3)
    r2 = dx*dx + dy*dy + dz*dz
    rinv = 1.0d0/sqrt(r2)
    power = nint(kdata(1))
    select case (power)
    case (1)
      val = rinv
    case (3)
      val = rinv*rinv*rinv
    case (5)
      val = rinv*rinv*rinv*rinv*rinv
    case default
      val = 0.0d0
    end select
  end subroutine invr_kernel_r64

  subroutine invr_kernel_r128(r_s, tau_s, r0j, kdata, val)
    real(r128), intent(in) :: r_s(3), tau_s(3), r0j(3), kdata(3)
    real(r128), intent(inout) :: val
    real(r128) :: dx, dy, dz, r2, rinv
    integer(8) :: power

    dx = r_s(1) - r0j(1)
    dy = r_s(2) - r0j(2)
    dz = r_s(3) - r0j(3)
    r2 = dx*dx + dy*dy + dz*dz
    rinv = 1.0_r128/sqrt(r2)
    power = nint(real(kdata(1), 8))
    select case (power)
    case (1)
      val = rinv
    case (3)
      val = rinv**3
    case (5)
      val = rinv**5
    case default
      val = 0.0_r128
    end select
  end subroutine invr_kernel_r128
end module test_interpolation_callbacks


program test_interpolation
  use, intrinsic :: iso_c_binding, only: c_funloc, c_funptr
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use linequaaadrature_mod, only: r128, bclaginterpweights_r128
  use lq_kernel_mod, only: build_target_nearroot_weights_local_r64, build_target_nearroot_weights_r64, &
                            build_target_nearroot_weights_r128
  use test_interpolation_callbacks, only: invr_kernel_r64, invr_kernel_r128
  implicit none

  ! ----- problem fixture (matches test_interpolation9.m) -----
  integer(8), parameter :: n = 16_8, ngrid = 150_8, npts = ngrid*ngrid
  real(8),    parameter :: rho = 1.4d0
  real(8),    parameter :: k_par = 0.25d0
  integer(8), parameter :: power_int = 1_8
  real(8),    parameter :: eps_machine = epsilon(1.0d0)

  ! ----- r64 GL setup -----
  real(8) :: tgl(n), wgl(n), dgl(n, n), legmat(n, n), vtmp(n, n), agl(n, n)

  ! ----- source curve -----
  real(8) :: xj(n), yj(n), zj(n)
  real(8) :: dxj(n), dyj(n), dzj(n), spj(n), stauj(3, n)
  real(8) :: xhat(n), yhat(n), zhat(n)

  ! ----- targets -----
  real(8) :: xt(npts), yt(npts), zt(npts)

  ! ----- r128 promotions of the r64 setup -----
  real(r128) :: tgl128(n), wgl128(n), dgl128(n, n), w_bclag128(n)
  real(r128) :: xj128(n), yj128(n), zj128(n), spj128(n)
  real(r128) :: xhat128(n), yhat128(n), zhat128(n)
  real(r128) :: xt128(npts), yt128(npts), zt128(npts)

  ! ----- per-target accumulators -----
  real(8) :: I_all(npts), absdiff_all(npts), reldiff_all(npts)
  real(8) :: reldiff128_all(npts), reldiff_lq64_r128_all(npts), reldiff_local_lq64_all(npts)

  ! ----- per-target scratch (r64 local path) -----
  type(c_funptr) :: invr_fun_r64
  real(8) :: kdata_r64(3)
  real(8) :: funvals0(n), weights(n)
  complex(8) :: troot
  logical    :: accepted
  real(8)    :: I_local, I_compress

  ! ----- per-target scratch (flat r64 LQ path) -----
  complex(8) :: troot64
  logical    :: accepted64
  real(8)    :: funvals64(n), weights64(n)
  real(8)    :: I_lq64

  ! ----- per-target scratch (r128 path) -----
  complex(16) :: troot128
  real(r128)  :: rho128
  logical     :: accepted128
  real(r128)  :: kdata_r128(3), funvals128(n), weights128(n)
  real(r128)  :: I_compress128_r128
  real(8)     :: I_compress128

  integer(8) :: p, kk, ii
  real(8)    :: max_abs_diff, max_rel_diff, max_rel_diff_vs_r128
  real(8)    :: max_rel_lq64_vs_r128, max_rel_local_vs_lq64
  real(8)    :: nan_r64

  ! ----- legacy interface blocks (gauss / legeexps from liblinequad.a) -----
  interface
    subroutine gauss(nn, x, w, d)
      integer(8), intent(in)    :: nn
      real(8),    intent(inout) :: x(nn), w(nn), d(nn, nn)
    end subroutine gauss

    subroutine legeexps(itype, nn, x, u, v, whts)
      integer(8), intent(in)    :: itype, nn
      real(8),    intent(inout) :: x(nn), u(nn, nn), v(nn, nn), whts(nn)
    end subroutine legeexps
  end interface

  ! ============================================================
  ! 1. r64 GL setup, source curve, targets
  ! ============================================================

  call gauss(n, tgl, wgl, dgl)
  call legeexps(2_8, n, tgl, legmat, vtmp, wgl)

  ! Vandermonde-like power-of-tgl matrix; matches MATLAB
  ! linequad_legacy_setup_mex output `agl` and is used inside the
  ! near-root weight builder as part of root-related machinery.
  agl(1, :) = 1.0d0
  do ii = 2, n
    agl(ii, :) = agl(ii - 1, :) * tgl(:)
  end do

  ! Source curve r(t) = (t, k_par t^2, t), evaluated at GL nodes.
  xj = tgl
  yj = k_par*tgl**2
  zj = tgl

  ! Analytic derivative r'(t) = (1, 2 k_par t, 1) -- matches dcurve in
  ! the MATLAB script (which also uses the analytic form rather than
  ! dgl*xj).
  dxj = 1.0d0
  dyj = 2.0d0*k_par*tgl
  dzj = 1.0d0
  spj = sqrt(dxj**2 + dyj**2 + dzj**2)
  do kk = 1, n
    stauj(1, kk) = dxj(kk)/spj(kk)
    stauj(2, kk) = dyj(kk)/spj(kk)
    stauj(3, kk) = dzj(kk)/spj(kk)
  end do

  ! Legendre coefficients used by the root refinement.
  call legendre_expand(n, tgl, wgl, xj, xhat)
  call legendre_expand(n, tgl, wgl, yj, yhat)
  call legendre_expand(n, tgl, wgl, zj, zhat)

  ! ngrid x ngrid target sheet, X(t,d) = r(t) + d*(0,1,-2 k_par t).
  call build_targets(xt, yt, zt)

  ! ============================================================
  ! 2. r128 promotions + edge structures (built once)
  ! ============================================================

  tgl128 = real(tgl, r128)
  wgl128 = real(wgl, r128)
  dgl128 = real(dgl, r128)
  xj128 = real(xj, r128)
  yj128 = real(yj, r128)
  zj128 = real(zj, r128)
  spj128 = real(spj, r128)
  xhat128 = real(xhat, r128)
  yhat128 = real(yhat, r128)
  zhat128 = real(zhat, r128)
  xt128 = real(xt, r128)
  yt128 = real(yt, r128)
  zt128 = real(zt, r128)

  rho128 = real(rho, r128)
  call bclaginterpweights_r128(n, tgl128, w_bclag128)

  ! ============================================================
  ! 3. per-target loop
  ! ============================================================

  invr_fun_r64 = c_funloc(invr_kernel_r64)
  kdata_r64 = 0.0d0
  kdata_r64(1) = real(power_int, 8)

  nan_r64 = ieee_value(1.0d0, ieee_quiet_nan)
  I_all = 0.0d0
  absdiff_all = nan_r64                 ! "nan" in MATLAB output
  reldiff_all = nan_r64
  reldiff128_all = nan_r64
  reldiff_lq64_r128_all = nan_r64
  reldiff_local_lq64_all = nan_r64

  do p = 1, npts
    ! ----- r64 local path (single composite call) -----
    funvals0 = 0.0d0
    weights = 0.0d0
    troot = (0.0d0, 0.0d0)
    accepted = .false.
    I_local = 0.0d0

    call build_target_nearroot_weights_local_r64( &
            n, tgl, wgl, legmat, &
            xj, yj, zj, spj, stauj, &
            xhat, yhat, zhat, &
            rho, xt(p), yt(p), zt(p), &
            invr_fun_r64, kdata_r64, &
            funvals0, weights, &
            troot, accepted, I_local)

    I_compress = 0.0d0
    do kk = 1, n
      I_compress = I_compress + funvals0(kk)*weights(kk)
    end do
    I_all(p) = I_compress

    if (accepted) then
      absdiff_all(p) = abs(I_compress - I_local)
      reldiff_all(p) = absdiff_all(p)/max(abs(I_local), eps_machine)
    end if

    ! ----- flat r64 LQ path -----
    funvals64 = 0.0d0
    weights64 = 0.0d0
    troot64 = (0.0d0, 0.0d0)
    accepted64 = .false.
    I_lq64 = 0.0d0
    call build_target_nearroot_weights_r64( &
            n, tgl, wgl, legmat, &
            xj, yj, zj, spj, stauj, &
            xhat, yhat, zhat, &
            rho, xt(p), yt(p), zt(p), &
            invr_fun_r64, kdata_r64, &
            funvals64, weights64, &
            troot64, accepted64, I_lq64, &
            adaptive_fallback=.false.)

    ! ----- r128 path -----
    kdata_r128 = 0.0_r128
    kdata_r128(1) = real(power_int, r128)
    funvals128 = 0.0_r128
    weights128 = 0.0_r128
    troot128 = cmplx(0.0_r128, 0.0_r128, kind=16)
    accepted128 = .false.
    I_compress128_r128 = 0.0_r128
    call build_target_nearroot_weights_r128( &
            n, tgl128, wgl128, dgl128, w_bclag128, dgl128, &
            xj128, yj128, zj128, spj128, real(stauj, r128), &
            xhat128, yhat128, zhat128, &
            rho128, xt128(p), yt128(p), zt128(p), &
            invr_kernel_r128, kdata_r128, &
            funvals128, weights128, &
            troot128, accepted128, I_compress128_r128)

    I_compress128 = real(I_compress128_r128, 8)
    reldiff128_all(p) = abs(I_compress - I_compress128) / &
                       max(abs(I_compress128), eps_machine)
    reldiff_lq64_r128_all(p) = abs(I_lq64 - I_compress128) / &
                               max(abs(I_compress128), eps_machine)
    reldiff_local_lq64_all(p) = abs(I_compress - I_lq64) / &
                                max(abs(I_lq64), eps_machine)
  end do

  ! ============================================================
  ! 4. Error summary (matches the MATLAB fprintf block)
  ! ============================================================

  max_abs_diff = nan_aware_max(absdiff_all)
  max_rel_diff = nan_aware_max(reldiff_all)
  max_rel_diff_vs_r128 = nan_aware_max(reldiff128_all)
  max_rel_lq64_vs_r128 = nan_aware_max(reldiff_lq64_r128_all)
  max_rel_local_vs_lq64 = nan_aware_max(reldiff_local_lq64_all)

  write (*, '(A)')
  write (*, '(A,I0)') 'Error summary, power = ', power_int
  write (*, '(A,ES28.18)') 'max abs local diff = ', max_abs_diff
  write (*, '(A,ES28.18)') 'max rel local diff = ', max_rel_diff
  write (*, '(A,ES28.18)') 'max rel local-vs-r128 = ', max_rel_diff_vs_r128
  write (*, '(A,ES28.18)') 'max rel local-vs-lq64 = ', max_rel_local_vs_lq64
  write (*, '(A,ES28.18)') 'max rel lq64-vs-r128 = ', max_rel_lq64_vs_r128

contains

  subroutine build_targets(x, y, z)
    real(8), intent(out) :: x(npts), y(npts), z(npts)
    real(8) :: tt, dd
    integer(8) :: jj, kk_local, idx

    idx = 0
    do jj = 1, ngrid
      dd = -2.0d0 + 4.0d0*real(jj - 1, 8)/real(ngrid - 1, 8)
      do kk_local = 1, ngrid
        tt = -2.0d0 + 4.0d0*real(kk_local - 1, 8)/real(ngrid - 1, 8)
        idx = idx + 1
        x(idx) = tt
        y(idx) = k_par*tt*tt + dd
        z(idx) = tt - dd*(2.0d0*k_par*tt)
      end do
    end do
  end subroutine build_targets

  subroutine legendre_expand(nn, tnodes, wnodes, values, coeffs)
    integer(8), intent(in)  :: nn
    real(8),    intent(in)  :: tnodes(nn), wnodes(nn), values(nn)
    real(8),    intent(out) :: coeffs(nn)
    real(8) :: pols(nn)
    integer(8) :: ell

    do ell = 0, nn - 1
      call legendre_values(nn, tnodes, ell, pols)
      coeffs(ell + 1) = 0.5d0*(2*ell + 1)*sum(wnodes*values*pols)
    end do
  end subroutine legendre_expand

  subroutine legendre_values(nn, x, deg, vals)
    integer(8), intent(in)  :: nn, deg
    real(8),    intent(in)  :: x(nn)
    real(8),    intent(out) :: vals(nn)
    real(8) :: p0(nn), p1(nn), p2(nn)
    integer(8) :: ell

    p0 = 1.0d0
    if (deg == 0) then
      vals = p0
      return
    end if
    p1 = x
    if (deg == 1) then
      vals = p1
      return
    end if
    do ell = 1, deg - 1
      p2 = ((2.0d0*ell + 1.0d0)*x*p1 - ell*p0)/real(ell + 1, 8)
      p0 = p1
      p1 = p2
    end do
    vals = p1
  end subroutine legendre_values

  real(r128) function bernstein_radius128(re_t, im_t)
    real(r128), intent(in) :: re_t, im_t
    complex(16) :: t
    t = cmplx(re_t, im_t, kind=16)
    bernstein_radius128 = abs(t + sqrt(t - 1.0_r128)*sqrt(t + 1.0_r128))
  end function bernstein_radius128

  ! Max over an array, skipping NaN entries (matches MATLAB max(...) over
  ! arrays with NaN, which ignores NaN).
  real(8) function nan_aware_max(arr)
    real(8), intent(in) :: arr(npts)
    real(8) :: m
    integer(8) :: ip
    logical :: have_any
    m = 0.0d0
    have_any = .false.
    do ip = 1, npts
      if (.not. (arr(ip) /= arr(ip))) then  ! true if not NaN
        if (.not. have_any) then
          m = arr(ip)
          have_any = .true.
        else if (arr(ip) > m) then
          m = arr(ip)
        end if
      end if
    end do
    if (.not. have_any) then
      nan_aware_max = nan_r64
    else
      nan_aware_max = m
    end if
  end function nan_aware_max

end program test_interpolation
