! lq_kernel_mex.f90
! Standalone MEX-facing wrapper for line_kernel_eval_r64.
! Accepts integer(8) fptr_int from MATLAB, constructs c_funptr locally,
! then calls lq_kernel_mod :: line_kernel_eval_r64.

subroutine line_kernel_eval_mex(m, r0, nbd, sbdnp, nquad, &
                                  sxbd, sxpbd, stangbd,    &
                                  fptr_int, kdata, funvals)
  use lq_kernel_mod, only: lke => line_kernel_eval_r64
  use iso_c_binding, only: c_funptr, c_null_funptr
  implicit none
  integer(8), intent(in)    :: m, nbd, sbdnp, nquad
  real(8),    intent(in)    :: r0(3,m)
  real(8),    intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
  integer(8), intent(in)    :: fptr_int
  real(8),    intent(in)    :: kdata(3,m)
  real(8),    intent(inout) :: funvals(nquad*sbdnp,m)   ! 2D from mwrap; same memory as (nquad,sbdnp,m)

  type(c_funptr) :: cfptr
  cfptr = transfer(fptr_int, c_null_funptr)
  ! Fortran allows passing 2D array where 3D is expected when sizes match (legacy)
  call lke(m, r0, nbd, sbdnp, nquad, sxbd, sxpbd, stangbd, cfptr, kdata, funvals)
end subroutine line_kernel_eval_mex

! ------------------------------------------------------------------
! line_quad_compress_mex
! funvals and sxbdw arrive as (nquad*sbdnp, m) from mwrap (2D).
! Fortran passes them straight to line_quad_compress_r64 which
! declares them as (nquad,sbdnp,m) — same contiguous memory.
! ------------------------------------------------------------------
subroutine line_quad_compress_mex(m, r0, nbd, sbdnp, nquad,      &
                                    sxbd, sxpbd, stangbd, sspbd,   &
                                    tgl, wgl, Dgl, w_bclag,        &
                                    Legmat, bclagmatlr,             &
                                    fptr_int, kdata, funvals, sxbdw)
  use lq_kernel_mod, only: lqc => line_quad_compress_r64
  use iso_c_binding, only: c_funptr, c_null_funptr
  implicit none
  integer(8), intent(in)    :: m, nbd, sbdnp, nquad
  real(8),    intent(in)    :: r0(3,m)
  real(8),    intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
  real(8),    intent(in)    :: sspbd(nbd)
  real(8),    intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
  real(8),    intent(in)    :: w_bclag(nquad)
  real(8),    intent(in)    :: Legmat(nquad,nquad), bclagmatlr(nquad,2)
  integer(8), intent(in)    :: fptr_int
  real(8),    intent(in)    :: kdata(3,m)
  real(8),    intent(inout) :: funvals(nquad*sbdnp,m)
  real(8),    intent(inout) :: sxbdw(nquad*sbdnp,m)

  type(c_funptr) :: cfptr
  cfptr = transfer(fptr_int, c_null_funptr)
  call lqc(m, r0, nbd, sbdnp, nquad, sxbd, sxpbd, stangbd, sspbd, &
            tgl, wgl, Dgl, w_bclag, Legmat, bclagmatlr,             &
            cfptr, kdata, funvals, sxbdw)
end subroutine line_quad_compress_mex

! ------------------------------------------------------------------
! line_quad_compress_nearroot_mex
! Like line_quad_compress_mex but wraps line_quad_compress_nearroot_r64.
! Extra args: root_re(m), root_im(m)  — real part / imag part of nearest
!             singularity in parameter space for each target.
!             root_ok_d(m)            — nonzero = use nearroot path,
!                                       0.0     = fall back to adaptive.
! root_ok_d arrives as double from MATLAB; converted to logical here.
! funvals/sxbdw: (nquad*sbdnp, m) from mwrap; same memory as (nquad,sbdnp,m).
! ------------------------------------------------------------------
! ------------------------------------------------------------------
! evaluate_solid_angle_integral_mex
! Thin wrapper: accepts use_nearroot as real(8) (0=false, nonzero=true),
! converts to logical, then calls evaluate_solid_angle_integral_r64.
! ------------------------------------------------------------------
subroutine evaluate_solid_angle_integral_mex(m, tx, n, sx, snx, sw, r_vert, nbd, &
                                               sxbd_in, use_nearroot_d, IalphaAsvestas)
  use solidangle_mod, only: esa => evaluate_solid_angle_integral_r64
  implicit none
  integer(8), intent(in)    :: m, n, nbd
  real(8),    intent(in)    :: tx(3,m), sx(3,n), snx(3,n), sw(n)
  real(8),    intent(in)    :: r_vert(3,3), sxbd_in(3,nbd)
  real(8),    intent(in)    :: use_nearroot_d
  real(8),    intent(inout) :: IalphaAsvestas(m)

  logical :: use_nearroot
  use_nearroot = (use_nearroot_d /= 0.0d0)
  call esa(m, tx, n, sx, snx, sw, r_vert, nbd, sxbd_in, use_nearroot, IalphaAsvestas)
end subroutine evaluate_solid_angle_integral_mex

subroutine line_quad_compress_nearroot_mex(m, r0, nbd, sbdnp, nquad,    &
                                             sxbd, sxpbd, stangbd, sspbd, &
                                             tgl, wgl, Dgl, w_bclag,      &
                                             Legmat, bclagmatlr,           &
                                             fptr_int, kdata, funvals, sxbdw, &
                                             root_re, root_im, root_ok_d)
  use lq_kernel_mod, only: lqn => line_quad_compress_nearroot_r64
  use iso_c_binding, only: c_funptr, c_null_funptr
  implicit none
  integer(8), intent(in)    :: m, nbd, sbdnp, nquad
  real(8),    intent(in)    :: r0(3,m)
  real(8),    intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
  real(8),    intent(in)    :: sspbd(nbd)
  real(8),    intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
  real(8),    intent(in)    :: w_bclag(nquad)
  real(8),    intent(in)    :: Legmat(nquad,nquad), bclagmatlr(nquad,2)
  integer(8), intent(in)    :: fptr_int
  real(8),    intent(in)    :: kdata(3,m)
  real(8),    intent(inout) :: funvals(nquad*sbdnp,m)
  real(8),    intent(inout) :: sxbdw(nquad*sbdnp,m)
  real(8),    intent(in)    :: root_re(m), root_im(m), root_ok_d(m)

  type(c_funptr) :: cfptr
  logical         :: root_ok(m)
  integer(8)      :: j

  cfptr = transfer(fptr_int, c_null_funptr)
  do j = 1, m
    root_ok(j) = (root_ok_d(j) /= 0.0d0)
  end do
  call lqn(m, r0, nbd, sbdnp, nquad, sxbd, sxpbd, stangbd, sspbd, &
            tgl, wgl, Dgl, w_bclag, Legmat, bclagmatlr,             &
            cfptr, kdata, funvals, sxbdw,                           &
            root_re, root_im, root_ok)
end subroutine line_quad_compress_nearroot_mex
