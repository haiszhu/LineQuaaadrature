! lq_kernel_mex.f90
! Standalone MEX-facing wrapper for line_kernel_eval_r64.
! Accepts integer(8) fptr_int from MATLAB, constructs c_funptr locally,
! then calls lq_kernel_mod :: line_kernel_eval_r64.

subroutine lqk_build_ssq_weights_mex(m, r0, rho_ssq, nbd, sxbd, sbdnp, &
    nquad, tgl, wgl, Legmat, w1, w3, w5, troot, xroot, &
    yroot, zroot, rfc, rfc_ssq)
  use lq_kernel_mod, only: build_ssq_weights_r64
  implicit none
  integer(8), intent(in) :: m, nbd, sbdnp, nquad
  real(8), intent(in) :: r0(3,m), rho_ssq, sxbd(3,nbd)
  real(8), intent(in) :: tgl(nquad), wgl(nquad), Legmat(nquad,nquad)
  real(8), intent(inout) :: w1(nquad,sbdnp,m)
  real(8), intent(inout) :: w3(nquad,sbdnp,m)
  real(8), intent(inout) :: w5(nquad,sbdnp,m)
  complex(8), intent(inout) :: troot(m,sbdnp), xroot(m,sbdnp)
  complex(8), intent(inout) :: yroot(m,sbdnp), zroot(m,sbdnp)
  integer(8), intent(inout) :: rfc(m,sbdnp), rfc_ssq(m)
  call build_ssq_weights_r64(m, r0, rho_ssq, nbd, sxbd, sbdnp, nquad, &
      tgl, wgl, Legmat, w1, w3, w5, troot, xroot, yroot, zroot, rfc, &
      rfc_ssq)
end subroutine lqk_build_ssq_weights_mex

subroutine lqk_estimate_nearroot_lengths_mex(t_root, root_imag_abs, max_len_each_side, len, lenl, lenr)
  use lq_kernel_mod, only: estimate_nearroot_lengths_r64
  implicit none
  real(8),    intent(in)  :: t_root, root_imag_abs
  integer(8), intent(in)  :: max_len_each_side
  integer(8), intent(out) :: len, lenl, lenr

  call estimate_nearroot_lengths_r64(t_root, root_imag_abs, max_len_each_side, len, lenl, lenr)
end subroutine lqk_estimate_nearroot_lengths_mex

subroutine lqk_build_nearroot_nodes_mex(t_root, nquad, tgl, wgl, len, lenl, lenr, t_up, w_up)
  use lq_kernel_mod, only: build_nearroot_nodes_r64
  implicit none
  real(8),    intent(in)    :: t_root
  integer(8), intent(in)    :: nquad, len, lenl, lenr
  real(8),    intent(in)    :: tgl(nquad), wgl(nquad)
  real(8),    intent(inout) :: t_up(nquad*(len-1)), w_up(nquad*(len-1))

  call build_nearroot_nodes_r64(t_root, nquad, tgl, wgl, len, lenl, lenr, t_up, w_up)
end subroutine lqk_build_nearroot_nodes_mex


subroutine lqk_build_nearroot_panels_local_mex( t0, nquad, npan, lenl, lenr, tgl, wgl, &
                                            xjhat, yjhat, zjhat, rbase, &
                                            tpan, upan, wpan, &
                                            xpan, ypan, zpan, &
                                            xdisp, ydisp, zdisp, &
                                            stangpan, sppan, dswpan)
  use lq_kernel_mod, only: build_nearroot_panels_local
  implicit none
  real(8),    intent(in)    :: t0
  integer(8), intent(in)    :: nquad, npan, lenl, lenr
  real(8),    intent(in)    :: tgl(nquad), wgl(nquad)
  real(8),    intent(in)    :: xjhat(nquad), yjhat(nquad), zjhat(nquad)
  real(8),    intent(in)    :: rbase(3)

  real(8),    intent(inout) :: tpan(nquad,npan), upan(nquad,npan), wpan(nquad,npan)
  real(8),    intent(inout) :: xpan(nquad,npan), ypan(nquad,npan), zpan(nquad,npan)
  real(8),    intent(inout) :: xdisp(nquad,npan), ydisp(nquad,npan), zdisp(nquad,npan)
  real(8),    intent(inout) :: stangpan(3,nquad,npan), sppan(nquad,npan), dswpan(nquad,npan)
  call build_nearroot_panels_local( t0, nquad, npan, lenl, lenr, tgl, wgl, & 
                                    xjhat, yjhat, zjhat, rbase, &
                                    tpan, upan, wpan, &
                                    xpan, ypan, zpan, &
                                    xdisp, ydisp, zdisp, &
                                    stangpan, sppan, dswpan)
end subroutine lqk_build_nearroot_panels_local_mex

subroutine lqk_line_kernel_eval_local_mex(nquad, npan, &
                                          xdisp, ydisp, zdisp, stangpan, sppan, dswpan, &
                                          target_loc, fptr_int, kdata, &
                                          funvals_local, integrand0_up, I_local, kval)
  use iso_c_binding, only: c_funptr, c_null_funptr
  use lq_kernel_mod, only: line_kernel_eval_local_r64
  implicit none
  integer(8), intent(in)    :: nquad, npan
  real(8),    intent(in)    :: xdisp(nquad,npan), ydisp(nquad,npan), zdisp(nquad,npan)
  real(8),    intent(in)    :: stangpan(3,nquad,npan)
  real(8),    intent(in)    :: sppan(nquad,npan), dswpan(nquad,npan)
  real(8),    intent(in)    :: target_loc(3)
  integer(8), intent(in)    :: fptr_int
  real(8),    intent(in)    :: kdata(3)

  real(8),    intent(inout) :: funvals_local(nquad,npan)
  real(8),    intent(inout) :: integrand0_up(nquad*npan)
  real(8),    intent(inout) :: I_local
  real(8),    intent(inout) :: kval(nquad,npan)

  type(c_funptr) :: cfptr
  cfptr = transfer(fptr_int, c_null_funptr)
  call line_kernel_eval_local_r64(nquad, npan, &
                                  xdisp, ydisp, zdisp, stangpan, sppan, dswpan, &
                                  target_loc, cfptr, kdata, &
                                  funvals_local, integrand0_up, I_local, kval)
end subroutine lqk_line_kernel_eval_local_mex

subroutine lqk_line_quad_compress_nearroot_local_mex(root_ok_i, nquad, npan, &
                                                    tgl, wgl, legmat, &
                                                    xhat, yhat, zhat, &
                                                    t0, tpan, wpan, &
                                                    target_loc, fptr_int, kdata, &
                                                    integrand0_up, &
                                                    funvals0, sxbdw)
  use iso_c_binding, only: c_funptr, c_null_funptr
  use lq_kernel_mod, only: line_quad_compress_nearroot_local_r64
  implicit none

  integer(8), intent(in) :: root_ok_i, nquad, npan
  real(8), intent(in) :: tgl(nquad), wgl(nquad), legmat(nquad,nquad)
  real(8), intent(in) :: xhat(nquad), yhat(nquad), zhat(nquad)
  real(8), intent(in) :: t0
  real(8), intent(in) :: tpan(nquad,npan), wpan(nquad,npan)
  real(8), intent(in) :: target_loc(3)
  integer(8), intent(in) :: fptr_int
  real(8), intent(in) :: kdata(3)
  real(8), intent(in) :: integrand0_up(nquad*npan)
  real(8), intent(inout) :: funvals0(nquad), sxbdw(nquad)

  logical :: root_ok
  type(c_funptr) :: cfptr

  root_ok = (root_ok_i /= 0_8)
  cfptr = transfer(fptr_int, c_null_funptr)

  call line_quad_compress_nearroot_local_r64(root_ok, nquad, npan, &
                                             tgl, wgl, legmat, &
                                             xhat, yhat, zhat, &
                                             t0, tpan, wpan, &
                                             target_loc, cfptr, kdata, &
                                             integrand0_up, &
                                             funvals0, sxbdw)
end subroutine lqk_line_quad_compress_nearroot_local_mex




subroutine lqk_eval_mex(m, r0, nbd, sbdnp, nquad, &
                                  sxbd, sxpbd, stangbd,    &
                                  fptr_int, kdata, funvals)
  use lq_kernel_mod, only: lqke => line_kernel_eval_r64
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
  call lqke(m, r0, nbd, sbdnp, nquad, sxbd, sxpbd, stangbd, cfptr, kdata, funvals)
end subroutine lqk_eval_mex

subroutine lqk_eval_r128_mex(m, r0, nbd, sbdnp, nquad, &
	                                  sxbd, sxpbd, stangbd,    &
	                                  fptr_int, kdata, funvals)
  use lq_kernel_mod, only: lqke_r128 => line_kernel_eval_r128
  use solidangle_mod, only: invr_kernel_r128
  use iso_c_binding, only: c_char, c_float128, c_funptr, c_int, c_long_long, c_null_char, c_null_funptr
  implicit none
  interface
    function hdf5_write_real128_matrix(file, name, n1, n2, vals) bind(C)
      import c_char, c_float128, c_int, c_long_long
      character(kind=c_char), intent(in) :: file(*), name(*)
      integer(c_long_long), value :: n1, n2
      real(c_float128), intent(in) :: vals(*)
      integer(c_int) :: hdf5_write_real128_matrix
    end function hdf5_write_real128_matrix
  end interface
  integer(8), intent(in)    :: m, nbd, sbdnp, nquad
  real(8),    intent(in)    :: r0(3,m)
  real(8),    intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
  integer(8), intent(in)    :: fptr_int
  real(8),    intent(in)    :: kdata(3,m)
  real(8),    intent(inout) :: funvals(nquad*sbdnp,m)   ! 2D from mwrap; same memory as (nquad,sbdnp,m)

  integer(c_int) :: h5_ok
  type(c_funptr) :: cfptr
  real(16) :: r0_128(3,m), sxbd_128(3,nbd), sxpbd_128(3,nbd), stangbd_128(3,nbd)
  real(16) :: kdata_128(3,m), funvals_128(nquad*sbdnp,m)

  cfptr = transfer(fptr_int, c_null_funptr)
  r0_128 = real(r0, 16)
  sxbd_128 = real(sxbd, 16)
  sxpbd_128 = real(sxpbd, 16)
  stangbd_128 = real(stangbd, 16)
  kdata_128 = real(kdata, 16)
  funvals_128 = real(funvals, 16)
  call lqke_r128(m, r0_128, nbd, sbdnp, nquad, sxbd_128, sxpbd_128, stangbd_128, &
                 invr_kernel_r128, kdata_128, funvals_128)
  h5_ok = hdf5_write_real128_matrix('lqk_eval_r128.h5'//c_null_char, &
                                    '/funvals'//c_null_char, nquad*sbdnp, m, funvals_128)
  funvals = real(funvals_128, 8)
end subroutine lqk_eval_r128_mex

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

