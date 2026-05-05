! lq_kernel_mex.f90
! Standalone MEX-facing wrapper for line_kernel_eval_r64.
! Accepts integer(8) fptr_int from MATLAB, constructs c_funptr locally,
! then calls lq_kernel_mod :: line_kernel_eval_r64.

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

subroutine lqk_line_quad_BrF_mex(nquad, n_up, ncol, &
                                 r_ell, rp_ell, &
                                 tgl, wgl, w_bclag, &
                                 t_up, w_up, &
                                 r0j, kdata, kernel_id, &
                                 Br, integrand0_up)
  use lq_kernel_mod, only: line_quad_BrF_r64
  implicit none
  integer(8), intent(in)    :: nquad, n_up, ncol, kernel_id
  real(8),    intent(in)    :: r_ell(3,nquad), rp_ell(3,nquad)
  real(8),    intent(in)    :: tgl(nquad), wgl(nquad), w_bclag(nquad)
  real(8),    intent(in)    :: t_up(n_up), w_up(n_up)
  real(8),    intent(in)    :: r0j(3), kdata(3)
  real(8),    intent(inout) :: Br(nquad, n_up)
  real(8),    intent(inout) :: integrand0_up(n_up, ncol)

  call line_quad_BrF_r64(nquad, n_up, ncol, &
                         r_ell, rp_ell, &
                         tgl, wgl, w_bclag, &
                         t_up, w_up, &
                         r0j, kdata, kernel_id, &
                         Br, integrand0_up)
end subroutine lqk_line_quad_BrF_mex

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

subroutine lqk_build_target_nearroot_weights_local_mex(nquad, tgl, wgl, legmat, &
                                                       xj, yj, zj, spj, stauj, &
                                                       xjhat, yjhat, zjhat, &
                                                       rho, xtk, ytk, ztk, &
                                                       fptr_int, kdata, &
                                                       funvals0, weights, &
                                                       troot, accepted_i, I_local)
  use iso_c_binding, only: c_funptr, c_null_funptr
  use lq_kernel_mod, only: build_target_nearroot_weights_local_r64
  implicit none
  integer(8), intent(in) :: nquad
  real(8), intent(in) :: tgl(nquad), wgl(nquad), legmat(nquad,nquad)
  real(8), intent(in) :: xj(nquad), yj(nquad), zj(nquad), spj(nquad)
  real(8), intent(in) :: stauj(3,nquad)
  real(8), intent(in) :: xjhat(nquad), yjhat(nquad), zjhat(nquad)
  real(8), intent(in) :: rho, xtk, ytk, ztk
  integer(8), intent(in) :: fptr_int
  real(8), intent(in) :: kdata(3)
  real(8), intent(inout) :: funvals0(nquad), weights(nquad)
  complex(8), intent(inout) :: troot
  integer(8), intent(inout) :: accepted_i
  real(8), intent(inout) :: I_local
  logical :: accepted
  type(c_funptr) :: cfptr
  accepted = (accepted_i /= 0_8)
  cfptr = transfer(fptr_int, c_null_funptr)
  call build_target_nearroot_weights_local_r64(nquad, tgl, wgl, legmat, &
                                               xj, yj, zj, spj, stauj, &
                                               xjhat, yjhat, zjhat, &
                                               rho, xtk, ytk, ztk, &
                                               cfptr, kdata, &
                                               funvals0, weights, &
                                               troot, accepted, I_local)
  accepted_i = merge(1_8, 0_8, accepted)
end subroutine lqk_build_target_nearroot_weights_local_mex

subroutine lqk_build_target_nearroot_weights_mex(nquad, tgl, wgl, legmat, &
                                                 xj, yj, zj, spj, stauj, &
                                                 xjhat, yjhat, zjhat, &
                                                 rho, xtk, ytk, ztk, &
                                                 fptr_int, kdata, &
                                                 funvals0, weights, &
                                                 troot, accepted_i, I_local, &
                                                 kernel_id, adaptive_fallback_i)
  use iso_c_binding, only: c_funptr, c_null_funptr
  use lq_kernel_mod, only: build_target_nearroot_weights_r64
  implicit none
  integer(8), intent(in) :: nquad
  real(8), intent(in) :: tgl(nquad), wgl(nquad), legmat(nquad,nquad)
  real(8), intent(in) :: xj(nquad), yj(nquad), zj(nquad), spj(nquad)
  real(8), intent(in) :: stauj(3,nquad)
  real(8), intent(in) :: xjhat(nquad), yjhat(nquad), zjhat(nquad)
  real(8), intent(in) :: rho, xtk, ytk, ztk
  integer(8), intent(in) :: fptr_int
  real(8), intent(in) :: kdata(3)
  real(8), intent(inout) :: funvals0(nquad), weights(nquad)
  complex(8), intent(inout) :: troot
  integer(8), intent(inout) :: accepted_i
  real(8), intent(inout) :: I_local
  integer(8), intent(in) :: kernel_id, adaptive_fallback_i
  logical :: accepted, adaptive_fallback
  type(c_funptr) :: cfptr

  accepted = (accepted_i /= 0_8)
  adaptive_fallback = (adaptive_fallback_i /= 0_8)
  cfptr = transfer(fptr_int, c_null_funptr)

  call build_target_nearroot_weights_r64(nquad, tgl, wgl, legmat, &
                                         xj, yj, zj, spj, stauj, &
                                         xjhat, yjhat, zjhat, &
                                         rho, xtk, ytk, ztk, &
                                         cfptr, kdata, &
                                         funvals0, weights, &
                                         troot, accepted, I_local, &
                                         kernel_id=kernel_id, &
                                         adaptive_fallback=adaptive_fallback)
  accepted_i = merge(1_8, 0_8, accepted)
end subroutine lqk_build_target_nearroot_weights_mex


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

subroutine lqk_line_quad_compress_nearroot_mex(m, r0, nbd, sbdnp, nquad,  &
                                     sxbd, sxpbd, stangbd, sspbd,       &
                                     tgl, wgl, Dgl, w_bclag,            &
                                     Legmat, bclagmatlr,                 &
                                     fptr_int, kdata, funvals, sxbdw,         &
                                     root_re, root_im, root_ok_i)
  use lq_kernel_mod, only: lqqcn => line_quad_compress_nearroot_r64 
  use iso_c_binding, only: c_funptr, c_null_funptr                                
  integer(8),   intent(in)    :: m, nbd, sbdnp, nquad
  real(8),      intent(in)    :: r0(3,m)
  real(8),      intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
  real(8),      intent(in)    :: sspbd(nbd)
  real(8),      intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
  real(8),      intent(in)    :: w_bclag(nquad)
  real(8),      intent(in)    :: Legmat(nquad,nquad), bclagmatlr(nquad,2)
  integer(8),   intent(in)    :: fptr_int
  real(8),      intent(in)    :: kdata(3,m)
  real(8),      intent(inout) :: funvals(nquad,sbdnp,m)
  real(8),      intent(inout) :: sxbdw(nquad,sbdnp,m)
  real(8),      intent(in)    :: root_re(m), root_im(m)
  integer(8),   intent(in)    :: root_ok_i(m)
  
  type(c_funptr) :: cfptr
  logical         :: root_ok(m)
  integer(8)      :: j

  cfptr = transfer(fptr_int, c_null_funptr)
  do j = 1, m
    root_ok(j) = (root_ok_i(j) /= 0_8)
  end do
  call lqqcn(m, r0, nbd, sbdnp, nquad, & 
            sxbd, sxpbd, stangbd, sspbd, &
            tgl, wgl, Dgl, w_bclag, &
            Legmat, bclagmatlr, &
            cfptr, kdata, funvals, sxbdw, &
            root_re, root_im, root_ok)
end subroutine lqk_line_quad_compress_nearroot_mex

subroutine lqk_line_quad_compress_nearroot_r128_mex(m, r0, nbd, sbdnp, nquad,  &
                                     sxbd, sxpbd, stangbd, sspbd,       &
                                     tgl, wgl, Dgl, w_bclag,            &
                                     Legmat, bclagmatlr,                 &
                                     fptr_int, kdata, funvals, sxbdw,         &
                                     root_re, root_im, root_ok_i)
  use lq_kernel_mod, only: lqqcn_r128 => line_quad_compress_nearroot_r128
  use solidangle_mod, only: invr_kernel_r128
  use iso_c_binding, only: c_char, c_float128, c_funptr, c_int, c_long_long, c_null_char, c_null_funptr
  implicit none
  interface
    function hdf5_write_two_real128_matrices(file, name1, vals1, name2, vals2, n1, n2) bind(C)
      import c_char, c_float128, c_int, c_long_long
      character(kind=c_char), intent(in) :: file(*), name1(*), name2(*)
      real(c_float128), intent(in) :: vals1(*), vals2(*)
      integer(c_long_long), value :: n1, n2
      integer(c_int) :: hdf5_write_two_real128_matrices
    end function hdf5_write_two_real128_matrices
  end interface
  integer(8),   intent(in)    :: m, nbd, sbdnp, nquad
  real(8),      intent(in)    :: r0(3,m)
  real(8),      intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
  real(8),      intent(in)    :: sspbd(nbd)
  real(8),      intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
  real(8),      intent(in)    :: w_bclag(nquad)
  real(8),      intent(in)    :: Legmat(nquad,nquad), bclagmatlr(nquad,2)
  integer(8),   intent(in)    :: fptr_int
  real(8),      intent(in)    :: kdata(3,m)
  real(8),      intent(inout) :: funvals(nquad,sbdnp,m)
  real(8),      intent(inout) :: sxbdw(nquad,sbdnp,m)
  real(8),      intent(in)    :: root_re(m), root_im(m)
  integer(8),   intent(in)    :: root_ok_i(m)
  
  type(c_funptr) :: cfptr
  logical         :: root_ok(m)
  integer(8)      :: j
  integer(c_int)   :: h5_ok

  real(16) :: r0_r128(3,m), kdata_r128(3,m), funvals_r128(nquad,sbdnp,m)
  real(16) :: sxbdw_r128(nquad,sbdnp,m)
  real(16) :: root_re_r128(m), root_im_r128(m)
  real(16) :: sxbd128(3,nbd), sxpbd128(3,nbd), stangbd128(3,nbd), sspbd128(nbd)
  real(16) :: tgl_r128(nquad), wgl_r128(nquad), Dgl_r128(nquad,nquad)
  real(16) :: w_bclag_r128(nquad)
  real(16) :: Legmat_r128(nquad,nquad), bclagmatlr_r128(nquad,2)
  logical   :: root_ok_r128(m)
  cfptr = transfer(fptr_int, c_null_funptr)
  do j = 1, m
    root_ok(j) = (root_ok_i(j) /= 0_8)
  end do

  r0_r128 = real(r0, 16)
  sxbd128 = real(sxbd, 16)
  sxpbd128 = real(sxpbd, 16)
  stangbd128 = real(stangbd, 16)
  sspbd128 = real(sspbd, 16)
  tgl_r128 = real(tgl, 16)
  wgl_r128 = real(wgl, 16)
  Dgl_r128 = real(Dgl, 16)
  w_bclag_r128 = real(w_bclag, 16)
  Legmat_r128 = real(Legmat, 16)
  bclagmatlr_r128 = real(bclagmatlr, 16)
  kdata_r128 = real(kdata, 16)
  funvals_r128 = real(funvals, 16)
  sxbdw_r128 = 0.0_16
  root_re_r128 = real(root_re, 16)
  root_im_r128 = real(root_im, 16)
  root_ok_r128 = root_ok
  call lqqcn_r128(m, r0_r128, nbd, sbdnp, nquad, sxbd128, sxpbd128, stangbd128, sspbd128, &
              tgl_r128, wgl_r128, Dgl_r128, w_bclag_r128, Legmat_r128, bclagmatlr_r128, &
              invr_kernel_r128, kdata_r128, funvals_r128, sxbdw_r128, &
              root_re_r128, root_im_r128, root_ok_r128)
  h5_ok = hdf5_write_two_real128_matrices('lqqcn_r128.h5'//c_null_char, &
                                          '/funvals'//c_null_char, funvals_r128, &
                                          '/sxbdw'//c_null_char, sxbdw_r128, &
                                          nquad*sbdnp, m)
  funvals = real(funvals_r128, 8)
  sxbdw = real(sxbdw_r128, 8)
end subroutine lqk_line_quad_compress_nearroot_r128_mex


subroutine lqk_eval_compress_nearroot_r128_mex(m, r0, nbd, sbdnp, nquad,  &
                                     sxbd, sxpbd, stangbd, sspbd,       &
                                     tgl, wgl, Dgl, w_bclag,            &
                                     Legmat, bclagmatlr,                 &
                                     fptr_int, kdata, funvals, sxbdw,         &
                                     root_re, root_im, root_ok_i)
  use lq_kernel_mod, only: lqke_r128 => line_kernel_eval_r128, &
                           lqqcn_r128 => line_quad_compress_nearroot_r128
  use linequaaadrature_mod, only: bclaginterpweights_r128, legeexps_r128
  use lq_adaptive_mod, only: line_quad_root_initial_guess_r128, line_quad_root_refine_r128
  use solidangle_mod, only: invr_kernel_r128
  use iso_c_binding, only: c_char, c_float128, c_funptr, c_int, c_long_long, c_null_char, c_null_funptr
  implicit none
  interface
    function hdf5_write_two_real128_matrices(file, name1, vals1, name2, vals2, n1, n2) bind(C)
      import c_char, c_float128, c_int, c_long_long
      character(kind=c_char), intent(in) :: file(*), name1(*), name2(*)
      real(c_float128), intent(in) :: vals1(*), vals2(*)
      integer(c_long_long), value :: n1, n2
      integer(c_int) :: hdf5_write_two_real128_matrices
    end function hdf5_write_two_real128_matrices
  end interface
  integer(8),   intent(in)    :: m, nbd, sbdnp, nquad
  real(8),      intent(in)    :: r0(3,m)
  real(8),      intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
  real(8),      intent(in)    :: sspbd(nbd)
  real(8),      intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
  real(8),      intent(in)    :: w_bclag(nquad)
  real(8),      intent(in)    :: Legmat(nquad,nquad), bclagmatlr(nquad,2)
  integer(8),   intent(in)    :: fptr_int
  real(8),      intent(in)    :: kdata(3,m)
  real(8),      intent(inout) :: funvals(nquad,sbdnp,m)
  real(8),      intent(inout) :: sxbdw(nquad,sbdnp,m)
  real(8),      intent(in)    :: root_re(m), root_im(m)
  integer(8),   intent(in)    :: root_ok_i(m)

  integer(8)     :: ell, j, k
  integer(c_int) :: h5_ok
  type(c_funptr) :: cfptr
  logical        :: root_ok_r128(m)
  real(16), parameter :: k_par_r128 = 0.25_16
  real(16) :: denoml, denomr, templ, tempr
  real(8)  :: pols(nquad), p0(nquad), p1(nquad), p2(nquad)
  real(8)  :: xhat(nquad), yhat(nquad), zhat(nquad)
  real(16) :: xj128(nquad), yj128(nquad), zj128(nquad)
  real(16) :: xhat128(nquad), yhat128(nquad), zhat128(nquad)
  real(16) :: r0_r128(3,m), kdata_r128(3,m), funvals_r128(nquad,sbdnp,m)
  real(16) :: sxbdw_r128(nquad,sbdnp,m)
  real(16) :: root_re_r128(m), root_im_r128(m)
  complex(16) :: tinit_r128, troot_r128
  integer(8) :: ifconv_r128
  real(16) :: sxbd128(3,nbd), sxpbd128(3,nbd), stangbd128(3,nbd), sspbd128(nbd)
  real(16) :: tgl_r128(nquad), wgl_r128(nquad), Dgl_r128(nquad,nquad)
  real(16) :: w_bclag_r128(nquad)
  real(16) :: Legmat_r128(nquad,nquad), vtmp_r128(nquad,nquad), bclagmatlr_r128(nquad,2)

  cfptr = transfer(fptr_int, c_null_funptr)
  do j = 1, m
    root_ok_r128(j) = (root_ok_i(j) /= 0_8)
  end do

  r0_r128 = real(r0, 16)
  kdata_r128 = real(kdata, 16)
  tgl_r128 = real(tgl, 16)
  wgl_r128 = real(wgl, 16)
  Dgl_r128 = real(Dgl, 16)
  sxbd128 = real(sxbd, 16)
  sspbd128 = real(sspbd, 16)
  root_re_r128 = real(root_re, 16)
  root_im_r128 = real(root_im, 16)
  sxbdw_r128 = 0.0_16

  xj128 = sxbd128(1, 1:nquad)
  yj128 = sxbd128(2, 1:nquad)
  zj128 = sxbd128(3, 1:nquad)

  do ell = 0, nquad - 1
    p0 = 1.0d0
    if (ell == 0) then
      pols = p0
    else
      p1 = tgl
      if (ell == 1) then
        pols = p1
      else
        do k = 1, ell - 1
          p2 = ((2.0d0*k + 1.0d0)*tgl*p1 - k*p0)/real(k + 1, 8)
          p0 = p1
          p1 = p2
        end do
        pols = p1
      end if
    end if
    xhat(ell + 1) = 0.5d0*(2*ell + 1)*sum(wgl*sxbd(1,1:nquad)*pols)
    yhat(ell + 1) = 0.5d0*(2*ell + 1)*sum(wgl*sxbd(2,1:nquad)*pols)
    zhat(ell + 1) = 0.5d0*(2*ell + 1)*sum(wgl*sxbd(3,1:nquad)*pols)
  end do
  xhat128 = real(xhat, 16)
  yhat128 = real(yhat, 16)
  zhat128 = real(zhat, 16)

  sxpbd128(1, :) = 1.0_16
  sxpbd128(2, :) = 2.0_16*k_par_r128*tgl_r128
  sxpbd128(3, :) = 1.0_16
  do k = 1, nquad
    stangbd128(:, k) = sxpbd128(:, k)/sspbd128(k)
  end do

  call bclaginterpweights_r128(nquad, tgl_r128, w_bclag_r128)
  call legeexps_r128(2_8, nquad, tgl_r128, Legmat_r128, vtmp_r128, wgl_r128)
  denoml = 0.0_16
  denomr = 0.0_16
  do k = 1, nquad
    templ = w_bclag_r128(k)/(-1.0_16 - tgl_r128(k))
    tempr = w_bclag_r128(k)/( 1.0_16 - tgl_r128(k))
    bclagmatlr_r128(k, 1) = templ
    bclagmatlr_r128(k, 2) = tempr
    denoml = denoml + templ
    denomr = denomr + tempr
  end do
  bclagmatlr_r128(:, 1) = bclagmatlr_r128(:, 1)/denoml
  bclagmatlr_r128(:, 2) = bclagmatlr_r128(:, 2)/denomr

  do j = 1, m
    call line_quad_root_initial_guess_r128(tgl_r128, xj128, yj128, zj128, nquad, &
                                           r0_r128(1,j), r0_r128(2,j), r0_r128(3,j), tinit_r128)
    call line_quad_root_refine_r128(xhat128, yhat128, zhat128, nquad, &
                                    r0_r128(1,j), r0_r128(2,j), r0_r128(3,j), &
                                    tinit_r128, troot_r128, ifconv_r128)
    root_re_r128(j) = real(troot_r128, 16)
    root_im_r128(j) = aimag(troot_r128)
    root_ok_r128(j) = root_ok_r128(j) .and. (ifconv_r128 == 1_8)
  end do

  call lqke_r128(m, r0_r128, nbd, sbdnp, nquad, sxbd128, sxpbd128, stangbd128, &
                 invr_kernel_r128, kdata_r128, funvals_r128)
  call lqqcn_r128(m, r0_r128, nbd, sbdnp, nquad, sxbd128, sxpbd128, stangbd128, sspbd128, &
                  tgl_r128, wgl_r128, Dgl_r128, w_bclag_r128, Legmat_r128, bclagmatlr_r128, &
                  invr_kernel_r128, kdata_r128, funvals_r128, sxbdw_r128, &
                  root_re_r128, root_im_r128, root_ok_r128)
  h5_ok = hdf5_write_two_real128_matrices('lqkecn_r128.h5'//c_null_char, &
                                          '/funvals'//c_null_char, funvals_r128, &
                                          '/sxbdw'//c_null_char, sxbdw_r128, &
                                          nquad*sbdnp, m)
  funvals = real(funvals_r128, 8)
  sxbdw = real(sxbdw_r128, 8)
end subroutine lqk_eval_compress_nearroot_r128_mex

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
                                             root_re, root_im, root_ok_i)
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
  real(8),    intent(in)    :: root_re(m), root_im(m)
  integer(8), intent(in)    :: root_ok_i(m)

  type(c_funptr) :: cfptr
  logical         :: root_ok(m)
  integer(8)      :: j

  cfptr = transfer(fptr_int, c_null_funptr)
  do j = 1, m
    root_ok(j) = (root_ok_i(j) /= 0_8)
  end do
  call lqn(m, r0, nbd, sbdnp, nquad, sxbd, sxpbd, stangbd, sspbd, &
            tgl, wgl, Dgl, w_bclag, Legmat, bclagmatlr,             &
            cfptr, kdata, funvals, sxbdw,                           &
            root_re, root_im, root_ok)
end subroutine line_quad_compress_nearroot_mex
