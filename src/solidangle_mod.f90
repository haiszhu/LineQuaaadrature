#ifdef BIESOLVER_WASM_SCALAR_ONLY
#define cpu_time lqs_scalar_noop_cpu_time
#endif
module solidangle_mod
  use linequaaadrature_mod
  use iso_c_binding, only: c_long_long, c_funloc, c_funptr, c_f_procpointer
#ifndef BIESOLVER_R64_ONLY
  use iso_c_binding, only: c_float128
#endif
  implicit none

#ifdef BIESOLVER_WASM_SCALAR_ONLY
  integer(8), parameter :: LQS_PROF = 0_8
#else
  integer(8), save :: LQS_PROF = 0_8
#endif
  integer(8), save :: lqs_rfc_hist(0:8) = 0_8
  integer(8), save :: lqs_len_hist(0:64) = 0_8
  integer(8), save :: lqs_nodes_adap = 0_8
  real(r64),  save :: lqs_t_adap = 0.0_r64
  real(r64),  save :: lqs_t_unif = 0.0_r64

contains

#ifdef BIESOLVER_WASM_SCALAR_ONLY
  subroutine lqs_scalar_noop_cpu_time(value)
    real(r64), intent(out) :: value
    value = 0.0_r64
  end subroutine lqs_scalar_noop_cpu_time
#endif

  ! ------------------------------------------------------------------
  ! evaluate_solid_angle_integral_r64
  ! Inputs:
  !   m           : number of target points
  !   tx(3,m)     : target positions
  !   n           : number of source quadrature nodes (VR grid size)
  !   sx(3,n)     : source positions
  !   snx(3,n)    : source normals (used for qhat)
  !   sw(n)       : source weights
  !   r_vert(3,3) : triangle vertices for circumcircle transform
  !   nbd         : number of boundary quad nodes (= 3*nquad)
  !   sxbd_in(3,nbd) : analytic boundary positions at GL nodes
  !   IalphaAsvestas(m) : output solid angle values
  ! ------------------------------------------------------------------
  subroutine evaluate_solid_angle_integral_r64(m, tx, n, sx, snx, sw, r_vert, nbd, sxbd_in, use_nearroot, IalphaAsvestas)
    use koorn_geom_mod,        only: circumcircle_transform_3d
    use lq_kernel_mod,         only: line_kernel_eval_r64, line_quad_compress_r64, &
                                     line_quad_compress_nearroot_r64, KERNEL_ASVESTAS
    use iso_c_binding,         only: c_funptr, c_funloc
    integer(8),  intent(in)    :: m, n, nbd
    real(r64),   intent(in)    :: tx(3,m), sx(3,n), snx(3,n), sw(n)
    real(r64),   intent(in)    :: r_vert(3,3), sxbd_in(3,nbd)
    logical,     intent(in)    :: use_nearroot
    real(r64),   intent(inout) :: IalphaAsvestas(m)

    integer(8) :: nquad, sbdnp
    integer(8) :: j, k, ell, idx_start, idx_end
    real(r64)  :: alpha, qhat(3), qnrm, R(3,3), c(3)
    real(r64)  :: templ, tempr, denoml, denomr, tgll, tglr
    real(r64)  :: sxp(3)
    type(c_funptr) :: cfptr

    real(r64), allocatable :: tgl(:), wgl(:), Dgl(:,:)
    real(r64), allocatable :: w_bclag(:)
    real(r64), allocatable :: Legmat(:,:), vtmp(:,:)
    real(r64), allocatable :: txnew(:,:), snxnew(:,:), sxpbd(:,:)
    real(r64), allocatable :: sxbd(:,:), stangbd(:,:), sspbd(:)
    real(r64), allocatable :: kdata(:,:)
    real(r64), allocatable :: funvals(:,:,:), sxbdw(:,:,:)
    real(r64), allocatable :: funvals_nr(:,:,:,:), sxbdw_nr(:,:,:,:)
    real(r64), allocatable :: bclagmatlr(:,:)
    real(r64), allocatable :: root_re(:,:), root_im(:,:)
    integer(8), allocatable :: root_ok(:,:)

    sbdnp = 3_8
    nquad = nbd / sbdnp

    allocate(tgl(nquad), wgl(nquad), Dgl(nquad,nquad))
    allocate(w_bclag(nquad))
    allocate(Legmat(nquad,nquad), vtmp(nquad,nquad))
    allocate(txnew(3,1), snxnew(3,n), sxpbd(3,nbd))
    allocate(sxbd(3,nbd), stangbd(3,nbd), sspbd(nbd))
    allocate(kdata(3,1))
    allocate(funvals(nquad,sbdnp,1), sxbdw(nquad,sbdnp,1))
    allocate(bclagmatlr(nquad,2))

    ! --- GL quadrature ---
    call gauss_r64(nquad, tgl, wgl, Dgl)
    call bclaginterpweights_r64(nquad, tgl, w_bclag)
    call legeexps_r64(2_8, nquad, tgl, Legmat, vtmp, wgl)

    ! --- bclagmatlr: barycentric interp weights at endpoints -1, +1 ---
    tgll = -1.0_r64;  tglr = 1.0_r64
    denoml = 0.0_r64;  denomr = 0.0_r64
    do k = 1, nquad
      templ = w_bclag(k) / (tgll - tgl(k))
      tempr = w_bclag(k) / (tglr - tgl(k))
      bclagmatlr(k,1) = templ
      bclagmatlr(k,2) = tempr
      denoml = denoml + templ
      denomr = denomr + tempr
    end do
    bclagmatlr(:,1) = bclagmatlr(:,1) / denoml
    bclagmatlr(:,2) = bclagmatlr(:,2) / denomr

    ! --- copy analytic sxbd, compute stangbd = Dgl*sxbd_ell / |Dgl*sxbd_ell| ---
    sxbd = sxbd_in
    do ell = 1, sbdnp
      idx_start = (ell-1)*nquad + 1
      idx_end   = ell*nquad
      do k = 1, nquad
        sxp(1) = sum(Dgl(k,:) * sxbd(1, idx_start:idx_end))
        sxp(2) = sum(Dgl(k,:) * sxbd(2, idx_start:idx_end))
        sxp(3) = sum(Dgl(k,:) * sxbd(3, idx_start:idx_end))
        sspbd(idx_start+k-1) = sqrt(sxp(1)**2 + sxp(2)**2 + sxp(3)**2)
        stangbd(:, idx_start+k-1) = sxp / sspbd(idx_start+k-1)
      end do
    end do

    ! --- circumcircle transform ---
    R = 0.0_r64;  c = 0.0_r64;  alpha = 0.0_r64
    call circumcircle_transform_3d(r_vert, R, c, alpha)

    do j = 1, n
      snxnew(:,j) = matmul(R, snx(:,j))
    end do
    do k = 1, nbd
      stangbd(:,k) = matmul(R, stangbd(:,k))
    end do
    txnew = 0.0_r64

    ! --- qhat = mean(snxnew) / norm (iside=0) ---
    qhat = 0.0_r64
    do j = 1, n
      qhat = qhat + snxnew(:,j)
    end do
    qhat = qhat / real(n, r64)
    qnrm = sqrt(qhat(1)**2 + qhat(2)**2 + qhat(3)**2)
    qhat = qhat / qnrm

    cfptr = c_funloc(asvestas_kernel_r64)
    kdata(:,1) = qhat

    allocate(funvals_nr(nquad,sbdnp,1,1), sxbdw_nr(nquad,sbdnp,1,1))
    allocate(root_re(1,sbdnp), root_im(1,sbdnp), root_ok(1,sbdnp))

    do j = 1, m

      do k = 1, nbd
        sxbd(:,k) = alpha * matmul(R, sxbd_in(:,k) - tx(:,j))
      end do

      sxpbd = 0.0_r64
      do ell = 1, sbdnp
        idx_start = (ell-1)*nquad + 1
        idx_end   = ell*nquad
        sxpbd(1,idx_start:idx_end) = matmul(Dgl, sxbd(1,idx_start:idx_end))
        sxpbd(2,idx_start:idx_end) = matmul(Dgl, sxbd(2,idx_start:idx_end))
        sxpbd(3,idx_start:idx_end) = matmul(Dgl, sxbd(3,idx_start:idx_end))
      end do

      funvals = 0.0_r64
      call line_kernel_eval_r64(1_8, txnew, nbd, sbdnp, nquad, &
                                sxbd, sxpbd, stangbd, cfptr, kdata, funvals)

      if (use_nearroot) then
        funvals_nr(:,:,:,1) = funvals
        sxbdw_nr = 0.0_r64
        root_re = 0.0_r64;  root_im = 0.0_r64;  root_ok = 0_8
        call line_quad_compress_nearroot_r64(1_8, txnew, nbd, sbdnp, nquad, 1_8, &
                                             sxbd, sxpbd, stangbd, sspbd,        &
                                             tgl, wgl, Dgl, w_bclag,             &
                                             Legmat, bclagmatlr,                 &
                                             cfptr, kdata, funvals_nr, sxbdw_nr, &
                                             root_re, root_im, root_ok, KERNEL_ASVESTAS)
        funvals = funvals_nr(:,:,:,1)
        sxbdw   = sxbdw_nr(:,:,:,1)
      else
        sxbdw = 0.0_r64
        call line_quad_compress_r64(1_8, txnew, nbd, sbdnp, nquad, &
                                    sxbd, sxpbd, stangbd, sspbd,   &
                                    tgl, wgl, Dgl, w_bclag,        &
                                    Legmat, bclagmatlr,            &
                                    cfptr, kdata, funvals, sxbdw)
      end if

      IalphaAsvestas(j) = sum(funvals(:,:,1) * sxbdw(:,:,1))
    end do

    deallocate(funvals_nr, sxbdw_nr, root_re, root_im, root_ok)

    deallocate(tgl, wgl, Dgl, w_bclag, Legmat, vtmp)
    deallocate(sxbd, stangbd, sspbd, txnew, snxnew, sxpbd)
    deallocate(kdata, funvals, sxbdw, bclagmatlr)

  end subroutine evaluate_solid_angle_integral_r64

#ifndef BIESOLVER_R64_ONLY
  subroutine evaluate_solid_angle_integral_r128(m, tx, n, sx, snx, sw, r_vert, nbd, sxbd_in, use_nearroot, IalphaAsvestas)
    use koorn_geom_mod, only: circumcircle_transform_3d_r128
    use lq_kernel_mod,  only: line_kernel_eval_r128, line_quad_compress_r128, &
                              line_quad_compress_nearroot_r128, KERNEL_ASVESTAS
    integer(8),  intent(in)    :: m, n, nbd
    real(r128),  intent(in)    :: tx(3,m), sx(3,n), snx(3,n), sw(n)
    real(r128),  intent(in)    :: r_vert(3,3), sxbd_in(3,nbd)
    logical,     intent(in)    :: use_nearroot
    real(r128),  intent(inout) :: IalphaAsvestas(m)

    integer(8) :: nquad, sbdnp
    integer(8) :: j, k, ell, idx_start, idx_end
    real(r128) :: alpha, qhat(3), qnrm, R(3,3), c(3)
    real(r128) :: templ, tempr, denoml, denomr, tgll, tglr
    real(r128) :: sxp(3)

    real(r128), allocatable :: tgl(:), wgl(:), Dgl(:,:)
    real(r128), allocatable :: w_bclag(:)
    real(r128), allocatable :: Legmat(:,:), vtmp(:,:)
    real(r128), allocatable :: txnew(:,:), snxnew(:,:), sxpbd(:,:)
    real(r128), allocatable :: sxbd(:,:), stangbd(:,:), sspbd(:)
    real(r128), allocatable :: kdata(:,:)
    real(r128), allocatable :: funvals(:,:,:), sxbdw(:,:,:)
    real(r128), allocatable :: funvals_nr(:,:,:,:), sxbdw_nr(:,:,:,:)
    real(r128), allocatable :: bclagmatlr(:,:)
    real(r128), allocatable :: root_re(:,:), root_im(:,:)
    integer(8), allocatable :: root_ok(:,:)

    sbdnp = 3_8
    nquad = nbd / sbdnp

    allocate(tgl(nquad), wgl(nquad), Dgl(nquad,nquad))
    allocate(w_bclag(nquad))
    allocate(Legmat(nquad,nquad), vtmp(nquad,nquad))
    allocate(txnew(3,m), snxnew(3,n), sxpbd(3,nbd))
    allocate(sxbd(3,nbd), stangbd(3,nbd), sspbd(nbd))
    allocate(kdata(3,m))
    allocate(funvals(nquad,sbdnp,m), sxbdw(nquad,sbdnp,m))
    allocate(bclagmatlr(nquad,2))

    ! --- GL quadrature ---
    call gauss_r128(nquad, tgl, wgl, Dgl)
    call bclaginterpweights_r128(nquad, tgl, w_bclag)
    call legeexps_r128(2_8, nquad, tgl, Legmat, vtmp, wgl)

    ! --- bclagmatlr: barycentric interp weights at endpoints -1, +1 ---
    tgll = -1.0_r128;  tglr = 1.0_r128
    denoml = 0.0_r128;  denomr = 0.0_r128
    do k = 1, nquad
      templ = w_bclag(k) / (tgll - tgl(k))
      tempr = w_bclag(k) / (tglr - tgl(k))
      bclagmatlr(k,1) = templ
      bclagmatlr(k,2) = tempr
      denoml = denoml + templ
      denomr = denomr + tempr
    end do
    bclagmatlr(:,1) = bclagmatlr(:,1) / denoml
    bclagmatlr(:,2) = bclagmatlr(:,2) / denomr

    ! --- copy analytic sxbd, compute stangbd = Dgl*sxbd_ell / |Dgl*sxbd_ell| ---
    sxbd = sxbd_in
    do ell = 1, sbdnp
      idx_start = (ell-1)*nquad + 1
      idx_end   = ell*nquad
      do k = 1, nquad
        sxp(1) = sum(Dgl(k,:) * sxbd(1, idx_start:idx_end))
        sxp(2) = sum(Dgl(k,:) * sxbd(2, idx_start:idx_end))
        sxp(3) = sum(Dgl(k,:) * sxbd(3, idx_start:idx_end))
        sspbd(idx_start+k-1) = sqrt(sxp(1)**2 + sxp(2)**2 + sxp(3)**2)
        stangbd(:, idx_start+k-1) = sxp / sspbd(idx_start+k-1)
      end do
    end do

    ! --- circumcircle transform ---
    R = 0.0_r128;  c = 0.0_r128;  alpha = 0.0_r128
    call circumcircle_transform_3d_r128(r_vert, R, c, alpha)

    ! --- apply transform ---
    do j = 1, m
      txnew(:,j) = alpha * matmul(R, tx(:,j) - c)
    end do
    do j = 1, n
      snxnew(:,j) = matmul(R, snx(:,j))
    end do
    do k = 1, nbd
      sxbd(:,k)    = alpha * matmul(R, sxbd(:,k) - c)
      stangbd(:,k) = matmul(R, stangbd(:,k))
    end do

    ! --- qhat = mean(snxnew) / norm ---
    qhat = 0.0_r128
    do j = 1, n
      qhat = qhat + snxnew(:,j)
    end do
    qhat = qhat / real(n, r128)
    qnrm = sqrt(qhat(1)**2 + qhat(2)**2 + qhat(3)**2)
    qhat = qhat / qnrm

    ! --- sxpbd = Dgl * sxbd per panel ---
    sxpbd = 0.0_r128
    do ell = 1, sbdnp
      idx_start = (ell-1)*nquad + 1
      idx_end   = ell*nquad
      sxpbd(1,idx_start:idx_end) = matmul(Dgl, sxbd(1,idx_start:idx_end))
      sxpbd(2,idx_start:idx_end) = matmul(Dgl, sxbd(2,idx_start:idx_end))
      sxpbd(3,idx_start:idx_end) = matmul(Dgl, sxbd(3,idx_start:idx_end))
    end do

    ! --- kdata, kernel eval, compression, accumulate ---
    do j = 1, m
      kdata(:,j) = qhat
    end do

    if (use_nearroot) then
      ! --- nearroot path: kernel eval + nearroot compress (mirrors r64) ---
      funvals = 0.0_r128
      call line_kernel_eval_r128(m, txnew, nbd, sbdnp, nquad, &
                                 sxbd, sxpbd, stangbd, asvestas_kernel_r128, kdata, funvals)
      allocate(funvals_nr(nquad,sbdnp,m,1), sxbdw_nr(nquad,sbdnp,m,1))
      allocate(root_re(m,sbdnp), root_im(m,sbdnp), root_ok(m,sbdnp))
      funvals_nr(:,:,:,1) = funvals
      sxbdw_nr = 0.0_r128
      root_re  = 0.0_r128
      root_im  = 0.0_r128
      root_ok  = 0_8
      call line_quad_compress_nearroot_r128(m, txnew, nbd, sbdnp, nquad, 1_8, &
                                            sxbd, sxpbd, stangbd, sspbd,       &
                                            tgl, wgl, Dgl, w_bclag,            &
                                            Legmat, bclagmatlr,                &
                                            asvestas_kernel_r128, kdata,       &
                                            funvals_nr, sxbdw_nr,              &
                                            root_re, root_im, root_ok,         &
                                            KERNEL_ASVESTAS)
      funvals = funvals_nr(:,:,:,1)
      sxbdw   = sxbdw_nr(:,:,:,1)
      deallocate(funvals_nr, sxbdw_nr, root_re, root_im, root_ok)

    else
      ! --- adaptive bisection path (all panels at once) ---
      funvals = 0.0_r128
      call line_kernel_eval_r128(m, txnew, nbd, sbdnp, nquad, &
                                  sxbd, sxpbd, stangbd, asvestas_kernel_r128, kdata, funvals)

      sxbdw = 0.0_r128
      call line_quad_compress_r128(m, txnew, nbd, sbdnp, nquad,     &
                                    sxbd, sxpbd, stangbd, sspbd,     &
                                    tgl, wgl, Dgl, w_bclag,          &
                                    Legmat, bclagmatlr,               &
                                    asvestas_kernel_r128, kdata, funvals, sxbdw)
    end if

    do j = 1, m
      IalphaAsvestas(j) = sum(funvals(:,:,j) * sxbdw(:,:,j))
    end do

    deallocate(tgl, wgl, Dgl, w_bclag, Legmat, vtmp)
    deallocate(sxbd, stangbd, sspbd, txnew, snxnew, sxpbd)
    deallocate(kdata, funvals, sxbdw, bclagmatlr)

  end subroutine evaluate_solid_angle_integral_r128

  ! could evaluate line integral along sbdnp panels with each panel nquad nodes 
  subroutine evaluate_line_integral_r128(m, r0, nbd, sbdnp, nquad,  &
                                     sxbd, sxpbd, stangbd, sspbd,   &
                                     fun, kdata, q_lq128)
    use lq_kernel_mod, only: kernel_iface_r128, line_kernel_eval_r128, &
                              line_quad_compress_nearroot_r128, KERNEL_INVR
    integer(8),     intent(in)     :: m, nbd, sbdnp, nquad
    real(r128),      intent(in)    :: r0(3,m)
    real(r128),      intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
    real(r128),      intent(in)    :: sspbd(nbd)
    procedure(kernel_iface_r128)   :: fun
    real(r128),  intent(in)        :: kdata(3,m)
    real(r128),  intent(inout)     :: q_lq128(m)

    integer(8) :: ell, j, k, idx_start, idx_end
    real(r128), allocatable :: tgl(:), wgl(:), Dgl(:,:), w_bclag(:)
    real(r128), allocatable :: Legmat(:,:), vtmp(:,:), bclagmatlr(:,:)
    real(r128), allocatable :: funvals_ell(:,:,:,:), sxbdw_ell(:,:,:,:)
    real(r128), allocatable :: root_re(:,:), root_im(:,:)
    integer(8), allocatable :: root_ok(:,:)
    real(r128) :: templ, tempr, denoml, denomr

    allocate(tgl(nquad), wgl(nquad), Dgl(nquad,nquad), w_bclag(nquad))
    allocate(Legmat(nquad,nquad), vtmp(nquad,nquad), bclagmatlr(nquad,2))
    allocate(root_re(m,1), root_im(m,1), root_ok(m,1))
    allocate(funvals_ell(nquad,1,m,1), sxbdw_ell(nquad,1,m,1))

    call gauss_r128(nquad, tgl, wgl, Dgl)
    call bclaginterpweights_r128(nquad, tgl, w_bclag)
    call legeexps_r128(2_8, nquad, tgl, Legmat, vtmp, wgl)

    denoml = 0.0_r128
    denomr = 0.0_r128
    do k = 1, nquad
      templ = w_bclag(k) / (-1.0_r128 - tgl(k))
      tempr = w_bclag(k) / ( 1.0_r128 - tgl(k))
      bclagmatlr(k,1) = templ
      bclagmatlr(k,2) = tempr
      denoml = denoml + templ
      denomr = denomr + tempr
    end do
    bclagmatlr(:,1) = bclagmatlr(:,1) / denoml
    bclagmatlr(:,2) = bclagmatlr(:,2) / denomr

    q_lq128 = 0.0_r128
    do ell = 1, sbdnp
      idx_start = (ell-1_8)*nquad + 1_8
      idx_end = ell*nquad

      funvals_ell = 0.0_r128
      call line_kernel_eval_r128(m, r0, nquad, 1_8, nquad, &
                                 sxbd(:,idx_start:idx_end), &
                                 sxpbd(:,idx_start:idx_end), &
                                 stangbd(:,idx_start:idx_end), &
                                 fun, kdata, funvals_ell(:,:,:,1))

      sxbdw_ell = 0.0_r128
      root_re = 0.0_r128
      root_im = 0.0_r128
      root_ok = 0_8
      call line_quad_compress_nearroot_r128(m, r0, nquad, 1_8, nquad, 1_8, &
                                            sxbd(:,idx_start:idx_end), &
                                            sxpbd(:,idx_start:idx_end), &
                                            stangbd(:,idx_start:idx_end), &
                                            sspbd(idx_start:idx_end), &
                                            tgl, wgl, Dgl, w_bclag, &
                                            Legmat, bclagmatlr, &
                                            fun, kdata, funvals_ell, sxbdw_ell, &
                                            root_re, root_im, root_ok, KERNEL_INVR)

      do j = 1, m
        if (root_ok(j,1) /= 0_8) then
          q_lq128(j) = q_lq128(j) + sum(funvals_ell(:,1,j,1) * sxbdw_ell(:,1,j,1))
        end if
      end do
    end do

    deallocate(tgl, wgl, Dgl, w_bclag, Legmat, vtmp, bclagmatlr)
    deallocate(root_re, root_im, root_ok, funvals_ell, sxbdw_ell)

  end subroutine evaluate_line_integral_r128
#endif


    ! ------------------------------------------------------------------
  ! eval_moments_funvals_r64
  ! Strict Fortran port of qotential/demo1.m :: momentsalladapt.
  !
  ! Output: funvals_pre(nbd, 2*(order+1), m) — packed [N | M] moments
  ! at the original source quadrature nodes for each target.
  !   - root accepted: near-root upsample + line_quad_BrF_r64 with
  !     KERNEL_MOMENTS_MN; funvals_pre_panel = Br * f_up.
  !   - root not accepted: moments_kernel_r64 direct at each panel
  !     source node (no bisection fallback).
  ! ------------------------------------------------------------------
  subroutine eval_moments_funvals_r64(m, tx, nbd, sxbd, nquad, order, &
                                      funvals_pre)
    use lq_kernel_mod,   only: line_quad_BrF_r64, KERNEL_MOMENTS_MN, &
                                estimate_nearroot_lengths_r64,        &
                                build_nearroot_nodes_r64
    use lq_adaptive_mod, only: line_quad_root_initial_guess_r64,      &
                                line_quad_root_refine_r64

    integer(8), intent(in)    :: m, nbd, nquad, order
    real(r64),  intent(in)    :: tx(3, m)
    real(r64),  intent(in)    :: sxbd(3, nbd)
    real(r64),  intent(inout) :: funvals_pre(nbd, 2_8*(order+1_8), m)

    ! ---- constants (mirroring momentsalladapt / solid-angle near-root) ----
    integer(8), parameter :: max_len_each_side = 12_8
    integer(8), parameter :: maxpan            = 128_8
    real(r64),  parameter :: rho               = 4.0_r64

    integer(8) :: ncol, npan, n_expa
    integer(8) :: ell, j, k, idx_start, idx_end, i
    integer(8) :: len, lenl, lenr, n_up, ifconv
    real(r64)  :: br_init, br_root
    complex(8) :: tinit, troot
    logical    :: accepted

    ! ---- GL quadrature + interp setup (allocated on entry, freed on exit) --
    real(r64), allocatable :: tgl(:), wgl(:), Dgl(:,:)
    real(r64), allocatable :: w_bclag(:)
    real(r64), allocatable :: Legmat(:,:), vtmp(:,:)

    ! ---- per-panel scratch ----
    real(r64), allocatable :: r_ell(:,:), rp_ell(:,:)            ! (3, nquad)
    real(r64), allocatable :: xj(:), yj(:), zj(:)                ! (nquad)
    real(r64), allocatable :: xyz_hat(:,:)                       ! (n_expa, 3)

    ! ---- near-root upsample scratch ----
    real(r64), allocatable :: t_up(:), w_up(:)                   ! (maxpan*nquad)
    real(r64), allocatable :: Br(:,:), f_up(:,:)                 ! (nquad, n_up), (n_up, ncol)

    ! ---- naive-branch scratch ----
    real(r64) :: r_s(3), tau_s(3), kdata_zero(3)
    real(r64), allocatable :: val_naive(:)                        ! (ncol)

    ! =================================================================
    ! Setup
    ! =================================================================
    npan   = nbd / nquad
    ncol   = 2_8 * (order + 1_8)
    n_expa = min(16_8, nquad)

    allocate(tgl(nquad), wgl(nquad), Dgl(nquad,nquad))
    allocate(w_bclag(nquad))
    allocate(Legmat(nquad,nquad), vtmp(nquad,nquad))

    call gauss_r64(nquad, tgl, wgl, Dgl)
    call bclaginterpweights_r64(nquad, tgl, w_bclag)
    call legeexps_r64(2_8, nquad, tgl, Legmat, vtmp, wgl)

    allocate(r_ell(3,nquad), rp_ell(3,nquad))
    allocate(xj(nquad), yj(nquad), zj(nquad))
    allocate(xyz_hat(n_expa, 3))

    allocate(t_up(maxpan*nquad), w_up(maxpan*nquad))
    allocate(Br(nquad, maxpan*nquad), f_up(maxpan*nquad, ncol))
    allocate(val_naive(ncol))

    funvals_pre = 0.0_r64
    kdata_zero  = 0.0_r64

    ! =================================================================
    ! Per-panel × per-target loop
    ! =================================================================
    do ell = 1, npan
      idx_start = (ell-1_8)*nquad + 1_8
      idx_end   = ell*nquad

      ! ---- per-panel xyz coords + derivatives + Legendre projection ----
      r_ell(1,:) = sxbd(1, idx_start:idx_end)
      r_ell(2,:) = sxbd(2, idx_start:idx_end)
      r_ell(3,:) = sxbd(3, idx_start:idx_end)
      xj = r_ell(1,:);  yj = r_ell(2,:);  zj = r_ell(3,:)

      rp_ell(1,:) = matmul(Dgl, xj)
      rp_ell(2,:) = matmul(Dgl, yj)
      rp_ell(3,:) = matmul(Dgl, zj)

      xyz_hat(:,1) = matmul(Legmat(1:n_expa,:), xj)
      xyz_hat(:,2) = matmul(Legmat(1:n_expa,:), yj)
      xyz_hat(:,3) = matmul(Legmat(1:n_expa,:), zj)

      do j = 1, m
        ! -------- root finding --------
        accepted = .false.
        call line_quad_root_initial_guess_r64(                         &
             tgl, xj, yj, zj, nquad,                                   &
             tx(1,j), tx(2,j), tx(3,j), tinit)
        ! bernstein_radius(re, im) = |z + sqrt(z-1)*sqrt(z+1)|, z = re + i*im
        br_init = abs(tinit + sqrt(tinit - 1.0_r64)*sqrt(tinit + 1.0_r64))

        if (br_init < 1.5_r64 * rho) then
          troot  = cmplx(0.0_r64, 0.0_r64, kind=r64)
          ifconv = 0
          call line_quad_root_refine_r64(                              &
               xyz_hat(:,1), xyz_hat(:,2), xyz_hat(:,3), n_expa,       &
               tx(1,j), tx(2,j), tx(3,j),                              &
               tinit, troot, ifconv)
          br_root = abs(troot + sqrt(troot - 1.0_r64)*sqrt(troot + 1.0_r64))
          if (ifconv == 1 .and. br_root < rho) accepted = .true.
        end if

        if (accepted) then
          ! ============ accepted: near-root BrF compress ============
          call estimate_nearroot_lengths_r64(                          &
               real(troot, r64), abs(aimag(troot)),                    &
               max_len_each_side, len, lenl, lenr)
          n_up = nquad * (len - 1_8)
          if (n_up > maxpan*nquad) error stop                          &
            'eval_moments_funvals_r64: n_up exceeds maxpan*nquad'

          call build_nearroot_nodes_r64(                               &
               real(troot, r64), nquad, tgl, wgl, len, lenl, lenr,     &
               t_up(1:n_up), w_up(1:n_up))

          call line_quad_BrF_r64(                                      &
               nquad, n_up, ncol,                                      &
               r_ell, rp_ell,                                          &
               tgl, wgl, w_bclag,                                      &
               t_up(1:n_up), w_up(1:n_up),                             &
               tx(:, j), real(troot, r64), aimag(troot),               &
               kdata_zero, KERNEL_MOMENTS_MN,                          &
               Br(:, 1:n_up), f_up(1:n_up, 1:ncol))

          ! Both halves at once (matches funvals_pre packed [N | M]).
          funvals_pre(idx_start:idx_end, 1:ncol, j) =                  &
               matmul(Br(:, 1:n_up), f_up(1:n_up, 1:ncol))

        else
          ! ============ not accepted: naive direct moments eval ============
          do k = 1, nquad
            i = idx_start + k - 1_8
            r_s   = sxbd(:, i)
            tau_s = rp_ell(:, k)        ! ignored inside moments_kernel_r64
            val_naive = 0.0_r64
            call moments_kernel_r64(r_s, tau_s, tx(:,j), kdata_zero,   &
                                    ncol, val_naive)
            funvals_pre(i, 1:ncol, j) = val_naive
          end do
        end if
      end do
    end do

    deallocate(tgl, wgl, Dgl, w_bclag, Legmat, vtmp)
    deallocate(r_ell, rp_ell, xj, yj, zj, xyz_hat)
    deallocate(t_up, w_up, Br, f_up, val_naive)

  end subroutine eval_moments_funvals_r64

  ! ------------------------------------------------------------------
  ! asvestas_kernel_r64  (bind(C) enables c_funloc + dlsym lookup)
  ! kdata(1:3) = qhat
  ! val = -[tau_s . (qhat x rhat)] / [|r| * (1 - qhat.rhat)]
  ! ------------------------------------------------------------------
  subroutine asvestas_kernel_r64(r_s, tau_s, r0j, kdata, val) bind(C)
    real(r64), intent(in)    :: r_s(3), tau_s(3), r0j(3), kdata(3)
    real(r64), intent(inout) :: val

    real(r64) :: rvec(3), rdist, rhat(3), qhat(3), qcrossrhat(3), qdotrhat

    qhat  = kdata(1:3)
    rvec  = r_s - r0j
    rdist = sqrt(rvec(1)**2 + rvec(2)**2 + rvec(3)**2)
    rhat  = rvec / rdist

    qcrossrhat(1) = qhat(2)*rhat(3) - qhat(3)*rhat(2)
    qcrossrhat(2) = qhat(3)*rhat(1) - qhat(1)*rhat(3)
    qcrossrhat(3) = qhat(1)*rhat(2) - qhat(2)*rhat(1)

    qdotrhat = qhat(1)*rhat(1) + qhat(2)*rhat(2) + qhat(3)*rhat(3)

    val = -(tau_s(1)*qcrossrhat(1) + tau_s(2)*qcrossrhat(2) + tau_s(3)*qcrossrhat(3)) &
           / (rdist * (1.0_r64 - qdotrhat))

  end subroutine asvestas_kernel_r64

  ! ------------------------------------------------------------------
  ! invr_kernel_r64  (bind(C) enables c_funloc + dlsym lookup)
  ! kdata(1) = power in {1,3,5}
  ! val = 1 / |r|^power
  ! ------------------------------------------------------------------
  subroutine invr_kernel_r64(r_s, tau_s, r0j, kdata, val) bind(C)
    real(r64), intent(in)    :: r_s(3), tau_s(3), r0j(3), kdata(3)
    real(r64), intent(inout) :: val

    real(r64) :: dx, dy, dz, r2, rinv
    integer :: power

    dx = r_s(1) - r0j(1)
    dy = r_s(2) - r0j(2)
    dz = r_s(3) - r0j(3)
    r2 = dx*dx + dy*dy + dz*dz
    rinv = 1.0_r64 / sqrt(r2)
    power = nint(kdata(1))

    select case (power)
    case (1)
      val = rinv
    case (3)
      val = rinv*rinv*rinv
    case (5)
      val = rinv*rinv*rinv*rinv*rinv
    case default
      val = 0.0_r64
    end select
  end subroutine invr_kernel_r64

  ! ------------------------------------------------------------------
  ! moments_kernel_r64  (vector bind(C) kernel)
  ! kdata(2) = flag
  ! dim = 2*(order+1)
  ! val = [N_0..N_order, M_0..M_order] at one source point
  ! ------------------------------------------------------------------
  subroutine moments_kernel_r64(r_s, tau_s, r0j, kdata, dim, val) bind(C)
    real(r64), intent(in)    :: r_s(3), tau_s(3), r0j(3), kdata(3)
    integer(c_long_long), value :: dim
    real(r64), intent(inout) :: val(dim)

    integer(8) :: order, flag, k
    real(r64) :: r0norm, r0norm_inv
    real(r64) :: r0dotr, rnorm, rnorm_inv, rnorm2_inv
    real(r64) :: r0dotr_over_rnorm2, r0norm2_over_rnorm2
    real(r64) :: r0mr_vec(3), r0mr, r0mr_inv, r0mr_over_rnorm2
    real(r64) :: r0dotr0mr_over_r0mr, rdotr0mr_over_r0mr
    real(r64) :: denominator1, denominator2, LMNcommon
    real(r64) :: LMcommon, r0normplusrnorm

    order = dim/2_8 - 1_8
    flag = nint(kdata(2), 8)
    val = 0.0_r64

    r0dotr = r0j(1)*r_s(1) + r0j(2)*r_s(2) + r0j(3)*r_s(3)
    rnorm = sqrt(r_s(1)**2 + r_s(2)**2 + r_s(3)**2)
    r0norm = sqrt(r0j(1)**2 + r0j(2)**2 + r0j(3)**2)
    rnorm_inv = 1.0_r64/rnorm
    rnorm2_inv = rnorm_inv**2
    r0norm_inv = 1.0_r64/r0norm
    r0dotr_over_rnorm2 = r0dotr*rnorm2_inv
    r0norm2_over_rnorm2 = r0norm**2*rnorm2_inv
    r0normplusrnorm = r0norm+rnorm

    r0mr_vec = r0j - r_s
    r0mr = sqrt(r0mr_vec(1)**2 + r0mr_vec(2)**2 + r0mr_vec(3)**2)
    r0mr_over_rnorm2 = r0mr*rnorm2_inv
    r0mr_inv = 1.0_r64/r0mr

    r0dotr0mr_over_r0mr = (r0j(1)*r0mr_vec(1) + r0j(2)*r0mr_vec(2) + &
                           r0j(3)*r0mr_vec(3))*r0mr_inv
    denominator1 = r0norm + r0dotr0mr_over_r0mr
    rdotr0mr_over_r0mr = (r_s(1)*r0mr_vec(1) + r_s(2)*r0mr_vec(2) + &
                          r_s(3)*r0mr_vec(3))*r0mr_inv
    denominator2 = rnorm + rdotr0mr_over_r0mr
    LMNcommon = r0normplusrnorm/(denominator1+denominator2)

    val(1) = log((r0normplusrnorm+r0mr)*LMNcommon*r0mr_inv)*rnorm_inv
    if (order >= 1_8) then
      val(2) = val(1)*r0dotr_over_rnorm2 + (r0mr-r0norm)*rnorm2_inv
    end if
    do k=2_8,order
      val(k+1_8) = real(2_8*k-1_8,r64)/real(k,r64)*r0dotr_over_rnorm2*val(k) - &
                   real(k-1_8,r64)/real(k,r64)*r0norm2_over_rnorm2*val(k-1_8) + &
                   1.0_r64/real(k,r64)*r0mr_over_rnorm2
    end do

    LMcommon = 1.0_r64/(r0norm*rnorm+r0dotr)
    val(order+2_8) = LMNcommon*LMcommon*(((r0normplusrnorm)*r0mr_inv+rnorm*r0norm_inv)* &
                       r0mr_inv-r0norm_inv)
    if (order >= 1_8) then
      val(order+3_8) = r0norm*val(order+2_8)/(r0norm+r0mr)
    end if
    do k=2_8,order
      val(order+2_8+k) = (r0dotr*val(order+1_8+k) + real(k-1_8,r64)*val(k-1_8) - &
                          r0mr_inv)*rnorm2_inv
    end do
  end subroutine moments_kernel_r64

#ifndef BIESOLVER_R64_ONLY
  subroutine invr_kernel_r128(r_s, tau_s, r0j, kdata, val) bind(C)
    real(c_float128), intent(in)    :: r_s(3), tau_s(3), r0j(3), kdata(3)
    real(c_float128), intent(inout) :: val

    real(r128) :: dx, dy, dz, r2, rinv
    integer(8) :: power

    dx = r_s(1) - r0j(1)
    dy = r_s(2) - r0j(2)
    dz = r_s(3) - r0j(3)
    r2 = dx*dx + dy*dy + dz*dz
    rinv = 1.0_r128 / sqrt(r2)
    power = nint(kdata(1), 8)

    select case (power)
    case (1_8)
      val = rinv
    case (3_8)
      val = rinv*rinv*rinv
    case (5_8)
      val = rinv*rinv*rinv*rinv*rinv
    case default
      val = 0.0_r128
    end select
  end subroutine invr_kernel_r128

  subroutine asvestas_kernel_r128(r_s, tau_s, r0j, kdata, val)
    real(r128), intent(in)    :: r_s(3), tau_s(3), r0j(3), kdata(3)
    real(r128), intent(inout) :: val

    real(r128) :: rvec(3), rdist, rhat(3), qhat(3), qcrossrhat(3), qdotrhat

    qhat  = kdata(1:3)
    rvec  = r_s - r0j
    rdist = sqrt(rvec(1)**2 + rvec(2)**2 + rvec(3)**2)
    rhat  = rvec / rdist

    qcrossrhat(1) = qhat(2)*rhat(3) - qhat(3)*rhat(2)
    qcrossrhat(2) = qhat(3)*rhat(1) - qhat(1)*rhat(3)
    qcrossrhat(3) = qhat(1)*rhat(2) - qhat(2)*rhat(1)

    qdotrhat = qhat(1)*rhat(1) + qhat(2)*rhat(2) + qhat(3)*rhat(3)

    val = -(tau_s(1)*qcrossrhat(1) + tau_s(2)*qcrossrhat(2) + tau_s(3)*qcrossrhat(3)) &
           / (rdist * (1.0_r128 - qdotrhat))

  end subroutine asvestas_kernel_r128

  subroutine moments_kernel_r128(r_s, tau_s, r0j, kdata, dim, val) bind(C)
    real(c_float128), intent(in)    :: r_s(3), tau_s(3), r0j(3), kdata(3)
    integer(c_long_long), value     :: dim
    real(c_float128), intent(inout) :: val(dim)

    integer(8) :: order, flag, k
    real(r128) :: r0norm, r0norm_inv
    real(r128) :: r0dotr, rnorm, rnorm_inv, rnorm2_inv
    real(r128) :: r0dotr_over_rnorm2, r0norm2_over_rnorm2
    real(r128) :: r0mr_vec(3), r0mr, r0mr_inv, r0mr_over_rnorm2
    real(r128) :: r0dotr0mr_over_r0mr, rdotr0mr_over_r0mr
    real(r128) :: denominator1, denominator2, LMNcommon
    real(r128) :: LMcommon, r0normplusrnorm

    order = dim/2_8 - 1_8
    flag = nint(kdata(2), 8)
    val = 0.0_r128

    r0dotr = r0j(1)*r_s(1) + r0j(2)*r_s(2) + r0j(3)*r_s(3)
    rnorm = sqrt(r_s(1)**2 + r_s(2)**2 + r_s(3)**2)
    r0norm = sqrt(r0j(1)**2 + r0j(2)**2 + r0j(3)**2)
    rnorm_inv = 1.0_r128/rnorm
    rnorm2_inv = rnorm_inv**2
    r0norm_inv = 1.0_r128/r0norm
    r0dotr_over_rnorm2 = r0dotr*rnorm2_inv
    r0norm2_over_rnorm2 = r0norm**2*rnorm2_inv
    r0normplusrnorm = r0norm+rnorm

    r0mr_vec = r0j - r_s
    r0mr = sqrt(r0mr_vec(1)**2 + r0mr_vec(2)**2 + r0mr_vec(3)**2)
    r0mr_over_rnorm2 = r0mr*rnorm2_inv
    r0mr_inv = 1.0_r128/r0mr

    r0dotr0mr_over_r0mr = (r0j(1)*r0mr_vec(1) + r0j(2)*r0mr_vec(2) + &
                           r0j(3)*r0mr_vec(3))*r0mr_inv
    denominator1 = r0norm + r0dotr0mr_over_r0mr
    rdotr0mr_over_r0mr = (r_s(1)*r0mr_vec(1) + r_s(2)*r0mr_vec(2) + &
                          r_s(3)*r0mr_vec(3))*r0mr_inv
    denominator2 = rnorm + rdotr0mr_over_r0mr
    LMNcommon = r0normplusrnorm/(denominator1+denominator2)

    val(1) = log((r0normplusrnorm+r0mr)*LMNcommon*r0mr_inv)*rnorm_inv
    if (order >= 1_8) then
      val(2) = val(1)*r0dotr_over_rnorm2 + (r0mr-r0norm)*rnorm2_inv
    end if
    do k=2_8,order
      val(k+1_8) = real(2_8*k-1_8,r128)/real(k,r128)*r0dotr_over_rnorm2*val(k) - &
                   real(k-1_8,r128)/real(k,r128)*r0norm2_over_rnorm2*val(k-1_8) + &
                   1.0_r128/real(k,r128)*r0mr_over_rnorm2
    end do

    LMcommon = 1.0_r128/(r0norm*rnorm+r0dotr)
    val(order+2_8) = LMNcommon*LMcommon*(((r0normplusrnorm)*r0mr_inv+rnorm*r0norm_inv)* &
                       r0mr_inv-r0norm_inv)
    if (order >= 1_8) then
      val(order+3_8) = r0norm*val(order+2_8)/(r0norm+r0mr)
    end if
    do k=2_8,order
      val(order+2_8+k) = (r0dotr*val(order+1_8+k) + real(k-1_8,r128)*val(k-1_8) - &
                          r0mr_inv)*rnorm2_inv
    end do
  end subroutine moments_kernel_r128
#endif

  subroutine evaluate_solid_angle_integral_fast_r64(m, r0, nbd, sbdnp, nquad, &
                                           sxbd, stangbd, sspbd, &
                                           len1, sxbd1, stangbd1, swbd1, &
                                           len2, sxbd2, stangbd2, swbd2, &
                                           len3, sxbd3, stangbd3, swbd3, &
                                           qhat, tgl, wgl, Dgl, w_bclag, bclagmatlr, &
                                           troot, xroot, yroot, zroot, rfc, &
                                           IalphaAsvestas, rho_in, &
                                           sxbd_raw, tx_raw, Rfr, alpha_fr, Legmat)
    use lq_kernel_mod,   only: estimate_nearroot_lengths_r64, &
                               build_nearroot_nodes_r64, bary_row_r64, &
                               update_refinement_codes_r64
    use lq_adaptive_mod, only: line_quad_root_initial_guess_r64, &
                               line_quad_root_refine_r64
    integer(8),   intent(in)    :: m, nquad, nbd, sbdnp
    real(r64),    intent(in)    :: r0(3,m), qhat(3)
    real(r64),    intent(in)    :: sxbd(3,nbd), stangbd(3,nbd), sspbd(nbd)
    integer(8),   intent(in)    :: len1, len2, len3
    real(r64),    intent(in)    :: sxbd1(3,len1*nbd), sxbd2(3,len2*nbd), sxbd3(3,len3*nbd)
    real(r64),    intent(in)    :: stangbd1(3,len1*nbd), stangbd2(3,len2*nbd)
    real(r64),    intent(in)    :: stangbd3(3,len3*nbd)
    real(r64),    intent(in)    :: swbd1(len1*nbd), swbd2(len2*nbd), swbd3(len3*nbd)
    real(r64),    intent(in)    :: tgl(nquad), wgl(nquad), w_bclag(nquad)
    real(r64),    intent(in)    :: Dgl(nquad,nquad), bclagmatlr(nquad,2)
    complex(r64), intent(in)    :: troot(m,sbdnp)
    real(r64),    intent(in)    :: xroot(m,sbdnp), yroot(m,sbdnp), zroot(m,sbdnp)
    integer(8),   intent(inout) :: rfc(m,sbdnp)
    real(r64),    intent(inout) :: IalphaAsvestas(m)
    real(r64),    intent(in), optional :: rho_in
    real(r64),    intent(in), optional :: sxbd_raw(3,nbd), tx_raw(3,m)
    real(r64),    intent(in), optional :: Rfr(3,3), alpha_fr
    real(r64),    intent(in), optional :: Legmat(nquad,nquad)

    real(r64), parameter :: FAC = 3.0_r64, COEF = 1.0_r64
    integer(8) :: ell, j, idx_ell_start, idx_ell_end, k
    real(r64) :: r_ell(3,nquad), tau_ell(3,nquad)
    real(r64) :: x_ell(nquad), y_ell(nquad), z_ell(nquad), w_ell(nquad)
    real(r64) :: r0j(3), r_root(3), t_rootjr
    integer(8) :: rfcj
    real(r64) :: rmr0(3,nquad), rmr0norm(nquad), rmr0hat(3,nquad)
    real(r64) :: qhatxrmr0hat1(nquad), qhatxrmr0hat2(nquad), qhatxrmr0hat3(nquad)
    real(r64) :: numerator0(nquad), denominator0(nquad), integrand0(nquad)
    real(r64) :: rlr(3,2), rl(3), rr(3), rp(3,nquad), DglT(nquad,nquad)
    real(r64) :: pan_len, sqn_dist, dvec_min
    integer(8) :: len, lenl, lenr, nquad_up
    real(r64), allocatable :: r_up(:,:), rp_up(:,:), t_up(:), w_ref(:), Bup(:,:)
    real(r64) :: spk, t1, t2, t3, d1, d2, d3, dn, h1, h2, h3
    real(r64) :: c1, c2, c3, num, den, acc
    real(r64) :: brow(nquad)
    real(r64) :: r_ellL(3,nquad), rpL(3,nquad), rlL(3), rrL(3), rlrL(3,2)
    real(r64) :: lqs_t0, lqs_t1
    complex(r64) :: browc(nquad), dvecc(nquad), rc(3), rpc(3), Fc, Fpc, dtc, tr
    integer(8) :: ifconv, it, kmin
    logical    :: shifted
    real(r64) :: r_ell1(3,len1*nquad), w_ell1(len1*nquad), tau_ell1(3,len1*nquad)
    real(r64) :: r_ell2(3,len2*nquad), w_ell2(len2*nquad), tau_ell2(3,len2*nquad)
    real(r64) :: r_ell3(3,len3*nquad), w_ell3(len3*nquad), tau_ell3(3,len3*nquad)
    integer(8) :: tmp_vec1(len1*nquad), idx_ell1(len1*nquad), len1nquad
    integer(8) :: tmp_vec2(len2*nquad), idx_ell2(len2*nquad), len2nquad
    integer(8) :: tmp_vec3(len3*nquad), idx_ell3(len3*nquad), len3nquad
    real(r64) :: rmr01(3,len1*nquad), rmr0norm1(len1*nquad), rmr0hat1(3,len1*nquad)
    real(r64) :: qhatxrmr0hat11(len1*nquad), qhatxrmr0hat21(len1*nquad)
    real(r64) :: qhatxrmr0hat31(len1*nquad)
    real(r64) :: numerator01(len1*nquad), denominator01(len1*nquad), integrand01(len1*nquad)
    real(r64) :: rmr02(3,len2*nquad), rmr0norm2(len2*nquad), rmr0hat2(3,len2*nquad)
    real(r64) :: qhatxrmr0hat12(len2*nquad), qhatxrmr0hat22(len2*nquad)
    real(r64) :: qhatxrmr0hat32(len2*nquad)
    real(r64) :: numerator02(len2*nquad), denominator02(len2*nquad), integrand02(len2*nquad)
    real(r64) :: rmr03(3,len3*nquad), rmr0norm3(len3*nquad), rmr0hat3(3,len3*nquad)
    real(r64) :: qhatxrmr0hat13(len3*nquad), qhatxrmr0hat23(len3*nquad)
    real(r64) :: qhatxrmr0hat33(len3*nquad)
    real(r64) :: numerator03(len3*nquad), denominator03(len3*nquad), integrand03(len3*nquad)

    DglT = transpose(Dgl)
    if (present(rho_in)) then
      call update_refinement_codes_r64(m, sbdnp, nquad, [len1, len2, len3], troot, rfc, rho_in)
    else
      call update_refinement_codes_r64(m, sbdnp, nquad, [len1, len2, len3], troot, rfc)
    end if
    len1nquad = len1*nquad
    tmp_vec1 = [(k, k = 1_8, len1nquad, 1_8)]
    len2nquad = len2*nquad
    tmp_vec2 = [(k, k = 1_8, len2nquad, 1_8)]
    len3nquad = len3*nquad
    tmp_vec3 = [(k, k = 1_8, len3nquad, 1_8)]
    allocate(t_up(198_8*nquad), w_ref(198_8*nquad))
    allocate(r_up(3,198_8*nquad), rp_up(3,198_8*nquad), Bup(nquad,198_8*nquad))

    do ell = 1, sbdnp
      idx_ell_start = (ell-1_8)*nquad + 1_8
      idx_ell_end   = ell*nquad
      r_ell   = sxbd(:, idx_ell_start:idx_ell_end)
      tau_ell = stangbd(:, idx_ell_start:idx_ell_end)
      w_ell   = sspbd(idx_ell_start:idx_ell_end)*wgl
      x_ell = r_ell(1,:)
      y_ell = r_ell(2,:)
      z_ell = r_ell(3,:)

      idx_ell1 = (ell-1_8)*len1nquad + tmp_vec1
      r_ell1   = sxbd1(:, idx_ell1)
      tau_ell1 = stangbd1(:, idx_ell1)
      w_ell1   = swbd1(idx_ell1)

      idx_ell2 = (ell-1_8)*len2nquad + tmp_vec2
      r_ell2   = sxbd2(:, idx_ell2)
      tau_ell2 = stangbd2(:, idx_ell2)
      w_ell2   = swbd2(idx_ell2)

      idx_ell3 = (ell-1_8)*len3nquad + tmp_vec3
      r_ell3   = sxbd3(:, idx_ell3)
      tau_ell3 = stangbd3(:, idx_ell3)
      w_ell3   = swbd3(idx_ell3)

      rp  = matmul(r_ell, DglT)
      rlr = matmul(r_ell, bclagmatlr)
      rl  = rlr(:,1)
      rr  = rlr(:,2)
      pan_len = sum(w_ell)

      do j = 1, m
        r0j = r0(:, j)
        t_rootjr  = real(troot(j,ell), r64)
        r_root(1) = xroot(j,ell)
        r_root(2) = yroot(j,ell)
        r_root(3) = zroot(j,ell)
        rfcj = rfc(j,ell)
        if (LQS_PROF == 1_8) then
          lqs_rfc_hist(min(rfcj, 8_8)) = lqs_rfc_hist(min(rfcj, 8_8)) + 1_8
          call cpu_time(lqs_t0)
        end if

        if (rfcj >= 4_8) then
          r_ellL = r_ell;  rpL = rp;  rlL = rl;  rrL = rr
          shifted = .false.
          if (present(sxbd_raw) .and. present(tx_raw) .and. present(Rfr) .and. &
              present(alpha_fr) .and. present(Legmat)) then
            do k = 1, nquad
              d1 = sxbd_raw(1, idx_ell_start+k-1_8) - tx_raw(1,j)
              d2 = sxbd_raw(2, idx_ell_start+k-1_8) - tx_raw(2,j)
              d3 = sxbd_raw(3, idx_ell_start+k-1_8) - tx_raw(3,j)
              r_ellL(1,k) = alpha_fr*(Rfr(1,1)*d1 + Rfr(1,2)*d2 + Rfr(1,3)*d3)
              r_ellL(2,k) = alpha_fr*(Rfr(2,1)*d1 + Rfr(2,2)*d2 + Rfr(2,3)*d3)
              r_ellL(3,k) = alpha_fr*(Rfr(3,1)*d1 + Rfr(3,2)*d2 + Rfr(3,3)*d3)
            end do
            rpL = matmul(r_ellL, DglT)
            tr = troot(j,ell)
            ifconv = 0_8
            do it = 1, 8
              dvecc = tr - tgl
              kmin = 1_8
              dvec_min = abs(dvecc(1))
              do k = 2, nquad
                if (abs(dvecc(k)) < dvec_min) then
                  kmin = k
                  dvec_min = abs(dvecc(k))
                end if
              end do
              if (abs(dvecc(kmin)) < 1.0e-14_r64) then
                rc  = r_ellL(:,kmin)
                rpc = rpL(:,kmin)
              else
                browc = w_bclag/dvecc
                browc = browc/sum(browc)
                rc(1)  = sum(r_ellL(1,:)*browc)
                rc(2)  = sum(r_ellL(2,:)*browc)
                rc(3)  = sum(r_ellL(3,:)*browc)
                rpc(1) = sum(rpL(1,:)*browc)
                rpc(2) = sum(rpL(2,:)*browc)
                rpc(3) = sum(rpL(3,:)*browc)
              end if
              Fc  = rc(1)*rc(1) + rc(2)*rc(2) + rc(3)*rc(3)
              Fpc = 2.0_r64*(rc(1)*rpc(1) + rc(2)*rpc(2) + rc(3)*rpc(3))
              dtc = -Fc/Fpc
              tr  = tr + dtc
              if (abs(dtc) < 1.0e-15_r64) then
                ifconv = 1_8
                exit
              end if
            end do
            if (ifconv == 1_8) then
              rlrL = matmul(r_ellL, bclagmatlr)
              rlL  = rlrL(:,1);  rrL = rlrL(:,2)
              t_rootjr = real(tr, r64)
              call bary_row_r64(nquad, tgl, w_bclag, t_rootjr, brow)
              r_root(1) = dot_product(r_ellL(1,:), brow)
              r_root(2) = dot_product(r_ellL(2,:), brow)
              r_root(3) = dot_product(r_ellL(3,:), brow)
              r0j = 0.0_r64
              shifted = .true.
            end if
          end if
          if (.not. shifted) then
            r_ellL = r_ell;  rpL = rp;  rlL = rl;  rrL = rr
          end if

          if (t_rootjr >= 1.0_r64) then
            sqn_dist = sqrt(dot_product(r0j - rrL, r0j - rrL))
          else if (t_rootjr <= -1.0_r64) then
            sqn_dist = sqrt(dot_product(r0j - rlL, r0j - rlL))
          else
            sqn_dist = sqrt(dot_product(r0j - r_root, r0j - r_root))
          end if
          call estimate_nearroot_lengths_r64(t_rootjr, &
                 2.0_r64*sqn_dist/pan_len, 99_8, len, lenl, lenr, FAC, COEF)
          nquad_up = nquad*(len - 1_8)
          if (LQS_PROF == 1_8) then
            lqs_len_hist(min(len, 64_8)) = lqs_len_hist(min(len, 64_8)) + 1_8
            lqs_nodes_adap = lqs_nodes_adap + nquad_up
          end if
          call build_nearroot_nodes_r64(t_rootjr, nquad, tgl, wgl, len, lenl, lenr, &
                                        t_up(1:nquad_up), w_ref(1:nquad_up))
          do k = 1, nquad_up
            call bary_row_r64(nquad, tgl, w_bclag, t_up(k), Bup(:,k))
          end do
          r_up(:,1:nquad_up)  = matmul(r_ellL, Bup(:,1:nquad_up))
          rp_up(:,1:nquad_up) = matmul(rpL, Bup(:,1:nquad_up))
          acc = 0.0_r64
          do k = 1, nquad_up
            spk = sqrt(rp_up(1,k)**2 + rp_up(2,k)**2 + rp_up(3,k)**2)
            t1  = rp_up(1,k)/spk
            t2  = rp_up(2,k)/spk
            t3  = rp_up(3,k)/spk
            d1  = r_up(1,k) - r0j(1)
            d2  = r_up(2,k) - r0j(2)
            d3  = r_up(3,k) - r0j(3)
            dn  = sqrt(d1*d1 + d2*d2 + d3*d3)
            h1  = d1/dn
            h2  = d2/dn
            h3  = d3/dn
            c1  = qhat(2)*h3 - qhat(3)*h2
            c2  = qhat(3)*h1 - qhat(1)*h3
            c3  = qhat(1)*h2 - qhat(2)*h1
            num = t1*c1 + t2*c2 + t3*c3
            den = dn*(1.0_r64 - (qhat(1)*h1 + qhat(2)*h2 + qhat(3)*h3))
            acc = acc - (num/den)*(spk*w_ref(k))
          end do
          IalphaAsvestas(j) = IalphaAsvestas(j) + acc
        else if (rfcj >= 3_8) then
          rmr03(1,:) = r_ell3(1,:) - r0j(1)
          rmr03(2,:) = r_ell3(2,:) - r0j(2)
          rmr03(3,:) = r_ell3(3,:) - r0j(3)
          rmr0norm3 = sqrt(rmr03(1,:)**2 + rmr03(2,:)**2 + rmr03(3,:)**2)
          rmr0hat3(1,:) = rmr03(1,:)/rmr0norm3
          rmr0hat3(2,:) = rmr03(2,:)/rmr0norm3
          rmr0hat3(3,:) = rmr03(3,:)/rmr0norm3
          qhatxrmr0hat13 = qhat(2)*rmr0hat3(3,:) - qhat(3)*rmr0hat3(2,:)
          qhatxrmr0hat23 = qhat(3)*rmr0hat3(1,:) - qhat(1)*rmr0hat3(3,:)
          qhatxrmr0hat33 = qhat(1)*rmr0hat3(2,:) - qhat(2)*rmr0hat3(1,:)
          numerator03 = tau_ell3(1,:)*qhatxrmr0hat13 &
                      + tau_ell3(2,:)*qhatxrmr0hat23 &
                      + tau_ell3(3,:)*qhatxrmr0hat33
          denominator03 = rmr0norm3*(1.0_r64 - (qhat(1)*rmr0hat3(1,:) &
                        + qhat(2)*rmr0hat3(2,:) + qhat(3)*rmr0hat3(3,:)))
          integrand03 = -numerator03/denominator03
          IalphaAsvestas(j) = IalphaAsvestas(j) + sum(integrand03*w_ell3)
        else if (rfcj >= 2_8) then
          rmr02(1,:) = r_ell2(1,:) - r0j(1)
          rmr02(2,:) = r_ell2(2,:) - r0j(2)
          rmr02(3,:) = r_ell2(3,:) - r0j(3)
          rmr0norm2 = sqrt(rmr02(1,:)**2 + rmr02(2,:)**2 + rmr02(3,:)**2)
          rmr0hat2(1,:) = rmr02(1,:)/rmr0norm2
          rmr0hat2(2,:) = rmr02(2,:)/rmr0norm2
          rmr0hat2(3,:) = rmr02(3,:)/rmr0norm2
          qhatxrmr0hat12 = qhat(2)*rmr0hat2(3,:) - qhat(3)*rmr0hat2(2,:)
          qhatxrmr0hat22 = qhat(3)*rmr0hat2(1,:) - qhat(1)*rmr0hat2(3,:)
          qhatxrmr0hat32 = qhat(1)*rmr0hat2(2,:) - qhat(2)*rmr0hat2(1,:)
          numerator02 = tau_ell2(1,:)*qhatxrmr0hat12 &
                      + tau_ell2(2,:)*qhatxrmr0hat22 &
                      + tau_ell2(3,:)*qhatxrmr0hat32
          denominator02 = rmr0norm2*(1.0_r64 - (qhat(1)*rmr0hat2(1,:) &
                        + qhat(2)*rmr0hat2(2,:) + qhat(3)*rmr0hat2(3,:)))
          integrand02 = -numerator02/denominator02
          IalphaAsvestas(j) = IalphaAsvestas(j) + sum(integrand02*w_ell2)
        else if (rfcj >= 1_8) then
          rmr01(1,:) = r_ell1(1,:) - r0j(1)
          rmr01(2,:) = r_ell1(2,:) - r0j(2)
          rmr01(3,:) = r_ell1(3,:) - r0j(3)
          rmr0norm1 = sqrt(rmr01(1,:)**2 + rmr01(2,:)**2 + rmr01(3,:)**2)
          rmr0hat1(1,:) = rmr01(1,:)/rmr0norm1
          rmr0hat1(2,:) = rmr01(2,:)/rmr0norm1
          rmr0hat1(3,:) = rmr01(3,:)/rmr0norm1
          qhatxrmr0hat11 = qhat(2)*rmr0hat1(3,:) - qhat(3)*rmr0hat1(2,:)
          qhatxrmr0hat21 = qhat(3)*rmr0hat1(1,:) - qhat(1)*rmr0hat1(3,:)
          qhatxrmr0hat31 = qhat(1)*rmr0hat1(2,:) - qhat(2)*rmr0hat1(1,:)
          numerator01 = tau_ell1(1,:)*qhatxrmr0hat11 &
                      + tau_ell1(2,:)*qhatxrmr0hat21 &
                      + tau_ell1(3,:)*qhatxrmr0hat31
          denominator01 = rmr0norm1*(1.0_r64 - (qhat(1)*rmr0hat1(1,:) &
                        + qhat(2)*rmr0hat1(2,:) + qhat(3)*rmr0hat1(3,:)))
          integrand01 = -numerator01/denominator01
          IalphaAsvestas(j) = IalphaAsvestas(j) + sum(integrand01*w_ell1)
        else
          rmr0(1,:) = x_ell - r0j(1)
          rmr0(2,:) = y_ell - r0j(2)
          rmr0(3,:) = z_ell - r0j(3)
          rmr0norm = sqrt(rmr0(1,:)**2 + rmr0(2,:)**2 + rmr0(3,:)**2)
          rmr0hat(1,:) = rmr0(1,:)/rmr0norm
          rmr0hat(2,:) = rmr0(2,:)/rmr0norm
          rmr0hat(3,:) = rmr0(3,:)/rmr0norm
          qhatxrmr0hat1 = qhat(2)*rmr0hat(3,:) - qhat(3)*rmr0hat(2,:)
          qhatxrmr0hat2 = qhat(3)*rmr0hat(1,:) - qhat(1)*rmr0hat(3,:)
          qhatxrmr0hat3 = qhat(1)*rmr0hat(2,:) - qhat(2)*rmr0hat(1,:)
          numerator0 = tau_ell(1,:)*qhatxrmr0hat1 &
                     + tau_ell(2,:)*qhatxrmr0hat2 &
                     + tau_ell(3,:)*qhatxrmr0hat3
          denominator0 = rmr0norm*(1.0_r64 - (qhat(1)*rmr0hat(1,:) &
                       + qhat(2)*rmr0hat(2,:) + qhat(3)*rmr0hat(3,:)))
          integrand0 = -numerator0/denominator0
          IalphaAsvestas(j) = IalphaAsvestas(j) + sum(integrand0*w_ell)
        end if
        if (LQS_PROF == 1_8) then
          call cpu_time(lqs_t1)
          if (rfcj >= 4_8) then
            lqs_t_adap = lqs_t_adap + (lqs_t1 - lqs_t0)
          else
            lqs_t_unif = lqs_t_unif + (lqs_t1 - lqs_t0)
          end if
        end if
      end do
    end do
    deallocate(t_up, w_ref, r_up, rp_up, Bup)
  end subroutine evaluate_solid_angle_integral_fast_r64


  subroutine evaluate_solid_angle_integral_fast_driver_r64(m, tx, n, sx, snx, sw, &
                                                           r_vert, nbd, sxbd_in, &
                                                           IalphaAsvestas)
    use koorn_geom_mod, only: circumcircle_transform_3d
    use lq_kernel_mod,  only: gauss_r64, bclaginterpweights_r64, bary_row_r64
    use lq_adaptive_mod, only: line_quad_root_initial_guess_r64, &
                               line_quad_root_refine_r64
    integer(8), intent(in)    :: m, n, nbd
    real(r64),  intent(in)    :: tx(3,m), sx(3,n), snx(3,n), sw(n)
    real(r64),  intent(in)    :: r_vert(3,3), sxbd_in(3,nbd)
    real(r64),  intent(inout) :: IalphaAsvestas(m)

    integer(8), parameter :: sbdnp = 3_8
    integer(8), parameter :: len1 = 4_8, len2 = 8_8, len3 = 16_8
    integer(8) :: nquad

    real(r64), allocatable :: tgl(:), wgl(:), Dgl(:,:), w_bclag(:)
    real(r64), allocatable :: Legmat(:,:), vtmp(:,:), bclagmatlr(:,:)
    real(r64), allocatable :: sxbd(:,:), stangbd(:,:), sspbd(:)
    real(r64), allocatable :: sxbd1(:,:), stangbd1(:,:), swbd1(:)
    real(r64), allocatable :: sxbd2(:,:), stangbd2(:,:), swbd2(:)
    real(r64), allocatable :: sxbd3(:,:), stangbd3(:,:), swbd3(:)
    real(r64), allocatable :: txn(:,:), xroot(:,:), yroot(:,:), zroot(:,:)
    real(r64), allocatable :: brow(:), xh(:), yh(:), zh(:)
    complex(r64), allocatable :: troot(:,:)
    integer(8),   allocatable :: rfc(:,:)

    complex(r64) :: tinit
    real(r64) :: R(3,3), cc(3), alpha, qhat(3), sxp(3)
    integer(8) :: ell, j, ii, i0, ifconv

    nquad = nbd/sbdnp

    allocate(tgl(nquad), wgl(nquad), Dgl(nquad,nquad), w_bclag(nquad))
    allocate(Legmat(nquad,nquad), vtmp(nquad,nquad), bclagmatlr(nquad,2))
    allocate(sxbd(3,nbd), stangbd(3,nbd), sspbd(nbd))
    allocate(sxbd1(3,len1*nbd), stangbd1(3,len1*nbd), swbd1(len1*nbd))
    allocate(sxbd2(3,len2*nbd), stangbd2(3,len2*nbd), swbd2(len2*nbd))
    allocate(sxbd3(3,len3*nbd), stangbd3(3,len3*nbd), swbd3(len3*nbd))
    allocate(txn(3,m), troot(m,nbd), rfc(m,nbd))
    allocate(xroot(m,nbd), yroot(m,nbd), zroot(m,nbd))
    allocate(brow(nquad), xh(nquad), yh(nquad), zh(nquad))

    call gauss_r64(nquad, tgl, wgl, Dgl)
    call bclaginterpweights_r64(nquad, tgl, w_bclag)
    call legeexps_r64(2_8, nquad, tgl, Legmat, vtmp, wgl)
    call endpoint_bary(nquad, tgl, w_bclag, bclagmatlr)

    sxbd = sxbd_in
    do ell = 1, sbdnp
      i0 = (ell-1_8)*nquad
      do ii = 1, nquad
        sxp(1) = sum(Dgl(ii,:)*sxbd(1, i0+1:i0+nquad))
        sxp(2) = sum(Dgl(ii,:)*sxbd(2, i0+1:i0+nquad))
        sxp(3) = sum(Dgl(ii,:)*sxbd(3, i0+1:i0+nquad))
        sspbd(i0+ii) = sqrt(sxp(1)**2 + sxp(2)**2 + sxp(3)**2)
        stangbd(:, i0+ii) = sxp/sspbd(i0+ii)
      end do
    end do

    call circumcircle_transform_3d(r_vert, R, cc, alpha)
    do j = 1, m
      txn(:,j) = alpha*matmul(R, tx(:,j) - cc)
    end do
    do ii = 1, nbd
      sxbd(:,ii)    = alpha*matmul(R, sxbd(:,ii) - cc)
      stangbd(:,ii) = matmul(R, stangbd(:,ii))
    end do
    sspbd = alpha*sspbd
    qhat = 0.0_r64
    do j = 1, n
      qhat = qhat + matmul(R, snx(:,j))
    end do
    qhat = qhat/sqrt(dot_product(qhat, qhat))

    call sub_level(len1, sxbd1, stangbd1, swbd1, len1*nbd)
    call sub_level(len2, sxbd2, stangbd2, swbd2, len2*nbd)
    call sub_level(len3, sxbd3, stangbd3, swbd3, len3*nbd)

    troot = cmplx(0.0_r64, 0.0_r64, kind=r64)
    xroot = 0.0_r64;  yroot = 0.0_r64;  zroot = 0.0_r64
    rfc = 1_8
    do ell = 1, sbdnp
      i0 = (ell-1_8)*nquad
      xh = matmul(Legmat, sxbd(1, i0+1:i0+nquad))
      yh = matmul(Legmat, sxbd(2, i0+1:i0+nquad))
      zh = matmul(Legmat, sxbd(3, i0+1:i0+nquad))
      do j = 1, m
        call line_quad_root_initial_guess_r64(tgl, sxbd(1,i0+1:i0+nquad), &
             sxbd(2,i0+1:i0+nquad), sxbd(3,i0+1:i0+nquad), nquad, &
             txn(1,j), txn(2,j), txn(3,j), tinit)
        ifconv = 0_8
        call line_quad_root_refine_r64(xh, yh, zh, nquad, &
             txn(1,j), txn(2,j), txn(3,j), tinit, troot(j,ell), ifconv)
        if (ifconv /= 1_8) then
          troot(j,ell) = cmplx(0.0_r64, 1.0e6_r64, kind=r64)
          rfc(j,ell) = 0_8
        end if
        call bary_row_r64(nquad, tgl, w_bclag, real(troot(j,ell), r64), brow)
        xroot(j,ell) = dot_product(sxbd(1, i0+1:i0+nquad), brow)
        yroot(j,ell) = dot_product(sxbd(2, i0+1:i0+nquad), brow)
        zroot(j,ell) = dot_product(sxbd(3, i0+1:i0+nquad), brow)
      end do
    end do

    call evaluate_solid_angle_integral_fast_r64(m, txn, nbd, sbdnp, nquad, &
         sxbd, stangbd, sspbd, &
         len1, sxbd1, stangbd1, swbd1, &
         len2, sxbd2, stangbd2, swbd2, &
         len3, sxbd3, stangbd3, swbd3, &
         qhat, tgl, wgl, Dgl, w_bclag, bclagmatlr, &
         troot, xroot, yroot, zroot, rfc, IalphaAsvestas, &
         rho_in=8.0_r64**(8.0_r64/real(nquad, r64)), &
         sxbd_raw=sxbd_in, tx_raw=tx, Rfr=R, alpha_fr=alpha, Legmat=Legmat)

  contains

    subroutine endpoint_bary(nq, t, wb, mat)
      integer(8), intent(in)  :: nq
      real(r64),  intent(in)  :: t(nq), wb(nq)
      real(r64),  intent(out) :: mat(nq,2)
      real(r64) :: d(nq)
      integer(8) :: q
      do q = 1, 2
        d = wb/((-1.0_r64 + 2.0_r64*real(q-1_8, r64)) - t)
        mat(:,q) = d/sum(d)
      end do
    end subroutine endpoint_bary

    subroutine sub_level(nlev, sxs, sts, sws, ntot)
      integer(8), intent(in)    :: nlev, ntot
      real(r64),  intent(inout) :: sxs(3,ntot), sts(3,ntot), sws(ntot)
      real(r64) :: hw, mc, tt, row(nquad), rp(3), sp, dx(nquad,3)
      integer(8) :: e2, p, q, idx, b0
      hw = 1.0_r64/real(nlev, r64)
      do e2 = 1, sbdnp
        b0 = (e2-1_8)*nquad
        dx(:,1) = matmul(Dgl, sxbd(1, b0+1:b0+nquad))
        dx(:,2) = matmul(Dgl, sxbd(2, b0+1:b0+nquad))
        dx(:,3) = matmul(Dgl, sxbd(3, b0+1:b0+nquad))
        do p = 1, nlev
          mc = -1.0_r64 + (2.0_r64*real(p, r64) - 1.0_r64)*hw
          do q = 1, nquad
            idx = (e2-1_8)*nlev*nquad + (p-1_8)*nquad + q
            tt  = mc + tgl(q)*hw
            call bary_row_r64(nquad, tgl, w_bclag, tt, row)
            sxs(1,idx) = dot_product(sxbd(1, b0+1:b0+nquad), row)
            sxs(2,idx) = dot_product(sxbd(2, b0+1:b0+nquad), row)
            sxs(3,idx) = dot_product(sxbd(3, b0+1:b0+nquad), row)
            rp(1) = dot_product(dx(:,1), row)
            rp(2) = dot_product(dx(:,2), row)
            rp(3) = dot_product(dx(:,3), row)
            sp = sqrt(dot_product(rp, rp))
            sts(:,idx) = rp/sp
            sws(idx)   = wgl(q)*hw*sp
          end do
        end do
      end do
    end subroutine sub_level

  end subroutine evaluate_solid_angle_integral_fast_driver_r64

end module solidangle_mod
