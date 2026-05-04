module solidangle_mod
  use linequaaadrature_mod
  implicit none

contains

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
                                     build_target_nearroot_weights_r64, KERNEL_ASVESTAS
    use lq_adaptive_mod, only: line_quad_root_initial_guess_r64, line_quad_root_refine_r64
    use iso_c_binding,         only: c_funptr, c_funloc
    integer(8),  intent(in)    :: m, n, nbd
    real(r64),   intent(in)    :: tx(3,m), sx(3,n), snx(3,n), sw(n)
    real(r64),   intent(in)    :: r_vert(3,3), sxbd_in(3,nbd)
    logical,     intent(in)    :: use_nearroot
    real(r64),   intent(inout) :: IalphaAsvestas(m)

    integer(8) :: nquad, sbdnp
    integer(8) :: j, k, ell, idx_start, idx_end
    integer(8) :: n_expa, ifconv
    real(r64)  :: alpha, qhat(3), qnrm, R(3,3), c(3)
    real(r64)  :: templ, tempr, denoml, denomr, tgll, tglr
    real(r64)  :: sxp(3), rho, br, I_local
    complex(8) :: tinit, troot
    type(c_funptr) :: cfptr

    real(r64), allocatable :: tgl(:), wgl(:), Dgl(:,:)
    real(r64), allocatable :: w_bclag(:)
    real(r64), allocatable :: Legmat(:,:), vtmp(:,:)
    real(r64), allocatable :: txnew(:,:), snxnew(:,:), sxpbd(:,:)
    real(r64), allocatable :: sxbd(:,:), stangbd(:,:), sspbd(:)
    real(r64), allocatable :: kdata(:,:)
    real(r64), allocatable :: funvals(:,:,:), sxbdw(:,:,:)
    real(r64), allocatable :: bclagmatlr(:,:)
    real(r64), allocatable :: root_re(:), root_im(:), xyz_hat(:,:)
    real(r64), allocatable :: funvals_ell(:,:,:), sxbdw_ell(:,:,:)
    real(r64), allocatable :: funvals0(:), weights(:)
    logical,   allocatable :: root_ok(:)
    logical :: accepted

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

    ! --- apply transform ---
    do j = 1, m
      txnew(:,j) = alpha * matmul(R, tx(:,j) - c)
    end do
    do j = 1, n
      snxnew(:,j) = matmul(R, snx(:,j))
    end do
    do k = 1, nbd
      sxbd(:,k)    = alpha * matmul(R, sxbd(:,k) - c)
      stangbd(:,k) = alpha * matmul(R, stangbd(:,k))
    end do

    ! --- qhat = mean(snxnew) / norm (iside=0) ---
    qhat = 0.0_r64
    do j = 1, n
      qhat = qhat + snxnew(:,j)
    end do
    qhat = qhat / real(n, r64)
    qnrm = sqrt(qhat(1)**2 + qhat(2)**2 + qhat(3)**2)
    qhat = qhat / qnrm

    ! --- sxpbd = Dgl * sxbd per panel ---
    sxpbd = 0.0_r64
    do ell = 1, sbdnp
      idx_start = (ell-1)*nquad + 1
      idx_end   = ell*nquad
      sxpbd(1,idx_start:idx_end) = matmul(Dgl, sxbd(1,idx_start:idx_end))
      sxpbd(2,idx_start:idx_end) = matmul(Dgl, sxbd(2,idx_start:idx_end))
      sxpbd(3,idx_start:idx_end) = matmul(Dgl, sxbd(3,idx_start:idx_end))
    end do

    ! --- kdata, cfptr ---
    do j = 1, m
      kdata(:,j) = qhat
    end do
    cfptr = c_funloc(asvestas_kernel_r64)

    if (use_nearroot) then
      ! --- nearroot path: per-panel root finding + compress ---
      n_expa = min(16_8, nquad)
      rho    = 4.0_r64**(16.0_r64 / real(nquad, r64))
      allocate(xyz_hat(n_expa, 3))
      allocate(funvals0(nquad), weights(nquad))
      funvals = 0.0_r64
      sxbdw   = 0.0_r64

      do ell = 1, sbdnp
        idx_start = (ell-1)*nquad + 1
        idx_end   = ell*nquad

        ! Legendre projection of panel coords (first n_expa modes)
        xyz_hat(:,1) = matmul(Legmat(1:n_expa,:), sxbd(1, idx_start:idx_end))
        xyz_hat(:,2) = matmul(Legmat(1:n_expa,:), sxbd(2, idx_start:idx_end))
        xyz_hat(:,3) = matmul(Legmat(1:n_expa,:), sxbd(3, idx_start:idx_end))

        do j = 1, m
          funvals0 = 0.0_r64
          weights  = 0.0_r64
          troot = cmplx(0.0_r64, 0.0_r64, kind=r64)
          accepted = .false.
          I_local = 0.0_r64
          ! n_expa_in=n_expa truncates the rootfinder polynomial to the
          ! first n_expa = min(16, nquad) Legendre coefficients of the
          ! curve, mirroring the un-wrapped path (and the r128 sibling).
          call build_target_nearroot_weights_r64(nquad, tgl, wgl, Legmat, &
                                                 sxbd(1,idx_start:idx_end), &
                                                 sxbd(2,idx_start:idx_end), &
                                                 sxbd(3,idx_start:idx_end), &
                                                 sspbd(idx_start:idx_end), &
                                                 stangbd(:,idx_start:idx_end), &
                                                 xyz_hat(:,1), xyz_hat(:,2), xyz_hat(:,3), &
                                                 rho, txnew(1,j), txnew(2,j), txnew(3,j), &
                                                 cfptr, kdata(:,j), &
                                                 funvals0, weights, &
                                                 troot, accepted, I_local, &
                                                 kernel_id=KERNEL_ASVESTAS, &
                                                 dgl_in=Dgl, w_bclag_in=w_bclag, &
                                                 sxpbd_in=sxpbd(:,idx_start:idx_end), &
                                                 n_expa_in=n_expa, &
                                                 adaptive_fallback=.true.)
          funvals(:,ell,j) = funvals0
          sxbdw(:,ell,j)   = weights
        end do
      end do

      deallocate(xyz_hat, funvals0, weights)

    else
      ! --- adaptive bisection path (all panels at once) ---
      funvals = 0.0_r64
      call line_kernel_eval_r64(m, txnew, nbd, sbdnp, nquad, &
                                 sxbd, sxpbd, stangbd, cfptr, kdata, funvals)
      sxbdw = 0.0_r64
      call line_quad_compress_r64(m, txnew, nbd, sbdnp, nquad,  &
                                   sxbd, sxpbd, stangbd, sspbd,  &
                                   tgl, wgl, Dgl, w_bclag,       &
                                   Legmat, bclagmatlr,            &
                                   cfptr, kdata, funvals, sxbdw)
    end if

    do j = 1, m
      IalphaAsvestas(j) = sum(funvals(:,:,j) * sxbdw(:,:,j))
    end do

    deallocate(tgl, wgl, Dgl, w_bclag, Legmat, vtmp)
    deallocate(sxbd, stangbd, sspbd, txnew, snxnew, sxpbd)
    deallocate(kdata, funvals, sxbdw, bclagmatlr)

  end subroutine evaluate_solid_angle_integral_r64

  subroutine evaluate_solid_angle_integral_r128(m, tx, n, sx, snx, sw, r_vert, nbd, sxbd_in, use_nearroot, IalphaAsvestas)
    use koorn_geom_mod, only: circumcircle_transform_3d_r128
    use lq_kernel_mod,  only: line_kernel_eval_r128, line_quad_compress_r128, &
                              line_quad_compress_nearroot_r128
    use lq_adaptive_mod, only: line_quad_root_initial_guess_r128, line_quad_root_refine_r128
    integer(8),  intent(in)    :: m, n, nbd
    real(r128),  intent(in)    :: tx(3,m), sx(3,n), snx(3,n), sw(n)
    real(r128),  intent(in)    :: r_vert(3,3), sxbd_in(3,nbd)
    logical,     intent(in)    :: use_nearroot
    real(r128),  intent(inout) :: IalphaAsvestas(m)

    integer(8) :: nquad, sbdnp
    integer(8) :: j, k, ell, idx_start, idx_end
    integer(8) :: n_expa, ifconv
    real(r128) :: alpha, qhat(3), qnrm, R(3,3), c(3)
    real(r128) :: templ, tempr, denoml, denomr, tgll, tglr
    real(r128) :: sxp(3)
    real(r128) :: rho, br
    complex(16) :: tinit, troot

    real(r128), allocatable :: tgl(:), wgl(:), Dgl(:,:)
    real(r128), allocatable :: w_bclag(:)
    real(r128), allocatable :: Legmat(:,:), vtmp(:,:)
    real(r128), allocatable :: txnew(:,:), snxnew(:,:), sxpbd(:,:)
    real(r128), allocatable :: sxbd(:,:), stangbd(:,:), sspbd(:)
    real(r128), allocatable :: kdata(:,:)
    real(r128), allocatable :: funvals(:,:,:), sxbdw(:,:,:)
    real(r128), allocatable :: bclagmatlr(:,:)
    real(r128), allocatable :: root_re(:), root_im(:), xyz_hat(:,:)
    real(r128), allocatable :: funvals_ell(:,:,:), sxbdw_ell(:,:,:)
    logical,    allocatable :: root_ok(:)

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
      stangbd(:,k) = alpha * matmul(R, stangbd(:,k))
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
      ! --- nearroot path: per-panel root finding + compress ---
      n_expa = min(16_8, nquad)
      rho    = 4.0_r128**(16.0_r128 / real(nquad, r128))
      allocate(root_re(m), root_im(m), root_ok(m))
      allocate(xyz_hat(n_expa, 3))
      allocate(funvals_ell(nquad, 1, m), sxbdw_ell(nquad, 1, m))
      funvals = 0.0_r128
      sxbdw   = 0.0_r128

      do ell = 1, sbdnp
        idx_start = (ell-1)*nquad + 1
        idx_end   = ell*nquad

        ! Legendre projection of panel coords (first n_expa modes)
        xyz_hat(:,1) = matmul(Legmat(1:n_expa,:), sxbd(1, idx_start:idx_end))
        xyz_hat(:,2) = matmul(Legmat(1:n_expa,:), sxbd(2, idx_start:idx_end))
        xyz_hat(:,3) = matmul(Legmat(1:n_expa,:), sxbd(3, idx_start:idx_end))

        ! Per-target root finding
        do j = 1, m
          call line_quad_root_initial_guess_r128(tgl, sxbd(1,idx_start:idx_end), &
                                          sxbd(2,idx_start:idx_end), &
                                          sxbd(3,idx_start:idx_end), &
                                          nquad, txnew(1,j), txnew(2,j), txnew(3,j), tinit)
          root_ok(j) = .false.
          root_re(j) = 0.0_r128
          root_im(j) = 0.0_r128
          br = abs(tinit + sqrt(tinit - 1.0_r128)*sqrt(tinit + 1.0_r128))
          if (br < 1.75_r128*rho) then
            call line_quad_root_refine_r128(xyz_hat(:,1), xyz_hat(:,2), xyz_hat(:,3), &
                                     n_expa, txnew(1,j), txnew(2,j), txnew(3,j), &
                                     tinit, troot, ifconv)
            br = abs(troot + sqrt(troot - 1.0_r128)*sqrt(troot + 1.0_r128))
            if (ifconv == 1_8 .and. br < rho) then
              root_ok(j) = .true.
              root_re(j) = real(troot, r128)
              root_im(j) = aimag(troot)
            end if
          end if
        end do

        ! Kernel eval for this panel (sbdnp=1)
        funvals_ell = 0.0_r128
        call line_kernel_eval_r128(m, txnew, nquad, 1_8, nquad, &
                                    sxbd(:, idx_start:idx_end), &
                                    sxpbd(:, idx_start:idx_end), &
                                    stangbd(:, idx_start:idx_end), &
                                    asvestas_kernel_r128, kdata, funvals_ell)

        ! Nearroot compress for this panel (sbdnp=1)
        sxbdw_ell = 0.0_r128
        call line_quad_compress_nearroot_r128(m, txnew, nquad, 1_8, nquad, &
                                               sxbd(:, idx_start:idx_end), &
                                               sxpbd(:, idx_start:idx_end), &
                                               stangbd(:, idx_start:idx_end), &
                                               sspbd(idx_start:idx_end), &
                                               tgl, wgl, Dgl, w_bclag, &
                                               Legmat, bclagmatlr, &
                                               asvestas_kernel_r128, kdata, funvals_ell, sxbdw_ell, &
                                               root_re, root_im, root_ok)
        funvals(:, ell, :) = funvals_ell(:, 1, :)
        sxbdw(:, ell, :)   = sxbdw_ell(:, 1, :)
      end do

      deallocate(root_re, root_im, root_ok, xyz_hat, funvals_ell, sxbdw_ell)

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
    use lq_kernel_mod, only: kernel_iface_r128, line_kernel_eval_r128, line_quad_compress_nearroot_r128
    use lq_adaptive_mod, only: line_quad_root_initial_guess_r128, line_quad_root_refine_r128
    integer(8),     intent(in)     :: m, nbd, sbdnp, nquad
    real(r128),      intent(in)    :: r0(3,m)
    real(r128),      intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
    real(r128),      intent(in)    :: sspbd(nbd)
    procedure(kernel_iface_r128)   :: fun
    real(r128),  intent(in)        :: kdata(3,m)
    real(r128),  intent(inout)     :: q_lq128(m)

    integer(8) :: ell, j, k, idx_start, idx_end, n_expa, ifconv
    real(r128) :: rho, br
    complex(16) :: tinit, troot
    real(r128), allocatable :: tgl(:), wgl(:), Dgl(:,:), w_bclag(:)
    real(r128), allocatable :: Legmat(:,:), vtmp(:,:), bclagmatlr(:,:)
    real(r128), allocatable :: root_re(:), root_im(:), xyz_hat(:,:)
    real(r128), allocatable :: funvals_ell(:,:,:), sxbdw_ell(:,:,:)
    logical, allocatable :: root_ok(:)
    real(r128) :: templ, tempr, denoml, denomr

    allocate(tgl(nquad), wgl(nquad), Dgl(nquad,nquad), w_bclag(nquad))
    allocate(Legmat(nquad,nquad), vtmp(nquad,nquad), bclagmatlr(nquad,2))
    allocate(root_re(m), root_im(m), root_ok(m))
    allocate(funvals_ell(nquad,1,m), sxbdw_ell(nquad,1,m))

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

    n_expa = min(16_8, nquad)
    rho = 4.0_r128**(16.0_r128 / real(nquad, r128))
    allocate(xyz_hat(n_expa,3))

    q_lq128 = 0.0_r128
    do ell = 1, sbdnp
      idx_start = (ell-1_8)*nquad + 1_8
      idx_end = ell*nquad

      xyz_hat(:,1) = matmul(Legmat(1:n_expa,:), sxbd(1,idx_start:idx_end))
      xyz_hat(:,2) = matmul(Legmat(1:n_expa,:), sxbd(2,idx_start:idx_end))
      xyz_hat(:,3) = matmul(Legmat(1:n_expa,:), sxbd(3,idx_start:idx_end))

      do j = 1, m
        call line_quad_root_initial_guess_r128(tgl, sxbd(1,idx_start:idx_end), &
                                               sxbd(2,idx_start:idx_end), &
                                               sxbd(3,idx_start:idx_end), &
                                               nquad, r0(1,j), r0(2,j), r0(3,j), tinit)
        root_ok(j) = .false.
        root_re(j) = 0.0_r128
        root_im(j) = 0.0_r128
        br = abs(tinit + sqrt(tinit - 1.0_r128)*sqrt(tinit + 1.0_r128))
        if (br < 1.75_r128*rho) then
          call line_quad_root_refine_r128(xyz_hat(:,1), xyz_hat(:,2), xyz_hat(:,3), &
                                          n_expa, r0(1,j), r0(2,j), r0(3,j), &
                                          tinit, troot, ifconv)
          br = abs(troot + sqrt(troot - 1.0_r128)*sqrt(troot + 1.0_r128))
          if (ifconv == 1_8 .and. br < rho) then
            root_ok(j) = .true.
            root_re(j) = real(troot, r128)
            root_im(j) = aimag(troot)
          end if
        end if
      end do

      funvals_ell = 0.0_r128
      call line_kernel_eval_r128(m, r0, nquad, 1_8, nquad, &
                                 sxbd(:,idx_start:idx_end), &
                                 sxpbd(:,idx_start:idx_end), &
                                 stangbd(:,idx_start:idx_end), &
                                 fun, kdata, funvals_ell)

      sxbdw_ell = 0.0_r128
      call line_quad_compress_nearroot_r128(m, r0, nquad, 1_8, nquad, &
                                            sxbd(:,idx_start:idx_end), &
                                            sxpbd(:,idx_start:idx_end), &
                                            stangbd(:,idx_start:idx_end), &
                                            sspbd(idx_start:idx_end), &
                                            tgl, wgl, Dgl, w_bclag, &
                                            Legmat, bclagmatlr, &
                                            fun, kdata, funvals_ell, sxbdw_ell, &
                                            root_re, root_im, root_ok)

      do j = 1, m
        if (root_ok(j)) then
          q_lq128(j) = q_lq128(j) + sum(funvals_ell(:,1,j) * sxbdw_ell(:,1,j))
        end if
      end do
    end do

    deallocate(tgl, wgl, Dgl, w_bclag, Legmat, vtmp, bclagmatlr)
    deallocate(root_re, root_im, root_ok, xyz_hat, funvals_ell, sxbdw_ell)
    
  end subroutine evaluate_line_integral_r128

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

  subroutine invr_kernel_r128(r_s, tau_s, r0j, kdata, val)
    real(r128), intent(in)    :: r_s(3), tau_s(3), r0j(3), kdata(3)
    real(r128), intent(inout) :: val

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

end module solidangle_mod
