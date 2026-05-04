! lq_kernel.f90
! Module: line_kernel_eval_r64 only.
! fun is type(c_funptr),value — declared before kdata (matches dummy-list order).
! Caller constructs c_funptr from integer(8) handle before calling.

module lq_kernel_mod
  use iso_c_binding, only: c_funptr, c_f_procpointer, c_double
  use linequaaadrature_mod, only: gauss_r64, gauss_r128, bclaginterpweights_r64, bclaginterpweights_r128
  implicit none

  integer, parameter :: r64  = 8
  integer, parameter :: r128 = 16

  ! ----------------------------------------------------------------
  ! Kernel-id selectors used by line_quad_BrF_r64 to pick the
  ! r128-internal kernel formula (BrF does its barycentric Lagrange
  ! interpolation and kernel evaluation in r128, so it can't go
  ! through the user-provided r64 callback).
  !
  ! To extend BrF to a new line kernel:
  !   1. Add a new KERNEL_<NAME> integer(8), parameter below.
  !   2. Add a `case (KERNEL_<NAME>)` branch inside line_quad_BrF_r64
  !      with the inline r128 formula.
  !   3. Have the caller of line_quad_compress_nearroot_r64 pass
  !      kernel_id=KERNEL_<NAME>.
  ! ----------------------------------------------------------------
  integer(8), parameter :: KERNEL_INVR     = 0_8   ! 1/|r-r0|^p, kdata(1)=p in {1,3,5}
  integer(8), parameter :: KERNEL_ASVESTAS = 1_8   ! Asvestas solid-angle kernel, kdata(1:3)=qhat

  logical :: lq_profile_enabled_r64 = .false.
  integer(8) :: lq_profile_kernel_eval_calls_r64 = 0_8
  integer(8) :: lq_profile_compress_calls_r64 = 0_8
  integer(8) :: lq_profile_compress_targets_r64 = 0_8
  integer(8) :: lq_profile_panel_checks_r64 = 0_8
  integer(8) :: lq_profile_total_nup_r64 = 0_8
  integer(8) :: lq_profile_max_nup_r64 = 0_8
  real(r64) :: lq_profile_kernel_eval_sec_r64 = 0.0_r64
  real(r64) :: lq_profile_setup_sec_r64 = 0.0_r64
  real(r64) :: lq_profile_geom_sec_r64 = 0.0_r64
  real(r64) :: lq_profile_bisect_sec_r64 = 0.0_r64
  real(r64) :: lq_profile_panel_map_sec_r64 = 0.0_r64
  real(r64) :: lq_profile_interp_sec_r64 = 0.0_r64
  real(r64) :: lq_profile_callback_sec_r64 = 0.0_r64
  real(r64) :: lq_profile_panel_sum_sec_r64 = 0.0_r64
  real(r64) :: lq_profile_panel_search_sec_r64 = 0.0_r64
  real(r64) :: lq_profile_pack_sec_r64 = 0.0_r64
  real(r64) :: lq_profile_br_sec_r64 = 0.0_r64
  real(r64) :: lq_profile_compress_sec_r64 = 0.0_r64
  real(r64) :: lq_profile_weight_sec_r64 = 0.0_r64

  abstract interface
    subroutine kernel_iface_r64(r_s, tau_s, r0j, kdata3, val) bind(C)
      import c_double
      real(c_double), intent(in)    :: r_s(3), tau_s(3), r0j(3), kdata3(3)
      real(c_double), intent(inout) :: val
    end subroutine kernel_iface_r64
  end interface

  abstract interface
    subroutine kernel_iface_r128(r_s, tau_s, r0j, kdata3, val)
      integer, parameter :: rr = 16
      real(rr), intent(in)    :: r_s(3), tau_s(3), r0j(3), kdata3(3)
      real(rr), intent(inout) :: val
    end subroutine kernel_iface_r128
  end interface

contains

  subroutine lq_profile_reset_r64()
    lq_profile_kernel_eval_calls_r64 = 0_8
    lq_profile_compress_calls_r64 = 0_8
    lq_profile_compress_targets_r64 = 0_8
    lq_profile_panel_checks_r64 = 0_8
    lq_profile_total_nup_r64 = 0_8
    lq_profile_max_nup_r64 = 0_8
    lq_profile_kernel_eval_sec_r64 = 0.0_r64
    lq_profile_setup_sec_r64 = 0.0_r64
    lq_profile_geom_sec_r64 = 0.0_r64
    lq_profile_bisect_sec_r64 = 0.0_r64
    lq_profile_panel_map_sec_r64 = 0.0_r64
    lq_profile_interp_sec_r64 = 0.0_r64
    lq_profile_callback_sec_r64 = 0.0_r64
    lq_profile_panel_sum_sec_r64 = 0.0_r64
    lq_profile_panel_search_sec_r64 = 0.0_r64
    lq_profile_pack_sec_r64 = 0.0_r64
    lq_profile_br_sec_r64 = 0.0_r64
    lq_profile_compress_sec_r64 = 0.0_r64
    lq_profile_weight_sec_r64 = 0.0_r64
  end subroutine lq_profile_reset_r64

  ! ----------------------------------------------------------------
  ! line_kernel_eval_r64
  ! Evaluate kernel at every (boundary quad node, target) pair.
  ! fun: c_funptr to subroutine(r_s(3), tau_s(3), r0j(3), kdata3(3), val)
  ! kdata(3,m): per-target kernel data (e.g. qhat for solid angle).
  ! funvals(nquad, sbdnp, m): output — kernel * speed at each node.
  !
  ! CRITICAL: fun declared before kdata to match dummy-list order.
  ! ----------------------------------------------------------------
  subroutine line_kernel_eval_r64(m, r0, nbd, sbdnp, nquad, &
                                   sxbd, sxpbd, stangbd,    &
                                   fun, kdata, funvals)
    integer(8),     intent(in)    :: m, nbd, sbdnp, nquad
    real(r64),      intent(in)    :: r0(3,m)
    real(r64),      intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
    type(c_funptr), value         :: fun       ! ← before kdata; matches dummy list
    real(r64),      intent(in)    :: kdata(3,m)
    real(r64),      intent(inout) :: funvals(nquad,sbdnp,m)

    procedure(kernel_iface_r64), pointer :: fptr
    integer(8) :: ell, j, q, idx
    integer :: c0, c1, rate
    real(r64)  :: r_s(3), tau_s(3), sp, val

    if (lq_profile_enabled_r64) call system_clock(c0, rate)
    call c_f_procpointer(fun, fptr)

    do ell = 1, sbdnp
      do q = 1, nquad
        idx   = (ell-1)*nquad + q
        r_s   = sxbd(:, idx)
        tau_s = stangbd(:, idx)
        sp    = sqrt(sxpbd(1,idx)**2 + sxpbd(2,idx)**2 + sxpbd(3,idx)**2)
        do j = 1, m
          val = 0.0_r64
          call fptr(r_s, tau_s, r0(:,j), kdata(:,j), val)
          funvals(q, ell, j) = val * sp
        end do
      end do
    end do

    if (lq_profile_enabled_r64) then
      call system_clock(c1)
      lq_profile_kernel_eval_calls_r64 = lq_profile_kernel_eval_calls_r64 + 1_8
      lq_profile_kernel_eval_sec_r64 = lq_profile_kernel_eval_sec_r64 + &
        real(c1 - c0, r64)/real(rate, r64)
    end if

  end subroutine line_kernel_eval_r64

  subroutine line_kernel_eval_local_r64(nquad, npan, &
                                      xdisp, ydisp, zdisp, stangpan, sppan, dswpan, &
                                      target_loc, fun, kdata, &
                                      funvals_local, integrand0_up, I_local, kval)
    integer(8),     intent(in)    :: nquad, npan
    real(r64),      intent(in)    :: xdisp(nquad,npan), ydisp(nquad,npan), zdisp(nquad,npan)
    real(r64),      intent(in)    :: stangpan(3,nquad,npan)
    real(r64),      intent(in)    :: sppan(nquad,npan), dswpan(nquad,npan)
    real(r64),      intent(in)    :: target_loc(3)
    type(c_funptr), value         :: fun
    real(r64),      intent(in)    :: kdata(3)
    real(r64),      intent(inout) :: funvals_local(nquad,npan)
    real(r64),      intent(inout) :: integrand0_up(nquad*npan)
    real(r64),      intent(inout) :: I_local
    real(r64),      intent(inout) :: kval(nquad,npan)

    procedure(kernel_iface_r64), pointer :: fptr
    integer(8) :: ipan, q, idx
    real(r64)  :: r_s(3), tau_s(3), val

    call c_f_procpointer(fun, fptr)

    I_local = 0.0_r64

    do ipan = 1, npan
      do q = 1, nquad
        idx = (ipan - 1_8)*nquad + q

        r_s(1) = xdisp(q,ipan)
        r_s(2) = ydisp(q,ipan)
        r_s(3) = zdisp(q,ipan)

        tau_s = stangpan(:,q,ipan)

        val = 0.0_r64
        call fptr(r_s, tau_s, target_loc, kdata, val)

        kval(q,ipan) = val
        funvals_local(q,ipan) = val*sppan(q,ipan)
        integrand0_up(idx) = val
        I_local = I_local + val*dswpan(q,ipan)
      end do
    end do
  end subroutine line_kernel_eval_local_r64

  ! ----------------------------------------------------------------
  ! bary_row_r64  (private)
  ! Barycentric type-II Lagrange basis values at xi:
  !   row(k) = L_k^{tgl}(xi)
  ! ----------------------------------------------------------------
  subroutine bary_row_r64(nquad, tgl, w_bclag, xi, row)
    integer(8), intent(in)  :: nquad
    real(r64),  intent(in)  :: tgl(nquad), w_bclag(nquad), xi
    real(r64),  intent(out) :: row(nquad)

    integer(8) :: k
    real(r64)  :: eps_hit, denom

    eps_hit = 100.0_r64 * epsilon(1.0_r64)
    row = 0.0_r64
    do k = 1, nquad
      if (abs(xi - tgl(k)) <= eps_hit * max(1.0_r64, abs(tgl(k)))) then
        row(k) = 1.0_r64
        return
      end if
    end do
    do k = 1, nquad
      row(k) = w_bclag(k) / (xi - tgl(k))
    end do
    denom = sum(row)
    row   = row / denom

  end subroutine bary_row_r64

  subroutine build_nearroot_panel_ends_r64(t_root, len, lenl, lenr, pan_t_end)
    real(r64),  intent(in)  :: t_root
    integer(8), intent(in)  :: len, lenl, lenr
    real(r64),  intent(out) :: pan_t_end(len)

    real(r64), parameter :: factor = 3.0_r64
    real(r64) :: pan_t_end1(lenl+1), pan_t_end2(lenr+1), factor_inv
    integer(8) :: k

    factor_inv = 1.0_r64/factor
    pan_t_end = 0.0_r64
    pan_t_end1 = 0.0_r64
    pan_t_end2 = 0.0_r64

    if (t_root >= 1.0_r64) then
      pan_t_end1(lenl+1) = 1.0_r64
      do k = lenl, 2, -1
        pan_t_end1(k) = factor_inv*pan_t_end1(k+1)
      end do
      do k = 1, len
        pan_t_end(k) = 1.0_r64 - 2.0_r64*pan_t_end1(lenl+2-k)
      end do
    else if (t_root <= -1.0_r64) then
      pan_t_end2(lenr+1) = 1.0_r64
      do k = lenr, 2, -1
        pan_t_end2(k) = factor_inv*pan_t_end2(k+1)
      end do
      do k = 1, len
        pan_t_end(k) = 2.0_r64*pan_t_end2(k) - 1.0_r64
      end do
    else
      pan_t_end1(lenl+1) = 1.0_r64
      do k = lenl, 2, -1
        pan_t_end1(k) = factor_inv*pan_t_end1(k+1)
      end do
      pan_t_end2(lenr+1) = 1.0_r64
      do k = lenr, 2, -1
        pan_t_end2(k) = factor_inv*pan_t_end2(k+1)
      end do
      do k = lenl + 1, 1, -1
        pan_t_end(k) = -1.0_r64 + (t_root + 1.0_r64)*(1.0_r64 - pan_t_end1(lenl+2-k))
      end do
      do k = 1, lenr + 1
        pan_t_end(lenl+k) = 1.0_r64 + (1.0_r64 - t_root)*(-1.0_r64 + pan_t_end2(k))
      end do
    end if
  end subroutine build_nearroot_panel_ends_r64

  subroutine build_nearroot_nodes_r64(t_root, nquad, tgl, wgl, len, lenl, lenr, &
                                      t_up, w_ref)
    real(r64),  intent(in)  :: t_root
    integer(8), intent(in)  :: nquad, len, lenl, lenr
    real(r64),  intent(in)  :: tgl(nquad), wgl(nquad)
    real(r64),  intent(out) :: t_up(nquad*(len-1)), w_ref(nquad*(len-1))

    real(r64) :: pan_t_end(len), pan_t_mid, pan_t_len
    integer(8) :: j, k, idx0

    call build_nearroot_panel_ends_r64(t_root, len, lenl, lenr, pan_t_end)
    do j = 1, len - 1
      pan_t_mid = 0.5_r64*(pan_t_end(j) + pan_t_end(j+1))
      pan_t_len = 0.5_r64*(pan_t_end(j+1) - pan_t_end(j))
      idx0 = (j - 1_8)*nquad
      do k = 1, nquad
        t_up(idx0+k) = pan_t_mid + tgl(k)*pan_t_len
        w_ref(idx0+k) = wgl(k)*pan_t_len
      end do
    end do
  end subroutine build_nearroot_nodes_r64

  subroutine build_nearroot_panels_local( t0, nquad, npan, lenl, lenr, tgl, wgl, &
                                    xjhat, yjhat, zjhat, rbase, &
                                    tpan, upan, wpan, &
                                    xpan, ypan, zpan, &
                                    xdisp, ydisp, zdisp, &
                                    stangpan, sppan, dswpan)
    real(r64),  intent(in)    :: t0
    integer(8), intent(in)    :: nquad, npan, lenl, lenr
    real(r64),  intent(in)    :: tgl(nquad), wgl(nquad)
    real(r64),  intent(in)    :: xjhat(nquad), yjhat(nquad), zjhat(nquad)
    real(r64),  intent(in)    :: rbase(3)
    real(r64),  intent(inout) :: tpan(nquad,npan), upan(nquad,npan), wpan(nquad,npan)
    real(r64),  intent(inout) :: xpan(nquad,npan), ypan(nquad,npan), zpan(nquad,npan)
    real(r64),  intent(inout) :: xdisp(nquad,npan), ydisp(nquad,npan), zdisp(nquad,npan)
    real(r64),  intent(inout) :: stangpan(3,nquad,npan), sppan(nquad,npan), dswpan(nquad,npan)
    real(r64), parameter :: factor = 2.0_r64
    integer(8) :: ipan, k
    real(r64)  :: rhoL, rhoR, ua, ub, h, c
    real(r64)  :: ubreak(npan+1)
    real(r64)  :: p0(nquad)
    real(r64)  :: tt(nquad), uu(nquad)
    real(r64)  :: ppan(nquad,nquad), dppan(nquad,nquad), qpan(nquad,nquad)
    real(r64)  :: dxpan(nquad), dypan(nquad), dzpan(nquad)
    real(r64)  :: rk, rkp1
    if (npan /= lenl + lenr) error stop 'build_nearroot_panels_local: npan /= lenl + lenr'
    if (npan < 1_8) error stop 'build_nearroot_panels_local: npan < 1'
    if (t0 < -1.0_r64 .or. t0 > 1.0_r64) error stop 'build_nearroot_panels_local: t0 outside [-1,1]'
    rhoL = t0 + 1.0_r64
    rhoR = 1.0_r64 - t0
    ubreak = 0.0_r64
    if (lenl == 0_8) then
      ubreak(1) = 0.0_r64
      do k = 1, lenr
        ubreak(k+1) = rhoR * factor**(-(lenr-k))
      end do
    else if (lenr == 0_8) then
      do k = 1, lenl
        ubreak(k) = -rhoL * factor**(-(k-1_8))
      end do
      ubreak(npan+1) = 0.0_r64
    else
      do k = 1, lenl
        ubreak(k) = -rhoL * factor**(-(k-1_8))
      end do
      ubreak(lenl+1) = 0.0_r64
      do k = 1, lenr
        ubreak(lenl+1+k) = rhoR * factor**(-(lenr-k))
      end do
    end if
    p0 = 0.0_r64
    p0(1) = 1.0_r64
    if (nquad >= 2_8) p0(2) = t0
    do k = 1, nquad-2
      rk = real(k, r64)
      rkp1 = real(k+1_8, r64)
      p0(k+2) = ((2.0_r64*rk + 1.0_r64)*t0*p0(k+1) - rk*p0(k))/rkp1
    end do
    do ipan = 1, npan
      ua = ubreak(ipan)
      ub = ubreak(ipan+1)
      h = 0.5_r64*(ub - ua)
      c = 0.5_r64*(ub + ua)
      upan(:,ipan) = c + h*tgl
      tpan(:,ipan) = t0 + upan(:,ipan)
      wpan(:,ipan) = h*wgl
      uu = upan(:,ipan)
      tt = tpan(:,ipan)
      ppan = 0.0_r64
      dppan = 0.0_r64
      qpan = 0.0_r64
      ppan(:,1) = 1.0_r64
      dppan(:,1) = 0.0_r64
      qpan(:,1) = 0.0_r64
      if (nquad >= 2_8) then
        ppan(:,2) = tt
        dppan(:,2) = 1.0_r64
        qpan(:,2) = uu
      end if
      do k = 1, nquad-2
        rk = real(k, r64)
        rkp1 = real(k+1_8, r64)
        ppan(:,k+2) = ((2.0_r64*rk + 1.0_r64)*tt*ppan(:,k+1) - rk*ppan(:,k))/rkp1
        dppan(:,k+2) = ((2.0_r64*rk + 1.0_r64)*(ppan(:,k+1) + tt*dppan(:,k+1)) - rk*dppan(:,k))/rkp1
        qpan(:,k+2) = ((2.0_r64*rk + 1.0_r64)*(tt*qpan(:,k+1) + uu*p0(k+1)) - rk*qpan(:,k))/rkp1
      end do
      xdisp(:,ipan) = matmul(qpan, xjhat)
      ydisp(:,ipan) = matmul(qpan, yjhat)
      zdisp(:,ipan) = matmul(qpan, zjhat)
      xpan(:,ipan) = rbase(1) + xdisp(:,ipan)
      ypan(:,ipan) = rbase(2) + ydisp(:,ipan)
      zpan(:,ipan) = rbase(3) + zdisp(:,ipan)
      dxpan = matmul(dppan, xjhat)
      dypan = matmul(dppan, yjhat)
      dzpan = matmul(dppan, zjhat)
      sppan(:,ipan) = sqrt(dxpan**2 + dypan**2 + dzpan**2)
      dswpan(:,ipan) = wpan(:,ipan) * sppan(:,ipan)
      stangpan(1,:,ipan) = dxpan / sppan(:,ipan)
      stangpan(2,:,ipan) = dypan / sppan(:,ipan)
      stangpan(3,:,ipan) = dzpan / sppan(:,ipan)
    end do
  end subroutine build_nearroot_panels_local

  subroutine check_nearroot_grid_selftest_r64()
    integer(8), parameter :: nquad_test = 4_8, len = 5_8, lenl = 2_8, lenr = 2_8
    real(r64) :: tgl_test(nquad_test), wgl_test(nquad_test)
    real(r64) :: t_up(nquad_test*(len-1)), w_ref(nquad_test*(len-1))
    integer(8) :: i

    tgl_test = [-0.8611363115940526_r64, -0.3399810435848563_r64, &
                 0.3399810435848563_r64,  0.8611363115940526_r64]
    wgl_test = [0.3478548451374539_r64, 0.6521451548625461_r64, &
                0.6521451548625461_r64, 0.3478548451374539_r64]
    call build_nearroot_nodes_r64(0.0_r64, nquad_test, tgl_test, wgl_test, &
                                  len, lenl, lenr, t_up, w_ref)
    do i = 2, size(t_up)
      if (t_up(i) <= t_up(i-1)) error stop 'nearroot grid selftest: non-increasing t_up'
    end do
    if (minval(w_ref) <= 0.0_r64) error stop 'nearroot grid selftest: non-positive weights'
  end subroutine check_nearroot_grid_selftest_r64

  subroutine estimate_nearroot_lengths_r64(t_root, root_imag_abs, max_len_each_side, len, lenl, lenr)
    real(r64),  intent(in)  :: t_root, root_imag_abs
    integer(8), intent(in)  :: max_len_each_side
    integer(8), intent(out) :: len, lenl, lenr

    ! real(r64), parameter :: factor = 3.0_r64
    real(r64), parameter :: factor = 2.0_r64
    real(r64) :: dist, target_width
    integer(8) :: levels

    ! Distance from the complex root to the real panel [-1,1].
    !
    ! Interior root projection:
    !   closest real point is t_root, so distance is |Im root|.
    !
    ! Exterior root projection:
    !   closest real point is the nearest endpoint, so include real offset.

    if (t_root >= 1.0_r64) then
      dist = sqrt((t_root - 1.0_r64)**2 + root_imag_abs**2)
    else if (t_root <= -1.0_r64) then
      dist = sqrt((t_root + 1.0_r64)**2 + root_imag_abs**2)
    else
      dist = root_imag_abs
    end if

    dist = max(dist, tiny(1.0_r64))
    target_width = min(2.0_r64, 2.0_r64*dist)

    levels = ceiling(log(2.0_r64/target_width)/log(factor)) + 1_8
    levels = max(1_8, min(levels, max_len_each_side))

    if (t_root >= 1.0_r64) then
      lenl = levels
      lenr = 0_8
      len = lenl + 1_8
    else if (t_root <= -1.0_r64) then
      lenl = 0_8
      lenr = levels
      len = lenr + 1_8
    else
      lenl = levels
      lenr = levels
      len = lenl + lenr + 1_8
    end if
  end subroutine estimate_nearroot_lengths_r64

  ! ----------------------------------------------------------------
  ! eval_integrand_vec_r64  (private)
  ! Evaluate kernel(t)*speed(t) at each t in t_pts via bary interp.
  ! fun declared before kdata (matches dummy-list order).
  ! ----------------------------------------------------------------
  subroutine eval_integrand_vec_r64(t_pts, np, r_ell, rp_ell, nquad, &
                                     tgl, w_bclag, fptr, r0j, kdata, vals)
    integer(8),     intent(in)  :: np, nquad
    real(r64),      intent(in)  :: t_pts(np)
    real(r64),      intent(in)  :: r_ell(3,nquad), rp_ell(3,nquad)
    real(r64),      intent(in)  :: tgl(nquad), w_bclag(nquad)
    procedure(kernel_iface_r64), pointer :: fptr
    real(r64),      intent(in)  :: r0j(3), kdata(3)
    real(r64),      intent(out) :: vals(np)

    integer(8) :: i, k, hit_idx
    integer :: c0, c1, rate
    real(r64)  :: row(nquad), r_t(3), rp_t(3), sp_t, tau_t(3), val
    real(r64)  :: xi, eps_hit, denom, wk

    eps_hit = 100.0_r64 * epsilon(1.0_r64)
    do i = 1, np
      if (lq_profile_enabled_r64) call system_clock(c0, rate)
      xi = t_pts(i)
      hit_idx = 0_8
      do k = 1, nquad
        if (abs(xi - tgl(k)) <= eps_hit * max(1.0_r64, abs(tgl(k)))) then
          hit_idx = k
          exit
        end if
      end do

      if (hit_idx > 0_8) then
        r_t = r_ell(:, hit_idx)
        rp_t = rp_ell(:, hit_idx)
      else
        denom = 0.0_r64
        do k = 1, nquad
          row(k) = w_bclag(k)/(xi - tgl(k))
          denom = denom + row(k)
        end do
        r_t = 0.0_r64
        rp_t = 0.0_r64
        do k = 1, nquad
          wk = row(k)/denom
          r_t(1) = r_t(1) + r_ell(1,k)*wk
          r_t(2) = r_t(2) + r_ell(2,k)*wk
          r_t(3) = r_t(3) + r_ell(3,k)*wk
          rp_t(1) = rp_t(1) + rp_ell(1,k)*wk
          rp_t(2) = rp_t(2) + rp_ell(2,k)*wk
          rp_t(3) = rp_t(3) + rp_ell(3,k)*wk
        end do
      end if
      sp_t  = sqrt(rp_t(1)**2 + rp_t(2)**2 + rp_t(3)**2)
      tau_t = rp_t / sp_t
      if (lq_profile_enabled_r64) then
        call system_clock(c1)
        lq_profile_interp_sec_r64 = lq_profile_interp_sec_r64 + &
          real(c1 - c0, r64)/real(rate, r64)
      end if
      val   = 0.0_r64
      if (lq_profile_enabled_r64) call system_clock(c0)
      call fptr(r_t, tau_t, r0j, kdata, val)
      if (lq_profile_enabled_r64) then
        call system_clock(c1)
        lq_profile_callback_sec_r64 = lq_profile_callback_sec_r64 + &
          real(c1 - c0, r64)/real(rate, r64)
      end if
      vals(i) = val * sp_t
    end do

  end subroutine eval_integrand_vec_r64

  subroutine eval_integrand_original_nodes_r64(r_ell, rp_ell, nquad, fptr, r0j, kdata, vals)
    integer(8),     intent(in)  :: nquad
    real(r64),      intent(in)  :: r_ell(3,nquad), rp_ell(3,nquad)
    procedure(kernel_iface_r64), pointer :: fptr
    real(r64),      intent(in)  :: r0j(3), kdata(3)
    real(r64),      intent(out) :: vals(nquad)

    integer(8) :: k
    integer :: c0, c1, rate
    real(r64) :: r_t(3), rp_t(3), sp_t, tau_t(3), val

    do k = 1, nquad
      if (lq_profile_enabled_r64) call system_clock(c0, rate)
      r_t = r_ell(:, k)
      rp_t = rp_ell(:, k)
      sp_t = sqrt(rp_t(1)**2 + rp_t(2)**2 + rp_t(3)**2)
      tau_t = rp_t/sp_t
      if (lq_profile_enabled_r64) then
        call system_clock(c1)
        lq_profile_interp_sec_r64 = lq_profile_interp_sec_r64 + &
          real(c1 - c0, r64)/real(rate, r64)
      end if

      val = 0.0_r64
      if (lq_profile_enabled_r64) call system_clock(c0)
      call fptr(r_t, tau_t, r0j, kdata, val)
      if (lq_profile_enabled_r64) then
        call system_clock(c1)
        lq_profile_callback_sec_r64 = lq_profile_callback_sec_r64 + &
          real(c1 - c0, r64)/real(rate, r64)
      end if
      vals(k) = val*sp_t
    end do
  end subroutine eval_integrand_original_nodes_r64

  ! ----------------------------------------------------------------
  ! panel_int_err_r64  (private)
  ! Fine (nquad) and coarse (nquad2) GL integrals on [aa,bb].
  ! err = |ih - il|.
  ! fhi: fine integrand values at mapped nodes (nquad).
  ! ----------------------------------------------------------------
  subroutine panel_int_err_r64(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                                nquad2, tgl2, wgl2, fptr, r0j, kdata,    &
                                aa, bb, ih, err, fhi)
    integer(8),     intent(in)  :: nquad, nquad2
    real(r64),      intent(in)  :: r_ell(3,nquad), rp_ell(3,nquad)
    real(r64),      intent(in)  :: tgl(nquad), wgl(nquad), w_bclag(nquad)
    real(r64),      intent(in)  :: tgl2(nquad2), wgl2(nquad2)
    procedure(kernel_iface_r64), pointer :: fptr
    real(r64),      intent(in)  :: r0j(3), kdata(3), aa, bb
    real(r64),      intent(out) :: ih, err, fhi(nquad)

    real(r64) :: c, h
    real(r64) :: t_fine(nquad),   v_fine(nquad)
    real(r64) :: t_coarse(nquad2), v_coarse(nquad2)
    real(r64) :: il
    integer(8) :: k
    integer :: c0, c1, rate

    if (lq_profile_enabled_r64) call system_clock(c0, rate)
    c = 0.5_r64 * (aa + bb)
    h = 0.5_r64 * (bb - aa)

    do k = 1, nquad
      t_fine(k) = c + h * tgl(k)
    end do
    if (lq_profile_enabled_r64) then
      call system_clock(c1)
      lq_profile_panel_checks_r64 = lq_profile_panel_checks_r64 + 1_8
      lq_profile_panel_map_sec_r64 = lq_profile_panel_map_sec_r64 + &
        real(c1 - c0, r64)/real(rate, r64)
    end if

    if (aa == -1.0_r64 .and. bb == 1.0_r64) then
      call eval_integrand_original_nodes_r64(r_ell, rp_ell, nquad, fptr, r0j, kdata, v_fine)
    else
      call eval_integrand_vec_r64(t_fine, nquad, r_ell, rp_ell, nquad, &
                                   tgl, w_bclag, fptr, r0j, kdata, v_fine)
    end if
    if (lq_profile_enabled_r64) call system_clock(c0)
    ih  = h * sum(wgl * v_fine)
    fhi = v_fine
    if (lq_profile_enabled_r64) then
      call system_clock(c1)
      lq_profile_panel_sum_sec_r64 = lq_profile_panel_sum_sec_r64 + &
        real(c1 - c0, r64)/real(rate, r64)
    end if

    if (lq_profile_enabled_r64) call system_clock(c0)
    do k = 1, nquad2
      t_coarse(k) = c + h * tgl2(k)
    end do
    if (lq_profile_enabled_r64) then
      call system_clock(c1)
      lq_profile_panel_map_sec_r64 = lq_profile_panel_map_sec_r64 + &
        real(c1 - c0, r64)/real(rate, r64)
    end if

    call eval_integrand_vec_r64(t_coarse, nquad2, r_ell, rp_ell, nquad, &
                                 tgl, w_bclag, fptr, r0j, kdata, v_coarse)
    if (lq_profile_enabled_r64) call system_clock(c0)
    il  = h * sum(wgl2 * v_coarse)
    err = abs(ih - il)
    if (lq_profile_enabled_r64) then
      call system_clock(c1)
      lq_profile_panel_sum_sec_r64 = lq_profile_panel_sum_sec_r64 + &
        real(c1 - c0, r64)/real(rate, r64)
    end if

  end subroutine panel_int_err_r64

  ! ----------------------------------------------------------------
  ! run_bisect_panel_r64  (private)
  ! Error-driven GL bisection on [-1,1].
  ! Splits highest-error panel until err_tot < tol*max(1,|result|)
  ! or npanels >= maxpan.
  ! Outputs t_up(n_up), w_up(n_up), f_up(n_up); n_up = npanels*nquad.
  ! Caller must size t_up/w_up/f_up to at least maxpan*nquad.
  ! ----------------------------------------------------------------
  subroutine run_bisect_panel_r64(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                                   nquad2, tgl2, wgl2, fptr, r0j, kdata,   &
                                   tol, maxpan, t_up, w_up, f_up, n_up)
    integer(8),     intent(in)    :: nquad, nquad2, maxpan
    real(r64),      intent(in)    :: r_ell(3,nquad), rp_ell(3,nquad)
    real(r64),      intent(in)    :: tgl(nquad), wgl(nquad), w_bclag(nquad)
    real(r64),      intent(in)    :: tgl2(nquad2), wgl2(nquad2)
    procedure(kernel_iface_r64), pointer :: fptr
    real(r64),      intent(in)    :: r0j(3), kdata(3), tol
    real(r64),      intent(out)   :: t_up(maxpan*nquad), w_up(maxpan*nquad), f_up(maxpan*nquad)
    integer(8),     intent(out)   :: n_up

    integer(8) :: npanels, idx_max, k, j, idx
    real(r64)  :: a_pan(maxpan), b_pan(maxpan)
    real(r64)  :: intval(maxpan), errp(maxpan), fpan(nquad,maxpan)
    real(r64)  :: result, err_tot, scale, target_err, mid, c, h
    real(r64)  :: ih_l, err_l, fhi_l(nquad)
    real(r64)  :: ih_r, err_r, fhi_r(nquad)
    integer :: c0, c1, rate

    npanels  = 1_8
    a_pan(1) = -1.0_r64
    b_pan(1) =  1.0_r64
    call panel_int_err_r64(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                            nquad2, tgl2, wgl2, fptr, r0j, kdata,    &
                            a_pan(1), b_pan(1), intval(1), errp(1), fpan(:,1))

    do
      if (lq_profile_enabled_r64) call system_clock(c0, rate)
      result     = sum(intval(1:npanels))
      err_tot    = sum(errp(1:npanels))
      scale      = max(1.0_r64, abs(result))
      target_err = max(tol, 100.0_r64 * epsilon(1.0_r64)) * scale
      if (err_tot <= target_err .or. npanels >= maxpan) exit

      idx_max = 1_8
      do k = 2, npanels
        if (errp(k) > errp(idx_max)) idx_max = k
      end do

      mid = 0.5_r64 * (a_pan(idx_max) + b_pan(idx_max))
      if (lq_profile_enabled_r64) then
        call system_clock(c1)
        lq_profile_panel_search_sec_r64 = lq_profile_panel_search_sec_r64 + &
          real(c1 - c0, r64)/real(rate, r64)
      end if
      call panel_int_err_r64(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                              nquad2, tgl2, wgl2, fptr, r0j, kdata,    &
                              a_pan(idx_max), mid, ih_l, err_l, fhi_l)
      call panel_int_err_r64(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                              nquad2, tgl2, wgl2, fptr, r0j, kdata,    &
                              mid, b_pan(idx_max), ih_r, err_r, fhi_r)

      npanels = npanels + 1_8
      a_pan(npanels)   = mid;    b_pan(npanels)   = b_pan(idx_max)
      intval(npanels)  = ih_r;   errp(npanels)    = err_r
      fpan(:,npanels)  = fhi_r
      b_pan(idx_max)   = mid;    intval(idx_max)  = ih_l
      errp(idx_max)    = err_l;  fpan(:,idx_max)  = fhi_l
    end do

    if (lq_profile_enabled_r64) call system_clock(c0, rate)
    n_up = npanels * nquad
    do j = 1, npanels
      c   = 0.5_r64 * (a_pan(j) + b_pan(j))
      h   = 0.5_r64 * (b_pan(j) - a_pan(j))
      idx = (j-1_8)*nquad
      do k = 1, nquad
        t_up(idx+k) = c + h * tgl(k)
        w_up(idx+k) = h * wgl(k)
        f_up(idx+k) = fpan(k,j)
      end do
    end do
    if (lq_profile_enabled_r64) then
      call system_clock(c1)
      lq_profile_pack_sec_r64 = lq_profile_pack_sec_r64 + &
        real(c1 - c0, r64)/real(rate, r64)
    end if

  end subroutine run_bisect_panel_r64

  ! ----------------------------------------------------------------
  ! line_quad_compress_r64  (public)
  ! Path B: error-driven GL bisection weight compression.
  ! For each (panel ell, target j):
  !   1. run_bisect_panel -> upsampled t_up, w_up, f_up
  !   2. build Br: Br(k,i) = (w_up(i)/wgl(k)) * L_k(t_up(i))
  !   3. compress: integrand0_compress = Br * (f_up/sp_up)
  !   4. sxbdw(k,ell,j) = integrand0_compress(k)*sp_ell(k)
  !                       / funvals(k,ell,j) * wgl(k)
  ! ----------------------------------------------------------------
  subroutine line_quad_compress_r64(m, r0, nbd, sbdnp, nquad,          &
                                     sxbd, sxpbd, stangbd, sspbd,       &
                                     tgl, wgl, Dgl, w_bclag,            &
                                     Legmat, bclagmatlr,                 &
                                     fun, kdata, funvals, sxbdw)
    integer(8),     intent(in)    :: m, nbd, sbdnp, nquad
    real(r64),      intent(in)    :: r0(3,m)
    real(r64),      intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
    real(r64),      intent(in)    :: sspbd(nbd)
    real(r64),      intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
    real(r64),      intent(in)    :: w_bclag(nquad)
    real(r64),      intent(in)    :: Legmat(nquad,nquad), bclagmatlr(nquad,2)
    type(c_funptr), value         :: fun
    real(r64),      intent(in)    :: kdata(3,m)
    real(r64),      intent(inout) :: funvals(nquad,sbdnp,m)
    real(r64),      intent(inout) :: sxbdw(nquad,sbdnp,m)

    integer(8), parameter :: maxpan = 128_8
    real(r64),  parameter :: tol    = 1.0e-14_r64

    integer(8) :: ell, j, i, k, idx_start, idx_end, nquad2, n_up
    integer :: c0, c1, rate
    real(r64)  :: r_ell(3,nquad), rp_ell(3,nquad), sp_ell(nquad)
    real(r64)  :: row(nquad), wgl_inv(nquad)
    real(r64)  :: t_up(maxpan*nquad), w_up(maxpan*nquad), f_up(maxpan*nquad)
    real(r64)  :: sp_up, Br(nquad,maxpan*nquad)
    real(r64)  :: integrand0_up(maxpan*nquad), integrand0_compress(nquad)
    real(r64), allocatable :: tgl2(:), wgl2(:), Dgl2(:,:)
    procedure(kernel_iface_r64), pointer :: fptr

    call c_f_procpointer(fun, fptr)
    if (lq_profile_enabled_r64) call system_clock(c0, rate)
    nquad2 = max(1_8, nquad / 2_8)
    allocate(tgl2(nquad2), wgl2(nquad2), Dgl2(nquad2,nquad2))
    ! use linequaaadrature_mod gauss for coarse rule — call directly via external
    ! (linequaaadrature_mod is compiled in the same library)
    call gauss_r64(nquad2, tgl2, wgl2, Dgl2)

    wgl_inv = 1.0_r64 / wgl
    if (lq_profile_enabled_r64) then
      call system_clock(c1)
      lq_profile_compress_calls_r64 = lq_profile_compress_calls_r64 + 1_8
      lq_profile_setup_sec_r64 = lq_profile_setup_sec_r64 + &
        real(c1 - c0, r64)/real(rate, r64)
    end if

    do ell = 1, sbdnp
      if (lq_profile_enabled_r64) call system_clock(c0)
      idx_start = (ell-1)*nquad + 1
      idx_end   = ell*nquad
      r_ell     = sxbd(:, idx_start:idx_end)
      ! rp_ell = Dgl * r_ell  (each row is derivative of one coord)
      rp_ell(1,:) = matmul(Dgl, r_ell(1,:))
      rp_ell(2,:) = matmul(Dgl, r_ell(2,:))
      rp_ell(3,:) = matmul(Dgl, r_ell(3,:))
      do k = 1, nquad
        sp_ell(k) = sqrt(sxpbd(1,idx_start+k-1)**2 + &
                         sxpbd(2,idx_start+k-1)**2 + &
                         sxpbd(3,idx_start+k-1)**2)
      end do
      if (lq_profile_enabled_r64) then
        call system_clock(c1)
        lq_profile_geom_sec_r64 = lq_profile_geom_sec_r64 + &
          real(c1 - c0, r64)/real(rate, r64)
      end if

      do j = 1, m
        if (lq_profile_enabled_r64) call system_clock(c0)
        call run_bisect_panel_r64(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                                   nquad2, tgl2, wgl2, fptr, r0(:,j), kdata(:,j), &
                                   tol, maxpan, t_up, w_up, f_up, n_up)
        if (lq_profile_enabled_r64) then
          call system_clock(c1)
          lq_profile_compress_targets_r64 = lq_profile_compress_targets_r64 + 1_8
          lq_profile_total_nup_r64 = lq_profile_total_nup_r64 + n_up
          lq_profile_max_nup_r64 = max(lq_profile_max_nup_r64, n_up)
          lq_profile_bisect_sec_r64 = lq_profile_bisect_sec_r64 + &
            real(c1 - c0, r64)/real(rate, r64)
        end if

        ! geometry at upsampled nodes (for sp_up)
        if (lq_profile_enabled_r64) call system_clock(c0)
        do i = 1, n_up
          call bary_row_r64(nquad, tgl, w_bclag, t_up(i), row)
          ! rp_up = rp_ell * row
          sp_up = 0.0_r64
          do k = 1, nquad
            sp_up = sp_up + (rp_ell(1,k)*row(k))**2 + &
                             (rp_ell(2,k)*row(k))**2 + &
                             (rp_ell(3,k)*row(k))**2
          end do
          sp_up = sqrt(sum((matmul(rp_ell, row))**2))
          integrand0_up(i) = f_up(i) / sp_up

          ! Br(:,i) = (w_up(i) * wgl_inv) * row
          Br(:,i) = w_up(i) * wgl_inv * row
        end do
        if (lq_profile_enabled_r64) then
          call system_clock(c1)
          lq_profile_br_sec_r64 = lq_profile_br_sec_r64 + &
            real(c1 - c0, r64)/real(rate, r64)
        end if

        ! integrand0_compress = Br(1:nquad, 1:n_up) * integrand0_up(1:n_up)
        if (lq_profile_enabled_r64) call system_clock(c0)
        do k = 1, nquad
          integrand0_compress(k) = sum(Br(k,1:n_up) * integrand0_up(1:n_up))
        end do
        if (lq_profile_enabled_r64) then
          call system_clock(c1)
          lq_profile_compress_sec_r64 = lq_profile_compress_sec_r64 + &
            real(c1 - c0, r64)/real(rate, r64)
        end if

        ! sxbdw(k,ell,j) = integrand0_compress(k)*sp_ell(k) / funvals(k,ell,j) * wgl(k)
        if (lq_profile_enabled_r64) call system_clock(c0)
        do k = 1, nquad
          if (abs(funvals(k,ell,j)) > 0.0_r64) then
            sxbdw(k,ell,j) = integrand0_compress(k) * sp_ell(k) / funvals(k,ell,j) * wgl(k)
          else
            sxbdw(k,ell,j) = wgl(k)
          end if
        end do
        if (lq_profile_enabled_r64) then
          call system_clock(c1)
          lq_profile_weight_sec_r64 = lq_profile_weight_sec_r64 + &
            real(c1 - c0, r64)/real(rate, r64)
        end if
      end do
    end do

    deallocate(tgl2, wgl2, Dgl2)

  end subroutine line_quad_compress_r64

  ! ----------------------------------------------------------------
  ! line_quad_BrF_r64
  !
  ! Build the per-target compression matrix row Br(:,i) and the
  ! reference integrand value integrand0_up(i) at each refined node
  ! t_up(i), with the barycentric Lagrange interpolation AND the
  ! kernel evaluation done in r128.
  !
  ! Why r128 internally: the geometry (r_t, rp_t) at t_up is built
  ! from r_ell/rp_ell via a barycentric Lagrange row, which is
  ! cancellation-prone near the GL nodes. Doing it in r128 stops
  ! that, and once we have r_t128 in r128 the kernel value is also
  ! cheap to keep in r128.
  !
  ! kernel_id selects which inline r128 kernel formula to use.
  ! BrF cannot go through the user-supplied r64 callback fptr,
  ! since that would silently downcast r_t128 -> r64 and lose the
  ! r128 advantage; instead, the supported kernels are reproduced
  ! in r128 in the dispatch below. To add a new kernel:
  !   1. add a KERNEL_<NAME> parameter near the top of lq_kernel_mod.
  !   2. add a `case (KERNEL_<NAME>)` branch in the dispatch below
  !      with the inline r128 formula.
  !   3. caller of compress_nearroot_r64 passes kernel_id=KERNEL_<NAME>.
  ! ----------------------------------------------------------------
  subroutine line_quad_BrF_r64(nquad, n_up, &
                               r_ell, rp_ell, &
                               tgl, wgl, w_bclag, &
                               t_up, w_up, &
                               r0j, kdata, kernel_id, &
                               Br, integrand0_up)
    integer(8), intent(in)  :: nquad, n_up
    real(r64),  intent(in)  :: r_ell(3,nquad), rp_ell(3,nquad)
    real(r64),  intent(in)  :: tgl(nquad), wgl(nquad), w_bclag(nquad)
    real(r64),  intent(in)  :: t_up(n_up), w_up(n_up)
    real(r64),  intent(in)  :: r0j(3), kdata(3)
    integer(8), intent(in)  :: kernel_id
    real(r64),  intent(out) :: Br(nquad, n_up)
    real(r64),  intent(out) :: integrand0_up(n_up)

    integer(8) :: i, k
    real(r128) :: r_t128(3), rp_t128(3), sp_t128, tau_t128(3)
    real(r128) :: rvec128(3), rdist128, rinv128, val128
    real(r128) :: rhat128(3), qhat128(3), qcrossrhat128(3), qdotrhat128
    real(r128) :: eps_hit128, denom128
    real(r128) :: row128(nquad)
    integer(8) :: power_int, hit_idx
    real(r128) :: r0j128(3), kdata128(3)
    real(r128) :: tgl128(nquad), wgl128(nquad), wgl_inv128(nquad), w_bclag128(nquad)
    real(r128) :: r_ell128(3,nquad), rp_ell128(3,nquad)
    real(r128) :: t_up128(n_up), w_up128(n_up)

    kdata128 = real(kdata, r128)
    r0j128 = real(r0j, r128)
    tgl128 = real(tgl, r128)
    wgl128 = real(wgl, r128)
    r_ell128 = real(r_ell, r128)
    rp_ell128 = real(rp_ell, r128)
    t_up128(1:n_up) = real(t_up(1:n_up), r128)
    w_up128(1:n_up) = real(w_up(1:n_up), r128)
    wgl_inv128 = 1.0_r128 / wgl128
    call bclaginterpweights_r128(nquad, tgl128, w_bclag128)

    ! Precompute kernel-specific loop-invariants outside the i loop.
    power_int = 0_8
    qhat128   = 0.0_r128
    select case (kernel_id)
    case (KERNEL_INVR)
      power_int = nint(real(kdata128(1), 8))
    case (KERNEL_ASVESTAS)
      qhat128 = kdata128(1:3)
    end select

    eps_hit128 = 10.0_r128 * epsilon(1.0_r128)
    do i = 1, n_up
      ! inlined bary_row_r128 fused with matmul(r_ell128,row128) and matmul(rp_ell128,row128)
      hit_idx = 0_8
      do k = 1, nquad
        if (abs(t_up128(i) - tgl128(k)) <= eps_hit128 * max(1.0_r128, abs(tgl128(k)))) then
          hit_idx = k
          exit
        end if
      end do
      if (hit_idx > 0_8) then
        row128 = 0.0_r128
        row128(hit_idx) = 1.0_r128
        r_t128  = r_ell128(:,  hit_idx)
        rp_t128 = rp_ell128(:, hit_idx)
      else
        denom128 = 0.0_r128
        r_t128   = 0.0_r128
        rp_t128  = 0.0_r128
        do k = 1, nquad
          row128(k)  = w_bclag128(k) / (t_up128(i) - tgl128(k))
          denom128   = denom128   + row128(k)
          r_t128(1)  = r_t128(1)  + r_ell128(1,k)  * row128(k)
          r_t128(2)  = r_t128(2)  + r_ell128(2,k)  * row128(k)
          r_t128(3)  = r_t128(3)  + r_ell128(3,k)  * row128(k)
          rp_t128(1) = rp_t128(1) + rp_ell128(1,k) * row128(k)
          rp_t128(2) = rp_t128(2) + rp_ell128(2,k) * row128(k)
          rp_t128(3) = rp_t128(3) + rp_ell128(3,k) * row128(k)
        end do
        row128  = row128  / denom128
        r_t128  = r_t128  / denom128
        rp_t128 = rp_t128 / denom128
      end if

      ! Geometry shared by every kernel.
      rvec128(1) = r_t128(1) - r0j128(1)
      rvec128(2) = r_t128(2) - r0j128(2)
      rvec128(3) = r_t128(3) - r0j128(3)
      rdist128   = sqrt(rvec128(1)**2 + rvec128(2)**2 + rvec128(3)**2)

      ! r128 kernel dispatch. See module-top comment for how to add a
      ! new kernel here.
      select case (kernel_id)
      case (KERNEL_INVR)
        rinv128 = 1.0_r128 / rdist128
        select case (power_int)
        case (1_8)
          val128 = rinv128
        case (3_8)
          val128 = rinv128**3
        case (5_8)
          val128 = rinv128**5
        case default
          val128 = 0.0_r128
        end select
      case (KERNEL_ASVESTAS)
        ! val = -[tau_s . (qhat x rhat)] / [|r| * (1 - qhat.rhat)]
        sp_t128  = sqrt(rp_t128(1)**2 + rp_t128(2)**2 + rp_t128(3)**2)
        tau_t128 = rp_t128 / sp_t128
        rhat128  = rvec128 / rdist128
        qcrossrhat128(1) = qhat128(2)*rhat128(3) - qhat128(3)*rhat128(2)
        qcrossrhat128(2) = qhat128(3)*rhat128(1) - qhat128(1)*rhat128(3)
        qcrossrhat128(3) = qhat128(1)*rhat128(2) - qhat128(2)*rhat128(1)
        qdotrhat128      = qhat128(1)*rhat128(1) + qhat128(2)*rhat128(2) + qhat128(3)*rhat128(3)
        val128 = -(tau_t128(1)*qcrossrhat128(1) + tau_t128(2)*qcrossrhat128(2) + tau_t128(3)*qcrossrhat128(3)) &
                 / (rdist128 * (1.0_r128 - qdotrhat128))
      case default
        val128 = 0.0_r128
      end select

      Br(:,i) = real(w_up128(i) * wgl_inv128 * row128, r64)
      integrand0_up(i) = real(val128, r64)
    end do

  end subroutine line_quad_BrF_r64

  subroutine line_quad_compress_nearroot_r64(m, r0, nbd, sbdnp, nquad,  &
                                     sxbd, sxpbd, stangbd, sspbd,       &
                                     tgl, wgl, Dgl, w_bclag,            &
                                     Legmat, bclagmatlr,                 &
                                     fun, kdata, funvals, sxbdw,         &
                                     root_re, root_im, root_ok, kernel_id)
    integer(8),     intent(in)    :: m, nbd, sbdnp, nquad
    real(r64),      intent(in)    :: r0(3,m)
    real(r64),      intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
    real(r64),      intent(in)    :: sspbd(nbd)
    real(r64),      intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
    real(r64),      intent(in)    :: w_bclag(nquad)
    real(r64),      intent(in)    :: Legmat(nquad,nquad), bclagmatlr(nquad,2)
    type(c_funptr), value         :: fun
    real(r64),      intent(in)    :: kdata(3,m)
    real(r64),      intent(inout) :: funvals(nquad,sbdnp,m)
    real(r64),      intent(inout) :: sxbdw(nquad,sbdnp,m)
    real(r64),      intent(in)    :: root_re(m), root_im(m)
    logical,        intent(in)    :: root_ok(m)
    ! Optional: which inline r128 kernel formula line_quad_BrF_r64 uses
    ! on the root_ok branch. Defaults to KERNEL_INVR (1/|r-r0|^p) so
    ! existing call sites that build on the inverse-power kernel
    ! (build_target_nearroot_weights_local_r64, the MEX wrappers, etc.)
    ! keep their old behavior with no source change. The solid-angle
    ! caller passes KERNEL_ASVESTAS.
    integer(8),     intent(in), optional :: kernel_id
    integer(8), parameter :: maxpan = 128_8, max_len_each_side = 12_8
    real(r64),  parameter :: tol = 1.0e-14_r64

    integer(8) :: ell, j, i, k, idx_start, idx_end, nquad2, n_up
    integer(8) :: len, lenl, lenr
    integer(8) :: k_id
    integer :: c0, c1, rate
    real(r64)  :: r_ell(3,nquad), rp_ell(3,nquad), sp_ell(nquad)
    real(r64)  :: row(nquad), wgl_inv(nquad)
    real(r64)  :: t_up(maxpan*nquad), w_up(maxpan*nquad), f_up(maxpan*nquad)
    real(r64)  :: sp_up, Br(nquad,maxpan*nquad)
    real(r64)  :: integrand0_up(maxpan*nquad), integrand0_compress(nquad)
    real(r64), allocatable :: tgl2(:), wgl2(:), Dgl2(:,:)
    procedure(kernel_iface_r64), pointer :: fptr

    if (present(kernel_id)) then
      k_id = kernel_id
    else
      k_id = KERNEL_INVR
    end if

    lq_profile_enabled_r64 = .false.  ! disable profiling for near-root code (too much overhead for small n_up)
    if (.not.lq_profile_enabled_r64) then

	    call c_f_procpointer(fun, fptr)
	    nquad2 = max(1_8, nquad / 2_8)
	    allocate(tgl2(nquad2), wgl2(nquad2), Dgl2(nquad2,nquad2))
	    call gauss_r64(nquad2, tgl2, wgl2, Dgl2)

	    wgl_inv = 1.0_r64 / wgl

    do ell = 1, sbdnp
      idx_start = (ell-1)*nquad + 1
      idx_end   = ell*nquad
      r_ell     = sxbd(:, idx_start:idx_end)
      rp_ell(1,:) = matmul(Dgl, r_ell(1,:))
      rp_ell(2,:) = matmul(Dgl, r_ell(2,:))
      rp_ell(3,:) = matmul(Dgl, r_ell(3,:))
      do k = 1, nquad
        sp_ell(k) = sqrt(sxpbd(1,idx_start+k-1)**2 + &
                         sxpbd(2,idx_start+k-1)**2 + &
                         sxpbd(3,idx_start+k-1)**2)
      end do

      do j = 1, m
        if (root_ok(j)) then

          call estimate_nearroot_lengths_r64(root_re(j), abs(root_im(j)), &
                                              max_len_each_side, len, lenl, lenr)
          n_up = nquad*(len - 1_8)
          if (n_up > maxpan*nquad) error stop 'line_quad_compress_nearroot_r128: n_up exceeds work array, you are asking too much, man/woman'
          call build_nearroot_nodes_r64(root_re(j), nquad, tgl, wgl, len, lenl, lenr, t_up(1:n_up), w_up(1:n_up))

          ! r128-internal barycentric interpolation + r128 kernel eval.
          ! BrF dispatches on k_id (KERNEL_INVR / KERNEL_ASVESTAS); the
          ! caller's r64 callback `fptr` is intentionally bypassed here
          ! because the cancellation-prone barycentric step is what we
          ! want in r128, and the matching r128 kernel formula is
          ! inlined inside BrF.
          call line_quad_BrF_r64(nquad, n_up, r_ell, rp_ell,         &
                                 tgl, wgl, w_bclag,                  &
                                 t_up(1:n_up), w_up(1:n_up),         &
                                 r0(:,j), kdata(:,j), k_id,          &
                                 Br(:,1:n_up), integrand0_up(1:n_up))
          ! block
          !   real(r128) :: r_t128(3), rp_t128(3), sp_t128
          !   real(r128) :: dx128, dy128, dz128, r2_128, rinv128, val128
          !   real(r128) :: eps_hit128, denom128
          !   real(r128) :: row128(nquad)
          !   integer(8) :: power_int, hit_idx
          !   real(r128) :: tgl128(nquad), wgl128(nquad), wgl_inv128(nquad), w_bclag128(nquad)
          !   ! real(r128) :: Dgl128(nquad,nquad)
          !   real(r128) :: r0j128(3), kdata128(3)
          !   real(r128) :: r_ell128(3,nquad), rp_ell128(3,nquad)
          !   real(r128) :: t_up128(maxpan*nquad), w_up128(maxpan*nquad)

          !   kdata128 = real(kdata(:,j), r128)
          !   r0j128 = real(r0(:,j), r128)
          !   tgl128 = real(tgl, r128)
          !   wgl128 = real(wgl, r128)
          !   ! Dgl128 = real(Dgl, r128)
          !   r_ell128 = real(r_ell, r128)
          !   rp_ell128 = real(sxpbd(:, idx_start:idx_end), r128)
          !   t_up128(1:n_up) = real(t_up(1:n_up), r128)
          !   w_up128(1:n_up) = real(w_up(1:n_up), r128)
          !   wgl_inv128 = 1.0_r128 / wgl128

          !   power_int = nint(real(kdata128(1), 8))
          !   eps_hit128 = 10.0_r128 * epsilon(1.0_r128)
          !   do i = 1, n_up
          !     ! inlined bary_row_r128 fused with matmul(r_ell128,row128) and matmul(rp_ell128,row128)
          !     hit_idx = 0_8
          !     do k = 1, nquad
          !       if (abs(t_up128(i) - tgl128(k)) <= eps_hit128 * max(1.0_r128, abs(tgl128(k)))) then
          !         hit_idx = k
          !         exit
          !       end if
          !     end do
          !     if (hit_idx > 0_8) then
          !       row128 = 0.0_r128
          !       row128(hit_idx) = 1.0_r128
          !       r_t128  = r_ell128(:,  hit_idx)
          !       rp_t128 = rp_ell128(:, hit_idx)
          !     else
          !       denom128 = 0.0_r128
          !       r_t128   = 0.0_r128
          !       rp_t128  = 0.0_r128
          !       do k = 1, nquad
          !         row128(k)  = w_bclag128(k) / (t_up128(i) - tgl128(k))
          !         denom128   = denom128   + row128(k)
          !         r_t128(1)  = r_t128(1)  + r_ell128(1,k)  * row128(k)
          !         r_t128(2)  = r_t128(2)  + r_ell128(2,k)  * row128(k)
          !         r_t128(3)  = r_t128(3)  + r_ell128(3,k)  * row128(k)
          !         rp_t128(1) = rp_t128(1) + rp_ell128(1,k) * row128(k)
          !         rp_t128(2) = rp_t128(2) + rp_ell128(2,k) * row128(k)
          !         rp_t128(3) = rp_t128(3) + rp_ell128(3,k) * row128(k)
          !       end do
          !       row128  = row128  / denom128
          !       r_t128  = r_t128  / denom128
          !       rp_t128 = rp_t128 / denom128
          !     end if
          !     ! sp_t128 = sqrt(rp_t128(1)**2 + rp_t128(2)**2 + rp_t128(3)**2)
          !     dx128 = r_t128(1) - r0j128(1)
          !     dy128 = r_t128(2) - r0j128(2)
          !     dz128 = r_t128(3) - r0j128(3)
          !     r2_128 = dx128*dx128 + dy128*dy128 + dz128*dz128
          !     rinv128 = 1.0_r128 / sqrt(r2_128)
          !     select case (power_int)
          !     case (1)
          !       val128 = rinv128
          !     case (3)
          !       val128 = rinv128**3
          !     case (5)
          !       val128 = rinv128**5
          !     case default
          !       val128 = 0.0_r128
          !     end select
          !     Br(:,i) = real(w_up128(i) * wgl_inv128 * row128, r64)
          !     integrand0_up(i) = real(val128, r64)
          !   end do
          ! end block

          do k = 1, nquad
            integrand0_compress(k) = sum(Br(k,1:n_up) * integrand0_up(1:n_up))
            if (abs(funvals(k,ell,j)) > 0.0_r64) then
              sxbdw(k,ell,j) = integrand0_compress(k) * sp_ell(k) / funvals(k,ell,j) * wgl(k)
            else
              sxbdw(k,ell,j) = wgl(k)
            end if
          end do
          
          ! ! mixed r128 implementation
          ! block
          !   real(r128) :: root_re128, root_im128
          !   real(r128) :: tgl128(nquad), wgl128(nquad)
          !   real(r128) :: t_up128(maxpan*nquad), w_up128(maxpan*nquad)

          !   root_re128 = real(root_re(j), r128)
          !   root_im128 = abs(real(root_im(j), r128))
          !   tgl128 = real(tgl, r128)
          !   wgl128 = real(wgl, r128)

          !   call estimate_nearroot_lengths_r128(root_re128, root_im128, &
          !                                       max_len_each_side, len, lenl, lenr)
          !   n_up = nquad*(len - 1_8)
          !   if (n_up > maxpan*nquad) error stop 'line_quad_compress_nearroot_r64: n_up exceeds work array'
          !   call build_nearroot_nodes_r128(root_re128, nquad, tgl128, wgl128, len, lenl, lenr, &
          !                                  t_up128(1:n_up), w_up128(1:n_up))
          !   t_up(1:n_up) = real(t_up128(1:n_up), r64)
          !   w_up(1:n_up) = real(w_up128(1:n_up), r64)
          ! end block

          
          ! block
          !   real(r128) :: t_up128(maxpan*nquad), r_ell128(3,nquad), rp_ell128(3,nquad)
          !   real(r128) :: tgl128(nquad), w_bclag128(nquad), r0j128(3), kdata128(3)
          !   real(r128) :: f_up128(maxpan*nquad)

          !   t_up128(1:n_up) = real(t_up(1:n_up), r128)
          !   r_ell128 = real(r_ell, r128)
          !   rp_ell128 = real(rp_ell, r128)
          !   tgl128 = real(tgl, r128)
          !   w_bclag128 = real(w_bclag, r128)
          !   r0j128 = real(r0(:,j), r128)
          !   kdata128 = real(kdata(:,j), r128)

          !   call eval_integrand_vec_r128(t_up128, n_up, r_ell128, rp_ell128, nquad, &
          !                                tgl128, w_bclag128, invr_kernel_r128_local, r0j128, kdata128, &
          !                                f_up128(1:n_up))
          !   f_up(1:n_up) = real(f_up128(1:n_up), r64)
          ! end block

          ! block
          !   real(r128) :: row128(nquad), rp_ell128(3,nquad)
          !   real(r128) :: wgl128(nquad), wgl_inv128(nquad), w_bclag128(nquad)
          !   real(r128) :: t_up128(maxpan*nquad), w_up128(maxpan*nquad), f_up128(maxpan*nquad)
          !   real(r128) :: sp_ell128(nquad), funvals128(nquad)
          !   real(r128) :: Br128(nquad,maxpan*nquad), integrand0_up128(maxpan*nquad)
          !   real(r128) :: integrand0_compress128(nquad), sp_up128

          !   rp_ell128 = real(rp_ell, r128)
          !   wgl128 = real(wgl, r128)
          !   wgl_inv128 = real(wgl_inv, r128)
          !   w_bclag128 = real(w_bclag, r128)
          !   t_up128(1:n_up) = real(t_up(1:n_up), r128)
          !   w_up128(1:n_up) = real(w_up(1:n_up), r128)
          !   f_up128(1:n_up) = real(f_up(1:n_up), r128)
          !   sp_ell128 = real(sp_ell, r128)
          !   funvals128 = real(funvals(:,ell,j), r128)

          !   do i = 1, n_up
          !     call bary_row_r128(nquad, real(tgl, r128), w_bclag128, t_up128(i), row128)
          !     sp_up128 = sqrt(sum((matmul(rp_ell128, row128))**2))
          !     integrand0_up128(i) = f_up128(i) / sp_up128
          !     Br128(:,i) = w_up128(i) * wgl_inv128 * row128
          !   end do

          !   do k = 1, nquad
          !     integrand0_compress128(k) = sum(Br128(k,1:n_up) * integrand0_up128(1:n_up))
          !     sxbdw(k,ell,j) = real(integrand0_compress128(k) * sp_ell128(k) / funvals128(k) * wgl128(k), r64)
          !   end do
          ! end block

          ! ! pure r64 implementation
          ! call estimate_nearroot_lengths_r64(root_re(j), abs(root_im(j)), &
          !                                   max_len_each_side, len, lenl, lenr)
          ! n_up = nquad*(len - 1_8)
          ! if (n_up > maxpan*nquad) error stop 'line_quad_compress_nearroot_r64: n_up exceeds work array'
          ! call build_nearroot_nodes_r64(root_re(j), nquad, tgl, wgl, len, lenl, lenr, &
          !                               t_up(1:n_up), w_up(1:n_up))
          
          ! call eval_integrand_vec_r64(t_up, n_up, r_ell, rp_ell, nquad, &
          !                             tgl, w_bclag, fptr, r0(:,j), kdata(:,j), &
          !                             f_up(1:n_up))
          
          ! do i = 1, n_up
          !   call bary_row_r64(nquad, tgl, w_bclag, t_up(i), row)
          !   sp_up = sqrt(sum((matmul(rp_ell, row))**2))
          !   integrand0_up(i) = f_up(i) / sp_up
          !   Br(:,i) = w_up(i) * wgl_inv * row
          ! end do
          
          ! do k = 1, nquad
          !   integrand0_compress(k) = sum(Br(k,1:n_up) * integrand0_up(1:n_up))
          !   sxbdw(k,ell,j) = integrand0_compress(k) * sp_ell(k) / funvals(k,ell,j) * wgl(k)
          ! end do

        else
          call run_bisect_panel_r64(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                                    nquad2, tgl2, wgl2, fptr, r0(:,j), kdata(:,j), &
                                    tol, maxpan, t_up, w_up, f_up, n_up)

          do i = 1, n_up
            call bary_row_r64(nquad, tgl, w_bclag, t_up(i), row)
            sp_up = sqrt(sum((matmul(rp_ell, row))**2))
            integrand0_up(i) = f_up(i) / sp_up
            Br(:,i) = w_up(i) * wgl_inv * row
          end do

          do k = 1, nquad
            integrand0_compress(k) = sum(Br(k,1:n_up) * integrand0_up(1:n_up))
            if (abs(funvals(k,ell,j)) > 0.0_r64) then
              sxbdw(k,ell,j) = integrand0_compress(k) * sp_ell(k) / funvals(k,ell,j) * wgl(k)
            else
              sxbdw(k,ell,j) = wgl(k)
            end if
          end do
        end if
      end do
    end do
    deallocate(tgl2, wgl2, Dgl2)
    end if

    if (lq_profile_enabled_r64) then
    call c_f_procpointer(fun, fptr)
    call system_clock(c0, rate)

    wgl_inv = 1.0_r64 / wgl
    call system_clock(c1)
    lq_profile_compress_calls_r64 = lq_profile_compress_calls_r64 + 1_8
    lq_profile_setup_sec_r64 = lq_profile_setup_sec_r64 + &
      real(c1 - c0, r64)/real(rate, r64)

    do ell = 1, sbdnp
      call system_clock(c0)
      idx_start = (ell-1)*nquad + 1
      idx_end   = ell*nquad
      r_ell     = sxbd(:, idx_start:idx_end)
      rp_ell(1,:) = matmul(Dgl, r_ell(1,:))
      rp_ell(2,:) = matmul(Dgl, r_ell(2,:))
      rp_ell(3,:) = matmul(Dgl, r_ell(3,:))
      do k = 1, nquad
        sp_ell(k) = sqrt(sxpbd(1,idx_start+k-1)**2 + &
                         sxpbd(2,idx_start+k-1)**2 + &
                         sxpbd(3,idx_start+k-1)**2)
      end do
      call system_clock(c1)
      lq_profile_geom_sec_r64 = lq_profile_geom_sec_r64 + &
        real(c1 - c0, r64)/real(rate, r64)
      
      do j = 1, m
        if (root_ok(j)) then
          call system_clock(c0)
          block
            real(r128) :: root_re128, root_im128
            real(r128) :: tgl128(nquad), wgl128(nquad)
            real(r128) :: t_up128(maxpan*nquad), w_up128(maxpan*nquad)

            root_re128 = real(root_re(j), r128)
            root_im128 = abs(real(root_im(j), r128))
            tgl128 = real(tgl, r128)
            wgl128 = real(wgl, r128)

            call estimate_nearroot_lengths_r128(root_re128, root_im128, &
                                                max_len_each_side, len, lenl, lenr)
            n_up = nquad*(len - 1_8)
            if (n_up > maxpan*nquad) error stop 'line_quad_compress_nearroot_r64: n_up exceeds work array'
            call build_nearroot_nodes_r128(root_re128, nquad, tgl128, wgl128, len, lenl, lenr, &
                                           t_up128(1:n_up), w_up128(1:n_up))
            t_up(1:n_up) = real(t_up128(1:n_up), r64)
            w_up(1:n_up) = real(w_up128(1:n_up), r64)
          end block
          ! call estimate_nearroot_lengths_r64(root_re(j), abs(root_im(j)), &
          !                                    max_len_each_side, len, lenl, lenr)
          ! n_up = nquad*(len - 1_8)
          ! if (n_up > maxpan*nquad) error stop 'line_quad_compress_nearroot_r64: n_up exceeds work array'
          ! call build_nearroot_nodes_r64(root_re(j), nquad, tgl, wgl, len, lenl, lenr, &
          !                               t_up(1:n_up), w_up(1:n_up))
          call eval_integrand_vec_r64(t_up, n_up, r_ell, rp_ell, nquad, &
                                      tgl, w_bclag, fptr, r0(:,j), kdata(:,j), &
                                      f_up(1:n_up))
          ! call eval_integrand_vec_r64(t_up, n_up, r_ell, rp_ell, nquad, &
          !                             tgl, w_bclag, fptr, r0(:,j), kdata(:,j), &
          !                             f_up(1:n_up))
          call system_clock(c1)
          lq_profile_compress_targets_r64 = lq_profile_compress_targets_r64 + 1_8
          lq_profile_total_nup_r64 = lq_profile_total_nup_r64 + n_up
          lq_profile_max_nup_r64 = max(lq_profile_max_nup_r64, n_up)
          lq_profile_bisect_sec_r64 = lq_profile_bisect_sec_r64 + &
            real(c1 - c0, r64)/real(rate, r64)
        else
          sxbdw(:,ell,j) = wgl
          cycle
        end if

        call system_clock(c0)
        do i = 1, n_up
          call bary_row_r64(nquad, tgl, w_bclag, t_up(i), row)
          sp_up = sqrt(sum((matmul(rp_ell, row))**2))
          integrand0_up(i) = f_up(i) / sp_up
          Br(:,i) = w_up(i) * wgl_inv * row
        end do
        call system_clock(c1)
        lq_profile_br_sec_r64 = lq_profile_br_sec_r64 + &
          real(c1 - c0, r64)/real(rate, r64)
        
        call system_clock(c0)
        do k = 1, nquad
          integrand0_compress(k) = sum(Br(k,1:n_up) * integrand0_up(1:n_up))
        end do
        call system_clock(c1)
        lq_profile_compress_sec_r64 = lq_profile_compress_sec_r64 + &
            real(c1 - c0, r64)/real(rate, r64)
        
        call system_clock(c0)
        do k = 1, nquad
          if (abs(funvals(k,ell,j)) > 0.0_r64) then
            sxbdw(k,ell,j) = integrand0_compress(k) * sp_ell(k) / funvals(k,ell,j) * wgl(k)
          else
            sxbdw(k,ell,j) = wgl(k)
          end if
        end do
        call system_clock(c1)
          lq_profile_weight_sec_r64 = lq_profile_weight_sec_r64 + &
            real(c1 - c0, r64)/real(rate, r64)
      end do
    end do
    end if 
  end subroutine line_quad_compress_nearroot_r64


  subroutine line_quad_compress_nearroot_local_r64(root_ok, nquad, npan, &
                                                  tgl, wgl, legmat, &
                                                  xhat, yhat, zhat, &
                                                  t0, tpan, wpan, &
                                                  target_loc, fun, kdata, &
                                                  integrand0_up, &
                                                  funvals0, sxbdw)
    logical,        intent(in)    :: root_ok
    integer(8),     intent(in)    :: nquad, npan
    real(r64),      intent(in)    :: tgl(nquad), wgl(nquad), legmat(nquad,nquad)
    real(r64),      intent(in)    :: xhat(nquad), yhat(nquad), zhat(nquad)
    real(r64),      intent(in)    :: t0
    real(r64),      intent(in)    :: tpan(nquad,npan), wpan(nquad,npan)
    real(r64),      intent(in)    :: target_loc(3)
    type(c_funptr), value         :: fun
    real(r64),      intent(in)    :: kdata(3)
    real(r64),      intent(in)    :: integrand0_up(nquad*npan)
    real(r64),      intent(inout) :: funvals0(nquad), sxbdw(nquad)
    procedure(kernel_iface_r64), pointer :: fptr
    integer(8) :: kk, k, q, ipan, idx, a
    real(r64) :: p0vec(nquad), u0(nquad)
    real(r64) :: pmat(nquad,nquad), dpmat(nquad,nquad), qmat(nquad,nquad)
    real(r64) :: xdisp0(nquad), ydisp0(nquad), zdisp0(nquad)
    real(r64) :: dx0(nquad), dy0(nquad), dz0(nquad), sp_orig(nquad)
    real(r64) :: r_s(3), tau_s(3), val
    real(r64) :: qrow(nquad), pvec(nquad), row(nquad)
    real(r64) :: integrand0_compress(nquad), wgl_inv(nquad)
    real(r64) :: u, rk, rkp1
    call c_f_procpointer(fun, fptr)
    p0vec = 0.0_r64
    p0vec(1) = 1.0_r64
    if (nquad >= 2_8) p0vec(2) = t0
    do kk = 1, nquad-2
      rk = real(kk,r64)
      rkp1 = real(kk+1_8,r64)
      p0vec(kk+2) = ((2.0_r64*rk + 1.0_r64)*t0*p0vec(kk+1) - rk*p0vec(kk))/rkp1
    end do
    u0 = tgl - t0
    pmat = 0.0_r64
    dpmat = 0.0_r64
    pmat(:,1) = 1.0_r64
    dpmat(:,1) = 0.0_r64
    if (nquad >= 2_8) then
      pmat(:,2) = tgl
      dpmat(:,2) = 1.0_r64
    end if
    do kk = 1, nquad-2
      rk = real(kk,r64)
      rkp1 = real(kk+1_8,r64)
      pmat(:,kk+2) = ((2.0_r64*rk + 1.0_r64)*tgl*pmat(:,kk+1) - rk*pmat(:,kk))/rkp1
      dpmat(:,kk+2) = ((2.0_r64*rk + 1.0_r64)*(pmat(:,kk+1) + tgl*dpmat(:,kk+1)) - rk*dpmat(:,kk))/rkp1
    end do
    dx0 = matmul(dpmat, xhat)
    dy0 = matmul(dpmat, yhat)
    dz0 = matmul(dpmat, zhat)
    sp_orig = sqrt(dx0**2 + dy0**2 + dz0**2)
    if (.not. root_ok) then
      sxbdw = wgl
      return
    end if
    wgl_inv = 1.0_r64/wgl
    integrand0_compress = 0.0_r64
    do ipan = 1, npan
      do q = 1, nquad
        idx = (ipan - 1_8)*nquad + q
        qrow = 0.0_r64
        u = tpan(q,ipan) - t0
        qrow(1) = 0.0_r64
        if (nquad >= 2_8) qrow(2) = u
        do kk = 1, nquad-2
          rk = real(kk,r64)
          rkp1 = real(kk+1_8,r64)
          qrow(kk+2) = ((2.0_r64*rk + 1.0_r64)*(tpan(q,ipan)*qrow(kk+1) + u*p0vec(kk+1)) - rk*qrow(kk))/rkp1
        end do
        pvec = p0vec + qrow
        do k = 1, nquad
          row(k) = 0.0_r64
          do a = 1, nquad
            row(k) = row(k) + pvec(a)*legmat(a,k)
          end do
          integrand0_compress(k) = integrand0_compress(k) + wpan(q,ipan)*wgl_inv(k)*row(k)*integrand0_up(idx)
        end do
      end do
    end do
    do k = 1, nquad
      if (abs(funvals0(k)) == 0.0_r64) error stop 'line_quad_compress_nearroot_local_r64: funvals0(k) is zero'
      sxbdw(k) = integrand0_compress(k)*sp_orig(k)/funvals0(k)*wgl(k)
    end do
  end subroutine line_quad_compress_nearroot_local_r64

  subroutine build_target_nearroot_weights_local_r64(nquad, tgl, wgl, legmat, &
                                                    xj, yj, zj, spj, stauj, &
                                                    xjhat, yjhat, zjhat, &
                                                    rho, xtk, ytk, ztk, &
                                                    fun, kdata, &
                                                    funvals0, weights, &
                                                    troot, accepted, I_local)
    use lq_adaptive_mod, only: line_quad_root_initial_guess_r64, line_quad_root_refine_r64
    implicit none
    integer(8), intent(in) :: nquad
    real(r64), intent(in) :: tgl(nquad), wgl(nquad), legmat(nquad,nquad)
    real(r64), intent(in) :: xj(nquad), yj(nquad), zj(nquad), spj(nquad)
    real(r64), intent(in) :: stauj(3,nquad)
    real(r64), intent(in) :: xjhat(nquad), yjhat(nquad), zjhat(nquad)
    real(r64), intent(in) :: rho, xtk, ytk, ztk
    type(c_funptr), value :: fun
    real(r64), intent(in) :: kdata(3)
    real(r64), intent(inout) :: funvals0(nquad), weights(nquad)
    complex(r64), intent(inout) :: troot
    logical, intent(inout) :: accepted
    real(r64), intent(inout) :: I_local
    integer(8), parameter :: max_len_each_side = 64_8
    procedure(kernel_iface_r64), pointer :: fptr
    complex(r64) :: tinit, zz
    real(r64) :: r0k(3), r_s(3), tau_s(3), val
    real(r64) :: tinit_re, tinit_im, troot_re, troot_im
    real(r64) :: br_init, br_root, t0
    real(r64) :: p0vec(nquad), rbase(3), target_loc(3)
    real(r64) :: rk, rkp1
    integer(8) :: converged, kk
    integer(8) :: len, lenl, lenr, npan
    r0k(1) = xtk
    r0k(2) = ytk
    r0k(3) = ztk
    accepted = .false.
    I_local = 0.0_r64
    troot = cmplx(0.0_r64, 0.0_r64, kind=r64)
    call line_quad_root_initial_guess_r64(tgl, xj, yj, zj, nquad, xtk, ytk, ztk, tinit)
    tinit_re = real(tinit, r64)
    tinit_im = aimag(tinit)
    br_init = bernstein_radius(tinit_re, tinit_im)
    if (br_init < 1.75_r64*rho) then
      converged = 0_8
      call line_quad_root_refine_r64(xjhat, yjhat, zjhat, nquad, xtk, ytk, ztk, tinit, troot, converged)
      troot_re = real(troot, r64)
      troot_im = aimag(troot)
      br_root = bernstein_radius(troot_re, troot_im)
      if (converged == 1_8 .and. br_root < rho) accepted = .true.
    end if
    call c_f_procpointer(fun, fptr)
    do kk = 1, nquad
      r_s(1) = xj(kk)
      r_s(2) = yj(kk)
      r_s(3) = zj(kk)
      tau_s = stauj(:,kk)
      val = 0.0_r64
      call fptr(r_s, tau_s, r0k, kdata, val)
      funvals0(kk) = val*spj(kk)
    end do
    if (.not. accepted) then
      weights = wgl
      return
    end if
    troot_re = real(troot, r64)
    troot_im = aimag(troot)
    if (troot_re < -1.0_r64) then
      t0 = -1.0_r64
    else if (troot_re > 1.0_r64) then
      t0 = 1.0_r64
    else
      t0 = troot_re
    end if
    p0vec = 0.0_r64
    p0vec(1) = 1.0_r64
    if (nquad >= 2_8) p0vec(2) = t0
    do kk = 1, nquad-2
      rk = real(kk, r64)
      rkp1 = real(kk+1_8, r64)
      p0vec(kk+2) = ((2.0_r64*rk + 1.0_r64)*t0*p0vec(kk+1) - rk*p0vec(kk))/rkp1
    end do
    rbase(1) = dot_product(p0vec, xjhat)
    rbase(2) = dot_product(p0vec, yjhat)
    rbase(3) = dot_product(p0vec, zjhat)
    call estimate_nearroot_lengths_r64(troot_re, abs(troot_im), max_len_each_side, len, lenl, lenr)
    npan = len - 1_8
    block
      real(r64) :: tpan(nquad,npan), upan(nquad,npan), wpan(nquad,npan)
      real(r64) :: xpan(nquad,npan), ypan(nquad,npan), zpan(nquad,npan)
      real(r64) :: xdisp(nquad,npan), ydisp(nquad,npan), zdisp(nquad,npan)
      real(r64) :: stangpan(3,nquad,npan), sppan(nquad,npan), dswpan(nquad,npan)
      real(r64) :: funvals_local(nquad,npan), integrand0_up(nquad*npan), kernval(nquad,npan)
      call build_nearroot_panels_local(t0, nquad, npan, lenl, lenr, tgl, wgl, &
                                      xjhat, yjhat, zjhat, rbase, &
                                      tpan, upan, wpan, &
                                      xpan, ypan, zpan, &
                                      xdisp, ydisp, zdisp, &
                                      stangpan, sppan, dswpan)
      target_loc(1) = xtk - rbase(1)
      target_loc(2) = ytk - rbase(2)
      target_loc(3) = ztk - rbase(3)
      funvals_local = 0.0_r64
      integrand0_up = 0.0_r64
      kernval = 0.0_r64
      I_local = 0.0_r64
      call line_kernel_eval_local_r64(nquad, npan, &
                                      xdisp, ydisp, zdisp, stangpan, sppan, dswpan, &
                                      target_loc, fun, kdata, &
                                      funvals_local, integrand0_up, I_local, kernval)
      call line_quad_compress_nearroot_local_r64(accepted, nquad, npan, &
                                                tgl, wgl, legmat, &
                                                xjhat, yjhat, zjhat, &
                                                t0, tpan, wpan, &
                                                target_loc, fun, kdata, &
                                                integrand0_up, &
                                                funvals0, weights)
    end block
  contains
    real(r64) function bernstein_radius(re_t, im_t)
      real(r64), intent(in) :: re_t, im_t
      complex(r64) :: zz
      zz = cmplx(re_t, im_t, kind=r64)
      bernstein_radius = abs(zz + sqrt(zz - 1.0_r64)*sqrt(zz + 1.0_r64))
    end function bernstein_radius
  end subroutine build_target_nearroot_weights_local_r64

  subroutine build_target_nearroot_weights_r64(nquad, tgl, wgl, legmat, &
                                               xj, yj, zj, spj, stauj, &
                                               xjhat, yjhat, zjhat, &
                                               rho, xtk, ytk, ztk, &
                                               fun, kdata, &
                                               funvals0, weights, &
                                               troot, accepted, I_local, &
                                               kernel_id, dgl_in, w_bclag_in, sxpbd_in, n_expa_in, adaptive_fallback)
    use lq_adaptive_mod, only: line_quad_root_initial_guess_r64, line_quad_root_refine_r64
    implicit none
    integer(8), intent(in) :: nquad
    real(r64), intent(in) :: tgl(nquad), wgl(nquad), legmat(nquad,nquad)
    real(r64), intent(in) :: xj(nquad), yj(nquad), zj(nquad), spj(nquad)
    real(r64), intent(in) :: stauj(3,nquad)
    ! Assumed-shape so the wrapper sees the caller's actual array
    ! length: solid-angle callers truncate the Legendre expansion of
    ! the curve to the first n_expa = min(16, nquad) coefficients
    ! before passing in xjhat/yjhat/zjhat (the high-degree tail is
    ! noisy and can mislead the rootfinder). With explicit-shape
    ! xjhat(nquad) the wrapper would over-read into adjacent memory
    ! when the caller's array is shorter than nquad.
    real(r64), intent(in) :: xjhat(:), yjhat(:), zjhat(:)
    real(r64), intent(in) :: rho, xtk, ytk, ztk
    type(c_funptr), value :: fun
    real(r64), intent(in) :: kdata(3)
    real(r64), intent(inout) :: funvals0(nquad), weights(nquad)
    complex(r64), intent(inout) :: troot
    logical, intent(inout) :: accepted
    real(r64), intent(inout) :: I_local
    integer(8), intent(in), optional :: kernel_id
    real(r64), intent(in), optional :: dgl_in(nquad,nquad), w_bclag_in(nquad)
    real(r64), intent(in), optional :: sxpbd_in(3,nquad)
    ! Polynomial size for the rootfinder. Defaults to size(xjhat) so
    ! the caller can simply pass arrays of the truncated length and
    ! skip this argument; pass it explicitly to truncate even further
    ! than the array's own length.
    integer(8), intent(in), optional :: n_expa_in
    ! If .true. and root finding is not accepted, still use the adaptive
    ! bisection fallback inside line_quad_compress_nearroot_r64. This is
    ! the default behavior for this general wrapper. Tests that compare
    ! against the plain-GL fallback can pass adaptive_fallback=.false.
    logical, intent(in), optional :: adaptive_fallback
    real(r64) :: sxbd(3,nquad), sxpbd(3,nquad), stangbd(3,nquad), sspbd(nquad)
    real(r64) :: dgl(nquad,nquad), w_bclag(nquad), bclagmatlr(nquad,2)
    real(r64) :: tgl_work(nquad), wgl_work(nquad)
    real(r64) :: r0(3,1), kdata1(3,1), funvals1(nquad,1,1), sxbdw1(nquad,1,1)
    real(r64) :: root_re(1), root_im(1)
    logical :: root_ok(1)
    complex(r64) :: tinit, zz
    real(r64) :: br_init, br_root
    integer(8) :: converged, k, k_id, n_expa
    logical :: use_adaptive_fallback

    if (present(kernel_id)) then
      k_id = kernel_id
    else
      k_id = KERNEL_INVR
    end if

    if (present(n_expa_in)) then
      n_expa = n_expa_in
    else
      n_expa = int(size(xjhat), kind=8)
    end if
    if (present(adaptive_fallback)) then
      use_adaptive_fallback = adaptive_fallback
    else
      use_adaptive_fallback = .true.
    end if

    accepted = .false.
    I_local = 0.0_r64
    troot = cmplx(0.0_r64, 0.0_r64, kind=r64)

    call line_quad_root_initial_guess_r64(tgl, xj, yj, zj, nquad, xtk, ytk, ztk, tinit)
    zz = cmplx(real(tinit, r64), aimag(tinit), kind=r64)
    br_init = abs(zz + sqrt(zz - 1.0_r64)*sqrt(zz + 1.0_r64))
    ! Bernstein gate: 1.75*rho matches the asvestas/solid-angle path
    ! (see evaluate_solid_angle_integral_r64 / _r128 and
    ! evaluate_line_integral_r128 in solidangle_mod.f90). The local
    ! variant build_target_nearroot_weights_local_r64 also uses 1.75.
    if (br_init < 1.75_r64*rho) then
      converged = 0_8
      call line_quad_root_refine_r64(xjhat(1:n_expa), yjhat(1:n_expa), zjhat(1:n_expa), &
                                     n_expa, xtk, ytk, ztk, tinit, troot, converged)
      zz = cmplx(real(troot, r64), aimag(troot), kind=r64)
      br_root = abs(zz + sqrt(zz - 1.0_r64)*sqrt(zz + 1.0_r64))
      if (converged == 1_8 .and. br_root < rho) accepted = .true.
    end if

    sxbd(1,:) = xj
    sxbd(2,:) = yj
    sxbd(3,:) = zj
    if (present(sxpbd_in)) then
      sxpbd = sxpbd_in
    else
      do k = 1, nquad
        sxpbd(:,k) = stauj(:,k)*spj(k)
      end do
    end if
    stangbd = stauj
    sspbd = spj

    if (present(dgl_in)) then
      dgl = dgl_in
    else
      tgl_work = tgl
      wgl_work = wgl
      call gauss_r64(nquad, tgl_work, wgl_work, dgl)
    end if
    if (present(w_bclag_in)) then
      w_bclag = w_bclag_in
    else
      call bclaginterpweights_r64(nquad, tgl, w_bclag)
    end if
    bclagmatlr = 0.0_r64

    r0(:,1) = [xtk, ytk, ztk]
    kdata1(:,1) = kdata
    funvals1 = 0.0_r64
    sxbdw1 = 0.0_r64
    call line_kernel_eval_r64(1_8, r0, nquad, 1_8, nquad, &
                              sxbd, sxpbd, stangbd, fun, kdata1, funvals1)
    if (.not. accepted .and. .not. use_adaptive_fallback) then
      funvals0 = funvals1(:,1,1)
      weights = wgl
      do k = 1, nquad
        I_local = I_local + funvals0(k)*weights(k)
      end do
      return
    end if

    root_re(1) = real(troot, r64)
    root_im(1) = aimag(troot)
    root_ok(1) = accepted
    call line_quad_compress_nearroot_r64(1_8, r0, nquad, 1_8, nquad, &
                                         sxbd, sxpbd, stangbd, sspbd, &
                                         tgl, wgl, dgl, w_bclag, &
                                         legmat, bclagmatlr, &
                                         fun, kdata1, funvals1, sxbdw1, &
                                         root_re, root_im, root_ok, &
                                         kernel_id=k_id)

    funvals0 = funvals1(:,1,1)
    weights = sxbdw1(:,1,1)
    do k = 1, nquad
      I_local = I_local + funvals0(k)*weights(k)
    end do
  end subroutine build_target_nearroot_weights_r64

  subroutine build_target_nearroot_weights_r128(nquad, tgl, wgl, dgl, w_bclag, &
                                                legmat, &
                                                xj, yj, zj, spj, stauj, &
                                                xjhat, yjhat, zjhat, &
                                                rho, xtk, ytk, ztk, &
                                                fun, kdata, &
                                                funvals0, weights, &
                                                troot, accepted, I_local)
    use lq_adaptive_mod, only: line_quad_root_initial_guess_r128, line_quad_root_refine_r128
    implicit none
    integer(8), intent(in) :: nquad
    real(r128), intent(in) :: tgl(nquad), wgl(nquad), dgl(nquad,nquad), w_bclag(nquad)
    real(r128), intent(in) :: legmat(nquad,nquad)
    real(r128), intent(in) :: xj(nquad), yj(nquad), zj(nquad), spj(nquad)
    real(r128), intent(in) :: stauj(3,nquad)
    real(r128), intent(in) :: xjhat(nquad), yjhat(nquad), zjhat(nquad)
    real(r128), intent(in) :: rho, xtk, ytk, ztk
    procedure(kernel_iface_r128) :: fun
    real(r128), intent(in) :: kdata(3)
    real(r128), intent(inout) :: funvals0(nquad), weights(nquad)
    complex(r128), intent(inout) :: troot
    logical, intent(inout) :: accepted
    real(r128), intent(inout) :: I_local
    real(r128) :: sxbd(3,nquad), sxpbd(3,nquad), stangbd(3,nquad), sspbd(nquad)
    real(r128) :: bclagmatlr(nquad,2)
    real(r128) :: r0(3,1), kdata1(3,1), funvals1(nquad,1,1), sxbdw1(nquad,1,1)
    real(r128) :: root_re(1), root_im(1)
    logical :: root_ok(1)
    complex(r128) :: tinit, zz
    real(r128) :: br_init, br_root
    integer(8) :: converged, k

    accepted = .false.
    I_local = 0.0_r128
    troot = cmplx(0.0_r128, 0.0_r128, kind=r128)

    call line_quad_root_initial_guess_r128(tgl, xj, yj, zj, nquad, xtk, ytk, ztk, tinit)
    zz = cmplx(real(tinit, r128), aimag(tinit), kind=r128)
    br_init = abs(zz + sqrt(zz - 1.0_r128)*sqrt(zz + 1.0_r128))
    if (br_init < 1.75_r128*rho) then
      converged = 0_8
      call line_quad_root_refine_r128(xjhat, yjhat, zjhat, nquad, xtk, ytk, ztk, tinit, troot, converged)
      zz = cmplx(real(troot, r128), aimag(troot), kind=r128)
      br_root = abs(zz + sqrt(zz - 1.0_r128)*sqrt(zz + 1.0_r128))
      if (converged == 1_8 .and. br_root < rho) accepted = .true.
    end if

    sxbd(1,:) = xj
    sxbd(2,:) = yj
    sxbd(3,:) = zj
    do k = 1, nquad
      sxpbd(:,k) = stauj(:,k)*spj(k)
    end do
    stangbd = stauj
    sspbd = spj

    bclagmatlr = 0.0_r128

    r0(:,1) = [xtk, ytk, ztk]
    kdata1(:,1) = kdata
    funvals1 = 0.0_r128
    sxbdw1 = 0.0_r128
    call line_kernel_eval_r128(1_8, r0, nquad, 1_8, nquad, &
                               sxbd, sxpbd, stangbd, fun, kdata1, funvals1)

    funvals0 = funvals1(:,1,1)
    if (accepted) then
      root_re(1) = real(troot, r128)
      root_im(1) = aimag(troot)
      root_ok(1) = .true.
      call line_quad_compress_nearroot_r128(1_8, r0, nquad, 1_8, nquad, &
                                            sxbd, sxpbd, stangbd, sspbd, &
                                            tgl, wgl, dgl, w_bclag, &
                                            legmat, bclagmatlr, &
                                            fun, kdata1, funvals1, sxbdw1, &
                                            root_re, root_im, root_ok)
      weights = sxbdw1(:,1,1)
    else
      weights = wgl
    end if

    do k = 1, nquad
      I_local = I_local + funvals0(k)*weights(k)
    end do
  end subroutine build_target_nearroot_weights_r128

  ! ================================================================
  ! r128 variants (real(16) quad precision)
  ! fun is a pure-Fortran procedure dummy (no bind(C)) throughout.
  ! ================================================================

  subroutine build_nearroot_panel_ends_r128(t_root, len, lenl, lenr, pan_t_end)
    real(r128), intent(in)  :: t_root
    integer(8), intent(in)  :: len, lenl, lenr
    real(r128), intent(out) :: pan_t_end(len)

    real(r128), parameter :: factor = 3.0_r128
    real(r128) :: pan_t_end1(lenl+1), pan_t_end2(lenr+1), factor_inv
    integer(8) :: k

    factor_inv = 1.0_r128/factor
    pan_t_end = 0.0_r128
    pan_t_end1 = 0.0_r128
    pan_t_end2 = 0.0_r128

    if (t_root >= 1.0_r128) then
      pan_t_end1(lenl+1) = 1.0_r128
      do k = lenl, 2, -1
        pan_t_end1(k) = factor_inv*pan_t_end1(k+1)
      end do
      do k = 1, len
        pan_t_end(k) = 1.0_r128 - 2.0_r128*pan_t_end1(lenl+2-k)
      end do
    else if (t_root <= -1.0_r128) then
      pan_t_end2(lenr+1) = 1.0_r128
      do k = lenr, 2, -1
        pan_t_end2(k) = factor_inv*pan_t_end2(k+1)
      end do
      do k = 1, len
        pan_t_end(k) = 2.0_r128*pan_t_end2(k) - 1.0_r128
      end do
    else
      pan_t_end1(lenl+1) = 1.0_r128
      do k = lenl, 2, -1
        pan_t_end1(k) = factor_inv*pan_t_end1(k+1)
      end do
      pan_t_end2(lenr+1) = 1.0_r128
      do k = lenr, 2, -1
        pan_t_end2(k) = factor_inv*pan_t_end2(k+1)
      end do
      do k = lenl + 1, 1, -1
        pan_t_end(k) = -1.0_r128 + (t_root + 1.0_r128)*(1.0_r128 - pan_t_end1(lenl+2-k))
      end do
      do k = 1, lenr + 1
        pan_t_end(lenl+k) = 1.0_r128 + (1.0_r128 - t_root)*(-1.0_r128 + pan_t_end2(k))
      end do
    end if
  end subroutine build_nearroot_panel_ends_r128

  subroutine build_nearroot_nodes_r128(t_root, nquad, tgl, wgl, len, lenl, lenr, &
                                       t_up, w_ref)
    real(r128), intent(in)  :: t_root
    integer(8), intent(in)  :: nquad, len, lenl, lenr
    real(r128), intent(in)  :: tgl(nquad), wgl(nquad)
    real(r128), intent(out) :: t_up(nquad*(len-1)), w_ref(nquad*(len-1))

    real(r128) :: pan_t_end(len), pan_t_mid, pan_t_len
    integer(8) :: j, k, idx0

    call build_nearroot_panel_ends_r128(t_root, len, lenl, lenr, pan_t_end)
    do j = 1, len - 1
      pan_t_mid = 0.5_r128*(pan_t_end(j) + pan_t_end(j+1))
      pan_t_len = 0.5_r128*(pan_t_end(j+1) - pan_t_end(j))
      idx0 = (j - 1_8)*nquad
      do k = 1, nquad
        t_up(idx0+k) = pan_t_mid + tgl(k)*pan_t_len
        w_ref(idx0+k) = wgl(k)*pan_t_len
      end do
    end do
  end subroutine build_nearroot_nodes_r128

  subroutine estimate_nearroot_lengths_r128(t_root, root_imag_abs, max_len_each_side, len, lenl, lenr)
    real(r128), intent(in)  :: t_root, root_imag_abs
    integer(8), intent(in)  :: max_len_each_side
    integer(8), intent(out) :: len, lenl, lenr

    real(r128), parameter :: factor = 3.0_r128
    real(r128) :: dist, target_width
    integer(8) :: levels

    dist = max(root_imag_abs, 1.0e-12_r128)
    target_width = min(2.0_r128, max(2.0_r128*dist, 1.0e-6_r128))
    levels = ceiling(log(2.0_r128/target_width)/log(factor)) + 1_8
    levels = max(1_8, min(levels, max_len_each_side))

    if (t_root >= 1.0_r128) then
      lenl = levels
      lenr = 0_8
      len = lenl + 1_8
    else if (t_root <= -1.0_r128) then
      lenl = 0_8
      lenr = levels
      len = lenr + 1_8
    else
      lenl = levels
      lenr = levels
      len = lenl + lenr + 1_8
    end if
  end subroutine estimate_nearroot_lengths_r128

  ! ----------------------------------------------------------------
  ! line_kernel_eval_r128
  ! ----------------------------------------------------------------
  subroutine line_kernel_eval_r128(m, r0, nbd, sbdnp, nquad, &
                                    sxbd, sxpbd, stangbd,    &
                                    fun, kdata, funvals)
    integer(8),  intent(in)    :: m, nbd, sbdnp, nquad
    real(r128),  intent(in)    :: r0(3,m)
    real(r128),  intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
    procedure(kernel_iface_r128)            :: fun
    real(r128),  intent(in)    :: kdata(3,m)
    real(r128),  intent(inout) :: funvals(nquad,sbdnp,m)

    integer(8) :: ell, j, q, idx
    real(r128) :: r_s(3), tau_s(3), sp, val

    do ell = 1, sbdnp
      do q = 1, nquad
        idx   = (ell-1)*nquad + q
        r_s   = sxbd(:, idx)
        tau_s = stangbd(:, idx)
        sp    = sqrt(sxpbd(1,idx)**2 + sxpbd(2,idx)**2 + sxpbd(3,idx)**2)
        do j = 1, m
          val = 0.0_r128
          call fun(r_s, tau_s, r0(:,j), kdata(:,j), val)
          funvals(q, ell, j) = val * sp
        end do
      end do
    end do

  end subroutine line_kernel_eval_r128

  ! ----------------------------------------------------------------
  ! bary_row_r128
  ! ----------------------------------------------------------------
  subroutine bary_row_r128(nquad, tgl, w_bclag, xi, row)
    integer(8),  intent(in)  :: nquad
    real(r128),  intent(in)  :: tgl(nquad), w_bclag(nquad), xi
    real(r128),  intent(out) :: row(nquad)

    integer(8) :: k
    real(r128) :: eps_hit, denom

    eps_hit = 100.0_r128 * epsilon(1.0_r128)
    row = 0.0_r128
    do k = 1, nquad
      if (abs(xi - tgl(k)) <= eps_hit * max(1.0_r128, abs(tgl(k)))) then
        row(k) = 1.0_r128
        return
      end if
    end do
    do k = 1, nquad
      row(k) = w_bclag(k) / (xi - tgl(k))
    end do
    denom = sum(row)
    row   = row / denom

  end subroutine bary_row_r128

  ! ----------------------------------------------------------------
  ! eval_integrand_vec_r128
  ! ----------------------------------------------------------------
  subroutine eval_integrand_vec_r128(t_pts, np, r_ell, rp_ell, nquad, &
                                      tgl, w_bclag, fun, r0j, kdata, vals)
    integer(8),  intent(in)  :: np, nquad
    real(r128),  intent(in)  :: t_pts(np)
    real(r128),  intent(in)  :: r_ell(3,nquad), rp_ell(3,nquad)
    real(r128),  intent(in)  :: tgl(nquad), w_bclag(nquad)
    procedure(kernel_iface_r128)           :: fun
    real(r128),  intent(in)  :: r0j(3), kdata(3)
    real(r128),  intent(out) :: vals(np)

    integer(8) :: i
    real(r128) :: row(nquad), r_t(3), rp_t(3), sp_t, tau_t(3), val

    do i = 1, np
      call bary_row_r128(nquad, tgl, w_bclag, t_pts(i), row)
      r_t   = matmul(r_ell,  row)
      rp_t  = matmul(rp_ell, row)
      sp_t  = sqrt(rp_t(1)**2 + rp_t(2)**2 + rp_t(3)**2)
      tau_t = rp_t / sp_t
      val   = 0.0_r128
      call fun(r_t, tau_t, r0j, kdata, val)
      vals(i) = val * sp_t
    end do

  end subroutine eval_integrand_vec_r128

  ! ----------------------------------------------------------------
  ! panel_int_err_r128
  ! ----------------------------------------------------------------
  subroutine panel_int_err_r128(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                                 nquad2, tgl2, wgl2, fun, r0j, kdata,     &
                                 aa, bb, ih, err, fhi)
    integer(8),  intent(in)  :: nquad, nquad2
    real(r128),  intent(in)  :: r_ell(3,nquad), rp_ell(3,nquad)
    real(r128),  intent(in)  :: tgl(nquad), wgl(nquad), w_bclag(nquad)
    real(r128),  intent(in)  :: tgl2(nquad2), wgl2(nquad2)
    procedure(kernel_iface_r128)          :: fun
    real(r128),  intent(in)  :: r0j(3), kdata(3), aa, bb
    real(r128),  intent(out) :: ih, err, fhi(nquad)

    real(r128) :: c, h
    real(r128) :: t_fine(nquad),    v_fine(nquad)
    real(r128) :: t_coarse(nquad2), v_coarse(nquad2)
    real(r128) :: il
    integer(8) :: k

    c = 0.5_r128 * (aa + bb)
    h = 0.5_r128 * (bb - aa)

    do k = 1, nquad
      t_fine(k) = c + h * tgl(k)
    end do
    call eval_integrand_vec_r128(t_fine, nquad, r_ell, rp_ell, nquad, &
                                  tgl, w_bclag, fun, r0j, kdata, v_fine)
    ih  = h * sum(wgl * v_fine)
    fhi = v_fine

    do k = 1, nquad2
      t_coarse(k) = c + h * tgl2(k)
    end do
    call eval_integrand_vec_r128(t_coarse, nquad2, r_ell, rp_ell, nquad, &
                                  tgl, w_bclag, fun, r0j, kdata, v_coarse)
    il  = h * sum(wgl2 * v_coarse)
    err = abs(ih - il)

  end subroutine panel_int_err_r128

  ! ----------------------------------------------------------------
  ! run_bisect_panel_r128
  ! ----------------------------------------------------------------
  subroutine run_bisect_panel_r128(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                                    nquad2, tgl2, wgl2, fun, r0j, kdata,    &
                                    tol, maxpan, t_up, w_up, f_up, n_up)
    integer(8),  intent(in)    :: nquad, nquad2, maxpan
    real(r128),  intent(in)    :: r_ell(3,nquad), rp_ell(3,nquad)
    real(r128),  intent(in)    :: tgl(nquad), wgl(nquad), w_bclag(nquad)
    real(r128),  intent(in)    :: tgl2(nquad2), wgl2(nquad2)
    procedure(kernel_iface_r128)              :: fun
    real(r128),  intent(in)    :: r0j(3), kdata(3), tol
    real(r128),  intent(out)   :: t_up(maxpan*nquad), w_up(maxpan*nquad), f_up(maxpan*nquad)
    integer(8),  intent(out)   :: n_up

    integer(8) :: npanels, idx_max, k, j, idx
    real(r128) :: a_pan(maxpan), b_pan(maxpan)
    real(r128) :: intval(maxpan), errp(maxpan), fpan(nquad,maxpan)
    real(r128) :: result, err_tot, scale, target_err, mid, c, h
    real(r128) :: ih_l, err_l, fhi_l(nquad)
    real(r128) :: ih_r, err_r, fhi_r(nquad)

    npanels  = 1_8
    a_pan(1) = -1.0_r128
    b_pan(1) =  1.0_r128
    call panel_int_err_r128(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                             nquad2, tgl2, wgl2, fun, r0j, kdata,     &
                             a_pan(1), b_pan(1), intval(1), errp(1), fpan(:,1))

    do
      result     = sum(intval(1:npanels))
      err_tot    = sum(errp(1:npanels))
      scale      = max(1.0_r128, abs(result))
      target_err = max(tol, 100.0_r128 * epsilon(1.0_r128)) * scale
      if (err_tot <= target_err .or. npanels >= maxpan) exit

      idx_max = 1_8
      do k = 2, npanels
        if (errp(k) > errp(idx_max)) idx_max = k
      end do

      mid = 0.5_r128 * (a_pan(idx_max) + b_pan(idx_max))
      call panel_int_err_r128(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                               nquad2, tgl2, wgl2, fun, r0j, kdata,     &
                               a_pan(idx_max), mid, ih_l, err_l, fhi_l)
      call panel_int_err_r128(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                               nquad2, tgl2, wgl2, fun, r0j, kdata,     &
                               mid, b_pan(idx_max), ih_r, err_r, fhi_r)

      npanels = npanels + 1_8
      a_pan(npanels)   = mid;    b_pan(npanels)   = b_pan(idx_max)
      intval(npanels)  = ih_r;   errp(npanels)    = err_r
      fpan(:,npanels)  = fhi_r
      b_pan(idx_max)   = mid;    intval(idx_max)  = ih_l
      errp(idx_max)    = err_l;  fpan(:,idx_max)  = fhi_l
    end do

    n_up = npanels * nquad
    do j = 1, npanels
      c   = 0.5_r128 * (a_pan(j) + b_pan(j))
      h   = 0.5_r128 * (b_pan(j) - a_pan(j))
      idx = (j-1_8)*nquad
      do k = 1, nquad
        t_up(idx+k) = c + h * tgl(k)
        w_up(idx+k) = h * wgl(k)
        f_up(idx+k) = fpan(k,j)
      end do
    end do

  end subroutine run_bisect_panel_r128

  ! ----------------------------------------------------------------
  ! line_quad_compress_r128
  ! ----------------------------------------------------------------
  subroutine line_quad_compress_r128(m, r0, nbd, sbdnp, nquad,          &
                                      sxbd, sxpbd, stangbd, sspbd,       &
                                      tgl, wgl, Dgl, w_bclag,            &
                                      Legmat, bclagmatlr,                 &
                                      fun, kdata, funvals, sxbdw)
    integer(8),  intent(in)    :: m, nbd, sbdnp, nquad
    real(r128),  intent(in)    :: r0(3,m)
    real(r128),  intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
    real(r128),  intent(in)    :: sspbd(nbd)
    real(r128),  intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
    real(r128),  intent(in)    :: w_bclag(nquad)
    real(r128),  intent(in)    :: Legmat(nquad,nquad), bclagmatlr(nquad,2)
    procedure(kernel_iface_r128)               :: fun
    real(r128),  intent(in)    :: kdata(3,m)
    real(r128),  intent(inout) :: funvals(nquad,sbdnp,m)
    real(r128),  intent(inout) :: sxbdw(nquad,sbdnp,m)

    integer(8), parameter :: maxpan = 128_8
    real(r128), parameter :: tol    = 1.0e-30_r128   ! ~quad machine eps * 1e5

    integer(8) :: ell, j, i, k, idx_start, idx_end, nquad2, n_up
    real(r128) :: r_ell(3,nquad), rp_ell(3,nquad), sp_ell(nquad)
    real(r128) :: row(nquad), wgl_inv(nquad)
    real(r128) :: t_up(maxpan*nquad), w_up(maxpan*nquad), f_up(maxpan*nquad)
    real(r128) :: sp_up, Br(nquad,maxpan*nquad)
    real(r128) :: integrand0_up(maxpan*nquad), integrand0_compress(nquad)
    real(r128), allocatable :: tgl2(:), wgl2(:), Dgl2(:,:)

    nquad2 = max(1_8, nquad / 2_8)
    allocate(tgl2(nquad2), wgl2(nquad2), Dgl2(nquad2,nquad2))
    call gauss_r128(nquad2, tgl2, wgl2, Dgl2)

    wgl_inv = 1.0_r128 / wgl

    do ell = 1, sbdnp
      idx_start = (ell-1)*nquad + 1
      idx_end   = ell*nquad
      r_ell     = sxbd(:, idx_start:idx_end)
      rp_ell(1,:) = matmul(Dgl, r_ell(1,:))
      rp_ell(2,:) = matmul(Dgl, r_ell(2,:))
      rp_ell(3,:) = matmul(Dgl, r_ell(3,:))
      do k = 1, nquad
        sp_ell(k) = sqrt(sxpbd(1,idx_start+k-1)**2 + &
                         sxpbd(2,idx_start+k-1)**2 + &
                         sxpbd(3,idx_start+k-1)**2)
      end do

      do j = 1, m
        call run_bisect_panel_r128(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                                    nquad2, tgl2, wgl2, fun, r0(:,j), kdata(:,j), &
                                    tol, maxpan, t_up, w_up, f_up, n_up)

        do i = 1, n_up
          call bary_row_r128(nquad, tgl, w_bclag, t_up(i), row)
          sp_up = sqrt(sum((matmul(rp_ell, row))**2))
          integrand0_up(i) = f_up(i) / sp_up
          Br(:,i) = w_up(i) * wgl_inv * row
        end do

        do k = 1, nquad
          integrand0_compress(k) = sum(Br(k,1:n_up) * integrand0_up(1:n_up))
        end do

        do k = 1, nquad
          if (abs(funvals(k,ell,j)) > 0.0_r128) then
            sxbdw(k,ell,j) = integrand0_compress(k) * sp_ell(k) / funvals(k,ell,j) * wgl(k)
          else
            sxbdw(k,ell,j) = wgl(k)
          end if
        end do
      end do
    end do

    deallocate(tgl2, wgl2, Dgl2)

  end subroutine line_quad_compress_r128

  ! ----------------------------------------------------------------
  ! line_quad_compress_nearroot_r128
  ! ----------------------------------------------------------------
  subroutine line_quad_compress_nearroot_r128(m, r0, nbd, sbdnp, nquad,  &
                                       sxbd, sxpbd, stangbd, sspbd,      &
                                       tgl, wgl, Dgl, w_bclag,           &
                                       Legmat, bclagmatlr,               &
                                       fun, kdata, funvals, sxbdw,        &
                                       root_re, root_im, root_ok)
    integer(8),  intent(in)    :: m, nbd, sbdnp, nquad
    real(r128),  intent(in)    :: r0(3,m)
    real(r128),  intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
    real(r128),  intent(in)    :: sspbd(nbd)
    real(r128),  intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
    real(r128),  intent(in)    :: w_bclag(nquad)
    real(r128),  intent(in)    :: Legmat(nquad,nquad), bclagmatlr(nquad,2)
    procedure(kernel_iface_r128) :: fun
    real(r128),  intent(in)    :: kdata(3,m)
    real(r128),  intent(inout) :: funvals(nquad,sbdnp,m)
    real(r128),  intent(inout) :: sxbdw(nquad,sbdnp,m)
    real(r128),  intent(in)    :: root_re(m), root_im(m)
    logical,     intent(in)    :: root_ok(m)

    integer(8), parameter :: maxpan = 128_8, max_len_each_side = 12_8
    real(r128), parameter :: tol    = 1.0e-30_r128

    integer(8) :: ell, j, i, k, idx_start, idx_end, nquad2, n_up
    integer(8) :: len, lenl, lenr
    real(r128) :: r_ell(3,nquad), rp_ell(3,nquad), sp_ell(nquad)
    real(r128) :: row(nquad), wgl_inv(nquad)
    real(r128) :: t_up(maxpan*nquad), w_up(maxpan*nquad), f_up(maxpan*nquad)
    real(r128) :: sp_up, Br(nquad,maxpan*nquad)
    real(r128) :: integrand0_up(maxpan*nquad), integrand0_compress(nquad)
    real(r128), allocatable :: tgl2(:), wgl2(:), Dgl2(:,:)

    nquad2 = max(1_8, nquad / 2_8)
    allocate(tgl2(nquad2), wgl2(nquad2), Dgl2(nquad2,nquad2))
    call gauss_r128(nquad2, tgl2, wgl2, Dgl2)

    wgl_inv = 1.0_r128 / wgl

    do ell = 1, sbdnp
      idx_start = (ell-1)*nquad + 1
      idx_end   = ell*nquad
      r_ell     = sxbd(:, idx_start:idx_end)
      rp_ell(1,:) = matmul(Dgl, r_ell(1,:))
      rp_ell(2,:) = matmul(Dgl, r_ell(2,:))
      rp_ell(3,:) = matmul(Dgl, r_ell(3,:))
      do k = 1, nquad
        sp_ell(k) = sqrt(sxpbd(1,idx_start+k-1)**2 + &
                         sxpbd(2,idx_start+k-1)**2 + &
                         sxpbd(3,idx_start+k-1)**2)
      end do

      do j = 1, m
        if (root_ok(j)) then
          call estimate_nearroot_lengths_r128(root_re(j), abs(root_im(j)), &
                                              max_len_each_side, len, lenl, lenr)
          n_up = nquad*(len - 1_8)
          if (n_up > maxpan*nquad) error stop 'line_quad_compress_nearroot_r128: n_up exceeds work array'
          call build_nearroot_nodes_r128(root_re(j), nquad, tgl, wgl, len, lenl, lenr, &
                                         t_up(1:n_up), w_up(1:n_up))
          call eval_integrand_vec_r128(t_up, n_up, r_ell, rp_ell, nquad, &
                                       tgl, w_bclag, fun, r0(:,j), kdata(:,j), &
                                       f_up(1:n_up))
        else
          call run_bisect_panel_r128(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                                     nquad2, tgl2, wgl2, fun, r0(:,j), kdata(:,j), &
                                     tol, maxpan, t_up, w_up, f_up, n_up)
        end if

        do i = 1, n_up
          call bary_row_r128(nquad, tgl, w_bclag, t_up(i), row)
          sp_up = sqrt(sum((matmul(rp_ell, row))**2))
          integrand0_up(i) = f_up(i) / sp_up
          Br(:,i) = w_up(i) * wgl_inv * row
        end do

        do k = 1, nquad
          integrand0_compress(k) = sum(Br(k,1:n_up) * integrand0_up(1:n_up))
        end do

        do k = 1, nquad
          if (abs(funvals(k,ell,j)) > 0.0_r128) then
            sxbdw(k,ell,j) = integrand0_compress(k) * sp_ell(k) / funvals(k,ell,j) * wgl(k)
          else
            sxbdw(k,ell,j) = wgl(k)
          end if
        end do
      end do
    end do

    deallocate(tgl2, wgl2, Dgl2)
  end subroutine line_quad_compress_nearroot_r128

end module lq_kernel_mod
