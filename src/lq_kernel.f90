! lq_kernel.f90
! Module: line_kernel_eval_r64 only.
! fun is type(c_funptr),value — declared before kdata (matches dummy-list order).
! Caller constructs c_funptr from integer(8) handle before calling.

module lq_kernel_mod
  use iso_c_binding, only: c_funptr, c_f_procpointer, c_double
  use linequaaadrature_mod, only: gauss_r64, gauss_r128
  implicit none

  integer, parameter :: r64  = 8
  integer, parameter :: r128 = 16

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

    real(r64), parameter :: factor = 3.0_r64
    real(r64) :: dist, target_width
    integer(8) :: levels

    dist = max(root_imag_abs, 1.0e-12_r64)
    target_width = min(2.0_r64, max(2.0_r64*dist, 1.0e-6_r64))
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

  subroutine line_quad_compress_nearroot_r64(m, r0, nbd, sbdnp, nquad,  &
                                     sxbd, sxpbd, stangbd, sspbd,       &
                                     tgl, wgl, Dgl, w_bclag,            &
                                     Legmat, bclagmatlr,                 &
                                     fun, kdata, funvals, sxbdw,         &
                                     root_re, root_im, root_ok)
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

    integer(8), parameter :: maxpan = 128_8, max_len_each_side = 12_8
    real(r64),  parameter :: tol    = 1.0e-14_r64

    integer(8) :: ell, j, i, k, idx_start, idx_end, nquad2, n_up
    integer(8) :: len, lenl, lenr
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
        if (root_ok(j)) then
          if (lq_profile_enabled_r64) call system_clock(c0)
          call estimate_nearroot_lengths_r64(root_re(j), abs(root_im(j)), &
                                             max_len_each_side, len, lenl, lenr)
          n_up = nquad*(len - 1_8)
          if (n_up > maxpan*nquad) error stop 'line_quad_compress_nearroot_r64: n_up exceeds work array'
          call build_nearroot_nodes_r64(root_re(j), nquad, tgl, wgl, len, lenl, lenr, &
                                        t_up(1:n_up), w_up(1:n_up))
          call eval_integrand_vec_r64(t_up, n_up, r_ell, rp_ell, nquad, &
                                      tgl, w_bclag, fptr, r0(:,j), kdata(:,j), &
                                      f_up(1:n_up))
          if (lq_profile_enabled_r64) then
            call system_clock(c1)
            lq_profile_compress_targets_r64 = lq_profile_compress_targets_r64 + 1_8
            lq_profile_total_nup_r64 = lq_profile_total_nup_r64 + n_up
            lq_profile_max_nup_r64 = max(lq_profile_max_nup_r64, n_up)
            lq_profile_bisect_sec_r64 = lq_profile_bisect_sec_r64 + &
              real(c1 - c0, r64)/real(rate, r64)
          end if
        else
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
        end if

        if (lq_profile_enabled_r64) call system_clock(c0)
        do i = 1, n_up
          call bary_row_r64(nquad, tgl, w_bclag, t_up(i), row)
          sp_up = sqrt(sum((matmul(rp_ell, row))**2))
          integrand0_up(i) = f_up(i) / sp_up
          Br(:,i) = w_up(i) * wgl_inv * row
        end do
        if (lq_profile_enabled_r64) then
          call system_clock(c1)
          lq_profile_br_sec_r64 = lq_profile_br_sec_r64 + &
            real(c1 - c0, r64)/real(rate, r64)
        end if

        if (lq_profile_enabled_r64) call system_clock(c0)
        do k = 1, nquad
          integrand0_compress(k) = sum(Br(k,1:n_up) * integrand0_up(1:n_up))
        end do
        if (lq_profile_enabled_r64) then
          call system_clock(c1)
          lq_profile_compress_sec_r64 = lq_profile_compress_sec_r64 + &
            real(c1 - c0, r64)/real(rate, r64)
        end if

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
  end subroutine line_quad_compress_nearroot_r64

  ! ================================================================
  ! r128 variants (real(16) quad precision)
  ! fun is a pure-Fortran procedure dummy (no bind(C)) throughout.
  ! ================================================================

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

end module lq_kernel_mod
