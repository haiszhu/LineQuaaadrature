! lq_kernel.f90
! Module: line_kernel_eval_r64 only.
! fun is type(c_funptr),value — declared before kdata (matches dummy-list order).
! Caller constructs c_funptr from integer(8) handle before calling.

module lq_kernel_mod
  use iso_c_binding, only: c_funptr, c_f_procpointer, c_double, c_long_long
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
  !   4. Scalar kernels write integrand0_up(i, 1); vector-output
  !      kernels expect ncol > 1 and write integrand0_up(i, 1:ncol).
  !      Existing INVR/ASVESTAS callers keep ncol=1.
  ! ----------------------------------------------------------------
  integer(8), parameter :: KERNEL_INVR        = 0_8   ! 1/|r-r0|^p, kdata(1)=p in {1,3,5}
  integer(8), parameter :: KERNEL_ASVESTAS    = 1_8   ! Asvestas solid-angle kernel, kdata(1:3)=qhat
  integer(8), parameter :: KERNEL_MOMENTS_MN  = 2_8   ! Line-segment moments; ncol=2*(order+1), packed [N|M]; kdata unused

  real(r64), parameter :: COEFF_I1(50) = [ &
    0.5_r64, -0.125_r64, 0.0625_r64, -0.0390625_r64, 0.02734375_r64, -0.0205078125_r64, &
    0.01611328125_r64, -0.013092041015625_r64, 0.0109100341796875_r64, &
    -0.009273529052734375_r64, 0.0080089569091796875_r64, &
    -0.0070078372955322265625_r64, 0.00619924068450927734375_r64, &
    -0.0055350363254547119140625_r64, 0.00498153269290924072265625_r64, &
    -0.0045145140029489994049072265625_r64, 0.00411617453210055828094482421875_r64, &
    -0.0037731599877588450908660888671875_r64, &
    0.0034752789360936731100082397460938_r64, &
    -0.0032146330158866476267576217651367_r64, 0.002985016371894744224846363067627_r64, &
    -0.0027814925283564662095159292221069_r64, &
    0.0026000908417245227610692381858826_r64, &
    -0.0024375851641167400885024107992649_r64, 0.002291330054269735683192266151309_r64, &
    -0.0021591379357541740091619431041181_r64, &
    0.0020391858282122754530973907094449_r64, &
    -0.0019299437302723321252528876357246_r64, &
    0.0018301190545685908084294624131871_r64, &
    -0.0017386131018401612680079892925278_r64, 0.001654486661428540561491473681599_r64, &
    -0.0015769325991740777226715608527741_r64, &
    0.0015052538446661650989137626321934_r64, &
    -0.0014388455868132460504322731043025_r64, &
    0.0013771807759498212196994613998324_r64, &
    -0.0013197982436185786688786505081727_r64, 0.001266292909417825479599786298382_r64, &
    -0.0012163076629934376317208473655511_r64, &
    0.0011695265990321515689623532361068_r64, &
    -0.0011256693515684458851262649897528_r64, &
    0.0010844863265110637185972552950058_r64, &
    -0.0010457546719928114429330676058984_r64, &
    0.0010092748578535273228307512940647_r64, &
    -0.00097486776042670252773424840903981_r64, &
    0.00094237216841247911014310679540515_r64, &
    -0.00091164264118163740002974461729411_r64, &
    0.00088254766327158514258198681035919_r64, &
    -0.00085496804879434810687629972253547_r64, &
    0.00082879555750472520564539258817214_r64, &
    -0.00080393169077958344947603081052697_r64 ]
  real(r64), parameter :: COEFF_I3(30) = [ &
    0.375_r64, -0.3125_r64, 0.2734375_r64, -0.24609375_r64, 0.2255859375_r64, &
    -0.20947265625_r64, 0.196380615234375_r64, -0.1854705810546875_r64, &
    0.17619705200195312_r64, -0.16818809509277344_r64, 0.16118025779724121_r64, &
    -0.15498101711273193_r64, 0.14944598078727722_r64, -0.14446444809436798_r64, &
    0.13994993409141898_r64, -0.13583375955931842_r64, 0.13206059957155958_r64, &
    -0.12858532063546591_r64, 0.12537068761957926_r64, -0.12238567124768451_r64, &
    0.11960417871932805_r64, -0.11700408787760352_r64, 0.11456650271348678_r64, &
    -0.11227517265921705_r64, 0.11011603472346287_r64, -0.1080768488952506_r64, &
    0.10614690516497827_r64, -0.10431678611040968_r64, 0.10257817300856951_r64, &
    -0.10092368634714097_r64 ]
  real(r64), parameter :: COEFF_I5(50) = [ &
    0.41666666666666667_r64, -0.546875_r64, 0.65625_r64, -0.751953125_r64, &
    0.837890625_r64, -0.91644287109375_r64, 0.98917643229166667_r64, &
    -1.0571823120117188_r64, 1.1212539672851562_r64, -1.1819885571797688_r64, &
    1.2398481369018555_r64, -1.2951985001564026_r64, 1.3483348488807678_r64, &
    -1.3994993409141898_r64, 1.4488934352993965_r64, -1.4966867951443419_r64, &
    1.5430238476255909_r64, -1.5880287098480039_r64, 1.6318089499691268_r64, &
    -1.6744585020705927_r64, 1.716059955538185_r64, -1.7566863749401307_r64, &
    1.7964027625474728_r64, -1.8352672453910479_r64, 1.8733320475176771_r64, &
    -1.9106442929696088_r64, 1.9472466740609806_r64, -1.9831780114990105_r64, &
    2.0184737269428195_r64, -2.053166244124649_r64, 2.0872853312704156_r64, &
    -2.1208583949627249_r64, 2.1539107335855205_r64, -2.1864657569281118_r64, &
    2.2185451773000304_r64, -2.2501691765378595_r64, 2.2813565525120505_r64, &
    -2.3121248481215879_r64, 2.3424904652638978_r64, -2.3724687658610248_r64, &
    2.4020741616913952_r64, -2.4313201945041962_r64, 2.4602196076688454_r64, &
    -2.4887844104258701_r64, 2.5170259356505609_r64, -2.5449548919111762_r64, &
    2.5725814104946672_r64, -2.5999150879811728_r64, 2.6269650248709331_r64, &
    -2.6537398607013483_r64 ]

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
    subroutine kernel_vec_iface_r64(r_s, tau_s, r0j, kdata3, dim, val) bind(C)
      import c_double, c_long_long
      real(c_double), intent(in)    :: r_s(3), tau_s(3), r0j(3), kdata3(3)
      integer(c_long_long), value   :: dim
      real(c_double), intent(inout) :: val(dim)
    end subroutine kernel_vec_iface_r64
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
  ! line_kernel_eval_vec_r64
  ! Vector-output sibling of line_kernel_eval_r64.
  ! fun: c_funptr to subroutine(r_s(3), tau_s(3), r0j(3), kdata3(3),
  !                              dim, val(dim))   bind(C)
  ! dim:  per-source-point output length (== ncol for moments).
  ! funvals(nquad, sbdnp, dim, m): pure kernel value at each node.
  !
  ! Convention note: unlike the scalar line_kernel_eval_r64 (which
  ! returns kernel*sp), this routine does NOT apply the curve-speed
  ! factor. Vector kernels (e.g. moments) are paired with downstream
  ! formulas that supply the line-element factor explicitly via spj,
  ! matching the q_moments_mex convention.
  ! ----------------------------------------------------------------
  subroutine line_kernel_eval_vec_r64(m, r0, nbd, sbdnp, nquad, dim, &
                                       sxbd, sxpbd, stangbd,         &
                                       fun, kdata, funvals)
    use iso_c_binding, only: c_long_long
    integer(8),     intent(in)    :: m, nbd, sbdnp, nquad, dim
    real(r64),      intent(in)    :: r0(3,m)
    real(r64),      intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
    type(c_funptr), value         :: fun
    real(r64),      intent(in)    :: kdata(3,m)
    real(r64),      intent(inout) :: funvals(nquad, sbdnp, dim, m)

    procedure(kernel_vec_iface_r64), pointer :: fptr
    integer(8) :: ell, j, q, idx
    real(r64)  :: r_s(3), tau_s(3)
    real(r64)  :: val(dim)
    integer(c_long_long) :: dim_c

    call c_f_procpointer(fun, fptr)
    dim_c = int(dim, c_long_long)

    do ell = 1, sbdnp
      do q = 1, nquad
        idx   = (ell-1)*nquad + q
        r_s   = sxbd(:, idx)
        tau_s = stangbd(:, idx)
        do j = 1, m
          val = 0.0_r64
          call fptr(r_s, tau_s, r0(:,j), kdata(:,j), dim_c, val)
          funvals(q, ell, :, j) = val
        end do
      end do
    end do
  end subroutine line_kernel_eval_vec_r64

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
                                    stangpan, sppan, dswpan, factor_in)
    real(r64),  intent(in)    :: t0
    integer(8), intent(in)    :: nquad, npan, lenl, lenr
    real(r64),  intent(in)    :: tgl(nquad), wgl(nquad)
    real(r64),  intent(in)    :: xjhat(nquad), yjhat(nquad), zjhat(nquad)
    real(r64),  intent(in)    :: rbase(3)
    real(r64),  intent(inout) :: tpan(nquad,npan), upan(nquad,npan), wpan(nquad,npan)
    real(r64),  intent(inout) :: xpan(nquad,npan), ypan(nquad,npan), zpan(nquad,npan)
    real(r64),  intent(inout) :: xdisp(nquad,npan), ydisp(nquad,npan), zdisp(nquad,npan)
    real(r64),  intent(inout) :: stangpan(3,nquad,npan), sppan(nquad,npan), dswpan(nquad,npan)
    real(r64),  intent(in), optional :: factor_in

    real(r64)  :: factor
    integer(8) :: ipan, k
    real(r64)  :: rhoL, rhoR, ua, ub, h, c
    real(r64)  :: ubreak(npan+1)
    real(r64)  :: p0(nquad)
    real(r64)  :: tt(nquad), uu(nquad)
    real(r64)  :: ppan(nquad,nquad), dppan(nquad,nquad), qpan(nquad,nquad)
    real(r64)  :: dxpan(nquad), dypan(nquad), dzpan(nquad)
    real(r64)  :: rk, rkp1
    factor = 2.0_r64
    if (present(factor_in)) factor = factor_in
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

  subroutine estimate_nearroot_lengths_r64(t_root, root_imag_abs, max_len_each_side, &
                                           len, lenl, lenr, factor_in, coeff_in)
    real(r64),  intent(in)  :: t_root, root_imag_abs
    integer(8), intent(in)  :: max_len_each_side
    integer(8), intent(out) :: len, lenl, lenr
    real(r64),  intent(in), optional :: factor_in, coeff_in

    real(r64) :: factor, coeff
    real(r64) :: dist, target_width
    integer(8) :: levels

    factor = 2.0_r64
    if (present(factor_in)) factor = factor_in
    coeff = 1.0_r64
    if (present(coeff_in)) coeff = coeff_in

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

    levels = ceiling(log(coeff*2.0_r64/target_width)/log(factor)) + 1_8
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
  ! t_up(i), with the barycentric Lagrange interpolation and the
  ! kernel evaluation done in r64.
  !
  ! KERNEL_INVR uses a local-coordinate shift around the real part
  ! of the singularity location, root_re. The bary inputs are
  ! pre-shifted (tgl_local = tgl - root_re, clamped at +/-1 if
  ! root_re is outside the panel), bary weights w_bclag_local are
  ! recomputed from tgl_local, and the geometry is accumulated as
  !   rvec = sum_k row(k) * (r_ell(:,k) - r0j(:))
  ! directly. This keeps every cancellation operand at the same
  ! small magnitude as the singularity and avoids the catastrophic
  ! R_t - r0j subtraction that the previous r128 path was paying
  ! to suppress. root_im is currently unused; reserved so callers
  ! can pass full complex root information without another API
  ! change later.
  !
  ! ASVESTAS / MOMENTS_MN / default branches keep the existing
  ! un-shifted bary-row and r_t/rp_t accumulation. Their precision
  ! profile is unchanged by the new INVR path.
  !
  ! kernel_id selects which inline r64 kernel formula to use. To
  ! add a new kernel:
  !   1. add a KERNEL_<NAME> parameter near the top of lq_kernel_mod.
  !   2. add a `case (KERNEL_<NAME>)` branch in the dispatch below.
  !   3. caller of compress_nearroot_r64 passes kernel_id=KERNEL_<NAME>.
  ! ----------------------------------------------------------------
  subroutine line_quad_BrF_r64(nquad, n_up, ncol, &
                               r_ell, rp_ell, &
                               tgl, wgl, w_bclag, &
                               t_up, w_up, &
                               r0j, root_re, root_im, kdata, kernel_id, &
                               Br, integrand0_up)
    integer(8), intent(in)  :: nquad, n_up, ncol
    real(r64),  intent(in)  :: r_ell(3,nquad), rp_ell(3,nquad)
    real(r64),  intent(in)  :: tgl(nquad), wgl(nquad), w_bclag(nquad)
    real(r64),  intent(in)  :: t_up(n_up), w_up(n_up)
    real(r64),  intent(in)  :: r0j(3), kdata(3)
    real(r64),  intent(in)  :: root_re, root_im
    integer(8), intent(in)  :: kernel_id
    real(r64),  intent(out) :: Br(nquad, n_up)
    real(r64),  intent(out) :: integrand0_up(n_up, ncol)

    integer(8) :: i, k
    real(r64) :: r_t128(3), rp_t128(3), sp_t128, tau_t128(3)
    real(r64) :: rvec128(3), rdist128, rinv128, val128
    real(r64) :: rhat128(3), qhat128(3), qcrossrhat128(3), qdotrhat128
    real(r64) :: eps_hit128, denom128
    real(r64) :: row128(nquad)
    integer(8) :: power_int, hit_idx
    real(r64) :: r0j128(3), kdata128(3)
    real(r64) :: r0norm128
    real(r64) :: tgl128(nquad), wgl128(nquad), wgl_inv128(nquad), w_bclag128(nquad)
    real(r64) :: r_ell128(3,nquad), rp_ell128(3,nquad)
    real(r64) :: t_up128(n_up), w_up128(n_up)
    ! INVR-only local-coordinate arrays (shifted to put the singularity at 0).
    real(r64) :: tgl_local(nquad), t_up_local(n_up), w_bclag_local(nquad)

    kdata128 = real(kdata, r64)
    r0j128 = real(r0j, r64)
    tgl128 = real(tgl, r64)
    wgl128 = real(wgl, r64)
    r_ell128 = real(r_ell, r64)
    rp_ell128 = real(rp_ell, r64)
    t_up128(1:n_up) = real(t_up(1:n_up), r64)
    w_up128(1:n_up) = real(w_up(1:n_up), r64)
    wgl_inv128 = 1.0_r64 / wgl128
    call bclaginterpweights_r64(nquad, tgl128, w_bclag128)
    r0norm128 = sqrt(r0j128(1)**2 + r0j128(2)**2 + r0j128(3)**2)

    ! INVR local-coordinate shift: bring the singularity to t=0 in
    ! the bary inputs. Clamp at +/-1 when root_re is outside [-1,1].
    ! root_im is intentionally unused here (reserved for API parity).
    if (root_re >= 1.0_r64) then
      tgl_local          = tgl128             - 1.0_r64
      t_up_local(1:n_up) = t_up128(1:n_up)    - 1.0_r64
    else if (root_re <= -1.0_r64) then
      tgl_local          = tgl128             + 1.0_r64
      t_up_local(1:n_up) = t_up128(1:n_up)    + 1.0_r64
    else
      tgl_local          = tgl128             - root_re
      t_up_local(1:n_up) = t_up128(1:n_up)    - root_re
    end if
    call bclaginterpweights_r64(nquad, tgl_local, w_bclag_local)

    ! Precompute kernel-specific loop-invariants outside the i loop.
    power_int = 0_8
    qhat128   = 0.0_r64
    select case (kernel_id)
    case (KERNEL_INVR)
      power_int = nint(real(kdata128(1), 8))
    case (KERNEL_ASVESTAS)
      qhat128 = kdata128(1:3)
    end select

    eps_hit128 = 10.0_r64 * epsilon(1.0_r64)
    do i = 1, n_up
      if (kernel_id == KERNEL_INVR) then
        ! ----- INVR: local-shifted bary, direct rvec accumulation -----
        ! Hit threshold matches utils/lqk_line_quad_BrF.m: bare epsilon
        ! is enough since tgl_local has the singularity at 0.
        hit_idx = 0_8
        do k = 1, nquad
          if (abs(t_up_local(i) - tgl_local(k)) <= epsilon(1.0_r64)) then
            hit_idx = k
            exit
          end if
        end do
        if (hit_idx > 0_8) then
          row128 = 0.0_r64
          row128(hit_idx) = 1.0_r64
          rvec128(1) = r_ell128(1, hit_idx) - r0j128(1)
          rvec128(2) = r_ell128(2, hit_idx) - r0j128(2)
          rvec128(3) = r_ell128(3, hit_idx) - r0j128(3)
        else
          denom128 = 0.0_r64
          rvec128  = 0.0_r64
          do k = 1, nquad
            row128(k)  = w_bclag_local(k) / (t_up_local(i) - tgl_local(k))
            denom128   = denom128   + row128(k)
            rvec128(1) = rvec128(1) + (r_ell128(1,k) - r0j128(1)) * row128(k)
            rvec128(2) = rvec128(2) + (r_ell128(2,k) - r0j128(2)) * row128(k)
            rvec128(3) = rvec128(3) + (r_ell128(3,k) - r0j128(3)) * row128(k)
          end do
          row128  = row128  / denom128
          rvec128 = rvec128 / denom128
        end if
        rdist128 = sqrt(rvec128(1)**2 + rvec128(2)**2 + rvec128(3)**2)

        rinv128 = 1.0_r64 / rdist128
        select case (power_int)
        case (1_8)
          val128 = rinv128
        case (3_8)
          val128 = rinv128**3
        case (5_8)
          val128 = rinv128**5
        case default
          val128 = 0.0_r64
        end select
        integrand0_up(i, 1) = real(val128, r64)
      else
      ! ----- Non-INVR (ASVESTAS / MOMENTS_MN / default): existing path -----
      ! inlined bary_row_r64 fused with matmul(r_ell128,row128) and matmul(rp_ell128,row128)
      hit_idx = 0_8
      do k = 1, nquad
        if (abs(t_up128(i) - tgl128(k)) <= eps_hit128 * max(1.0_r64, abs(tgl128(k)))) then
          hit_idx = k
          exit
        end if
      end do
      if (hit_idx > 0_8) then
        row128 = 0.0_r64
        row128(hit_idx) = 1.0_r64
        r_t128  = r_ell128(:,  hit_idx)
        rp_t128 = rp_ell128(:, hit_idx)
      else
        denom128 = 0.0_r64
        r_t128   = 0.0_r64
        rp_t128  = 0.0_r64
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

      ! Geometry shared by non-INVR kernels.
      rvec128(1) = r_t128(1) - r0j128(1)
      rvec128(2) = r_t128(2) - r0j128(2)
      rvec128(3) = r_t128(3) - r0j128(3)
      rdist128   = sqrt(rvec128(1)**2 + rvec128(2)**2 + rvec128(3)**2)

      ! r64 kernel dispatch. See module-top comment for how to add a
      ! new kernel here.
      select case (kernel_id)
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
                 / (rdist128 * (1.0_r64 - qdotrhat128))
        integrand0_up(i, 1) = real(val128, r64)
      case (KERNEL_MOMENTS_MN)
        ! Line-segment moments {N_0..N_order, M_0..M_order} packed
        ! [N|M] in val128(1..ncol), with ncol = 2*(order+1).
        ! Recurrence math: src/qotential_legacy_mod.f90 :: moments_r64.
        block
          integer(8) :: order_loc, kk
          real(r64) :: val128(ncol)        ! shadows outer scalar val128
          real(r64) :: rnorm128, rnorm_inv128, rnorm2_inv128
          real(r64) :: r0dotr128, r0dotr_over_rnorm2_128
          real(r64) :: r0norm2_over_rnorm2_128
          real(r64) :: r0mr_inv128, r0mr_over_rnorm2_128
          real(r64) :: r0normplusrnorm128
          real(r64) :: r0dotr0mr_over_r0mr128, rdotr0mr_over_r0mr128
          real(r64) :: denominator1_128, denominator2_128
          real(r64) :: LMNcommon128, LMcommon128

          order_loc = ncol/2_8 - 1_8

          rnorm128       = sqrt(r_t128(1)**2 + r_t128(2)**2 + r_t128(3)**2)
          rnorm_inv128   = 1.0_r64 / rnorm128
          rnorm2_inv128  = rnorm_inv128 * rnorm_inv128
          r0dotr128      = r0j128(1)*r_t128(1) + r0j128(2)*r_t128(2) &
                         + r0j128(3)*r_t128(3)
          r0dotr_over_rnorm2_128  = r0dotr128 * rnorm2_inv128
          r0norm2_over_rnorm2_128 = (r0norm128 * r0norm128) * rnorm2_inv128
          r0mr_inv128             = 1.0_r64 / rdist128            ! r0mr ≡ rdist128
          r0mr_over_rnorm2_128    = rdist128 * rnorm2_inv128
          r0normplusrnorm128      = r0norm128 + rnorm128

          r0dotr0mr_over_r0mr128 = (r0j128(1)*(r0j128(1)-r_t128(1)) + &
                                    r0j128(2)*(r0j128(2)-r_t128(2)) + &
                                    r0j128(3)*(r0j128(3)-r_t128(3))) * r0mr_inv128
          denominator1_128 = r0norm128 + r0dotr0mr_over_r0mr128
          rdotr0mr_over_r0mr128 = (r_t128(1)*(r0j128(1)-r_t128(1)) + &
                                   r_t128(2)*(r0j128(2)-r_t128(2)) + &
                                   r_t128(3)*(r0j128(3)-r_t128(3))) * r0mr_inv128
          denominator2_128 = rnorm128 + rdotr0mr_over_r0mr128
          LMNcommon128 = r0normplusrnorm128 / (denominator1_128 + denominator2_128)

          ! N0 -> val128(1)
          val128(1) = log((r0normplusrnorm128 + rdist128) * LMNcommon128 * r0mr_inv128) &
                      * rnorm_inv128
          ! N1 -> val128(2) (only if order >= 1)
          if (order_loc >= 1_8) then
            val128(2) = val128(1) * r0dotr_over_rnorm2_128 &
                      + (rdist128 - r0norm128) * rnorm2_inv128
          end if
          ! N recurrence k=2..order_loc -> val128(k+1)
          do kk = 2_8, order_loc
            val128(kk+1) = real(2_8*kk - 1_8, r64) / real(kk, r64)        &
                             * r0dotr_over_rnorm2_128 * val128(kk)          &
                         - real(kk - 1_8,    r64) / real(kk, r64)         &
                             * r0norm2_over_rnorm2_128 * val128(kk-1)       &
                         + 1.0_r64 / real(kk, r64) * r0mr_over_rnorm2_128
          end do

          ! M0 -> val128(order_loc+2)
          LMcommon128 = 1.0_r64 / (r0norm128 * rnorm128 + r0dotr128)
          val128(order_loc + 2_8) = LMNcommon128 * LMcommon128 *                            &
              (((r0normplusrnorm128) * r0mr_inv128 + rnorm128 / r0norm128) * r0mr_inv128 &
               - 1.0_r64 / r0norm128)
          ! M1 -> val128(order_loc+3) (only if order >= 1)
          if (order_loc >= 1_8) then
            val128(order_loc + 3_8) = r0norm128 * val128(order_loc + 2_8) &
                                      / (r0norm128 + rdist128)
          end if
          ! M recurrence k=2..order_loc -> val128(order_loc+2+k); reads N_{k-2} = val128(k-1)
          do kk = 2_8, order_loc
            val128(order_loc + 2_8 + kk) = (r0dotr128 * val128(order_loc + 1_8 + kk) +    &
                                            real(kk - 1_8, r64) * val128(kk - 1_8) -    &
                                            r0mr_inv128) * rnorm2_inv128
          end do

          integrand0_up(i, :) = real(val128, r64)
        end block
      case default
        val128 = 0.0_r64
        integrand0_up(i, 1) = real(val128, r64)
      end select
      end if

      Br(:,i) = real(w_up128(i) * wgl_inv128 * row128, r64)
    end do

  end subroutine line_quad_BrF_r64

  ! ----------------------------------------------------------------
  ! line_quad_BrF_r128 — r128 sibling of line_quad_BrF_r64.
  ! Verbatim port with kind/literal substitutions; no bind(C); kernel
  ! argument is a procedure(kernel_iface_r128) dummy. Same INVR
  ! local-coord shift, same ASVESTAS / MOMENTS_MN dispatch.
  ! ----------------------------------------------------------------
  subroutine line_quad_BrF_r128(nquad, n_up, ncol, &
                                r_ell, rp_ell, &
                                tgl, wgl, w_bclag, &
                                t_up, w_up, &
                                r0j, root_re, root_im, kdata, kernel_id, &
                                Br, integrand0_up)
    integer(8),  intent(in)  :: nquad, n_up, ncol
    real(r128),  intent(in)  :: r_ell(3,nquad), rp_ell(3,nquad)
    real(r128),  intent(in)  :: tgl(nquad), wgl(nquad), w_bclag(nquad)
    real(r128),  intent(in)  :: t_up(n_up), w_up(n_up)
    real(r128),  intent(in)  :: r0j(3), kdata(3)
    real(r128),  intent(in)  :: root_re, root_im
    integer(8),  intent(in)  :: kernel_id
    real(r128),  intent(out) :: Br(nquad, n_up)
    real(r128),  intent(out) :: integrand0_up(n_up, ncol)

    integer(8) :: i, k
    real(r128) :: r_t128(3), rp_t128(3), sp_t128, tau_t128(3)
    real(r128) :: rvec128(3), rdist128, rinv128, val128
    real(r128) :: rhat128(3), qhat128(3), qcrossrhat128(3), qdotrhat128
    real(r128) :: eps_hit128, denom128
    real(r128) :: row128(nquad)
    integer(8) :: power_int, hit_idx
    real(r128) :: r0j128(3), kdata128(3)
    real(r128) :: r0norm128
    real(r128) :: tgl128(nquad), wgl128(nquad), wgl_inv128(nquad), w_bclag128(nquad)
    real(r128) :: r_ell128(3,nquad), rp_ell128(3,nquad)
    real(r128) :: t_up128(n_up), w_up128(n_up)
    real(r128) :: tgl_local(nquad), t_up_local(n_up), w_bclag_local(nquad)

    kdata128 = kdata
    r0j128   = r0j
    tgl128   = tgl
    wgl128   = wgl
    r_ell128 = r_ell
    rp_ell128 = rp_ell
    t_up128(1:n_up) = t_up(1:n_up)
    w_up128(1:n_up) = w_up(1:n_up)
    wgl_inv128 = 1.0_r128 / wgl128
    call bclaginterpweights_r128(nquad, tgl128, w_bclag128)
    r0norm128 = sqrt(r0j128(1)**2 + r0j128(2)**2 + r0j128(3)**2)

    ! INVR local-coordinate shift (mirrors r64 verbatim).
    if (root_re >= 1.0_r128) then
      tgl_local          = tgl128       - 1.0_r128
      t_up_local(1:n_up) = t_up128(1:n_up) - 1.0_r128
    else if (root_re <= -1.0_r128) then
      tgl_local          = tgl128       + 1.0_r128
      t_up_local(1:n_up) = t_up128(1:n_up) + 1.0_r128
    else
      tgl_local          = tgl128       - root_re
      t_up_local(1:n_up) = t_up128(1:n_up) - root_re
    end if
    call bclaginterpweights_r128(nquad, tgl_local, w_bclag_local)

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
      if (kernel_id == KERNEL_INVR) then
        ! ----- INVR: local-shifted bary, direct rvec accumulation -----
        hit_idx = 0_8
        do k = 1, nquad
          if (abs(t_up_local(i) - tgl_local(k)) <= epsilon(1.0_r128)) then
            hit_idx = k
            exit
          end if
        end do
        if (hit_idx > 0_8) then
          row128 = 0.0_r128
          row128(hit_idx) = 1.0_r128
          rvec128(1) = r_ell128(1, hit_idx) - r0j128(1)
          rvec128(2) = r_ell128(2, hit_idx) - r0j128(2)
          rvec128(3) = r_ell128(3, hit_idx) - r0j128(3)
        else
          denom128 = 0.0_r128
          rvec128  = 0.0_r128
          do k = 1, nquad
            row128(k)  = w_bclag_local(k) / (t_up_local(i) - tgl_local(k))
            denom128   = denom128   + row128(k)
            rvec128(1) = rvec128(1) + (r_ell128(1,k) - r0j128(1)) * row128(k)
            rvec128(2) = rvec128(2) + (r_ell128(2,k) - r0j128(2)) * row128(k)
            rvec128(3) = rvec128(3) + (r_ell128(3,k) - r0j128(3)) * row128(k)
          end do
          row128  = row128  / denom128
          rvec128 = rvec128 / denom128
        end if
        rdist128 = sqrt(rvec128(1)**2 + rvec128(2)**2 + rvec128(3)**2)

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
        integrand0_up(i, 1) = val128
      else
      ! ----- Non-INVR (ASVESTAS / MOMENTS_MN / default): existing path -----
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

      rvec128(1) = r_t128(1) - r0j128(1)
      rvec128(2) = r_t128(2) - r0j128(2)
      rvec128(3) = r_t128(3) - r0j128(3)
      rdist128   = sqrt(rvec128(1)**2 + rvec128(2)**2 + rvec128(3)**2)

      select case (kernel_id)
      case (KERNEL_ASVESTAS)
        sp_t128  = sqrt(rp_t128(1)**2 + rp_t128(2)**2 + rp_t128(3)**2)
        tau_t128 = rp_t128 / sp_t128
        rhat128  = rvec128 / rdist128
        qcrossrhat128(1) = qhat128(2)*rhat128(3) - qhat128(3)*rhat128(2)
        qcrossrhat128(2) = qhat128(3)*rhat128(1) - qhat128(1)*rhat128(3)
        qcrossrhat128(3) = qhat128(1)*rhat128(2) - qhat128(2)*rhat128(1)
        qdotrhat128      = qhat128(1)*rhat128(1) + qhat128(2)*rhat128(2) + qhat128(3)*rhat128(3)
        val128 = -(tau_t128(1)*qcrossrhat128(1) + tau_t128(2)*qcrossrhat128(2) + tau_t128(3)*qcrossrhat128(3)) &
                 / (rdist128 * (1.0_r128 - qdotrhat128))
        integrand0_up(i, 1) = val128
      case (KERNEL_MOMENTS_MN)
        block
          integer(8) :: order_loc, kk
          real(r128) :: val128(ncol)
          real(r128) :: rnorm128, rnorm_inv128, rnorm2_inv128
          real(r128) :: r0dotr128, r0dotr_over_rnorm2_128
          real(r128) :: r0norm2_over_rnorm2_128
          real(r128) :: r0mr_inv128, r0mr_over_rnorm2_128
          real(r128) :: r0normplusrnorm128
          real(r128) :: r0dotr0mr_over_r0mr128, rdotr0mr_over_r0mr128
          real(r128) :: denominator1_128, denominator2_128
          real(r128) :: LMNcommon128, LMcommon128

          order_loc = ncol/2_8 - 1_8

          rnorm128       = sqrt(r_t128(1)**2 + r_t128(2)**2 + r_t128(3)**2)
          rnorm_inv128   = 1.0_r128 / rnorm128
          rnorm2_inv128  = rnorm_inv128 * rnorm_inv128
          r0dotr128      = r0j128(1)*r_t128(1) + r0j128(2)*r_t128(2) &
                         + r0j128(3)*r_t128(3)
          r0dotr_over_rnorm2_128  = r0dotr128 * rnorm2_inv128
          r0norm2_over_rnorm2_128 = (r0norm128 * r0norm128) * rnorm2_inv128
          r0mr_inv128             = 1.0_r128 / rdist128
          r0mr_over_rnorm2_128    = rdist128 * rnorm2_inv128
          r0normplusrnorm128      = r0norm128 + rnorm128

          r0dotr0mr_over_r0mr128 = (r0j128(1)*(r0j128(1)-r_t128(1)) + &
                                    r0j128(2)*(r0j128(2)-r_t128(2)) + &
                                    r0j128(3)*(r0j128(3)-r_t128(3))) * r0mr_inv128
          denominator1_128 = r0norm128 + r0dotr0mr_over_r0mr128
          rdotr0mr_over_r0mr128 = (r_t128(1)*(r0j128(1)-r_t128(1)) + &
                                   r_t128(2)*(r0j128(2)-r_t128(2)) + &
                                   r_t128(3)*(r0j128(3)-r_t128(3))) * r0mr_inv128
          denominator2_128 = rnorm128 + rdotr0mr_over_r0mr128
          LMNcommon128 = r0normplusrnorm128 / (denominator1_128 + denominator2_128)

          val128(1) = log((r0normplusrnorm128 + rdist128) * LMNcommon128 * r0mr_inv128) &
                      * rnorm_inv128
          if (order_loc >= 1_8) then
            val128(2) = val128(1) * r0dotr_over_rnorm2_128 &
                      + (rdist128 - r0norm128) * rnorm2_inv128
          end if
          do kk = 2_8, order_loc
            val128(kk+1) = real(2_8*kk - 1_8, r128) / real(kk, r128)        &
                             * r0dotr_over_rnorm2_128 * val128(kk)          &
                         - real(kk - 1_8,    r128) / real(kk, r128)         &
                             * r0norm2_over_rnorm2_128 * val128(kk-1)       &
                         + 1.0_r128 / real(kk, r128) * r0mr_over_rnorm2_128
          end do

          LMcommon128 = 1.0_r128 / (r0norm128 * rnorm128 + r0dotr128)
          val128(order_loc + 2_8) = LMNcommon128 * LMcommon128 *                            &
              (((r0normplusrnorm128) * r0mr_inv128 + rnorm128 / r0norm128) * r0mr_inv128 &
               - 1.0_r128 / r0norm128)
          if (order_loc >= 1_8) then
            val128(order_loc + 3_8) = r0norm128 * val128(order_loc + 2_8) &
                                      / (r0norm128 + rdist128)
          end if
          do kk = 2_8, order_loc
            val128(order_loc + 2_8 + kk) = (r0dotr128 * val128(order_loc + 1_8 + kk) +    &
                                            real(kk - 1_8, r128) * val128(kk - 1_8) -    &
                                            r0mr_inv128) * rnorm2_inv128
          end do

          integrand0_up(i, :) = val128
        end block
      case default
        val128 = 0.0_r128
        integrand0_up(i, 1) = val128
      end select
      end if

      Br(:,i) = w_up128(i) * wgl_inv128 * row128
    end do
  end subroutine line_quad_BrF_r128

  ! ----------------------------------------------------------------
  ! line_quad_compress_nearroot_r64
  !
  ! Fortran twin of utils/lqk_line_quad_compress_nearroot.m. For every
  ! (panel ell, target j) pair: run the rootfinder to locate the
  ! closest complex pre-image of r0 on the curve; if accepted, build
  ! near-root upsampled compression weights via line_quad_BrF_r64 and
  ! the column-loop weight formula; otherwise fall back to plain GL
  ! weights (weights(:,ell,j,q) = wgl).
  !
  ! funvals (intent in): per-source-quad-node kernel values
  !   (kernel * sp for scalar kernels, or pure moments for ncol > 1).
  !   Used as divisor in the weight formula. Caller is expected to
  !   pre-fill via line_kernel_eval_r64 / line_kernel_eval_vec_r64.
  ! weights (intent inout): output compression weights.
  ! root_re, root_im, root_ok (intent inout): rootfinder outputs.
  !   Reset to zero at entry, then populated per (j, ell).
  ! kernel_id: KERNEL_INVR / KERNEL_ASVESTAS / KERNEL_MOMENTS_MN; passed
  !   through to line_quad_BrF_r64 for its inline r64 kernel dispatch.
  !
  ! NOTE: stangbd, sspbd, bclagmatlr appear in the signature for parity
  ! with the .m wrapper (and downstream symmetry); they are not used
  ! by the body. Same parity reason explains keeping fun in the list.
  ! ----------------------------------------------------------------
  subroutine line_quad_compress_nearroot_r64(m, r0, nbd, sbdnp, nquad, ncol,    &
                                             sxbd, sxpbd, stangbd, sspbd,       &
                                             tgl, wgl, Dgl, w_bclag,            &
                                             Legmat, bclagmatlr,                &
                                             fun, kdata, funvals, weights,      &
                                             root_re, root_im, root_ok, kernel_id)
    use lq_adaptive_mod, only: line_quad_root_initial_guess_r64, line_quad_root_refine_r64
    integer(8),     intent(in)    :: m, nbd, sbdnp, nquad, ncol
    real(r64),      intent(in)    :: r0(3,m)
    real(r64),      intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
    real(r64),      intent(in)    :: sspbd(nbd)
    real(r64),      intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
    real(r64),      intent(in)    :: w_bclag(nquad)
    real(r64),      intent(in)    :: Legmat(nquad,nquad), bclagmatlr(nquad,2)
    type(c_funptr), value         :: fun
    real(r64),      intent(in)    :: kdata(3,m)
    real(r64),      intent(in)    :: funvals(nquad,sbdnp,m,ncol)
    real(r64),      intent(inout) :: weights(nquad,sbdnp,m,ncol)
    real(r64),      intent(inout) :: root_re(m,sbdnp), root_im(m,sbdnp)
    integer(8),     intent(inout) :: root_ok(m,sbdnp)
    integer(8),     intent(in)    :: kernel_id

    integer(8), parameter :: maxpan = 128_8
    integer(8), parameter :: max_len_each_side = 24_8

    real(r64)    :: rho
    integer(8)   :: n_expa
    real(r64)    :: r_ell(3,nquad), rp_ell(3,nquad), sp_ell(nquad)
    real(r64)    :: xyz_hat(nquad,3)
    real(r64)    :: t_up(maxpan*nquad), w_up(maxpan*nquad)
    real(r64)    :: Br(nquad, maxpan*nquad)
    real(r64)    :: integrand0_up(maxpan*nquad, ncol)
    real(r64)    :: integrand0_compress(nquad)
    integer(8)   :: ell, j, k, q, idx_start, idx_end
    integer(8)   :: len, lenl, lenr, n_up
    integer(8)   :: ifconv
    complex(r64) :: tinit, troot, zz
    real(r64)    :: bern

    ! root outputs reset on entry (mirrors lqk_line_quad_compress_nearroot.m)
    root_re = 0.0_r64
    root_im = 0.0_r64
    root_ok = 0_8

    n_expa = min(16_8, nquad)
    rho    = 8.0_r64**(16.0_r64 / real(nquad, r64))

    do ell = 1, sbdnp
      idx_start = (ell-1_8)*nquad + 1_8
      idx_end   = ell*nquad
      r_ell     = sxbd(:, idx_start:idx_end)
      rp_ell(1,:) = matmul(Dgl, r_ell(1,:))
      rp_ell(2,:) = matmul(Dgl, r_ell(2,:))
      rp_ell(3,:) = matmul(Dgl, r_ell(3,:))
      do k = 1, nquad
        sp_ell(k) = sqrt(sxpbd(1,idx_start+k-1_8)**2 + &
                         sxpbd(2,idx_start+k-1_8)**2 + &
                         sxpbd(3,idx_start+k-1_8)**2)
      end do

      xyz_hat(1:n_expa,1) = matmul(Legmat(1:n_expa,:), r_ell(1,:))
      xyz_hat(1:n_expa,2) = matmul(Legmat(1:n_expa,:), r_ell(2,:))
      xyz_hat(1:n_expa,3) = matmul(Legmat(1:n_expa,:), r_ell(3,:))

      do j = 1, m
        ! initial complex root guess
        tinit = cmplx(0.0_r64, 0.0_r64, kind=r64)
        call line_quad_root_initial_guess_r64(tgl, r_ell(1,:), r_ell(2,:), r_ell(3,:), &
                                              nquad, r0(1,j), r0(2,j), r0(3,j), tinit)
        zz = cmplx(real(tinit, r64), aimag(tinit), kind=r64)
        bern = abs(zz + sqrt(zz - 1.0_r64) * sqrt(zz + 1.0_r64))

        ! Bernstein gate -> refine -> acceptance
        if (bern < 1.75_r64 * rho) then
          troot  = cmplx(0.0_r64, 0.0_r64, kind=r64)
          ifconv = 0_8
          call line_quad_root_refine_r64(xyz_hat(1:n_expa,1), xyz_hat(1:n_expa,2), xyz_hat(1:n_expa,3), &
                                         n_expa, r0(1,j), r0(2,j), r0(3,j), &
                                         tinit, troot, ifconv)
          zz = cmplx(real(troot, r64), aimag(troot), kind=r64)
          bern = abs(zz + sqrt(zz - 1.0_r64) * sqrt(zz + 1.0_r64))
          if (ifconv == 1_8 .and. bern < rho) then
            root_re(j,ell) = real(troot, r64)
            root_im(j,ell) = aimag(troot)
            root_ok(j,ell) = 1_8
          end if
        end if

        if (root_ok(j,ell) /= 0_8) then
          call estimate_nearroot_lengths_r64(root_re(j,ell), abs(root_im(j,ell)), &
                                             max_len_each_side, len, lenl, lenr)
          n_up = nquad*(len - 1_8)
          if (n_up > maxpan*nquad) error stop 'line_quad_compress_nearroot_r64: n_up exceeds work array'
          call build_nearroot_nodes_r64(root_re(j,ell), nquad, tgl, wgl, len, lenl, lenr, &
                                        t_up(1:n_up), w_up(1:n_up))

          call line_quad_BrF_r64(nquad, n_up, ncol, r_ell, rp_ell,             &
                                 tgl, wgl, w_bclag,                            &
                                 t_up(1:n_up), w_up(1:n_up),                   &
                                 r0(:,j), root_re(j,ell), root_im(j,ell),      &
                                 kdata(:,j), kernel_id,                        &
                                 Br(:, 1:n_up), integrand0_up(1:n_up, 1:ncol))
          do q = 1, ncol
            integrand0_compress = matmul(Br(:, 1:n_up), integrand0_up(1:n_up, q))
            do k = 1, nquad
              weights(k,ell,j,q) = integrand0_compress(k) * sp_ell(k) / funvals(k,ell,j,q) * wgl(k)
            end do
          end do
        else
          do q = 1, ncol
            weights(:,ell,j,q) = wgl
          end do
        end if
      end do
    end do

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
  end subroutine build_target_nearroot_weights_local_r64

  real(r64) function bernstein_radius(re_t, im_t)
    real(r64), intent(in) :: re_t, im_t
    complex(r64) :: zz
    zz = cmplx(re_t, im_t, kind=r64)
    bernstein_radius = abs(zz + sqrt(zz - 1.0_r64)*sqrt(zz + 1.0_r64))
  end function bernstein_radius

  subroutine update_refinement_codes_r64(m, sbdnp, nquad, lens, troot, rfc, rho_in)
    integer(8),   intent(in)    :: m, sbdnp, nquad, lens(3)
    complex(r64), intent(in)    :: troot(m, sbdnp)
    integer(8),   intent(inout) :: rfc(m, sbdnp)
    real(r64),    intent(in), optional :: rho_in

    real(r64)    :: rho, panlen, panmid
    complex(r64) :: z
    integer(8)   :: ell, j, lev, i

    rho = 8.0_r64**(16.0_r64/real(nquad, r64))
    if (present(rho_in)) rho = rho_in

    do ell = 1, sbdnp
      do j = 1, m
        if (rfc(j,ell) > 0_8 .and. &
            bernstein_radius(real(troot(j,ell), r64), aimag(troot(j,ell))) < rho) then
          rfc(j,ell) = 1_8
        else
          rfc(j,ell) = 0_8
        end if
      end do
      do lev = 1, 3
        panlen = 2.0_r64/real(lens(lev), r64)
        panmid = -panlen/2.0_r64 - 1.0_r64
        do i = 1, lens(lev)
          panmid = panmid + panlen
          do j = 1, m
            if (rfc(j,ell) > lev - 1_8) then
              z = (troot(j,ell) - panmid)*2.0_r64/panlen
              if (bernstein_radius(real(z, r64), aimag(z)) < rho) rfc(j,ell) = lev + 1_8
            end if
          end do
        end do
      end do
    end do
  end subroutine update_refinement_codes_r64

  subroutine build_target_nearroot_weights_r64(nquad, ncol, tgl, wgl, legmat, &
                                               xj, yj, zj, spj, stauj, &
                                               xjhat, yjhat, zjhat, &
                                               rho, xtk, ytk, ztk, &
                                               fun, kdata, &
                                               funvals0, weights, &
                                               troot, accepted, I_local, &
                                               kernel_id, dgl_in, w_bclag_in, sxpbd_in, n_expa_in, adaptive_fallback)
    use lq_adaptive_mod, only: line_quad_root_initial_guess_r64, line_quad_root_refine_r64
    implicit none
    integer(8), intent(in) :: nquad, ncol
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
    real(r64), intent(inout) :: funvals0(nquad,ncol), weights(nquad,ncol)
    complex(r64), intent(inout) :: troot
    logical, intent(inout) :: accepted
    real(r64), intent(inout) :: I_local(ncol)
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
    integer(8), parameter :: maxpan = 128_8, max_len_each_side = 12_8
    real(r64) :: r0(3,1), kdata1(3,1)
    real(r64) :: funvals1(nquad,1,1,1), sxbdw1(nquad,1,1,1)
    real(r64) :: funvals_vec(nquad,1,ncol,1)
    real(r64) :: rp_ell(3,nquad)
    real(r64) :: t_up(maxpan*nquad), w_up(maxpan*nquad)
    real(r64) :: Br(nquad,maxpan*nquad), f_up(maxpan*nquad,ncol)
    real(r64) :: integrand0_compress(nquad)
    real(r64) :: root_re(1,1), root_im(1,1)
    integer(8) :: root_ok(1,1)
    complex(r64) :: tinit, zz
    real(r64) :: br_init, br_root
    integer(8) :: converged, k, q, k_id, n_expa
    integer(8) :: len, lenl, lenr, n_up
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

    if (k_id == KERNEL_MOMENTS_MN) then
      funvals_vec = 0.0_r64
      call line_kernel_eval_vec_r64(1_8, r0, nquad, 1_8, nquad, ncol, &
                                    sxbd, sxpbd, stangbd, fun, kdata1, funvals_vec)
      do q = 1, ncol
        funvals0(:,q) = funvals_vec(:,1,q,1)
      end do

      if (.not. accepted) then
        do q = 1, ncol
          weights(:,q) = wgl*spj
          I_local(q) = sum(funvals0(:,q)*weights(:,q))
        end do
        return
      end if

      call estimate_nearroot_lengths_r64(real(troot, r64), abs(aimag(troot)), &
                                         max_len_each_side, len, lenl, lenr)
      n_up = nquad*(len - 1_8)
      if (n_up > maxpan*nquad) error stop 'build_target_nearroot_weights_r64: n_up exceeds work array'
      call build_nearroot_nodes_r64(real(troot, r64), nquad, tgl, wgl, len, lenl, lenr, &
                                    t_up(1:n_up), w_up(1:n_up))
      rp_ell(1,:) = matmul(dgl, xj)
      rp_ell(2,:) = matmul(dgl, yj)
      rp_ell(3,:) = matmul(dgl, zj)
      call line_quad_BrF_r64(nquad, n_up, ncol, sxbd, rp_ell, &
                             tgl, wgl, w_bclag, t_up(1:n_up), w_up(1:n_up), &
                             r0(:,1), real(troot, r64), aimag(troot), kdata, k_id, &
                             Br(:,1:n_up), f_up(1:n_up,1:ncol))

      do q = 1, ncol
        integrand0_compress = matmul(Br(:,1:n_up), f_up(1:n_up,q))
        do k = 1, nquad
          weights(k,q) = integrand0_compress(k) * spj(k) / funvals0(k,q) * wgl(k)
        end do
        I_local(q) = sum(funvals0(:,q)*weights(:,q))
      end do
      return
    end if

    call line_kernel_eval_r64(1_8, r0, nquad, 1_8, nquad, &
                              sxbd, sxpbd, stangbd, fun, kdata1, funvals1(:,:,:,1))

    ! Three flat dispatch cases on the scalar path:
    !   1. plain-GL fallback   (root rejected, adaptive_fallback off)
    !   2. bary-Lagrange       (root accepted -> line_quad_compress_nearroot_r64)
    !   3. adaptive bisection  (root rejected, adaptive_fallback on -> line_quad_compress_r64)
    if (.not. accepted .and. .not. use_adaptive_fallback) then
      funvals0(:,1) = funvals1(:,1,1,1)
      weights(:,1)  = wgl
    else if (accepted) then
      root_re(1,1) = real(troot, r64)
      root_im(1,1) = aimag(troot)
      root_ok(1,1) = 1_8
      call line_quad_compress_nearroot_r64(1_8, r0, nquad, 1_8, nquad, 1_8, &
                                           sxbd, sxpbd, stangbd, sspbd, &
                                           tgl, wgl, dgl, w_bclag, &
                                           legmat, bclagmatlr, &
                                           fun, kdata1, funvals1, sxbdw1, &
                                           root_re, root_im, root_ok, k_id)
      funvals0(:,1) = funvals1(:,1,1,1)
      weights(:,1)  = sxbdw1(:,1,1,1)
    else
      call line_quad_compress_r64(1_8, r0, nquad, 1_8, nquad, &
                                  sxbd, sxpbd, stangbd, sspbd, &
                                  tgl, wgl, dgl, w_bclag, &
                                  legmat, bclagmatlr, &
                                  fun, kdata1, funvals1(:,:,:,1), sxbdw1(:,:,:,1))
      funvals0(:,1) = funvals1(:,1,1,1)
      weights(:,1)  = sxbdw1(:,1,1,1)
    end if

    do k = 1, nquad
      I_local(1) = I_local(1) + funvals0(k,1)*weights(k,1)
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
    real(r128) :: r0(3,1), kdata1(3,1)
    real(r128) :: funvals1(nquad,1,1,1), sxbdw1(nquad,1,1,1)
    real(r128) :: root_re(1,1), root_im(1,1)
    integer(8) :: root_ok(1,1)
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
                               sxbd, sxpbd, stangbd, fun, kdata1, funvals1(:,:,:,1))

    funvals0 = funvals1(:,1,1,1)
    if (accepted) then
      root_re(1,1) = real(troot, r128)
      root_im(1,1) = aimag(troot)
      root_ok(1,1) = 1_8
      call line_quad_compress_nearroot_r128(1_8, r0, nquad, 1_8, nquad, 1_8, &
                                            sxbd, sxpbd, stangbd, sspbd, &
                                            tgl, wgl, dgl, w_bclag, &
                                            legmat, bclagmatlr, &
                                            fun, kdata1, funvals1, sxbdw1, &
                                            root_re, root_im, root_ok, KERNEL_INVR)
      weights = sxbdw1(:,1,1,1)
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

    ! factor matches build_nearroot_panel_ends_r128's panel-shrink ratio (3.0)
    ! so the level estimate is consistent with the actual panel layout.
    real(r128), parameter :: factor = 3.0_r128
    real(r128) :: dist, target_width
    integer(8) :: levels

    ! Distance from the complex root to the real panel [-1,1].
    !
    ! Interior root projection:
    !   closest real point is t_root, so distance is |Im root|.
    !
    ! Exterior root projection:
    !   closest real point is the nearest endpoint, so include real offset.

    if (t_root >= 1.0_r128) then
      dist = sqrt((t_root - 1.0_r128)**2 + root_imag_abs**2)
    else if (t_root <= -1.0_r128) then
      dist = sqrt((t_root + 1.0_r128)**2 + root_imag_abs**2)
    else
      dist = root_imag_abs
    end if

    dist = max(dist, tiny(1.0_r128))
    target_width = min(2.0_r128, 2.0_r128*dist)

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
  ! ----------------------------------------------------------------
  ! line_quad_compress_nearroot_r128
  !
  ! r128 sibling of line_quad_compress_nearroot_r64 — verbatim port
  ! with kind/literal substitutions; no bind(C); kernel argument is
  ! a procedure(kernel_iface_r128) dummy. Same internal rootfinding,
  ! ncol-aware 4D shapes, and bary-Lagrange near-root weight formula
  ! as the r64 sibling. Calls line_quad_BrF_r128 instead of _r64.
  !
  ! See line_quad_compress_nearroot_r64 docblock for behavioural details.
  ! ----------------------------------------------------------------
  subroutine line_quad_compress_nearroot_r128(m, r0, nbd, sbdnp, nquad, ncol,    &
                                              sxbd, sxpbd, stangbd, sspbd,       &
                                              tgl, wgl, Dgl, w_bclag,            &
                                              Legmat, bclagmatlr,                &
                                              fun, kdata, funvals, weights,      &
                                              root_re, root_im, root_ok, kernel_id)
    use lq_adaptive_mod, only: line_quad_root_initial_guess_r128, line_quad_root_refine_r128
    integer(8),  intent(in)    :: m, nbd, sbdnp, nquad, ncol
    real(r128),  intent(in)    :: r0(3,m)
    real(r128),  intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
    real(r128),  intent(in)    :: sspbd(nbd)
    real(r128),  intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
    real(r128),  intent(in)    :: w_bclag(nquad)
    real(r128),  intent(in)    :: Legmat(nquad,nquad), bclagmatlr(nquad,2)
    procedure(kernel_iface_r128) :: fun
    real(r128),  intent(in)    :: kdata(3,m)
    real(r128),  intent(in)    :: funvals(nquad,sbdnp,m,ncol)
    real(r128),  intent(inout) :: weights(nquad,sbdnp,m,ncol)
    real(r128),  intent(inout) :: root_re(m,sbdnp), root_im(m,sbdnp)
    integer(8),  intent(inout) :: root_ok(m,sbdnp)
    integer(8),  intent(in)    :: kernel_id

    integer(8), parameter :: maxpan = 128_8
    integer(8), parameter :: max_len_each_side = 24_8
    real(r128), parameter :: tol               = 1.0e-30_r128

    real(r128)    :: rho
    integer(8)    :: n_expa
    real(r128)    :: r_ell(3,nquad), rp_ell(3,nquad), sp_ell(nquad)
    real(r128)    :: xyz_hat(nquad,3)
    real(r128)    :: t_up(maxpan*nquad), w_up(maxpan*nquad)
    real(r128)    :: Br(nquad, maxpan*nquad)
    real(r128)    :: integrand0_up(maxpan*nquad, ncol)
    real(r128)    :: integrand0_compress(nquad)
    integer(8)    :: ell, j, k, q, idx_start, idx_end, i
    integer(8)    :: len, lenl, lenr, n_up
    integer(8)    :: ifconv
    complex(16)   :: tinit, troot, zz
    real(r128)    :: bern

    ! adaptive-bisection fallback scratch (used in the else branch where
    ! the rootfinder rejects; mirrors line_quad_compress_r128's body).
    integer(8)              :: nquad2
    real(r128), allocatable :: tgl2(:), wgl2(:), Dgl2(:,:)
    real(r128)              :: wgl_inv(nquad)
    real(r128)              :: row(nquad)
    real(r128)              :: f_up(maxpan*nquad)
    real(r128)              :: integrand0_up_ad(maxpan*nquad)
    real(r128)              :: sp_up

    ! root outputs reset on entry (mirrors lqk_line_quad_compress_nearroot.m)
    root_re = 0.0_r128
    root_im = 0.0_r128
    root_ok = 0_8

    n_expa = min(16_8, nquad)
    ! heuristic for nearroot acceptance radius (~2 * r64 rho); validated
    ! against adaptive bisection in test/solid_angle/debug_solid_angle_r128.f90
    ! at r128 noise-floor agreement across the in-range epsilon sweep.
    rho    = 2.0_r128 * (4.0_r128**(16.0_r128 / real(nquad, r128)))

    ! one-time bisection-fallback setup (used in the else branch below)
    nquad2 = max(1_8, nquad / 2_8)
    allocate(tgl2(nquad2), wgl2(nquad2), Dgl2(nquad2,nquad2))
    call gauss_r128(nquad2, tgl2, wgl2, Dgl2)
    wgl_inv = 1.0_r128 / wgl

    do ell = 1, sbdnp
      idx_start = (ell-1_8)*nquad + 1_8
      idx_end   = ell*nquad
      r_ell     = sxbd(:, idx_start:idx_end)
      rp_ell(1,:) = matmul(Dgl, r_ell(1,:))
      rp_ell(2,:) = matmul(Dgl, r_ell(2,:))
      rp_ell(3,:) = matmul(Dgl, r_ell(3,:))
      do k = 1, nquad
        sp_ell(k) = sqrt(sxpbd(1,idx_start+k-1_8)**2 + &
                         sxpbd(2,idx_start+k-1_8)**2 + &
                         sxpbd(3,idx_start+k-1_8)**2)
      end do

      xyz_hat(1:n_expa,1) = matmul(Legmat(1:n_expa,:), r_ell(1,:))
      xyz_hat(1:n_expa,2) = matmul(Legmat(1:n_expa,:), r_ell(2,:))
      xyz_hat(1:n_expa,3) = matmul(Legmat(1:n_expa,:), r_ell(3,:))

      do j = 1, m
        ! initial complex root guess
        tinit = cmplx(0.0_r128, 0.0_r128, kind=r128)
        call line_quad_root_initial_guess_r128(tgl, r_ell(1,:), r_ell(2,:), r_ell(3,:), &
                                               nquad, r0(1,j), r0(2,j), r0(3,j), tinit)
        zz = cmplx(real(tinit, r128), aimag(tinit), kind=r128)
        bern = abs(zz + sqrt(zz - 1.0_r128) * sqrt(zz + 1.0_r128))

        ! Bernstein gate -> refine -> acceptance
        if (bern < 1.75_r128 * rho) then
          troot  = cmplx(0.0_r128, 0.0_r128, kind=r128)
          ifconv = 0_8
          call line_quad_root_refine_r128(xyz_hat(1:n_expa,1), xyz_hat(1:n_expa,2), xyz_hat(1:n_expa,3), &
                                          n_expa, r0(1,j), r0(2,j), r0(3,j), &
                                          tinit, troot, ifconv)
          zz = cmplx(real(troot, r128), aimag(troot), kind=r128)
          bern = abs(zz + sqrt(zz - 1.0_r128) * sqrt(zz + 1.0_r128))
          if (ifconv == 1_8 .and. bern < rho) then
            root_re(j,ell) = real(troot, r128)
            root_im(j,ell) = aimag(troot)
            root_ok(j,ell) = 1_8
          end if
        end if

        if (root_ok(j,ell) /= 0_8) then
          call estimate_nearroot_lengths_r128(root_re(j,ell), abs(root_im(j,ell)), &
                                              max_len_each_side, len, lenl, lenr)
          n_up = nquad*(len - 1_8)
          if (n_up > maxpan*nquad) error stop 'line_quad_compress_nearroot_r128: n_up exceeds work array'
          call build_nearroot_nodes_r128(root_re(j,ell), nquad, tgl, wgl, len, lenl, lenr, &
                                         t_up(1:n_up), w_up(1:n_up))

          call line_quad_BrF_r128(nquad, n_up, ncol, r_ell, rp_ell,             &
                                  tgl, wgl, w_bclag,                            &
                                  t_up(1:n_up), w_up(1:n_up),                   &
                                  r0(:,j), root_re(j,ell), root_im(j,ell),      &
                                  kdata(:,j), kernel_id,                        &
                                  Br(:, 1:n_up), integrand0_up(1:n_up, 1:ncol))

          do q = 1, ncol
            integrand0_compress = matmul(Br(:, 1:n_up), integrand0_up(1:n_up, q))
            do k = 1, nquad
              weights(k,ell,j,q) = integrand0_compress(k) * sp_ell(k) / funvals(k,ell,j,q) * wgl(k)
            end do
          end do
        else
          ! Rootfinder rejected the panel/target: fall back to adaptive
          ! bisection (mirrors line_quad_compress_r128's body) instead of
          ! the previous plain-GL (wgl) fallback. Validated against the
          ! in-range epsilon sweep in test/solid_angle/debug_solid_angle_r128.f90.
          call run_bisect_panel_r128(r_ell, rp_ell, nquad, tgl, wgl, w_bclag, &
                                     nquad2, tgl2, wgl2, fun, r0(:,j), kdata(:,j), &
                                     tol, maxpan, t_up, w_up, f_up, n_up)
          do i = 1_8, n_up
            call bary_row_r128(nquad, tgl, w_bclag, t_up(i), row)
            sp_up = sqrt(sum((matmul(rp_ell, row))**2))
            integrand0_up_ad(i) = f_up(i) / sp_up
            Br(:,i) = w_up(i) * wgl_inv * row
          end do
          do k = 1_8, nquad
            integrand0_compress(k) = sum(Br(k,1:n_up) * integrand0_up_ad(1:n_up))
          end do
          do q = 1_8, ncol
            do k = 1_8, nquad
              weights(k,ell,j,q) = integrand0_compress(k) * sp_ell(k) &
                                    / funvals(k,ell,j,q) * wgl(k)
            end do
          end do
        end if
      end do
    end do

    deallocate(tgl2, wgl2, Dgl2)

  end subroutine line_quad_compress_nearroot_r128


  real(r64) function legendre_eval_r64(n, chat, t) result(v)
    integer(8), intent(in) :: n
    real(r64),  intent(in) :: chat(n), t
    real(r64)  :: pkm2, pkm1, pk
    integer(8) :: k
    pkm2 = 1.0_r64
    v    = chat(1)*pkm2
    if (n > 1_8) then
      pkm1 = t
      v = v + chat(2)*pkm1
      do k = 3, n
        pk = ((2.0_r64*real(k-2,r64) + 1.0_r64)*t*pkm1 &
              - real(k-2,r64)*pkm2)/real(k-1,r64)
        v = v + chat(k)*pk
        pkm2 = pkm1;  pkm1 = pk
      end do
    end if
  end function legendre_eval_r64


  subroutine build_ssq_weights_r64(m, r0, rho_ssq, nbd, sxbd, sbdnp, nquad, &
                                   tgl, wgl, Legmat, w1, w3, w5, &
                                   troot, xroot, yroot, zroot, rfc, rfc_ssq)
    use lq_adaptive_mod, only: line_quad_root_initial_guess_r64, &
                               line_quad_root_refine_r64
    integer(8),   intent(in)    :: m, nbd, sbdnp, nquad
    real(r64),    intent(in)    :: r0(3,m), rho_ssq, sxbd(3,nbd)
    real(r64),    intent(in)    :: tgl(nquad), wgl(nquad), Legmat(nquad,nquad)
    real(r64),    intent(inout) :: w1(nquad,sbdnp,m)
    real(r64),    intent(inout) :: w3(nquad,sbdnp,m)
    real(r64),    intent(inout) :: w5(nquad,sbdnp,m)
    complex(r64), intent(inout) :: troot(m,sbdnp), xroot(m,sbdnp)
    complex(r64), intent(inout) :: yroot(m,sbdnp), zroot(m,sbdnp)
    integer(8),   intent(inout) :: rfc(m,sbdnp), rfc_ssq(m)

    real(r64)    :: r_ell(3,nquad), rho, tr
    real(r64)    :: xh(nquad), yh(nquad), zh(nquad)
    complex(r64) :: tinit, tk
    integer(8)   :: ell, j, i0, n_expa, converged

    rho    = 4.0_r64**(16.0_r64/real(nquad, r64))
    n_expa = min(16_8, nquad)

    do ell = 1, sbdnp
      i0    = (ell-1_8)*nquad
      r_ell = sxbd(:, i0+1:i0+nquad)

      xh = matmul(Legmat, r_ell(1,:))
      yh = matmul(Legmat, r_ell(2,:))
      zh = matmul(Legmat, r_ell(3,:))

      do j = 1, m
        w1(:,ell,j)  = 0.0_r64
        w3(:,ell,j)  = 0.0_r64
        w5(:,ell,j)  = 0.0_r64
        rfc(j,ell)   = 0_8
        troot(j,ell) = (0.0_r64, 0.0_r64)
        xroot(j,ell) = (0.0_r64, 0.0_r64)
        yroot(j,ell) = (0.0_r64, 0.0_r64)
        zroot(j,ell) = (0.0_r64, 0.0_r64)

        call line_quad_root_initial_guess_r64(tgl, r_ell(1,:), r_ell(2,:), &
             r_ell(3,:), nquad, r0(1,j), r0(2,j), r0(3,j), tinit)

        if (bernstein_radius(real(tinit,r64), aimag(tinit)) < 1.5_r64*rho) then
          converged = 0_8
          call line_quad_root_refine_r64(xh(1:n_expa), yh(1:n_expa), zh(1:n_expa), &
               n_expa, r0(1,j), r0(2,j), r0(3,j), tinit, tk, converged)
          troot(j,ell) = tk
          tr = real(tk, r64)
          xroot(j,ell) = cmplx(legendre_eval_r64(n_expa, xh(1:n_expa), tr), 0.0_r64, r64)
          yroot(j,ell) = cmplx(legendre_eval_r64(n_expa, yh(1:n_expa), tr), 0.0_r64, r64)
          zroot(j,ell) = cmplx(legendre_eval_r64(n_expa, zh(1:n_expa), tr), 0.0_r64, r64)
          if (converged == 1_8 .and. &
              bernstein_radius(real(tk,r64), aimag(tk)) < rho) rfc(j,ell) = converged
        end if

        if (rfc(j,ell) == 1_8) then
          block
            real(r64) :: p1(nquad), p3(nquad), p5(nquad)
            real(r64) :: zr, zi, b, c, d, u1, u2, ww, zi_over_w
            real(r64) :: f, f1, f2, arg1, arg2, s, x1, x2, Fs1, Fs2
            real(r64) :: bx2, bx2p, aa, aa2
            integer(8) :: Ns, i, kk, jj
            logical    :: in_cone, outside_interval, use_series

            zr = real(tk, r64)
            zi = aimag(tk)
            b  = -2.0_r64*zr
            c  = zr*zr + zi*zi
            d  = zi*zi
            u1 = sqrt((1.0_r64+zr)*(1.0_r64+zr) + zi*zi)
            u2 = sqrt((1.0_r64-zr)*(1.0_r64-zr) + zi*zi)

            if (4.0_r64*abs(zi) < 1.0_r64 - abs(zr)) then
              Ns = 11_8
              x1 = 1.0_r64 - abs(zr)
              f  = 0.0_r64
              bx2  = zi*zi/(x1*x1)
              bx2p = 1.0_r64
              do i = 1, Ns
                bx2p = bx2p*bx2
                f = f + COEFF_I1(i)*bx2p
              end do
              arg1 = (1.0_r64 - abs(zr))*f
              arg2 = 1.0_r64 + abs(zr) + sqrt((1.0_r64+abs(zr))*(1.0_r64+abs(zr)) + zi*zi)
              p1(1) = log(arg2) - log(arg1)
            else
              arg1 = -1.0_r64 + abs(zr) + sqrt((-1.0_r64+abs(zr))*(-1.0_r64+abs(zr)) + zi*zi)
              arg2 =  1.0_r64 + abs(zr) + sqrt(( 1.0_r64+abs(zr))*( 1.0_r64+abs(zr)) + zi*zi)
              p1(1) = log(arg2) - log(arg1)
            end if
            if (nquad > 1_8) then
              p1(2) = u2 - u1 - b/2.0_r64*p1(1)
              s = 1.0_r64
              do i = 2, nquad-1
                s = -s
                p1(i+1) = (u2 - s*u1 + 0.5_r64*(1.0_r64-2.0_r64*real(i,r64))*b*p1(i) &
                           - real(i-1,r64)*c*p1(i-1))/real(i,r64)
              end do
            end if

            ww = min(abs(1.0_r64+zr), abs(1.0_r64-zr))
            zi = abs(zi)
            zi_over_w = zi/ww
            outside_interval = (abs(zr) > 1.0_r64)

            in_cone    = (zi_over_w < 0.6_r64)
            use_series = (outside_interval .and. in_cone)
            if (.not. use_series) then
              p3(1) = (b+2.0_r64)/(2.0_r64*d*u2) - (b-2.0_r64)/(2.0_r64*d*u1)
            else
              if      (zi_over_w < 0.01_r64) then;  Ns = 4_8
              else if (zi_over_w < 0.1_r64)  then;  Ns = 10_8
              else if (zi_over_w < 0.2_r64)  then;  Ns = 15_8
              else;                                 Ns = 30_8
              end if
              x1 = -1.0_r64 - zr
              x2 =  1.0_r64 - zr
              f1 = 0.0_r64;  f2 = 0.0_r64
              bx2 = zi*zi/(x1*x1);  bx2p = 1.0_r64
              do i = 1, Ns
                bx2p = bx2p*bx2
                f1 = f1 + COEFF_I3(i)*bx2p
              end do
              bx2 = zi*zi/(x2*x2);  bx2p = 1.0_r64
              do i = 1, Ns
                bx2p = bx2p*bx2
                f2 = f2 + COEFF_I3(i)*bx2p
              end do
              Fs1 = abs(x1)/(x1*x1*x1)*(-0.5_r64 + f1)
              Fs2 = abs(x2)/(x2*x2*x2)*(-0.5_r64 + f2)
              p3(1) = Fs2 - Fs1
            end if
            if (nquad > 1_8) then
              p3(2) = 1.0_r64/u1 - 1.0_r64/u2 - b/2.0_r64*p3(1)
              do i = 2, nquad-1
                p3(i+1) = p1(i-1) - b*p3(i) - c*p3(i-1)
              end do
            end if

            in_cone    = (zi_over_w < 0.7_r64)
            use_series = (outside_interval .and. in_cone)
            if (.not. use_series) then
              p5(1) = (2.0_r64+b)/(6.0_r64*d*u2**3) - (-2.0_r64+b)/(6.0_r64*d*u1**3) &
                      + 2.0_r64/(3.0_r64*d)*p3(1)
            else
              if      (zi_over_w < 0.01_r64) then;  Ns = 4_8
              else if (zi_over_w < 0.2_r64)  then;  Ns = 10_8
              else if (zi_over_w < 0.5_r64)  then;  Ns = 24_8
              else if (zi_over_w < 0.6_r64)  then;  Ns = 35_8
              else;                                 Ns = 50_8
              end if
              x1 = -1.0_r64 - zr
              x2 =  1.0_r64 - zr
              f1 = 0.0_r64;  f2 = 0.0_r64
              bx2 = zi*zi/(x1*x1);  bx2p = 1.0_r64
              do i = 1, Ns
                bx2p = bx2p*bx2
                f1 = f1 + COEFF_I5(i)*bx2p
              end do
              bx2 = zi*zi/(x2*x2);  bx2p = 1.0_r64
              do i = 1, Ns
                bx2p = bx2p*bx2
                f2 = f2 + COEFF_I5(i)*bx2p
              end do
              Fs1 = 1.0_r64/(x1*x1*x1*abs(x1))*(-0.25_r64 + f1)
              Fs2 = 1.0_r64/(x2*x2*x2*abs(x2))*(-0.25_r64 + f2)
              p5(1) = Fs2 - Fs1
            end if
            if (nquad > 1_8) then
              p5(2) = 1.0_r64/(3.0_r64*u1*u1*u1) - 1.0_r64/(3.0_r64*u2*u2*u2) &
                      - 0.5_r64*b*p5(1)
              do i = 2, nquad-1
                p5(i+1) = p3(i-1) - b*p5(i) - c*p5(i-1)
              end do
            end if

            do kk = 1, nquad
              do jj = nquad, kk+1, -1
                p1(jj) = p1(jj) - tgl(kk)*p1(jj-1)
                p3(jj) = p3(jj) - tgl(kk)*p3(jj-1)
                p5(jj) = p5(jj) - tgl(kk)*p5(jj-1)
              end do
            end do
            do kk = nquad-1, 1, -1
              do jj = kk+1, nquad
                p1(jj) = p1(jj)/(tgl(jj) - tgl(jj-kk))
                p3(jj) = p3(jj)/(tgl(jj) - tgl(jj-kk))
                p5(jj) = p5(jj)/(tgl(jj) - tgl(jj-kk))
              end do
              do jj = kk, nquad-1
                p1(jj) = p1(jj) - p1(jj+1)
                p3(jj) = p3(jj) - p3(jj+1)
                p5(jj) = p5(jj) - p5(jj+1)
              end do
            end do

            do i = 1, nquad
              aa  = abs(cmplx(tgl(i), 0.0_r64, r64) - tk)
              aa2 = aa*aa
              w1(i,ell,j) = p1(i)*aa
              w3(i,ell,j) = p3(i)*aa*aa2
              w5(i,ell,j) = p5(i)*aa*aa2*aa2
            end do
          end block
        else
          w1(:,ell,j) = wgl
          w3(:,ell,j) = wgl
          w5(:,ell,j) = wgl
        end if

        if (rfc(j,ell) > 0_8 .and. abs(aimag(troot(j,ell))) < rho_ssq &
            .and. abs(real(troot(j,ell), r64)) < 1.0_r64 + rho_ssq) rfc_ssq(j) = 1_8
      end do
    end do
  end subroutine build_ssq_weights_r64

end module lq_kernel_mod
