subroutine lqe_line_integral_r128_mex(m, r0, nbd, sbdnp, nquad,  &
                                     sxbd, sxpbd, stangbd, sspbd,   &
                                     fptr_int, kdata, q_lq64)
  use solidangle_mod, only: lqeli_r128 => evaluate_line_integral_r128, invr_kernel_r128
  integer(8),     intent(in)     :: m, nbd, sbdnp, nquad
  real(8),      intent(in)    :: r0(3,m)
  real(8),      intent(in)    :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
  real(8),      intent(in)    :: sspbd(nbd)
  integer(8),   intent(in)    :: fptr_int
  real(8),  intent(in)        :: kdata(3,m)
  real(8),  intent(inout)     :: q_lq64(m)

  integer :: io_unit
  integer(8) :: j
  real(16) :: r0_128(3,m), sxbd_128(3,nbd), sxpbd_128(3,nbd), stangbd_128(3,nbd)
  real(16) :: sspbd_128(nbd), kdata_128(3,m), q_lq128(m)

  r0_128 = real(r0, 16)
  sxbd_128 = real(sxbd, 16)
  sxpbd_128 = real(sxpbd, 16)
  stangbd_128 = real(stangbd, 16)
  sspbd_128 = real(sspbd, 16)
  kdata_128 = real(kdata, 16)
  q_lq128 = real(q_lq64, 16)

  call lqeli_r128(m, r0_128, nbd, sbdnp, nquad,  &
                  sxbd_128, sxpbd_128, stangbd_128, sspbd_128,   &
                  invr_kernel_r128, kdata_128, q_lq128)
  open (newunit=io_unit, file='lqeli_r128.txt', status='replace', action='write')
  write (io_unit, '(A)') '% j q_lq128'
  do j = 1, m
    write (io_unit, '(I8,1X,ES46.36)') j, q_lq128(j)
  end do
  close (io_unit)
  q_lq64 = real(q_lq128, 8)
end subroutine lqe_line_integral_r128_mex

! ------------------------------------------------------------------
! lqs_eval_moments_funvals_mex
! Thin wrapper around solidangle_mod :: eval_moments_funvals_r64.
! Output funvals_pre(nbd, ncol, m) stores M moments only.
! Here ncol = 2*(order+1), matching qotential's momentsalladapt:
! compute full [N_0..N_{2*order+1} | M_0..M_{2*order+1}], return M.
subroutine lqs_evaluate_solid_angle_integral_fast_mex(m, tx, n, sx, snx, sw, r_vert, nbd, sxbd_in, IalphaAsvestas)
  use solidangle_mod, only: evaluate_solid_angle_integral_fast_driver_r64
  implicit none
  integer(8), intent(in)    :: m, n, nbd
  real(8),    intent(in)    :: tx(3,m), sx(3,n), snx(3,n), sw(n)
  real(8),    intent(in)    :: r_vert(3,3), sxbd_in(3,nbd)
  real(8),    intent(inout) :: IalphaAsvestas(m)

  call evaluate_solid_angle_integral_fast_driver_r64(m, tx, n, sx, snx, sw, r_vert, &
                                                     nbd, sxbd_in, IalphaAsvestas)
end subroutine lqs_evaluate_solid_angle_integral_fast_mex

! ------------------------------------------------------------------
subroutine lqs_eval_moments_funvals_mex(m, tx, nbd, sxbd, nquad, order, ncol, &
                                        funvals_pre)
  use solidangle_mod, only: eval_moments_funvals_r64
  implicit none
  integer(8), intent(in)    :: m, nbd, nquad, order, ncol
  real(8),    intent(in)    :: tx(3, m)
  real(8),    intent(in)    :: sxbd(3, nbd)
  real(8),    intent(inout) :: funvals_pre(nbd, ncol, m)

  real(8), allocatable :: funvals_full(:,:,:)
  integer(8) :: moment_order

  moment_order = 2_8*order + 1_8
  allocate(funvals_full(nbd, 2_8*ncol, m))

  call eval_moments_funvals_r64(m, tx, nbd, sxbd, nquad, moment_order, funvals_full)
  funvals_pre = funvals_full(:, ncol+1_8:2_8*ncol, :)

  deallocate(funvals_full)
end subroutine lqs_eval_moments_funvals_mex
