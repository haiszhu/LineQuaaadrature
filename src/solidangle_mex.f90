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
