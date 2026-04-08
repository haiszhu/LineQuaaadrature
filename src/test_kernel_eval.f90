! test_kernel_eval.f90
! Test: line_kernel_eval_r64 via lq_kernel_mod.
! Caller converts integer(8) handle -> c_funptr before calling module procedure.

subroutine my_kernel(r_s, tau_s, r0j, kdata3, val)
  implicit none
  real(8), intent(in)    :: r_s(3), tau_s(3), r0j(3), kdata3(3)
  real(8), intent(inout) :: val
  real(8) :: rvec(3), rdist
  rvec  = r_s - r0j
  rdist = sqrt(rvec(1)**2 + rvec(2)**2 + rvec(3)**2)
  val   = 1.0d0 / rdist
end subroutine my_kernel

program test_kernel_eval
  use lq_kernel_mod, only: line_kernel_eval_r64, gauss_r64
  use iso_c_binding, only: c_funloc, c_funptr, c_null_funptr
  implicit none

  external :: my_kernel

  integer(8), parameter :: nquad = 4, m = 1, nbd = 4, sbdnp = 1

  real(8) :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
  real(8) :: sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd)
  real(8) :: r0(3,m), kdata(3,m), funvals(nquad,sbdnp,m)
  type(c_funptr) :: cfptr
  integer(8) :: q

  call gauss_r64(nquad, tgl, wgl, Dgl)

  ! Straight line: x from -1 to 1, y=z=0
  do q = 1, nquad
    sxbd(1,q) = tgl(q); sxbd(2,q) = 0.0d0; sxbd(3,q) = 0.0d0
    sxpbd(1,q) = 1.0d0; sxpbd(2,q) = 0.0d0; sxpbd(3,q) = 0.0d0
    stangbd(1,q) = 1.0d0; stangbd(2,q) = 0.0d0; stangbd(3,q) = 0.0d0
  end do

  r0(1,1) = 0.0d0; r0(2,1) = 1.0d0; r0(3,1) = 0.0d0
  kdata(1,1) = 0.0d0; kdata(2,1) = 1.0d0; kdata(3,1) = 0.0d0
  funvals = 0.0d0

  ! Convert function address to c_funptr before passing to module procedure
  cfptr = c_funloc(my_kernel)
  call line_kernel_eval_r64(m, r0, nbd, sbdnp, nquad, &
                             sxbd, sxpbd, stangbd,     &
                             cfptr, kdata, funvals)

  print *, 'funvals(:,1,1) =', funvals(:,1,1)

end program test_kernel_eval
