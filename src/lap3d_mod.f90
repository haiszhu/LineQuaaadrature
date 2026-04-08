! lap3d_mod.f90
! Laplace 3D double-layer potential utilities.
!
! Kernel convention (matches FMM3D / Lap3dDLPmat):
!   K(r_t, r_s) = (r_t - r_s) . n_s / (4*pi * |r_t - r_s|^3)
!
! Public:
!   lap3ddlpmat_r64     — dense DLP matrix K(m,n), entry (i,j) = kernel * w(j)
!   lap3ddlp_direct_r64 — direct O(N^2) DLP with general density sigma
!   lap3ddlp_direct_r128 — same at real(16)

module lap3d_mod
  use linequaaadrature_mod, only: r64, r128
  implicit none
  private
  public :: lap3ddlpmat_r64, lap3ddlp_direct_r64, lap3ddlp_direct_r128, lap3ddlpmat_r128

contains

  ! ------------------------------------------------------------------
  ! lap3ddlpmat_r64
  ! Dense DLP matrix.
  !   A(i,j) = [(r0_i - r_j) . rn_j / (4*pi*|r0_i-r_j|^3)] * w(j)
  !
  ! Inputs:
  !   m         : number of targets
  !   r0(3,m)   : target positions
  !   n         : number of sources
  !   r(3,n)    : source positions
  !   rn(3,n)   : source outward normals
  !   w(n)      : source quadrature weights
  !
  ! Output:
  !   A(m,n)    : DLP matrix (includes w, NOT the density)
  ! ------------------------------------------------------------------
  subroutine lap3ddlpmat_r64(m, r0, n, r, rn, w, A)
    integer(8), intent(in)  :: m, n
    real(r64),  intent(in)  :: r0(3,m), r(3,n), rn(3,n), w(n)
    real(r64),  intent(out) :: A(m,n)

    real(r64), parameter :: pi4inv = 1.0_r64 / (4.0_r64 * &
        3.14159265358979323846264338327950288_r64)
    real(r64), parameter :: rr_min = 1.0e-300_r64

    integer(8) :: i, j
    real(r64)  :: dx, dy, dz, rr, rinv3, dotn

    do j = 1_8, n
      do i = 1_8, m
        dx = r0(1,i) - r(1,j)
        dy = r0(2,i) - r(2,j)
        dz = r0(3,i) - r(3,j)
        rr = dx*dx + dy*dy + dz*dz
        if (rr > rr_min) then
          rinv3    = 1.0_r64 / (sqrt(rr) * rr)
          dotn     = dx*rn(1,j) + dy*rn(2,j) + dz*rn(3,j)
          A(i,j)   = pi4inv * dotn * rinv3 * w(j)
        else
          A(i,j) = 0.0_r64
        end if
      end do
    end do

  end subroutine lap3ddlpmat_r64

  ! ------------------------------------------------------------------
  ! lap3ddlp_direct_r64
  ! Direct O(N^2) DLP evaluation with general density sigma:
  !   u(i) = sum_j K(r0_i, r_j) * sigma(j) * w(j)
  !        = sum_j [(r0_i-r_j).rn_j / (4*pi*|...|^3)] * sigma(j) * w(j)
  !
  ! Replaces Lap3dDLPfmm(t, s, sigma, eps) from MATLAB test.m.
  ! For the test.m use case pass sigma = -ones(n).
  !
  ! Inputs:
  !   m         : number of targets
  !   r0(3,m)   : target positions
  !   n         : number of source nodes
  !   r(3,n)    : source positions
  !   rn(3,n)   : source normals
  !   w(n)      : source quadrature weights
  !   sigma(n)  : density values
  !
  ! Output:
  !   u(m)      : DLP values
  ! ------------------------------------------------------------------
  subroutine lap3ddlp_direct_r64(m, r0, n, r, rn, w, sigma, u)
    integer(8), intent(in)  :: m, n
    real(r64),  intent(in)  :: r0(3,m), r(3,n), rn(3,n), w(n), sigma(n)
    real(r64),  intent(out) :: u(m)

    real(r64), parameter :: pi4inv = 1.0_r64 / (4.0_r64 * &
        3.14159265358979323846264338327950288_r64)
    real(r64), parameter :: rr_min = 1.0e-300_r64

    integer(8) :: i, j
    real(r64)  :: dx, dy, dz, rr, rinv3, dotn, acc

    u = 0.0_r64

!$omp parallel do default(shared) private(i,j,dx,dy,dz,rr,rinv3,dotn,acc) schedule(static)
    do i = 1_8, m
      acc = 0.0_r64
      do j = 1_8, n
        dx = r0(1,i) - r(1,j)
        dy = r0(2,i) - r(2,j)
        dz = r0(3,i) - r(3,j)
        rr = dx*dx + dy*dy + dz*dz
        if (rr > rr_min) then
          rinv3 = 1.0_r64 / (sqrt(rr) * rr)
          dotn  = dx*rn(1,j) + dy*rn(2,j) + dz*rn(3,j)
          acc   = acc + pi4inv * dotn * rinv3 * sigma(j) * w(j)
        end if
      end do
      u(i) = acc
    end do
!$omp end parallel do

  end subroutine lap3ddlp_direct_r64

  ! ------------------------------------------------------------------
  ! lap3ddlp_direct_r128
  ! Same as lap3ddlp_direct_r64 at real(16) precision.
  ! ------------------------------------------------------------------
  subroutine lap3ddlp_direct_r128(m, r0, n, r, rn, w, sigma, u)
    integer(8),  intent(in)  :: m, n
    real(r128),  intent(in)  :: r0(3,m), r(3,n), rn(3,n), w(n), sigma(n)
    real(r128),  intent(out) :: u(m)

    real(r128), parameter :: pi4inv = 1.0_r128 / (4.0_r128 * &
        3.14159265358979323846264338327950288419716939937510_r128)
    real(r128), parameter :: rr_min = 1.0e-300_r128

    integer(8)  :: i, j
    real(r128)  :: dx, dy, dz, rr, rinv3, dotn, acc

    u = 0.0_r128

!$omp parallel do default(shared) private(i,j,dx,dy,dz,rr,rinv3,dotn,acc) schedule(static)
    do i = 1_8, m
      acc = 0.0_r128
      do j = 1_8, n
        dx = r0(1,i) - r(1,j)
        dy = r0(2,i) - r(2,j)
        dz = r0(3,i) - r(3,j)
        rr = dx*dx + dy*dy + dz*dz
        if (rr > rr_min) then
          rinv3 = 1.0_r128 / (sqrt(rr) * rr)
          dotn  = dx*rn(1,j) + dy*rn(2,j) + dz*rn(3,j)
          acc   = acc + pi4inv * dotn * rinv3 * sigma(j) * w(j)
        end if
      end do
      u(i) = acc
    end do
!$omp end parallel do

  end subroutine lap3ddlp_direct_r128

  ! ------------------------------------------------------------------
  ! lap3ddlpmat_r128
  ! Dense DLP matrix at real(16) precision.
  ! ------------------------------------------------------------------
  subroutine lap3ddlpmat_r128(m, r0, n, r, rn, w, A)
    integer(8), intent(in)  :: m, n
    real(r128), intent(in)  :: r0(3,m), r(3,n), rn(3,n), w(n)
    real(r128), intent(out) :: A(m,n)

    real(r128), parameter :: pi4inv = 1.0_r128 / (4.0_r128 * &
        3.14159265358979323846264338327950288419716939937510_r128)
    real(r128), parameter :: rr_min = 1.0e-300_r128

    integer(8) :: i, j
    real(r128) :: dx, dy, dz, rr, rinv3, dotn

    do j = 1_8, n
      do i = 1_8, m
        dx = r0(1,i) - r(1,j)
        dy = r0(2,i) - r(2,j)
        dz = r0(3,i) - r(3,j)
        rr = dx*dx + dy*dy + dz*dz
        if (rr > rr_min) then
          rinv3  = 1.0_r128 / (sqrt(rr) * rr)
          dotn   = dx*rn(1,j) + dy*rn(2,j) + dz*rn(3,j)
          A(i,j) = pi4inv * dotn * rinv3 * w(j)
        else
          A(i,j) = 0.0_r128
        end if
      end do
    end do

  end subroutine lap3ddlpmat_r128

end module lap3d_mod
