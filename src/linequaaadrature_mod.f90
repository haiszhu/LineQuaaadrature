module linequaaadrature_mod
  implicit none

  integer, parameter :: r64  = 8
  integer, parameter :: r128 = 16

contains

  ! ------------------------------------------------------------------
  ! gauss
  ! Gauss-Legendre nodes tgl, weights wgl, and differentiation matrix Dgl
  ! on [-1,1].  Algorithm: Newton-Raphson on P_n (cosine initial guess),
  ! D-matrix via Fornberg formula.  Mirrors gauss.m (von Winckel / Bowei Wu).
  ! ------------------------------------------------------------------
  subroutine gauss_r64(n, tgl, wgl, Dgl)
    integer(8), intent(in)    :: n
    real(r64),  intent(inout) :: tgl(n), wgl(n), Dgl(n,n)

    integer(8) :: i, j, k, m
    real(r64)  :: z, z1, p0, p1, p2, pp, pi, tol
    real(r64)  :: a(n)

    pi  = acos(-1.0_r64)
    tol = epsilon(1.0_r64)
    m   = (n + 1_8) / 2_8

    ! ---- nodes and weights (Newton-Raphson, exploit symmetry) ----
    do i = 1, m
      z = cos(pi * (real(i,r64) - 0.25_r64) / (real(n,r64) + 0.5_r64))
      do
        p0 = 1.0_r64;  p1 = z
        do k = 2, n
          p2 = ((2.0_r64*real(k,r64) - 1.0_r64)*z*p1 &
               - (real(k,r64) - 1.0_r64)*p0) / real(k,r64)
          p0 = p1;  p1 = p2
        end do
        pp = real(n,r64) * (z*p1 - p0) / (z*z - 1.0_r64)
        z1 = z
        z  = z1 - p1/pp
        if (abs(z - z1) <= tol * max(1.0_r64, abs(z))) exit
      end do
      tgl(i)       = -z
      tgl(n+1-i)   =  z
      wgl(i)       = 2.0_r64 / ((1.0_r64 - z*z) * pp*pp)
      wgl(n+1-i)   = wgl(i)
    end do

    ! ---- differentiation matrix (Fornberg) ----
    do k = 1, n
      a(k) = 1.0_r64
      do j = 1, n
        if (j /= k) a(k) = a(k) * (tgl(k) - tgl(j))
      end do
    end do
    do k = 1, n
      do j = 1, n
        if (j /= k) then
          Dgl(j,k) = (a(j)/a(k)) / (tgl(j) - tgl(k))
        end if
      end do
      Dgl(k,k) = 0.0_r64
      do j = 1, n
        if (j /= k) Dgl(k,k) = Dgl(k,k) + 1.0_r64/(tgl(k) - tgl(j))
      end do
    end do

  end subroutine gauss_r64

  subroutine gauss_r128(n, tgl, wgl, Dgl)
    integer(8), intent(in)    :: n
    real(r128), intent(inout) :: tgl(n), wgl(n), Dgl(n,n)

    integer(8) :: i, j, k, m
    real(r128) :: z, z1, p0, p1, p2, pp, pi, tol
    real(r128) :: a(n)

    pi  = acos(-1.0_r128)
    tol = epsilon(1.0_r128)
    m   = (n + 1_8) / 2_8

    ! ---- nodes and weights ----
    do i = 1, m
      z = cos(pi * (real(i,r128) - 0.25_r128) / (real(n,r128) + 0.5_r128))
      do
        p0 = 1.0_r128;  p1 = z
        do k = 2, n
          p2 = ((2.0_r128*real(k,r128) - 1.0_r128)*z*p1 &
               - (real(k,r128) - 1.0_r128)*p0) / real(k,r128)
          p0 = p1;  p1 = p2
        end do
        pp = real(n,r128) * (z*p1 - p0) / (z*z - 1.0_r128)
        z1 = z
        z  = z1 - p1/pp
        if (abs(z - z1) <= tol * max(1.0_r128, abs(z))) exit
      end do
      tgl(i)       = -z
      tgl(n+1-i)   =  z
      wgl(i)       = 2.0_r128 / ((1.0_r128 - z*z) * pp*pp)
      wgl(n+1-i)   = wgl(i)
    end do

    ! ---- differentiation matrix (Fornberg) ----
    do k = 1, n
      a(k) = 1.0_r128
      do j = 1, n
        if (j /= k) a(k) = a(k) * (tgl(k) - tgl(j))
      end do
    end do
    do k = 1, n
      do j = 1, n
        if (j /= k) then
          Dgl(j,k) = (a(j)/a(k)) / (tgl(j) - tgl(k))
        end if
      end do
      Dgl(k,k) = 0.0_r128
      do j = 1, n
        if (j /= k) Dgl(k,k) = Dgl(k,k) + 1.0_r128/(tgl(k) - tgl(j))
      end do
    end do

  end subroutine gauss_r128

  ! ------------------------------------------------------------------
  ! bclaginterpweights
  ! Barycentric Lagrange interpolation weights for nodes tgl:
  !   w_bclag(k) = 1 / prod_{j /= k} (tgl(k) - tgl(j))
  ! Used in type-II barycentric formula:
  !   L_k(xi) = (w_bclag(k)/(xi-tgl(k))) / sum_j(w_bclag(j)/(xi-tgl(j)))
  ! ------------------------------------------------------------------
  subroutine bclaginterpweights_r64(n, tgl, w_bclag)
    integer(8), intent(in)    :: n
    real(r64),  intent(in)    :: tgl(n)
    real(r64),  intent(inout) :: w_bclag(n)

    integer(8) :: j, k
    real(r64)  :: prod

    do k = 1, n
      prod = 1.0_r64
      do j = 1, n
        if (j /= k) prod = prod * (tgl(k) - tgl(j))
      end do
      w_bclag(k) = 1.0_r64 / prod
    end do

  end subroutine bclaginterpweights_r64

  subroutine bclaginterpweights_r128(n, tgl, w_bclag)
    integer(8), intent(in)    :: n
    real(r128), intent(in)    :: tgl(n)
    real(r128), intent(inout) :: w_bclag(n)

    integer(8) :: j, k
    real(r128) :: prod

    do k = 1, n
      prod = 1.0_r128
      do j = 1, n
        if (j /= k) prod = prod * (tgl(k) - tgl(j))
      end do
      w_bclag(k) = 1.0_r128 / prod
    end do

  end subroutine bclaginterpweights_r128

  ! ------------------------------------------------------------------
  ! legeexps
  ! GL nodes tgl, weights wgl, and Legendre expansion matrices.
  !   itype=0 : nodes only
  !   itype=1 : nodes + weights
  !   itype=2 : nodes + weights + Legmat (u) + v
  !
  ! v(i,j)      = P_{j-1}(tgl(i))          coefs -> values  (Vandermonde)
  ! Legmat(i,j) = ((2i-1)/2)*P_{i-1}(tgl(j))*wgl(j)  values -> coefs (= u)
  !
  ! Mirrors legeexps.f (Rokhlin) but in F90 and precision-generic.
  ! ------------------------------------------------------------------
  subroutine legeexps_r64(itype, n, tgl, Legmat, v, wgl)
    integer(8), intent(in)    :: itype, n
    real(r64),  intent(inout) :: tgl(n), Legmat(n,n), v(n,n), wgl(n)

    integer(8) :: i, j
    real(r64)  :: pols(n), d, Dgl_tmp(n,n)

    call gauss_r64(n, tgl, wgl, Dgl_tmp)
    if (itype < 2) return

    ! build v: v(i,j) = P_{j-1}(tgl(i))
    do i = 1, n
      call legepols_r64(tgl(i), n, pols)
      do j = 1, n
        v(i,j) = pols(j)
      end do
    end do

    ! build Legmat (u): u(i,j) = (2i-1)/2 * v(j,i) * wgl(j)
    do i = 1, n
      d = real(2*i-1, r64) / 2.0_r64
      do j = 1, n
        Legmat(i,j) = d * v(j,i) * wgl(j)
      end do
    end do

  end subroutine legeexps_r64

  subroutine legeexps_r128(itype, n, tgl, Legmat, v, wgl)
    integer(8), intent(in)    :: itype, n
    real(r128), intent(inout) :: tgl(n), Legmat(n,n), v(n,n), wgl(n)

    integer(8) :: i, j
    real(r128) :: pols(n), d, Dgl_tmp(n,n)

    call gauss_r128(n, tgl, wgl, Dgl_tmp)
    if (itype < 2) return

    do i = 1, n
      call legepols_r128(tgl(i), n, pols)
      do j = 1, n
        v(i,j) = pols(j)
      end do
    end do

    do i = 1, n
      d = real(2*i-1, r128) / 2.0_r128
      do j = 1, n
        Legmat(i,j) = d * v(j,i) * wgl(j)
      end do
    end do

  end subroutine legeexps_r128

  ! ------------------------------------------------------------------
  ! legepols_r64 / legepols_r128
  ! 3-term recurrence: pols(k) = P_{k-1}(x), k=1..n.
  ! Mirrors legepols in legeexps.f.
  ! ------------------------------------------------------------------
  subroutine legepols_r64(x, n, pols)
    real(r64),  intent(in)  :: x
    integer(8), intent(in)  :: n
    real(r64),  intent(out) :: pols(n)

    integer(8) :: k
    real(r64)  :: pkm1, pk, pkp1

    pols(1) = 1.0_r64
    if (n == 1) return
    pols(2) = x
    if (n == 2) return
    pkm1 = 1.0_r64;  pk = x
    do k = 1, n-2
      pkp1 = ((2*k+1)*x*pk - k*pkm1) / real(k+1, r64)
      pols(k+2) = pkp1
      pkm1 = pk;  pk = pkp1
    end do

  end subroutine legepols_r64

  subroutine legepols_r128(x, n, pols)
    real(r128), intent(in)  :: x
    integer(8), intent(in)  :: n
    real(r128), intent(out) :: pols(n)

    integer(8) :: k
    real(r128) :: pkm1, pk, pkp1

    pols(1) = 1.0_r128
    if (n == 1) return
    pols(2) = x
    if (n == 2) return
    pkm1 = 1.0_r128;  pk = x
    do k = 1, n-2
      pkp1 = ((2*k+1)*x*pk - k*pkm1) / real(k+1, r128)
      pols(k+2) = pkp1
      pkm1 = pk;  pk = pkp1
    end do

  end subroutine legepols_r128

end module linequaaadrature_mod
