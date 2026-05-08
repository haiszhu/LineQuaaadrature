! koorn_geom.f90
! Module: Koornwinder polynomials, Vioreanu-Rokhlin nodes,
!         boundary panel generation, circumcircle transform.
!
! Ported from rrq-legacy/src/koornexps.f90 and linequad.f.
! No LAPACK dependency: uses a local LU solver for small matrices.
!
! Public:
!   koorn_vals2coefs_coefs2vals(korder, kpols, umatr, vmatr)
!   line3quadr_3dline(x_uvs, korder, kpols, umatr, nquad, tgl, wgl, Dgl,
!                     sbdnp, tpan, nbd, sxbd, swbd, stangbd, sspbd, r_vert0b)
!   circumcircle_transform_3d(r_vert, R, c, alpha)

module koorn_geom_mod
  implicit none
  integer, parameter :: r64  = 8
  integer, parameter :: r128 = 16
  private
  public :: koorn_vals2coefs_coefs2vals, line3quadr_3dline, &
            circumcircle_transform_3d, circumcircle_transform_3d_r128, &
            lqkg_setup_target_r128

contains

  ! ----------------------------------------------------------------
  ! koorn_pols  (private)
  ! Orthonormal Koornwinder polynomials on the simplex (0,0),(1,0),(0,1).
  ! Returns npols = (nmax+1)*(nmax+2)/2 values in pols.
  ! Ordering: (n,k) = (0,0),(1,0),(1,1),(2,0),(2,1),(2,2),...
  ! ----------------------------------------------------------------
  subroutine koorn_pols(uv, nmax, npols, pols)
    real(r64), intent(in)  :: uv(2)
    integer(8), intent(in) :: nmax
    integer(8), intent(out):: npols
    real(r64), intent(out) :: pols(*)

    real(r64) :: legpols(0:100), jacpols(0:100,0:100)
    real(r64) :: u, v, z, y, x, sc, an, bn, cn
    integer(8) :: k, n, iii

    u = uv(1);  v = uv(2)
    z = 2.0_r64*u + v - 1.0_r64
    y = 1.0_r64 - v

    ! P_k(z/y)*y^k  via recurrence (handles y=0 via the legpols array)
    legpols(0) = 1.0_r64
    legpols(1) = z
    do k = 1, nmax
      legpols(k+1) = ((2*k+1)*z*legpols(k) - k*legpols(k-1)*y*y) / real(k+1,r64)
    end do

    ! Jacobi P^{(0,2k+1)} recurrence for each k
    x = 1.0_r64 - 2.0_r64*v
    do k = 0, nmax
      jacpols(0,k) = 1.0_r64
      jacpols(1,k) = (-real(2*k+1,r64) + real(2+2*k+1,r64)*x) / 2.0_r64
      do n = 1, nmax-k-1
        an = real((2*n+2*k+1+1)*(2*n+2*k+1+2),r64) / real(2*(n+1)*(n+2*k+1+1),r64)
        bn = -real((2*k+1)**2,r64) * real(2*n+2*k+1+1,r64) &
             / real(2*(n+1)*(n+2*k+1+1)*(2*n+2*k+1),r64)
        cn = real(n*(n+2*k+1)*(2*n+2*k+1+2),r64) &
             / real((n+1)*(n+2*k+1+1)*(2*n+2*k+1),r64)
        jacpols(n+1,k) = (an*x + bn)*jacpols(n,k) - cn*jacpols(n-1,k)
      end do
    end do

    ! Assemble orthonormal Koornwinder polynomials
    iii = 0
    do n = 0, nmax
      do k = 0, n
        sc = sqrt(1.0_r64 / real((2*k+1)*(2*n+2),r64))
        iii = iii + 1
        pols(iii) = legpols(k) * jacpols(n-k,k) / sc
      end do
    end do
    npols = iii

  end subroutine koorn_pols

  ! ----------------------------------------------------------------
  ! get_vioreanu_nodes  (private)
  ! Load precomputed Vioreanu-Rokhlin nodes of order norder.
  ! ----------------------------------------------------------------
  subroutine get_vioreanu_nodes(norder, npols, uvs)
    integer(8), intent(in)  :: norder, npols
    real(r64),  intent(out) :: uvs(2,npols)
    INCLUDE 'koorn-uvs-dat.txt'
  end subroutine get_vioreanu_nodes

  ! ----------------------------------------------------------------
  ! koorn_coefs2vals  (private)
  ! Build Vandermonde matrix: amat(i,j) = K_j(uvs(:,i)).
  ! ----------------------------------------------------------------
  subroutine koorn_coefs2vals(nmax, npols, uvs, amat)
    integer(8), intent(in)  :: nmax, npols
    real(r64),  intent(in)  :: uvs(2,npols)
    real(r64),  intent(out) :: amat(npols,npols)

    real(r64)  :: pols(2000)
    integer(8) :: i, j, npols2

    do i = 1, npols
      call koorn_pols(uvs(:,i), nmax, npols2, pols)
      do j = 1, npols
        amat(i,j) = pols(j)
      end do
    end do

  end subroutine koorn_coefs2vals

  ! ----------------------------------------------------------------
  ! lu_solve  (private)
  ! Solve A*X = B in-place: on entry A(n,n) and B(n,k);
  ! on exit A is overwritten with LU, B with the solution.
  ! Partial pivoting.  For small n (~36).
  ! ----------------------------------------------------------------
  subroutine lu_solve(n, A, k, B)
    integer(8), intent(in)    :: n, k
    real(r64),  intent(inout) :: A(n,n), B(n,k)

    integer(8) :: i, j, p, col, pivot_row
    real(r64)  :: pivot, tmp, factor, rowA(n), rowB(k)

    do col = 1, n
      ! Find pivot
      pivot_row = col
      pivot     = abs(A(col,col))
      do i = col+1, n
        if (abs(A(i,col)) > pivot) then
          pivot = abs(A(i,col))
          pivot_row = i
        end if
      end do
      ! Swap rows in A and B
      if (pivot_row /= col) then
        rowA = A(col,:);    A(col,:)        = A(pivot_row,:);   A(pivot_row,:)  = rowA
        rowB = B(col,:);    B(col,:)        = B(pivot_row,:);   B(pivot_row,:) = rowB
      end if
      ! Eliminate
      do i = col+1, n
        factor = A(i,col) / A(col,col)
        A(i,col) = factor
        do j = col+1, n
          A(i,j) = A(i,j) - factor*A(col,j)
        end do
        do p = 1, k
          B(i,p) = B(i,p) - factor*B(col,p)
        end do
      end do
    end do

    ! Back substitution
    do col = n, 1, -1
      do p = 1, k
        tmp = B(col,p)
        do j = col+1, n
          tmp = tmp - A(col,j)*B(j,p)
        end do
        B(col,p) = tmp / A(col,col)
      end do
    end do

  end subroutine lu_solve

  ! ----------------------------------------------------------------
  ! koorn_vals2coefs_coefs2vals  (public)
  ! Compute vals->coefs matrix umatr and coefs->vals matrix vmatr
  ! for Koornwinder polynomials at Vioreanu-Rokhlin nodes of order korder.
  ! ----------------------------------------------------------------
  subroutine koorn_vals2coefs_coefs2vals(korder, kpols, umatr, vmatr)
    integer(8), intent(in)    :: korder, kpols
    real(r64),  intent(inout) :: umatr(kpols,kpols), vmatr(kpols,kpols)

    real(r64)  :: xys(2,kpols), vtmp(kpols,kpols)
    integer(8) :: i, info

    call get_vioreanu_nodes(korder, kpols, xys)
    call koorn_coefs2vals(korder, kpols, xys, vmatr)

    ! umatr = inv(vmatr): solve vmatr * umatr = I
    vtmp  = vmatr
    umatr = 0.0_r64
    do i = 1, kpols
      umatr(i,i) = 1.0_r64
    end do
    call lu_solve(kpols, vtmp, kpols, umatr)

  end subroutine koorn_vals2coefs_coefs2vals

  ! ----------------------------------------------------------------
  ! line3quadr_3dline  (public)
  ! Generate GL boundary quadrature nodes/weights/tangents for a
  ! curved triangle boundary, using Koornwinder parameterization.
  ! Ported from linequad.f:5683.
  !
  ! x_uvs(3,kpols): source points on VR grid (3D positions)
  ! umatr(kpols,kpols): vals->coefs matrix from koorn_vals2coefs_coefs2vals
  ! nquad: GL order per panel
  ! tgl(nquad), wgl(nquad), Dgl(nquad,nquad): GL nodes/weights/D-matrix
  ! sbdnp: number of panels (=3 for triangle boundary)
  ! tpan(sbdnp+1): panel breakpoints in [0,2pi]
  ! nbd = sbdnp*nquad
  ! Output: sxbd(3,nbd), swbd(nbd), stangbd(3,nbd), sspbd(nbd), r_vert0b(3,3)
  ! ----------------------------------------------------------------
  subroutine line3quadr_3dline(x_uvs, korder, kpols, umatr, nquad, &
                                tgl, wgl, Dgl, sbdnp, tpan, nbd,   &
                                sxbd, swbd, stangbd, sspbd, r_vert0b)
    integer(8), intent(in)    :: korder, kpols, nquad, sbdnp, nbd
    real(r64),  intent(in)    :: x_uvs(3,kpols), umatr(kpols,kpols)
    real(r64),  intent(in)    :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
    real(r64),  intent(in)    :: tpan(sbdnp+1)
    real(r64),  intent(inout) :: sxbd(3,nbd), swbd(nbd)
    real(r64),  intent(inout) :: stangbd(3,nbd), sspbd(nbd)
    real(r64),  intent(inout) :: r_vert0b(3,3)

    real(r64)  :: coefs_xyz(kpols,3)
    real(r64)  :: tlo(sbdnp), thi(sbdnp), pt(sbdnp)
    real(r64)  :: uvlo(2,sbdnp), uvhi(2,sbdnp)
    real(r64)  :: xlo(3,sbdnp), xhi(3,sbdnp)
    real(r64)  :: uv(2), pols(kpols), sxp(3,nbd), t(nbd)
    real(r64)  :: uvend(2,3), sx_ell(3,nquad)
    real(r64)  :: thold1, thold2, tholdinv1, pi
    integer(8) :: k, ell, i, ii, ii_start, ii_end, itmp

    pi       = 4.0_r64 * atan(1.0_r64)
    thold1   = 2.0_r64 * pi / 3.0_r64
    thold2   = 4.0_r64 * pi / 3.0_r64
    tholdinv1 = 1.0_r64 / thold1

    tlo = tpan(1:sbdnp)
    thi = tpan(2:sbdnp+1)
    pt  = thi - tlo

    ! Koornwinder coefs for xyz: coefs_xyz(kpols,3) = umatr * x_uvs^T
    coefs_xyz = matmul(umatr, transpose(x_uvs))

    ! uvlo, uvhi: reference (u,v) coords of panel endpoints
    uvlo = 0.0_r64;  uvhi = 0.0_r64
    do k = 1, sbdnp
      call tparam_to_uv(tlo(k), thold1, thold2, tholdinv1, uvlo(:,k))
      call tparam_to_uv(thi(k), thold1, thold2, tholdinv1, uvhi(:,k))
    end do

    ! Panel endpoint positions (for xlo, xhi — not used after vertex extraction)
    do k = 1, sbdnp
      call koorn_pols(uvlo(:,k), korder, ii, pols)
      xlo(:,k) = matmul(pols(1:kpols), coefs_xyz)
      call koorn_pols(uvhi(:,k), korder, ii, pols)
      xhi(:,k) = matmul(pols(1:kpols), coefs_xyz)
    end do

    ! Triangle vertices: uv = (0,0), (1,0), (0,1)
    uvend(:,1) = [0.0_r64, 0.0_r64]
    uvend(:,2) = [1.0_r64, 0.0_r64]
    uvend(:,3) = [0.0_r64, 1.0_r64]
    do k = 1, 3
      call koorn_pols(uvend(:,k), korder, ii, pols)
      r_vert0b(:,k) = matmul(pols(1:kpols), coefs_xyz)
    end do

    ! GL nodes mapped to each panel
    do ell = 1, sbdnp
      itmp = (ell-1)*nquad
      do i = 1, nquad
        ii = itmp + i
        t(ii)    = tlo(ell) + 0.5_r64*(1.0_r64 + tgl(i)) * pt(ell)
        swbd(ii) = 0.5_r64 * wgl(i) * pt(ell)
      end do
    end do

    ! Position sxbd and derivative sxp at all nodes
    sxbd = 0.0_r64;  sxp = 0.0_r64
    do k = 1, nbd
      call tparam_to_uv(t(k), thold1, thold2, tholdinv1, uv)
      call koorn_pols(uv, korder, ii, pols)
      sxbd(:,k) = matmul(pols(1:kpols), coefs_xyz)
    end do

    ! Derivative via differentiation matrix per panel
    do ell = 1, sbdnp
      ii_start = (ell-1)*nquad + 1
      ii_end   = ell*nquad
      sx_ell   = sxbd(:, ii_start:ii_end)
      ! sxp = (2/pt) * D * sx_ell^T  (D acts on GL index)
      sxp(1, ii_start:ii_end) = (2.0_r64/pt(ell)) * matmul(Dgl, sx_ell(1,:))
      sxp(2, ii_start:ii_end) = (2.0_r64/pt(ell)) * matmul(Dgl, sx_ell(2,:))
      sxp(3, ii_start:ii_end) = (2.0_r64/pt(ell)) * matmul(Dgl, sx_ell(3,:))
    end do

    sspbd = sqrt(sxp(1,:)**2 + sxp(2,:)**2 + sxp(3,:)**2)
    stangbd(1,:) = sxp(1,:) / sspbd
    stangbd(2,:) = sxp(2,:) / sspbd
    stangbd(3,:) = sxp(3,:) / sspbd
    swbd = swbd * sspbd

  end subroutine line3quadr_3dline

  ! helper: convert parameter t in [0,2pi) to (u,v) simplex coords
  subroutine tparam_to_uv(t, thold1, thold2, tholdinv1, uv)
    real(r64), intent(in)  :: t, thold1, thold2, tholdinv1
    real(r64), intent(out) :: uv(2)
    uv = 0.0_r64
    if (t < thold1) then
      uv(2) = 1.0_r64 - tholdinv1*t
    else if (t < thold2) then
      uv(1) = tholdinv1*(t - thold1)
    else
      uv(1) = 1.0_r64 - tholdinv1*(t - thold2)
      uv(2) = tholdinv1*(t - thold2)
    end if
  end subroutine tparam_to_uv

  ! ----------------------------------------------------------------
  ! circumcircle_transform_3d  (public)
  ! Similarity transform mapping a triangle's circumcircle to the
  ! unit circle in the xy-plane:  x_new = alpha * R * (x - c)
  !
  ! r_vert(3,3): triangle vertices (columns)
  ! R(3,3):      rotation matrix (output)
  ! c(3):        circumcenter (output)
  ! alpha:       1/circumradius (output)
  ! ----------------------------------------------------------------
  subroutine circumcircle_transform_3d(r_vert, R, c, alpha)
    real(r64), intent(in)    :: r_vert(3,3)
    real(r64), intent(inout) :: R(3,3), c(3), alpha

    real(r64) :: r1(3), r2(3), r3(3)
    real(r64) :: a(3), b(3), axb(3)
    real(r64) :: aa, bb, ab, denom, s, t
    real(r64) :: R_circ, e1(3), e2(3), e3(3), norm_e

    r1 = r_vert(:,1);  r2 = r_vert(:,2);  r3 = r_vert(:,3)
    a  = r2 - r1;      b  = r3 - r1

    aa = dot_product(a,a);  bb = dot_product(b,b);  ab = dot_product(a,b)

    ! a x b
    axb(1) = a(2)*b(3) - a(3)*b(2)
    axb(2) = a(3)*b(1) - a(1)*b(3)
    axb(3) = a(1)*b(2) - a(2)*b(1)

    denom = 2.0_r64 * dot_product(axb, axb)
    s = bb*(aa - ab) / denom
    t = aa*(bb - ab) / denom
    c = r1 + s*a + t*b

    R_circ = sqrt(dot_product(r1-c, r1-c))
    alpha  = 1.0_r64 / R_circ

    ! Rotation frame
    e3 = axb / sqrt(dot_product(axb,axb))        ! patch normal -> +z
    e1 = (r1 - c) / R_circ                        ! toward r1    -> +x
    ! e2 = e3 x e1
    e2(1) = e3(2)*e1(3) - e3(3)*e1(2)
    e2(2) = e3(3)*e1(1) - e3(1)*e1(3)
    e2(3) = e3(1)*e1(2) - e3(2)*e1(1)

    R(1,:) = e1;  R(2,:) = e2;  R(3,:) = e3    ! rows = new basis vectors

  end subroutine circumcircle_transform_3d

  ! ----------------------------------------------------------------
  ! circumcircle_transform_3d_r128  (real(16) port)
  ! ----------------------------------------------------------------
  subroutine circumcircle_transform_3d_r128(r_vert, R, c, alpha)
    real(r128), intent(in)    :: r_vert(3,3)
    real(r128), intent(inout) :: R(3,3), c(3), alpha

    real(r128) :: r1(3), r2(3), r3(3)
    real(r128) :: a(3), b(3), axb(3)
    real(r128) :: aa, bb, ab, denom, s, t
    real(r128) :: R_circ, e1(3), e2(3), e3(3)

    r1 = r_vert(:,1);  r2 = r_vert(:,2);  r3 = r_vert(:,3)
    a  = r2 - r1;      b  = r3 - r1

    aa = dot_product(a,a);  bb = dot_product(b,b);  ab = dot_product(a,b)

    axb(1) = a(2)*b(3) - a(3)*b(2)
    axb(2) = a(3)*b(1) - a(1)*b(3)
    axb(3) = a(1)*b(2) - a(2)*b(1)

    denom = 2.0_r128 * dot_product(axb, axb)
    s = bb*(aa - ab) / denom
    t = aa*(bb - ab) / denom
    c = r1 + s*a + t*b

    R_circ = sqrt(dot_product(r1-c, r1-c))
    alpha  = 1.0_r128 / R_circ

    e3 = axb / sqrt(dot_product(axb,axb))
    e1 = (r1 - c) / R_circ
    e2(1) = e3(2)*e1(3) - e3(3)*e1(2)
    e2(2) = e3(3)*e1(1) - e3(1)*e1(3)
    e2(3) = e3(1)*e1(2) - e3(2)*e1(1)

    R(1,:) = e1;  R(2,:) = e2;  R(3,:) = e3

  end subroutine circumcircle_transform_3d_r128

  ! ----------------------------------------------------------------
  ! lqkg_setup_target_r128
  ! Bundles the per-triangle r128 setup that the orchestration
  ! lqs_evaluate_solid_angle_integral_r128 used to do in MATLAB:
  !   1. sxbd = sxbd_in;
  !      per panel ell: stangbd, sspbd from Dgl * sxbd_in / |Dgl * sxbd_in|.
  !   2. circumcircle_transform_3d_r128(r_vert, R, c, alpha).
  !   3. txnew  = alpha * R * (tx  - c)
  !      snxnew =          R * snx
  !      sxbd   = alpha * R * (sxbd - c)
  !      stangbd =          R * stangbd
  !   4. qhat = mean(snxnew(:,j) over j) / |...|
  !   5. per panel ell: sxpbd = Dgl * sxbd.
  !   6. kdata(:,j) = qhat for j = 1..m   (Asvestas convention).
  ! All compute is r128; outputs are exposed for HDF5 write at the mex layer.
  ! ----------------------------------------------------------------
  subroutine lqkg_setup_target_r128(m, n, nbd, sbdnp, nquad,                  &
                                    r_vert, tx, snx, sxbd_in, Dgl,            &
                                    R, c, alpha,                              &
                                    sxbd, sxpbd, stangbd, sspbd,              &
                                    txnew, snxnew, qhat, kdata)
    integer(8), intent(in)    :: m, n, nbd, sbdnp, nquad
    real(r128), intent(in)    :: r_vert(3,3)
    real(r128), intent(in)    :: tx(3,m), snx(3,n)
    real(r128), intent(in)    :: sxbd_in(3,nbd)
    real(r128), intent(in)    :: Dgl(nquad,nquad)
    real(r128), intent(out)   :: R(3,3), c(3), alpha
    real(r128), intent(out)   :: sxbd(3,nbd), sxpbd(3,nbd)
    real(r128), intent(out)   :: stangbd(3,nbd), sspbd(nbd)
    real(r128), intent(out)   :: txnew(3,m), snxnew(3,n)
    real(r128), intent(out)   :: qhat(3), kdata(3,m)

    integer(8) :: ell, k, q, j, idx_start, idx_end, idx
    real(r128) :: sxp(3), qnrm

    ! 1. stangbd, sspbd from Dgl * sxbd_in (per panel)
    sxbd = sxbd_in
    do ell = 1_8, sbdnp
      idx_start = (ell - 1_8)*nquad + 1_8
      idx_end   = ell*nquad
      do k = 1_8, nquad
        sxp(1) = 0.0_r128;  sxp(2) = 0.0_r128;  sxp(3) = 0.0_r128
        do q = 1_8, nquad
          idx = idx_start + q - 1_8
          sxp(1) = sxp(1) + Dgl(k,q) * sxbd(1, idx)
          sxp(2) = sxp(2) + Dgl(k,q) * sxbd(2, idx)
          sxp(3) = sxp(3) + Dgl(k,q) * sxbd(3, idx)
        end do
        idx = idx_start + k - 1_8
        sspbd(idx)        = sqrt(sxp(1)**2 + sxp(2)**2 + sxp(3)**2)
        stangbd(:, idx)   = sxp / sspbd(idx)
      end do
    end do

    ! 2. circumcircle transform
    call circumcircle_transform_3d_r128(r_vert, R, c, alpha)

    ! 3. apply transform
    do j = 1_8, m
      txnew(:,j) = alpha * matmul(R, tx(:,j) - c)
    end do
    do j = 1_8, n
      snxnew(:,j) = matmul(R, snx(:,j))
    end do
    do k = 1_8, nbd
      sxbd(:,k)    = alpha * matmul(R, sxbd(:,k) - c)
      stangbd(:,k) = matmul(R, stangbd(:,k))
    end do

    ! 4. qhat = mean(snxnew, 2) / norm
    qhat = 0.0_r128
    do j = 1_8, n
      qhat = qhat + snxnew(:,j)
    end do
    qhat = qhat / real(n, r128)
    qnrm = sqrt(qhat(1)**2 + qhat(2)**2 + qhat(3)**2)
    qhat = qhat / qnrm

    ! 5. sxpbd = Dgl * sxbd per panel
    sxpbd = 0.0_r128
    do ell = 1_8, sbdnp
      idx_start = (ell - 1_8)*nquad + 1_8
      idx_end   = ell*nquad
      do k = 1_8, nquad
        do q = 1_8, nquad
          idx = idx_start + q - 1_8
          sxpbd(1, idx_start + k - 1_8) = sxpbd(1, idx_start + k - 1_8) + Dgl(k,q) * sxbd(1, idx)
          sxpbd(2, idx_start + k - 1_8) = sxpbd(2, idx_start + k - 1_8) + Dgl(k,q) * sxbd(2, idx)
          sxpbd(3, idx_start + k - 1_8) = sxpbd(3, idx_start + k - 1_8) + Dgl(k,q) * sxbd(3, idx)
        end do
      end do
    end do

    ! 6. kdata = qhat for all j (Asvestas)
    do j = 1_8, m
      kdata(:,j) = qhat
    end do

  end subroutine lqkg_setup_target_r128

end module koorn_geom_mod
