! ellipsoid_mesh_mod.f90
! Triangular mesh generator for an ellipsoid surface.
! Port of create_ellipsoid_tri_mesh.m.
!
! Public:
!   create_ellipsoid_tri_mesh_r64
!   areal_quad_triangle_r64   — conical product rule on 2D reference simplex (drop-in for VR)
!
! Private helpers:
!   get_vr_nodes_wts_r64  — VR nodes + weights (INCLUDE koorn-uvs/wts-dat.txt)
!   ellipsoidparam_r64    — ellipsoid surface point, normal, speed
!   face_map_nodes_r64    — map parameter-space points to 3D for one cube face
!   face_map_point_r64    — scalar version of face_map_nodes_r64
!   edge_nodes_r64        — GL-mapped points along a triangle edge

module ellipsoid_mesh_mod
  use linequaaadrature_mod, only: gauss_r64, r64, gauss_r128, r128
  implicit none
  private
  public :: create_ellipsoid_tri_mesh_r64, areal_quad_triangle_r64
  public :: create_ellipsoid_tri_mesh_r128, areal_quad_triangle_r128

contains

  ! ------------------------------------------------------------------
  ! create_ellipsoid_tri_mesh_r64
  ! Build a VR-quadrature triangular mesh on an ellipsoid surface.
  !
  ! Inputs:
  !   mp, np       : panel counts per cube face (ntri = 12*mp*np)
  !   p            : Vioreanu order (nvr = p*(p+1)/2 nodes per triangle)
  !   ratio        : ellipsoid axis ratio (x-axis = 1, y=z=ratio)
  !   nq           : if > 0, use areal quadrature (nvr = 3*nq^2);
  !                  if = 0, use Vioreanu-Rokhlin nodes (nvr = p*(p+1)/2)
  !   nquad_bdry   : GL nodes per triangle edge (nbd = 3*nquad_bdry)
  !   nvr          : caller must set: p*(p+1)/2 (VR) or 3*nq^2 (areal)
  !   ntri         : 12*mp*np   (caller must compute)
  !
  ! Outputs (all intent(inout), caller pre-allocates):
  !   x(3,nvr,ntri)          : quadrature positions
  !   nx(3,nvr,ntri)         : outward unit normals at VR nodes
  !   w(nvr,ntri)            : quadrature weights
  !   xbd(3,3*nquad_bdry,ntri) : boundary GL nodes (3 edges x nquad_bdry)
  !   tri2face(ntri)         : cube face index (1..6) for each triangle
  !   tri2cell(2,ntri)       : (i,j) panel index within face
  !   tri_vert(3,3,ntri)     : 3D coordinates of triangle vertices
  !   ptr(ntri+1)            : CSR pointer: ptr(k) = 1 + (k-1)*nvr
  ! ------------------------------------------------------------------
  subroutine create_ellipsoid_tri_mesh_r64(mp, np, p, nq, ratio, nquad_bdry, nvr, ntri, &
      x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr)
    integer(8), intent(in)    :: mp, np, p, nq, nquad_bdry, nvr, ntri
    real(r64),  intent(in)    :: ratio
    real(r64),  intent(inout) :: x(3,nvr,ntri), nx(3,nvr,ntri), w(nvr,ntri)
    real(r64),  intent(inout) :: xbd(3,3*nquad_bdry,ntri)
    integer(8), intent(inout) :: tri2face(ntri)
    integer(8), intent(inout) :: tri2cell(2,ntri)
    real(r64),  intent(inout) :: tri_vert(3,3,ntri)
    integer(8), intent(inout) :: ptr(ntri+1)

    real(r64) :: pi
    real(r64) :: rotation(6)
    real(r64) :: V(2,3)       ! parameter-space vertices of current triangle
    real(r64) :: vll(2), vlr(2), vur(2), vul(2)
    real(r64) :: rll(3), rlr(3), rur(3), rul(3)
    real(r64) :: r_v(3)
    real(r64) :: d1, d2, jac, wwt_k, det2x2
    integer(8) :: iface, i, j, kk, iv, k_idx, itri, nb3

    real(r64), allocatable :: uvs(:,:), vr_wts(:)
    real(r64), allocatable :: xg(:), wg(:), Dg(:,:)
    real(r64), allocatable :: tgl_a(:), wgl_a(:), Dgl_a(:,:)
    real(r64), allocatable :: x12(:,:)
    real(r64), allocatable :: xt(:,:), nxt(:,:), spt(:), tht(:)
    real(r64), allocatable :: e1(:,:), e2(:,:), e3(:,:)
    real(r64), allocatable :: xb1(:,:), nxb1(:,:), spb1(:), thb1(:)
    real(r64), allocatable :: xb2(:,:), nxb2(:,:), spb2(:), thb2(:)
    real(r64), allocatable :: xb3(:,:), nxb3(:,:), spb3(:), thb3(:)

    pi = 4.0_r64 * atan(1.0_r64)
    rotation(1) =  0.0_r64
    rotation(2) =  0.0_r64
    rotation(3) = -pi / 2.0_r64
    rotation(4) =  0.0_r64
    rotation(5) = -pi
    rotation(6) = -1.5_r64 * pi

    allocate(uvs(2,nvr), vr_wts(nvr))
    allocate(xg(nquad_bdry), wg(nquad_bdry), Dg(nquad_bdry,nquad_bdry))
    allocate(x12(2,nvr))
    allocate(xt(3,nvr), nxt(3,nvr), spt(nvr), tht(nvr))
    allocate(e1(2,nquad_bdry), e2(2,nquad_bdry), e3(2,nquad_bdry))
    allocate(xb1(3,nquad_bdry), nxb1(3,nquad_bdry), spb1(nquad_bdry), thb1(nquad_bdry))
    allocate(xb2(3,nquad_bdry), nxb2(3,nquad_bdry), spb2(nquad_bdry), thb2(nquad_bdry))
    allocate(xb3(3,nquad_bdry), nxb3(3,nquad_bdry), spb3(nquad_bdry), thb3(nquad_bdry))

    if (nq > 0_8) then
      ! Areal quadrature: use conical product rule (nvr = 3*nq^2)
      allocate(tgl_a(nq), wgl_a(nq), Dgl_a(nq,nq))
      call gauss_r64(nq, tgl_a, wgl_a, Dgl_a)
      call areal_quad_triangle_r64(nq, tgl_a, wgl_a, uvs, vr_wts)
      deallocate(tgl_a, wgl_a, Dgl_a)
    else
      ! Vioreanu-Rokhlin nodes (nvr = p*(p+1)/2)
      call get_vr_nodes_wts_r64(p - 1_8, nvr, uvs, vr_wts)
    end if
    call gauss_r64(nquad_bdry, xg, wg, Dg)

    itri = 0_8
    do iface = 1_8, 6_8
      do j = 1_8, np
        do i = 1_8, mp

          ! parameter panel corners in [-1,1]^2
          vll(1) = 2.0_r64*(j-1)/np - 1.0_r64;  vll(2) = 2.0_r64*(i-1)/mp - 1.0_r64
          vlr(1) = 2.0_r64*j/np     - 1.0_r64;  vlr(2) = 2.0_r64*(i-1)/mp - 1.0_r64
          vul(1) = 2.0_r64*(j-1)/np - 1.0_r64;  vul(2) = 2.0_r64*i/mp     - 1.0_r64
          vur(1) = 2.0_r64*j/np     - 1.0_r64;  vur(2) = 2.0_r64*i/mp     - 1.0_r64

          call face_map_point_r64(vll(1), vll(2), iface, ratio, rotation, rll)
          call face_map_point_r64(vur(1), vur(2), iface, ratio, rotation, rur)
          call face_map_point_r64(vul(1), vul(2), iface, ratio, rotation, rul)
          call face_map_point_r64(vlr(1), vlr(2), iface, ratio, rotation, rlr)

          d1 = sum((rll - rur)**2)
          d2 = sum((rul - rlr)**2)

          do kk = 1_8, 2_8
            itri = itri + 1_8

            ! Choose triangle vertices (2D parameter-space columns)
            if (d1 + 1.0e-13_r64 > d2) then
              ! split along rul-rlr diagonal
              if (kk == 1) then
                V(:,1) = vll;  V(:,2) = vlr;  V(:,3) = vul
              else
                V(:,1) = vur;  V(:,2) = vul;  V(:,3) = vlr
              end if
            else
              ! split along rll-rur diagonal
              if (kk == 1) then
                V(:,1) = vul;  V(:,2) = vll;  V(:,3) = vur
              else
                V(:,1) = vlr;  V(:,2) = vur;  V(:,3) = vll
              end if
            end if

            ! Map VR nodes to parameter space via barycentric coords:
            !   x12 = V1 + (V2-V1)*uvs(1,:) + (V3-V1)*uvs(2,:)
            do iv = 1_8, nvr
              x12(1,iv) = V(1,1) + (V(1,2)-V(1,1))*uvs(1,iv) + (V(1,3)-V(1,1))*uvs(2,iv)
              x12(2,iv) = V(2,1) + (V(2,2)-V(2,1))*uvs(1,iv) + (V(2,3)-V(2,1))*uvs(2,iv)
            end do

            call face_map_nodes_r64(x12(1,:), x12(2,:), nvr, iface, ratio, rotation, &
                                    xt, nxt, spt, tht)

            ! 2x2 Jacobian det of the parameter-space triangle map
            det2x2 = (V(1,2)-V(1,1))*(V(2,3)-V(2,1)) - (V(1,3)-V(1,1))*(V(2,2)-V(2,1))
            jac = abs(det2x2)

            ! Quadrature weights: w = sp * wwt * vr_wts * jac
            !   wwt = cos(th) / (1 + x12(1)^2 + x12(2)^2)
            !   th  = azimuth from cart2sph (first output of face_map_nodes)
            do iv = 1_8, nvr
              wwt_k = cos(tht(iv)) / (1.0_r64 + x12(1,iv)**2 + x12(2,iv)**2)
              w(iv, itri) = spt(iv) * wwt_k * vr_wts(iv) * jac
            end do

            x(:,:,itri)  = xt
            nx(:,:,itri) = nxt

            ! Edge boundary nodes (3 edges, each nquad_bdry GL points)
            call edge_nodes_r64(V(:,1), V(:,2), xg, nquad_bdry, e1)
            call edge_nodes_r64(V(:,2), V(:,3), xg, nquad_bdry, e2)
            call edge_nodes_r64(V(:,3), V(:,1), xg, nquad_bdry, e3)

            call face_map_nodes_r64(e1(1,:), e1(2,:), nquad_bdry, iface, ratio, rotation, &
                                    xb1, nxb1, spb1, thb1)
            call face_map_nodes_r64(e2(1,:), e2(2,:), nquad_bdry, iface, ratio, rotation, &
                                    xb2, nxb2, spb2, thb2)
            call face_map_nodes_r64(e3(1,:), e3(2,:), nquad_bdry, iface, ratio, rotation, &
                                    xb3, nxb3, spb3, thb3)

            nb3 = nquad_bdry
            xbd(:, 1:nb3,          itri) = xb1
            xbd(:, nb3+1:2*nb3,    itri) = xb2
            xbd(:, 2*nb3+1:3*nb3,  itri) = xb3

            ! Triangle vertices in 3D
            call face_map_point_r64(V(1,1), V(2,1), iface, ratio, rotation, r_v)
            tri_vert(:,1,itri) = r_v
            call face_map_point_r64(V(1,2), V(2,2), iface, ratio, rotation, r_v)
            tri_vert(:,2,itri) = r_v
            call face_map_point_r64(V(1,3), V(2,3), iface, ratio, rotation, r_v)
            tri_vert(:,3,itri) = r_v

            tri2face(itri)   = iface
            tri2cell(1,itri) = i
            tri2cell(2,itri) = j

          end do   ! kk
        end do   ! i
      end do   ! j
    end do   ! iface

    ! CSR pointer: ptr(k) = 1 + (k-1)*nvr
    do k_idx = 1_8, ntri + 1_8
      ptr(k_idx) = 1_8 + (k_idx - 1_8) * nvr
    end do

    deallocate(uvs, vr_wts)
    deallocate(xg, wg, Dg)
    deallocate(x12, xt, nxt, spt, tht)
    deallocate(e1, e2, e3)
    deallocate(xb1, nxb1, spb1, thb1)
    deallocate(xb2, nxb2, spb2, thb2)
    deallocate(xb3, nxb3, spb3, thb3)

  end subroutine create_ellipsoid_tri_mesh_r64

  ! ------------------------------------------------------------------
  ! get_vr_nodes_wts_r64  (private)
  ! Load precomputed Vioreanu-Rokhlin nodes and weights of order norder.
  ! uvs(2,npols): (u,v) coords on the simplex
  ! wts(npols)  : quadrature weights
  !
  ! Variable names norder/uvs/wts are required by the INCLUDE files.
  ! ------------------------------------------------------------------
  subroutine get_vr_nodes_wts_r64(norder, npols, uvs, wts)
    integer(8), intent(in)  :: norder, npols
    real(r64),  intent(out) :: uvs(2,npols), wts(npols)
    INCLUDE 'koorn-uvs-dat.txt'
    INCLUDE 'koorn-wts-dat.txt'
  end subroutine get_vr_nodes_wts_r64

  ! ------------------------------------------------------------------
  ! ellipsoidparam_r64  (private)
  ! Ellipsoid surface parametrisation at n points.
  !
  ! Inputs:
  !   p_ang(n) : azimuth angle phi  (first arg of ellipsoidparam in MATLAB)
  !   t_ang(n) : elevation angle theta (second arg)
  !   ratio    : ellipsoid axis ratio (y=z axis; x-axis = 1)
  !
  ! Outputs:
  !   x_out(3,n)  : surface position r = f*[cos(t)*cos(p); cos(t)*sin(p); sin(t)]
  !   nx_out(3,n) : outward unit normal = cross(rp, rt) / |cross(rp, rt)|
  !   sp(n)       : surface area element = |cross(rp, rt)|
  ! ------------------------------------------------------------------
  subroutine ellipsoidparam_r64(p_ang, t_ang, n, ratio, x_out, nx_out, sp)
    integer(8), intent(in)  :: n
    real(r64),  intent(in)  :: p_ang(n), t_ang(n), ratio
    real(r64),  intent(out) :: x_out(3,n), nx_out(3,n), sp(n)

    integer(8) :: k
    real(r64)  :: p, t, f, ft, fp
    real(r64)  :: rt(3), rp(3), nxk(3), spk

    do k = 1_8, n
      p = p_ang(k);  t = t_ang(k)

      f  = 1.0_r64 / sqrt(cos(p)**2*cos(t)**2 + ratio**2*sin(p)**2*cos(t)**2 &
                          + ratio**2*sin(t)**2)

      ft = -0.5_r64 * f**3 * ( -2.0_r64*cos(p)**2*cos(t)*sin(t) &
                                -2.0_r64*ratio**2*sin(p)**2*cos(t)*sin(t) &
                                +2.0_r64*ratio**2*sin(t)*cos(t) )

      fp = -0.5_r64 * f**3 * ( -2.0_r64*cos(p)*sin(p)*cos(t)**2 &
                                +2.0_r64*ratio**2*sin(p)*cos(p)*cos(t)**2 )

      x_out(1,k) = f * cos(t) * cos(p)
      x_out(2,k) = f * cos(t) * sin(p)
      x_out(3,k) = f * sin(t)

      rt(1) = ft*cos(t)*cos(p) - f*sin(t)*cos(p)
      rt(2) = ft*cos(t)*sin(p) - f*sin(t)*sin(p)
      rt(3) = ft*sin(t)        + f*cos(t)

      rp(1) = fp*cos(t)*cos(p) - f*cos(t)*sin(p)
      rp(2) = fp*cos(t)*sin(p) + f*cos(t)*cos(p)
      rp(3) = fp*sin(t)

      ! nx = cross(rp, rt)
      nxk(1) = rp(2)*rt(3) - rp(3)*rt(2)
      nxk(2) = rp(3)*rt(1) - rp(1)*rt(3)
      nxk(3) = rp(1)*rt(2) - rp(2)*rt(1)

      spk    = sqrt(nxk(1)**2 + nxk(2)**2 + nxk(3)**2)
      sp(k)  = spk
      nx_out(:,k) = nxk / spk
    end do

  end subroutine ellipsoidparam_r64

  ! ------------------------------------------------------------------
  ! face_map_nodes_r64  (private)
  ! Map n parameter-space points (x1,x2) in [-1,1]^2 to 3D ellipsoid
  ! positions for a given cube face.
  !
  ! Cube-sphere projection: q = [1; x1; x2] / sqrt(1+x1^2+x2^2)
  ! cart2sph(q(1), q(2), q(3)):
  !   azimuth   th  = atan2(q(2), q(1))
  !   elevation phi = atan2(q(3), sqrt(q(1)^2+q(2)^2))
  !
  ! Faces 1,4 : call ellipsoidparam(th + floor(face/3)*pi, phi)
  ! Faces 2,3,5,6: call ellipsoidparam(th + pi/2, phi), then Rx(rotation(face))
  !
  ! th(n): azimuth angle (first cart2sph output), needed for wwt in weights.
  ! ------------------------------------------------------------------
  subroutine face_map_nodes_r64(x1, x2, n, iface, ratio, rotation, x_out, nx_out, sp, th)
    integer(8), intent(in)  :: n, iface
    real(r64),  intent(in)  :: x1(n), x2(n), ratio, rotation(6)
    real(r64),  intent(out) :: x_out(3,n), nx_out(3,n), sp(n), th(n)

    real(r64) :: pi
    real(r64) :: qnrm, q1, q2, q3, azimuth, elevation
    real(r64) :: ang, cang, sang, tmp2, tmp3
    real(r64), allocatable :: p_arr(:), t_arr(:)
    integer(8) :: k

    pi = 4.0_r64 * atan(1.0_r64)

    allocate(p_arr(n), t_arr(n))

    do k = 1_8, n
      qnrm = sqrt(1.0_r64 + x1(k)**2 + x2(k)**2)
      q1   = 1.0_r64 / qnrm
      q2   = x1(k)   / qnrm
      q3   = x2(k)   / qnrm

      ! cart2sph: azimuth = atan2(q2,q1),  elevation = atan2(q3, sqrt(q1^2+q2^2))
      azimuth   = atan2(q2, q1)
      elevation = atan2(q3, sqrt(q1**2 + q2**2))

      th(k) = azimuth

      if (iface == 1_8 .or. iface == 4_8) then
        p_arr(k) = azimuth + floor(real(iface,r64) / 3.0_r64) * pi
      else
        p_arr(k) = azimuth + 0.5_r64 * pi
      end if
      t_arr(k) = elevation
    end do

    call ellipsoidparam_r64(p_arr, t_arr, n, ratio, x_out, nx_out, sp)

    ! Apply x-axis rotation for faces 2, 3, 5, 6
    ! Rx = [1 0 0; 0 cos(ang) sin(ang); 0 -sin(ang) cos(ang)]
    if (iface /= 1_8 .and. iface /= 4_8) then
      ang  = rotation(iface)
      cang = cos(ang);  sang = sin(ang)
      do k = 1_8, n
        tmp2 = cang * x_out(2,k) + sang * x_out(3,k)
        tmp3 =-sang * x_out(2,k) + cang * x_out(3,k)
        x_out(2,k) = tmp2;  x_out(3,k) = tmp3

        tmp2 = cang * nx_out(2,k) + sang * nx_out(3,k)
        tmp3 =-sang * nx_out(2,k) + cang * nx_out(3,k)
        nx_out(2,k) = tmp2;  nx_out(3,k) = tmp3
      end do
    end if

    deallocate(p_arr, t_arr)

  end subroutine face_map_nodes_r64

  ! ------------------------------------------------------------------
  ! face_map_point_r64  (private)
  ! Scalar version of face_map_nodes_r64: map one (x1,x2) -> r(3).
  ! ------------------------------------------------------------------
  subroutine face_map_point_r64(x1, x2, iface, ratio, rotation, r)
    integer(8), intent(in)  :: iface
    real(r64),  intent(in)  :: x1, x2, ratio, rotation(6)
    real(r64),  intent(out) :: r(3)
    real(r64) :: x1v(1), x2v(1), x_out(3,1), nx_out(3,1), sp(1), thv(1)
    x1v(1) = x1;  x2v(1) = x2
    call face_map_nodes_r64(x1v, x2v, 1_8, iface, ratio, rotation, x_out, nx_out, sp, thv)
    r = x_out(:,1)
  end subroutine face_map_point_r64

  ! ------------------------------------------------------------------
  ! edge_nodes_r64  (private)
  ! Compute nquad GL-mapped points along triangle edge from a to b
  ! in parameter space, using GL nodes s(nquad) in [-1,1].
  !
  ! e(:,i) = 0.5*(1-s(i))*a + 0.5*(1+s(i))*b
  ! ------------------------------------------------------------------
  subroutine edge_nodes_r64(a, b, s, nquad, e)
    integer(8), intent(in)  :: nquad
    real(r64),  intent(in)  :: a(2), b(2), s(nquad)
    real(r64),  intent(out) :: e(2,nquad)
    integer(8) :: i
    do i = 1_8, nquad
      e(:,i) = 0.5_r64*(1.0_r64 - s(i))*a + 0.5_r64*(1.0_r64 + s(i))*b
    end do
  end subroutine edge_nodes_r64

  ! ------------------------------------------------------------------
  ! areal_quad_triangle_r64
  ! Conical product rule quadrature on the 2D reference simplex
  !   {(u,v) : u >= 0, v >= 0, u+v <= 1}
  ! with vertices V1=(0,0), V2=(1,0), V3=(0,1).
  !
  ! Centroid decomposition into 3 sub-triangles:
  !   sub 1: {c, V1, V2},  sub 2: {c, V2, V3},  sub 3: {c, V3, V1}
  !   c = (1/3, 1/3)
  ! Each sub-triangle uses the conical map (t,s) -> c + t*(1-s)*ea + t*s*eb
  ! with 2D Jacobian t * |det([ea, eb])|  (constant in s, = 1/3 for all three).
  !
  ! Output is a DROP-IN REPLACEMENT for get_vr_nodes_wts_r64:
  !   uvs(2, 3*nq^2)  : (u,v) nodes on the reference simplex
  !   wts(3*nq^2)     : quadrature weights; sum(wts) = 1/2 (area of simplex)
  !
  ! Usage: identical to VR nodes — map through the same affine + surface
  ! Jacobian pipeline (face_map_nodes_r64, etc.) unchanged.
  !
  ! Inputs:
  !   nq              : GL order in both t and s; 3*nq^2 nodes total
  !   tgl(nq), wgl(nq): GL nodes/weights on [-1,1] from gauss_r64
  ! ------------------------------------------------------------------
  subroutine areal_quad_triangle_r64(nq, tgl, wgl, uvs, wts)
    integer(8), intent(in)  :: nq
    real(r64),  intent(in)  :: tgl(nq), wgl(nq)
    real(r64),  intent(out) :: uvs(2, 3*nq*nq)
    real(r64),  intent(out) :: wts(3*nq*nq)

    ! Reference simplex vertices and centroid
    real(r64), parameter :: V1(2) = [0.0_r64, 0.0_r64]
    real(r64), parameter :: V2(2) = [1.0_r64, 0.0_r64]
    real(r64), parameter :: V3(2) = [0.0_r64, 1.0_r64]

    real(r64) :: c(2), ea(2), eb(2), Jac, ti, wti, si, wsi
    integer(8) :: m, it, is, k

    c = (V1 + V2 + V3) / 3.0_r64      ! = (1/3, 1/3)

    do m = 1_8, 3_8
      select case (m)
        case (1); ea = V1 - c;  eb = V2 - c
        case (2); ea = V2 - c;  eb = V3 - c
        case (3); ea = V3 - c;  eb = V1 - c
      end select

      ! 2D Jacobian: |det([ea, eb])| = |ea(1)*eb(2) - ea(2)*eb(1)|
      ! For the reference simplex all three = 1/3
      Jac = abs(ea(1)*eb(2) - ea(2)*eb(1))

      do it = 1_8, nq
        ti  = 0.5_r64 * (tgl(it) + 1.0_r64)   ! map [-1,1] -> [0,1]
        wti = 0.5_r64 * wgl(it)
        do is = 1_8, nq
          si  = 0.5_r64 * (tgl(is) + 1.0_r64)
          wsi = 0.5_r64 * wgl(is)

          k = (m-1_8)*nq*nq + (it-1_8)*nq + is

          uvs(:,k) = c + ti*(1.0_r64-si)*ea + ti*si*eb
          wts(k)   = ti * wti * wsi * Jac
        end do
      end do
    end do

  end subroutine areal_quad_triangle_r64

  ! ================================================================
  ! r128 variants
  ! Note: VR nodes are tabulated r64 data, so r128 mesh builder
  ! requires nq > 0 (areal quadrature only).
  ! ================================================================

  ! ----------------------------------------------------------------
  ! areal_quad_triangle_r128
  ! ----------------------------------------------------------------
  subroutine areal_quad_triangle_r128(nq, tgl, wgl, uvs, wts)
    integer(8),  intent(in)  :: nq
    real(r128),  intent(in)  :: tgl(nq), wgl(nq)
    real(r128),  intent(out) :: uvs(2, 3*nq*nq)
    real(r128),  intent(out) :: wts(3*nq*nq)

    real(r128), parameter :: zero = 0.0_r128, one = 1.0_r128, third = one/3.0_r128
    real(r128) :: V1(2), V2(2), V3(2), c(2), ea(2), eb(2), Jac
    real(r128) :: ti, wti, si, wsi
    integer(8) :: m, it, is, k

    V1 = [zero, zero];  V2 = [one, zero];  V3 = [zero, one]
    c  = (V1 + V2 + V3) / 3.0_r128

    do m = 1_8, 3_8
      select case (m)
        case (1); ea = V1 - c;  eb = V2 - c
        case (2); ea = V2 - c;  eb = V3 - c
        case (3); ea = V3 - c;  eb = V1 - c
      end select
      Jac = abs(ea(1)*eb(2) - ea(2)*eb(1))

      do it = 1_8, nq
        ti  = 0.5_r128 * (tgl(it) + one)
        wti = 0.5_r128 * wgl(it)
        do is = 1_8, nq
          si  = 0.5_r128 * (tgl(is) + one)
          wsi = 0.5_r128 * wgl(is)
          k = (m-1_8)*nq*nq + (it-1_8)*nq + is
          uvs(:,k) = c + ti*(one-si)*ea + ti*si*eb
          wts(k)   = ti * wti * wsi * Jac
        end do
      end do
    end do

  end subroutine areal_quad_triangle_r128

  ! ----------------------------------------------------------------
  ! ellipsoidparam_r128  (private)
  ! ----------------------------------------------------------------
  subroutine ellipsoidparam_r128(p_ang, t_ang, n, ratio, x_out, nx_out, sp)
    integer(8),  intent(in)  :: n
    real(r128),  intent(in)  :: p_ang(n), t_ang(n), ratio
    real(r128),  intent(out) :: x_out(3,n), nx_out(3,n), sp(n)

    integer(8) :: k
    real(r128) :: p, t, f, ft, fp
    real(r128) :: rt(3), rp(3), nxk(3), spk

    do k = 1_8, n
      p = p_ang(k);  t = t_ang(k)
      f  = 1.0_r128 / sqrt(cos(p)**2*cos(t)**2 + ratio**2*sin(p)**2*cos(t)**2 &
                           + ratio**2*sin(t)**2)
      ft = -0.5_r128 * f**3 * ( -2.0_r128*cos(p)**2*cos(t)*sin(t) &
                                 -2.0_r128*ratio**2*sin(p)**2*cos(t)*sin(t) &
                                 +2.0_r128*ratio**2*sin(t)*cos(t) )
      fp = -0.5_r128 * f**3 * ( -2.0_r128*cos(p)*sin(p)*cos(t)**2 &
                                 +2.0_r128*ratio**2*sin(p)*cos(p)*cos(t)**2 )

      x_out(1,k) = f * cos(t) * cos(p)
      x_out(2,k) = f * cos(t) * sin(p)
      x_out(3,k) = f * sin(t)

      rt(1) = ft*cos(t)*cos(p) - f*sin(t)*cos(p)
      rt(2) = ft*cos(t)*sin(p) - f*sin(t)*sin(p)
      rt(3) = ft*sin(t)        + f*cos(t)

      rp(1) = fp*cos(t)*cos(p) - f*cos(t)*sin(p)
      rp(2) = fp*cos(t)*sin(p) + f*cos(t)*cos(p)
      rp(3) = fp*sin(t)

      nxk(1) = rp(2)*rt(3) - rp(3)*rt(2)
      nxk(2) = rp(3)*rt(1) - rp(1)*rt(3)
      nxk(3) = rp(1)*rt(2) - rp(2)*rt(1)

      spk = sqrt(nxk(1)**2 + nxk(2)**2 + nxk(3)**2)
      sp(k) = spk
      nx_out(:,k) = nxk / spk
    end do

  end subroutine ellipsoidparam_r128

  ! ----------------------------------------------------------------
  ! face_map_nodes_r128  (private)
  ! ----------------------------------------------------------------
  subroutine face_map_nodes_r128(x1, x2, n, iface, ratio, rotation, x_out, nx_out, sp, th)
    integer(8),  intent(in)  :: n, iface
    real(r128),  intent(in)  :: x1(n), x2(n), ratio, rotation(6)
    real(r128),  intent(out) :: x_out(3,n), nx_out(3,n), sp(n), th(n)

    real(r128) :: pi
    real(r128) :: qnrm, q1, q2, q3, azimuth, elevation
    real(r128) :: ang, cang, sang, tmp2, tmp3
    real(r128), allocatable :: p_arr(:), t_arr(:)
    integer(8) :: k

    pi = 4.0_r128 * atan(1.0_r128)
    allocate(p_arr(n), t_arr(n))

    do k = 1_8, n
      qnrm = sqrt(1.0_r128 + x1(k)**2 + x2(k)**2)
      q1   = 1.0_r128 / qnrm
      q2   = x1(k)   / qnrm
      q3   = x2(k)   / qnrm

      azimuth   = atan2(q2, q1)
      elevation = atan2(q3, sqrt(q1**2 + q2**2))
      th(k) = azimuth

      if (iface == 1_8 .or. iface == 4_8) then
        p_arr(k) = azimuth + floor(real(iface,r128) / 3.0_r128) * pi
      else
        p_arr(k) = azimuth + 0.5_r128 * pi
      end if
      t_arr(k) = elevation
    end do

    call ellipsoidparam_r128(p_arr, t_arr, n, ratio, x_out, nx_out, sp)

    if (iface /= 1_8 .and. iface /= 4_8) then
      ang  = rotation(iface)
      cang = cos(ang);  sang = sin(ang)
      do k = 1_8, n
        tmp2 = cang * x_out(2,k) + sang * x_out(3,k)
        tmp3 =-sang * x_out(2,k) + cang * x_out(3,k)
        x_out(2,k) = tmp2;  x_out(3,k) = tmp3

        tmp2 = cang * nx_out(2,k) + sang * nx_out(3,k)
        tmp3 =-sang * nx_out(2,k) + cang * nx_out(3,k)
        nx_out(2,k) = tmp2;  nx_out(3,k) = tmp3
      end do
    end if

    deallocate(p_arr, t_arr)

  end subroutine face_map_nodes_r128

  ! ----------------------------------------------------------------
  ! face_map_point_r128  (private)
  ! ----------------------------------------------------------------
  subroutine face_map_point_r128(x1, x2, iface, ratio, rotation, r)
    integer(8),  intent(in)  :: iface
    real(r128),  intent(in)  :: x1, x2, ratio, rotation(6)
    real(r128),  intent(out) :: r(3)
    real(r128) :: x1v(1), x2v(1), x_out(3,1), nx_out(3,1), spv(1), thv(1)
    x1v(1) = x1;  x2v(1) = x2
    call face_map_nodes_r128(x1v, x2v, 1_8, iface, ratio, rotation, x_out, nx_out, spv, thv)
    r = x_out(:,1)
  end subroutine face_map_point_r128

  ! ----------------------------------------------------------------
  ! edge_nodes_r128  (private)
  ! ----------------------------------------------------------------
  subroutine edge_nodes_r128(a, b, s, nquad, e)
    integer(8),  intent(in)  :: nquad
    real(r128),  intent(in)  :: a(2), b(2), s(nquad)
    real(r128),  intent(out) :: e(2,nquad)
    integer(8) :: i
    do i = 1_8, nquad
      e(:,i) = 0.5_r128*(1.0_r128 - s(i))*a + 0.5_r128*(1.0_r128 + s(i))*b
    end do
  end subroutine edge_nodes_r128

  ! ----------------------------------------------------------------
  ! create_ellipsoid_tri_mesh_r128
  ! Requires nq > 0 (areal quadrature; VR nodes have no r128 data).
  ! Interface identical to r64 version.
  ! ----------------------------------------------------------------
  subroutine create_ellipsoid_tri_mesh_r128(mp, np, p, nq, ratio, nquad_bdry, nvr, ntri, &
      x, nx, w, xbd, tri2face, tri2cell, tri_vert, ptr)
    integer(8),  intent(in)    :: mp, np, p, nq, nquad_bdry, nvr, ntri
    real(r128),  intent(in)    :: ratio
    real(r128),  intent(inout) :: x(3,nvr,ntri), nx(3,nvr,ntri), w(nvr,ntri)
    real(r128),  intent(inout) :: xbd(3,3*nquad_bdry,ntri)
    integer(8),  intent(inout) :: tri2face(ntri)
    integer(8),  intent(inout) :: tri2cell(2,ntri)
    real(r128),  intent(inout) :: tri_vert(3,3,ntri)
    integer(8),  intent(inout) :: ptr(ntri+1)

    real(r128) :: pi
    real(r128) :: rotation(6)
    real(r128) :: V(2,3)
    real(r128) :: vll(2), vlr(2), vur(2), vul(2)
    real(r128) :: rll(3), rlr(3), rur(3), rul(3)
    real(r128) :: r_v(3)
    real(r128) :: d1, d2, jac, wwt_k, det2x2
    integer(8) :: iface, i, j, kk, iv, k_idx, itri, nb3

    real(r128), allocatable :: uvs(:,:), vr_wts(:)
    real(r128), allocatable :: xg(:), wg(:), Dg(:,:)
    real(r128), allocatable :: tgl_a(:), wgl_a(:), Dgl_a(:,:)
    real(r128), allocatable :: x12(:,:)
    real(r128), allocatable :: xt(:,:), nxt(:,:), spt(:), tht(:)
    real(r128), allocatable :: e1(:,:), e2(:,:), e3(:,:)
    real(r128), allocatable :: xb1(:,:), nxb1(:,:), spb1(:), thb1(:)
    real(r128), allocatable :: xb2(:,:), nxb2(:,:), spb2(:), thb2(:)
    real(r128), allocatable :: xb3(:,:), nxb3(:,:), spb3(:), thb3(:)

    pi = 4.0_r128 * atan(1.0_r128)
    rotation(1) =  0.0_r128
    rotation(2) =  0.0_r128
    rotation(3) = -pi / 2.0_r128
    rotation(4) =  0.0_r128
    rotation(5) = -pi
    rotation(6) = -1.5_r128 * pi

    allocate(uvs(2,nvr), vr_wts(nvr))
    allocate(xg(nquad_bdry), wg(nquad_bdry), Dg(nquad_bdry,nquad_bdry))
    allocate(x12(2,nvr))
    allocate(xt(3,nvr), nxt(3,nvr), spt(nvr), tht(nvr))
    allocate(e1(2,nquad_bdry), e2(2,nquad_bdry), e3(2,nquad_bdry))
    allocate(xb1(3,nquad_bdry), nxb1(3,nquad_bdry), spb1(nquad_bdry), thb1(nquad_bdry))
    allocate(xb2(3,nquad_bdry), nxb2(3,nquad_bdry), spb2(nquad_bdry), thb2(nquad_bdry))
    allocate(xb3(3,nquad_bdry), nxb3(3,nquad_bdry), spb3(nquad_bdry), thb3(nquad_bdry))

    ! Areal quadrature (nq > 0 required; VR nodes have no r128 tabulated data)
    allocate(tgl_a(nq), wgl_a(nq), Dgl_a(nq,nq))
    call gauss_r128(nq, tgl_a, wgl_a, Dgl_a)
    call areal_quad_triangle_r128(nq, tgl_a, wgl_a, uvs, vr_wts)
    deallocate(tgl_a, wgl_a, Dgl_a)

    call gauss_r128(nquad_bdry, xg, wg, Dg)

    itri = 0_8
    do iface = 1_8, 6_8
      do j = 1_8, np
        do i = 1_8, mp

          vll(1) = 2.0_r128*(j-1)/np - 1.0_r128;  vll(2) = 2.0_r128*(i-1)/mp - 1.0_r128
          vlr(1) = 2.0_r128*j/np     - 1.0_r128;  vlr(2) = 2.0_r128*(i-1)/mp - 1.0_r128
          vul(1) = 2.0_r128*(j-1)/np - 1.0_r128;  vul(2) = 2.0_r128*i/mp     - 1.0_r128
          vur(1) = 2.0_r128*j/np     - 1.0_r128;  vur(2) = 2.0_r128*i/mp     - 1.0_r128

          call face_map_point_r128(vll(1), vll(2), iface, ratio, rotation, rll)
          call face_map_point_r128(vur(1), vur(2), iface, ratio, rotation, rur)
          call face_map_point_r128(vul(1), vul(2), iface, ratio, rotation, rul)
          call face_map_point_r128(vlr(1), vlr(2), iface, ratio, rotation, rlr)

          d1 = sum((rll - rur)**2)
          d2 = sum((rul - rlr)**2)

          do kk = 1_8, 2_8
            itri = itri + 1_8

            if (d1 + 1.0e-13_r128 > d2) then
              if (kk == 1) then
                V(:,1) = vll;  V(:,2) = vlr;  V(:,3) = vul
              else
                V(:,1) = vur;  V(:,2) = vul;  V(:,3) = vlr
              end if
            else
              if (kk == 1) then
                V(:,1) = vul;  V(:,2) = vll;  V(:,3) = vur
              else
                V(:,1) = vlr;  V(:,2) = vur;  V(:,3) = vll
              end if
            end if

            do iv = 1_8, nvr
              x12(1,iv) = V(1,1) + (V(1,2)-V(1,1))*uvs(1,iv) + (V(1,3)-V(1,1))*uvs(2,iv)
              x12(2,iv) = V(2,1) + (V(2,2)-V(2,1))*uvs(1,iv) + (V(2,3)-V(2,1))*uvs(2,iv)
            end do

            call face_map_nodes_r128(x12(1,:), x12(2,:), nvr, iface, ratio, rotation, &
                                     xt, nxt, spt, tht)

            det2x2 = (V(1,2)-V(1,1))*(V(2,3)-V(2,1)) - (V(1,3)-V(1,1))*(V(2,2)-V(2,1))
            jac = abs(det2x2)

            do iv = 1_8, nvr
              wwt_k = cos(tht(iv)) / (1.0_r128 + x12(1,iv)**2 + x12(2,iv)**2)
              w(iv, itri) = spt(iv) * wwt_k * vr_wts(iv) * jac
            end do

            x(:,:,itri)  = xt
            nx(:,:,itri) = nxt

            call edge_nodes_r128(V(:,1), V(:,2), xg, nquad_bdry, e1)
            call edge_nodes_r128(V(:,2), V(:,3), xg, nquad_bdry, e2)
            call edge_nodes_r128(V(:,3), V(:,1), xg, nquad_bdry, e3)

            call face_map_nodes_r128(e1(1,:), e1(2,:), nquad_bdry, iface, ratio, rotation, &
                                     xb1, nxb1, spb1, thb1)
            call face_map_nodes_r128(e2(1,:), e2(2,:), nquad_bdry, iface, ratio, rotation, &
                                     xb2, nxb2, spb2, thb2)
            call face_map_nodes_r128(e3(1,:), e3(2,:), nquad_bdry, iface, ratio, rotation, &
                                     xb3, nxb3, spb3, thb3)

            nb3 = nquad_bdry
            xbd(:, 1:nb3,         itri) = xb1
            xbd(:, nb3+1:2*nb3,   itri) = xb2
            xbd(:, 2*nb3+1:3*nb3, itri) = xb3

            call face_map_point_r128(V(1,1), V(2,1), iface, ratio, rotation, r_v)
            tri_vert(:,1,itri) = r_v
            call face_map_point_r128(V(1,2), V(2,2), iface, ratio, rotation, r_v)
            tri_vert(:,2,itri) = r_v
            call face_map_point_r128(V(1,3), V(2,3), iface, ratio, rotation, r_v)
            tri_vert(:,3,itri) = r_v

            tri2face(itri)   = iface
            tri2cell(1,itri) = i
            tri2cell(2,itri) = j

          end do   ! kk
        end do   ! i
      end do   ! j
    end do   ! iface

    do k_idx = 1_8, ntri + 1_8
      ptr(k_idx) = 1_8 + (k_idx - 1_8) * nvr
    end do

    deallocate(uvs, vr_wts)
    deallocate(xg, wg, Dg)
    deallocate(x12, xt, nxt, spt, tht)
    deallocate(e1, e2, e3)
    deallocate(xb1, nxb1, spb1, thb1)
    deallocate(xb2, nxb2, spb2, thb2)
    deallocate(xb3, nxb3, spb3, thb3)

  end subroutine create_ellipsoid_tri_mesh_r128

end module ellipsoid_mesh_mod
