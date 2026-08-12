module torus_mesh_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use linequaaadrature_mod, only: r64
  use stellarator_mesh_mod, only: get_vr_nodes_wts_r64
  implicit none
  private

  public :: create_torus_tri_mesh_r64

contains

  subroutine create_torus_tri_mesh_r64(mp, np, p, radii, scales, nosc, &
      orig, ntri, sx, snx, sw, rts, rps, ier)
    integer(8), intent(in) :: mp, np, p, nosc, orig, ntri
    real(r64), intent(in) :: radii(3), scales(3)
    real(r64), intent(inout) :: sx(3,ntri*p*(p+1_8)/2_8)
    real(r64), intent(inout) :: snx(3,ntri*p*(p+1_8)/2_8)
    real(r64), intent(inout) :: sw(ntri*p*(p+1_8)/2_8)
    real(r64), intent(inout) :: rts(3,ntri*p*(p+1_8)/2_8)
    real(r64), intent(inout) :: rps(3,ntri*p*(p+1_8)/2_8)
    integer(8), intent(inout) :: ier

    integer(8) :: hdim, i, j, side, k, itri, off
    real(r64) :: pi, ts1, ts2, phi0, phi1, theta0, theta1
    real(r64) :: d1, d2, phi, theta, speed
    real(r64) :: vll(2), vlr(2), vul(2), vur(2)
    real(r64) :: vertices(2,3,2), xll(3), xlr(3), xul(3), xur(3)
    real(r64) :: x(3), xphi(3), xtheta(3), normal(3)
    real(r64), allocatable :: uvs(:,:), wts(:)

    sx = 0.0_r64
    snx = 0.0_r64
    sw = 0.0_r64
    rts = 0.0_r64
    rps = 0.0_r64
    ier = 1_8

    if (mp < 1_8 .or. np < 1_8 .or. p < 1_8 .or. p > 20_8) return
    if (orig /= 0_8 .and. orig /= 1_8) return
    if (ntri /= 2_8*mp*np) return

    ier = 2_8
    if (nosc < 0_8) return
    if (.not. all(ieee_is_finite(radii))) return
    if (.not. all(ieee_is_finite(scales))) return
    if (radii(1) <= 0.0_r64 .or. radii(2) <= 0.0_r64) return
    if (any(scales <= 0.0_r64)) return

    hdim = p*(p+1_8)/2_8
    allocate(uvs(2,hdim), wts(hdim))
    call get_vr_nodes_wts_r64(p-1_8, hdim, uvs, wts)

    pi = acos(-1.0_r64)
    ts1 = pi/real(np,r64)
    ts2 = pi/real(mp,r64)
    itri = 0_8

    do j = 1_8, np
      phi0 = real(orig-2_8+2_8*j,r64)*ts1
      phi1 = real(orig+2_8*j,r64)*ts1
      do i = 1_8, mp
        theta0 = real(orig-2_8+2_8*i,r64)*ts2
        theta1 = real(orig+2_8*i,r64)*ts2

        vll = [phi0,theta0]
        vlr = [phi1,theta0]
        vul = [phi0,theta1]
        vur = [phi1,theta1]
        call torus_param_r64(radii,scales,nosc,phi0,theta0,xll,xphi,xtheta)
        call torus_param_r64(radii,scales,nosc,phi1,theta0,xlr,xphi,xtheta)
        call torus_param_r64(radii,scales,nosc,phi0,theta1,xul,xphi,xtheta)
        call torus_param_r64(radii,scales,nosc,phi1,theta1,xur,xphi,xtheta)
        d1 = sum((xll-xur)**2)
        d2 = sum((xul-xlr)**2)

        if (d1+1.0e-13_r64 > d2) then
          vertices(:,1,1) = vll
          vertices(:,2,1) = vlr
          vertices(:,3,1) = vul
          vertices(:,1,2) = vur
          vertices(:,2,2) = vul
          vertices(:,3,2) = vlr
        else
          vertices(:,1,1) = vul
          vertices(:,2,1) = vll
          vertices(:,3,1) = vur
          vertices(:,1,2) = vlr
          vertices(:,2,2) = vur
          vertices(:,3,2) = vll
        end if

        do side = 1_8, 2_8
          itri = itri+1_8
          off = (itri-1_8)*hdim
          do k = 1_8, hdim
            phi = vertices(1,1,side) &
                 +(vertices(1,2,side)-vertices(1,1,side))*uvs(1,k) &
                 +(vertices(1,3,side)-vertices(1,1,side))*uvs(2,k)
            theta = vertices(2,1,side) &
                   +(vertices(2,2,side)-vertices(2,1,side))*uvs(1,k) &
                   +(vertices(2,3,side)-vertices(2,1,side))*uvs(2,k)
            call torus_param_r64(radii,scales,nosc,phi,theta,x,xphi,xtheta)
            rts(:,off+k) = 2.0_r64*ts2*xtheta
            rps(:,off+k) = 2.0_r64*ts1*xphi
            normal(1) = rps(2,off+k)*rts(3,off+k) &
                       -rps(3,off+k)*rts(2,off+k)
            normal(2) = rps(3,off+k)*rts(1,off+k) &
                       -rps(1,off+k)*rts(3,off+k)
            normal(3) = rps(1,off+k)*rts(2,off+k) &
                       -rps(2,off+k)*rts(1,off+k)
            speed = sqrt(sum(normal*normal))
            if (.not. ieee_is_finite(speed) .or. speed <= 0.0_r64) then
              ier = 3_8
              return
            end if
            sx(:,off+k) = x
            snx(:,off+k) = normal/speed
            sw(off+k) = speed*wts(k)
          end do
        end do
      end do
    end do

    ier = 0_8
  end subroutine create_torus_tri_mesh_r64

  subroutine torus_param_r64(radii, scales, nosc, phi, theta, x, xphi, xtheta)
    real(r64), intent(in) :: radii(3), scales(3), phi, theta
    integer(8), intent(in) :: nosc
    real(r64), intent(out) :: x(3), xphi(3), xtheta(3)

    real(r64) :: radius, radius_phi, radius_theta, wave

    wave = real(nosc,r64)*phi
    radius = radii(2)+radii(1)*cos(theta)+radii(3)*cos(wave)
    radius_phi = -real(nosc,r64)*radii(3)*sin(wave)
    radius_theta = -radii(1)*sin(theta)

    x(1) = scales(1)*radius*cos(phi)
    x(2) = scales(2)*radius*sin(phi)
    x(3) = scales(3)*radii(1)*sin(theta)
    xphi(1) = scales(1)*(radius_phi*cos(phi)-radius*sin(phi))
    xphi(2) = scales(2)*(radius_phi*sin(phi)+radius*cos(phi))
    xphi(3) = 0.0_r64
    xtheta(1) = scales(1)*radius_theta*cos(phi)
    xtheta(2) = scales(2)*radius_theta*sin(phi)
    xtheta(3) = scales(3)*radii(1)*cos(theta)
  end subroutine torus_param_r64

end module torus_mesh_mod
