subroutine lqtm_create_torus_tri_mesh_mex(mp, np, p, radii, scales, nosc, &
    orig, ntri, sx, snx, sw, rts, rps, ier)
  use torus_mesh_mod, only: create_torus_tri_mesh_r64
  implicit none
  integer(8), intent(in) :: mp, np, p, nosc, orig, ntri
  real(8), intent(in) :: radii(3), scales(3)
  real(8), intent(inout) :: sx(3,ntri*p*(p+1_8)/2_8)
  real(8), intent(inout) :: snx(3,ntri*p*(p+1_8)/2_8)
  real(8), intent(inout) :: sw(ntri*p*(p+1_8)/2_8)
  real(8), intent(inout) :: rts(3,ntri*p*(p+1_8)/2_8)
  real(8), intent(inout) :: rps(3,ntri*p*(p+1_8)/2_8)
  integer(8), intent(inout) :: ier

  call create_torus_tri_mesh_r64(mp,np,p,radii,scales,nosc,orig,ntri, &
      sx,snx,sw,rts,rps,ier)
end subroutine lqtm_create_torus_tri_mesh_mex
