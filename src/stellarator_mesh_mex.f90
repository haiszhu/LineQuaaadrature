subroutine lqsm_stellarator_mesh_init_mex(mp, np, p, nfp, nmode, mn, rc, zs, &
    restol, cap, charts, nchart, ntri, ier)
  use stellarator_mesh_mod, only: stellarator_mesh_init_r64
  implicit none
  integer(8), intent(in) :: mp, np, p, nfp, nmode, cap
  integer(8), intent(in) :: mn(2,nmode)
  real(8), intent(in) :: rc(nmode), zs(nmode), restol
  real(8), intent(inout) :: charts(6,cap)
  integer(8), intent(inout) :: nchart, ntri, ier

  call stellarator_mesh_init_r64(mp, np, p, nfp, nmode, mn, rc, zs, &
       restol, cap, charts, nchart, ntri, ier)
end subroutine lqsm_stellarator_mesh_init_mex

subroutine lqsm_create_stellarator_tri_mesh_mex(mp, np, p, nfp, nmode, mn, &
    rc, zs, nchart, charts, ntri, sx, snx, sw, rts, rps, ier)
  use stellarator_mesh_mod, only: create_stellarator_tri_mesh_r64
  implicit none
  integer(8), intent(in) :: mp, np, p, nfp, nmode, nchart, ntri
  integer(8), intent(in) :: mn(2,nmode)
  real(8), intent(in) :: rc(nmode), zs(nmode), charts(6,nchart)
  real(8), intent(inout) :: sx(3,ntri*p*(p+1_8)/2_8)
  real(8), intent(inout) :: snx(3,ntri*p*(p+1_8)/2_8)
  real(8), intent(inout) :: sw(ntri*p*(p+1_8)/2_8)
  real(8), intent(inout) :: rts(3,ntri*p*(p+1_8)/2_8)
  real(8), intent(inout) :: rps(3,ntri*p*(p+1_8)/2_8)
  integer(8), intent(inout) :: ier

  call create_stellarator_tri_mesh_r64(mp, np, p, nfp, nmode, mn, rc, zs, &
       nchart, charts, ntri, sx, snx, sw, rts, rps, ier)
end subroutine lqsm_create_stellarator_tri_mesh_mex

subroutine lqsm_stellarator_tri_uv2x_mex(mp, np, p, nfp, nmode, mn, rc, zs, &
    nchart, charts, itri, nuv, uv, x, ier)
  use stellarator_mesh_mod, only: stellarator_tri_uv2x_r64
  implicit none
  integer(8), intent(in) :: mp, np, p, nfp, nmode, nchart, itri, nuv
  integer(8), intent(in) :: mn(2,nmode)
  real(8), intent(in) :: rc(nmode), zs(nmode), charts(6,nchart), uv(2,nuv)
  real(8), intent(inout) :: x(3,nuv)
  integer(8), intent(inout) :: ier

  call stellarator_tri_uv2x_r64(mp, np, p, nfp, nmode, mn, rc, zs, &
       nchart, charts, itri, nuv, uv, x, ier)
end subroutine lqsm_stellarator_tri_uv2x_mex
