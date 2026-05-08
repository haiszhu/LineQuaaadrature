! test_solid_angle.f90
! Fortran port of test.m: accuracy test for the near-field solid-angle
! correction to the Laplace double-layer potential on an ellipsoid.
!
! Two checks:
!   Part A — sanity: sum IalphaAsvestas over all triangles for 2 fixed targets.
!              exterior -> omega ≈ 0,  interior -> omega ≈ 4*pi.
!
!   Part B — full DLP accuracy (Option A, direct sum replacing Lap3dDLPfmm):
!              1. u_direct(i) = -sum_{all j} K(r_ti, r_sj) * w_j
!              2. Near-field correction loop (mirrors test.m):
!                   K_corr(near) += -IalphaAsv/(4pi) - K_ij_naive*ones
!              3. u = u_direct - K_corr
!              Expected: |u| ~ machine eps for ratio=1 (sphere).
!
! Compile via Makefile:  make test
! Run:                   ./build/test_solid_angle

program test_solid_angle
  use linequaaadrature_mod, only: r64
  use ellipsoid_mesh_mod,   only: create_ellipsoid_tri_mesh_r64
  use solidangle_mod,       only: evaluate_solid_angle_integral_r64
  use lap3d_mod,            only: lap3ddlpmat_r64, lap3ddlp_direct_r64
  use omp_lib,              only: omp_get_wtime
  implicit none

  ! ---- problem parameters ----
  integer(8), parameter :: order      = 14_8
  integer(8), parameter :: mp         = 8_8
  integer(8), parameter :: np         = 8_8
  real(r64),  parameter :: ratio      = 1.0_r64    ! unit sphere
  integer(8), parameter :: nquad_bdry = order + 8_8
  ! quad_mode: 0 = Vioreanu-Rokhlin,  nq > 0 = areal rule (nvr = 3*nq^2)
  integer(8), parameter :: nq         = 0_8
  integer(8), parameter :: nvr        = order*(order+1_8)/2_8  ! VR: set nq=0
  ! integer(8), parameter :: nq         = 12_8
  ! integer(8), parameter :: nvr        = 3_8*nq*nq             ! areal: uncomment both
  integer(8), parameter :: ntri       = 12_8*mp*np
  integer(8), parameter :: nplotpts   = 20_8        ! grid points per axis

  real(r64), parameter :: pi = 3.14159265358979323846264338327950288_r64

  ! ---- mesh ----
  real(r64),  allocatable :: x(:,:,:), nx_s(:,:,:), w(:,:)
  real(r64),  allocatable :: xbd(:,:,:), tri_vert(:,:,:)
  integer(8), allocatable :: tri2face(:), tri2cell(:,:), ptr(:)

  ! ---- all sources flattened ----
  integer(8)              :: N_src
  real(r64),  allocatable :: sx_all(:,:), snx_all(:,:), sw_all(:), sigma_all(:)

  ! ---- target points ----
  integer(8)              :: ntarget
  real(r64),  allocatable :: tx(:,:)

  ! ---- result vectors ----
  real(r64),  allocatable :: u(:), K_corr(:)

  ! ---- per-triangle work (all private in OMP loop) ----
  integer(8)              :: ntc
  integer(8), allocatable :: near_idx(:)
  real(r64),  allocatable :: tcj(:,:), IalphaAsv(:), Kmat(:,:)
  real(r64)               :: delta_k

  ! ---- Part A sanity check ----
  real(r64) :: tx2(3,2), omega_sum(2), ialpha2(2)

  ! ---- misc ----
  real(r64) :: domain(6), px, py, pz, qpoint(3), qradii
  real(r64) :: err_max, err_l2
  real(8)   :: t0, t1
  integer(8) :: k, i, ix, iy, iz

  logical :: use_nearroot
  character(len=32) :: mode
  integer :: nargs

  ! ==============================================================
  use_nearroot = .true.

  nargs = command_argument_count()
  if (nargs >= 1) then
    call get_command_argument(1, mode)
    select case (trim(adjustl(mode)))
    case ('1', 'true', 'TRUE', 'True', 'nearroot', 'NEARROOT', 'Nearroot')
      use_nearroot = .true.
    case ('0', 'false', 'FALSE', 'False', 'adaptive', 'ADAPTIVE', 'Adaptive')
      use_nearroot = .false.
    case default
      write(*,'(A,A)') 'Unknown mode: ', trim(mode)
      write(*,'(A)')   'Use: ./build/test_solid_angle adaptive'
      write(*,'(A)')   '  or ./build/test_solid_angle nearroot'
      stop 1
    end select
  end if
  write(*,'(/,A)') '=== test_solid_angle ==='
  write(*,'(A,I0,A,I0,A,I0,A,F4.1)') &
      'order=', order, '  mp=np=', mp, '  ntri=', ntri, '  ratio=', ratio
  write(*,'(A,L1)') 'use_nearroot = ', use_nearroot

  ! ==============================================================
  ! Build ellipsoid mesh
  ! ==============================================================
  allocate(x(3,nvr,ntri), nx_s(3,nvr,ntri), w(nvr,ntri))
  allocate(xbd(3,3*nquad_bdry,ntri), tri_vert(3,3,ntri))
  allocate(tri2face(ntri), tri2cell(2,ntri), ptr(ntri+1))

  t0 = omp_get_wtime()
  call create_ellipsoid_tri_mesh_r64(mp, np, order, nq, ratio, nquad_bdry, nvr, ntri, &
      x, nx_s, w, xbd, tri2face, tri2cell, tri_vert, ptr)
  t1 = omp_get_wtime()
  write(*,'(A,F7.2,A)') 'Mesh built in ', t1-t0, ' s'

  ! ==============================================================
  ! Part A: sanity check — sum solid angle over all triangles
  !   exterior point -> 0,  interior point -> 4*pi
  ! ==============================================================
  write(*,'(/,A)') '--- Part A: solid angle sum ---'

  ! Two exterior points: one far, one near the surface.
  ! Note: the Asvestas kernel has a 1/(1-qhat.rhat) singularity that
  ! makes it valid only for exterior targets; interior is not tested here.
  tx2(:,1) = [1.5_r64, 0.0_r64, 0.0_r64]   ! exterior, far
  tx2(:,2) = [1.1_r64, 0.0_r64, 0.0_r64]   ! exterior, near surface

  omega_sum = 0.0_r64
  do k = 1_8, ntri
    ialpha2 = 0.0_r64
    call evaluate_solid_angle_integral_r64(2_8, tx2, nvr,               &
        x(:,:,k), nx_s(:,:,k), w(:,k),                                  &
        tri_vert(:,:,k), 3_8*nquad_bdry, xbd(:,:,k), use_nearroot, ialpha2)
    omega_sum = omega_sum + ialpha2
  end do

  write(*,'(A,ES10.3)') 'far exterior  (1.5,0,0): |omega_sum| = ', abs(omega_sum(1))
  write(*,'(A,ES10.3)') 'near exterior (1.1,0,0): |omega_sum| = ', abs(omega_sum(2))

  ! ==============================================================
  ! Part B: full DLP accuracy (Option A)
  ! ==============================================================
  write(*,'(/,A)') '--- Part B: DLP accuracy (direct sum + near correction) ---'

  ! Flatten all sources
  N_src = ntri * nvr
  allocate(sx_all(3,N_src), snx_all(3,N_src), sw_all(N_src), sigma_all(N_src))
  sx_all    = reshape(x,    [3_8, N_src])
  snx_all   = reshape(nx_s, [3_8, N_src])
  sw_all    = reshape(w,    [N_src])
  sigma_all = -1.0_r64
  write(*,'(A,I0)') 'N_src = ', N_src

  ! Generate exterior target grid
  domain = [-0.5_r64, 2.0_r64, -0.5_r64, 2.0_r64, -1.25_r64, 1.25_r64]
  allocate(tx(3, nplotpts**3))
  ntarget = 0_8
  do iz = 1_8, nplotpts
    do iy = 1_8, nplotpts
      do ix = 1_8, nplotpts
        px = domain(1) + (domain(2)-domain(1)) * (ix-1) / real(nplotpts-1, r64)
        py = domain(3) + (domain(4)-domain(3)) * (iy-1) / real(nplotpts-1, r64)
        pz = domain(5) + (domain(6)-domain(5)) * (iz-1) / real(nplotpts-1, r64)
        if (px**2 + (py*ratio)**2 + (pz*ratio)**2 > 1.0_r64 + 1.0e-8_r64) then
          ntarget = ntarget + 1_8
          tx(1,ntarget) = px
          tx(2,ntarget) = py
          tx(3,ntarget) = pz
        end if
      end do
    end do
  end do
  write(*,'(A,I0,A)') 'ntarget = ', ntarget, ' exterior points'

  ! Direct DLP: u(i) = -sum_{all j} K(r_ti, r_sj) * w_j
  allocate(u(ntarget), K_corr(ntarget))
  K_corr = 0.0_r64

  write(*,'(A)') 'Direct DLP sum ...'
  t0 = omp_get_wtime()
  call lap3ddlp_direct_r64(ntarget, tx(:,1:ntarget), N_src, sx_all, snx_all, sw_all, sigma_all, u)
  t1 = omp_get_wtime()
  write(*,'(A,F7.2,A)') '  done in ', t1-t0, ' s'

  ! Near-field correction loop (parallelised over panels)
  ! near_idx is allocated per iteration inside the OMP loop (private to each thread).

  write(*,'(A)') 'Near-field correction ...'
  t0 = omp_get_wtime()

  !$omp parallel do schedule(dynamic) &
  !$omp   default(shared) &
  !$omp   private(ntc, near_idx, tcj, IalphaAsv, Kmat, qpoint, qradii, i, delta_k)
  do k = 1_8, ntri
    qpoint = sum(x(:,:,k), dim=2) / real(nvr, r64)
    qradii = 1.75_r64 * sqrt(sum(w(:,k)))

    allocate(near_idx(ntarget))
    ntc = 0_8
    do i = 1_8, ntarget
      if (  (tx(1,i)-qpoint(1))**2 &
          + (tx(2,i)-qpoint(2))**2 &
          + (tx(3,i)-qpoint(3))**2 < qradii**2 ) then
        ntc = ntc + 1_8
        near_idx(ntc) = i
      end if
    end do

    if (ntc > 0_8) then
      allocate(tcj(3,ntc), IalphaAsv(ntc), Kmat(ntc,nvr))

      do i = 1_8, ntc
        tcj(:,i) = tx(:, near_idx(i))
      end do

      IalphaAsv = 0.0_r64
      call evaluate_solid_angle_integral_r64(ntc, tcj, nvr,             &
          x(:,:,k), nx_s(:,:,k), w(:,k),                                &
          tri_vert(:,:,k), 3_8*nquad_bdry, xbd(:,:,k), use_nearroot, IalphaAsv)

      call lap3ddlpmat_r64(ntc, tcj, nvr, x(:,:,k), nx_s(:,:,k), w(:,k), Kmat)

      do i = 1_8, ntc
        delta_k = - IalphaAsv(i) / (4.0_r64*pi) - sum(Kmat(i,:))
        !$omp atomic
        K_corr(near_idx(i)) = K_corr(near_idx(i)) + delta_k
      end do

      deallocate(tcj, IalphaAsv, Kmat)
    end if

    deallocate(near_idx)
  end do
  !$omp end parallel do

  t1 = omp_get_wtime()
  write(*,'(A,F7.2,A)') '  done in ', t1-t0, ' s'

  ! Final:  u = u_direct + K_corr * (-1)   (mirrors test.m)
  do i = 1_8, ntarget
    u(i) = u(i) - K_corr(i)
  end do

  ! Report: DLP(1) = 0 outside
  err_max = maxval(abs(u(1:ntarget)))
  err_l2  = sqrt(sum(u(1:ntarget)**2) / real(ntarget, r64))
  write(*,'(/,A)')       '--- Result: DLP(1) = 0 outside ---'
  write(*,'(A,ES10.3)')  'max |u|  = ', err_max
  write(*,'(A,ES10.3)')  'rms |u|  = ', err_l2

  ! Cleanup
  deallocate(x, nx_s, w, xbd, tri_vert, tri2face, tri2cell, ptr)
  deallocate(sx_all, snx_all, sw_all, sigma_all)
  deallocate(tx, u, K_corr)

end program test_solid_angle
