! test_solid_angle_r128.f90
! Quad-precision (real(16)) port of test_solid_angle.f90.
! Uses areal quadrature (nq > 0) since VR nodes have no r128 tabulated data.
!
! Two checks:
!   Part A — solid angle sum for 2 exterior targets.
!   Part B — full DLP accuracy with near-field solid-angle correction.
!
! Compile via Makefile:  make test_r128
! Run:                   ./build/test_solid_angle_r128

program test_solid_angle_r128
  use linequaaadrature_mod, only: r128
  use ellipsoid_mesh_mod,   only: create_ellipsoid_tri_mesh_r128
  use solidangle_mod,       only: evaluate_solid_angle_integral_r128
  use lap3d_mod,            only: lap3ddlp_direct_r128, lap3ddlpmat_r128
  use omp_lib,              only: omp_get_wtime
  implicit none

  ! ---- problem parameters ----
  integer(8), parameter :: order      = 16_8
  integer(8), parameter :: mp         = 8_8
  integer(8), parameter :: np         = 8_8
  real(r128), parameter :: ratio      = 1.0_r128    ! unit sphere
  integer(8), parameter :: nquad_bdry = order + 16_8
  integer(8), parameter :: nq         = 16_8        ! areal rule: nvr = 3*nq^2
  integer(8), parameter :: nvr        = 3_8*nq*nq
  integer(8), parameter :: ntri       = 12_8*mp*np
  integer(8), parameter :: nplotpts   = 10_8

  real(r128), parameter :: pi = &
    3.14159265358979323846264338327950288419716939937510_r128

  ! ---- mesh ----
  real(r128),  allocatable :: x(:,:,:), nx_s(:,:,:), w(:,:)
  real(r128),  allocatable :: xbd(:,:,:), tri_vert(:,:,:)
  integer(8),  allocatable :: tri2face(:), tri2cell(:,:), ptr(:)

  ! ---- all sources flattened ----
  integer(8)              :: N_src
  real(r128),  allocatable :: sx_all(:,:), snx_all(:,:), sw_all(:), sigma_all(:)

  ! ---- target points ----
  integer(8)              :: ntarget
  real(r128),  allocatable :: tx(:,:)

  ! ---- result vectors ----
  real(r128),  allocatable :: u(:), K_corr(:)

  ! ---- per-triangle work ----
  integer(8)              :: ntc
  integer(8),  allocatable :: near_idx(:)
  real(r128),  allocatable :: tcj(:,:), IalphaAsv(:), Kmat(:,:)
  real(r128)               :: delta_k

  ! ---- Part A sanity check ----
  real(r128) :: tx2(3,4), omega_sum(4), ialpha2(4)

  ! ---- misc ----
  real(r128) :: domain(6), px, py, pz, qpoint(3), qradii
  real(r128) :: err_max, err_l2
  real(8)    :: t0, t1
  integer(8) :: k, i, ix, iy, iz

  logical :: use_nearroot
  character(len=32) :: mode
  integer :: nargs

  ! ==============================================================
  use_nearroot = .false.

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
      write(*,'(A)')   'Use: ./build/test_solid_angle_r128 adaptive'
      write(*,'(A)')   '  or ./build/test_solid_angle_r128 nearroot'
      stop 1
    end select
  end if

  write(*,'(/,A)') '=== test_solid_angle_r128 (quad precision) ==='
  write(*,'(A,I0,A,I0,A,I0,A,I0,A,F4.1)') &
      'order=', order, '  nq=', nq, '  nvr=', nvr, &
      '  ntri=', ntri, '  ratio=', real(ratio, 8)
  write(*,'(A,L1)') 'use_nearroot = ', use_nearroot

  ! ==============================================================
  ! Build ellipsoid mesh (r128, areal quadrature)
  ! ==============================================================
  allocate(x(3,nvr,ntri), nx_s(3,nvr,ntri), w(nvr,ntri))
  allocate(xbd(3,3*nquad_bdry,ntri), tri_vert(3,3,ntri))
  allocate(tri2face(ntri), tri2cell(2,ntri), ptr(ntri+1))

  t0 = omp_get_wtime()
  call create_ellipsoid_tri_mesh_r128(mp, np, order, nq, ratio, nquad_bdry, nvr, ntri, &
      x, nx_s, w, xbd, tri2face, tri2cell, tri_vert, ptr)
  t1 = omp_get_wtime()
  write(*,'(A,F7.2,A)') 'Mesh built in ', t1-t0, ' s'

  ! ==============================================================
  ! Part A: solid angle sum
  ! ==============================================================
  write(*,'(/,A)') '--- Part A: solid angle sum ---'

  tx2(:,1) = [1.5_r128, 0.0_r128, 0.0_r128]   ! exterior, far
  tx2(:,2) = [1.001_r128, 0.0_r128, 0.0_r128]   ! exterior, near surface
  tx2(:,3) = [1.000001_r128, 0.0_r128, 0.0_r128]
  tx2(:,4) = [1.000000001_r128, 0.0_r128, 0.0_r128] 

  omega_sum = 0.0_r128
  do k = 1_8, ntri
    ialpha2 = 0.0_r128
    call evaluate_solid_angle_integral_r128(4_8, tx2, nvr,               &
        x(:,:,k), nx_s(:,:,k), w(:,k),                                   &
        tri_vert(:,:,k), 3_8*nquad_bdry, xbd(:,:,k), use_nearroot, ialpha2)
    omega_sum = omega_sum + ialpha2
  end do

  write(*,'(A,ES14.6)') 'far exterior  (1.5,0,0): |omega_sum| = ', &
      real(abs(omega_sum(1)), 8)
  write(*,'(A,ES14.6)') 'near exterior (1.001,0,0): |omega_sum| = ', &
      real(abs(omega_sum(2)), 8)
  write(*,'(A,ES14.6)') 'near exterior (1.000001,0,0): |omega_sum| = ', &
      real(abs(omega_sum(3)), 8)
  write(*,'(A,ES14.6)') 'near exterior (1.000000001,0,0): |omega_sum| = ', &
      real(abs(omega_sum(4)), 8)

  ! ==============================================================
  ! Part B: full DLP accuracy
  ! ==============================================================
  write(*,'(/,A)') '--- Part B: DLP accuracy (direct sum + near correction) ---'

  N_src = ntri * nvr
  allocate(sx_all(3,N_src), snx_all(3,N_src), sw_all(N_src), sigma_all(N_src))
  sx_all    = reshape(x,    [3_8, N_src])
  snx_all   = reshape(nx_s, [3_8, N_src])
  sw_all    = reshape(w,    [N_src])
  sigma_all = -1.0_r128
  write(*,'(A,I0)') 'N_src = ', N_src

  ! Generate exterior target grid
  domain = [-0.5_r128, 2.0_r128, -0.5_r128, 2.0_r128, -1.25_r128, 1.25_r128]
  allocate(tx(3, nplotpts**3))
  ntarget = 0_8
  do iz = 1_8, nplotpts
    do iy = 1_8, nplotpts
      do ix = 1_8, nplotpts
        px = domain(1) + (domain(2)-domain(1)) * (ix-1) / real(nplotpts-1, r128)
        py = domain(3) + (domain(4)-domain(3)) * (iy-1) / real(nplotpts-1, r128)
        pz = domain(5) + (domain(6)-domain(5)) * (iz-1) / real(nplotpts-1, r128)
        if (px**2 + (py*ratio)**2 + (pz*ratio)**2 > 1.0_r128 + 1.0e-8_r128) then
          ntarget = ntarget + 1_8
          tx(1,ntarget) = px
          tx(2,ntarget) = py
          tx(3,ntarget) = pz
        end if
      end do
    end do
  end do
  write(*,'(A,I0,A)') 'ntarget = ', ntarget, ' exterior points'

  allocate(u(ntarget), K_corr(ntarget))
  K_corr = 0.0_r128

  write(*,'(A)') 'Direct DLP sum ...'
  t0 = omp_get_wtime()
  call lap3ddlp_direct_r128(ntarget, tx(:,1:ntarget), N_src, sx_all, snx_all, sw_all, sigma_all, u)
  t1 = omp_get_wtime()
  write(*,'(A,F7.2,A)') '  done in ', t1-t0, ' s'

  write(*,'(A)') 'Near-field correction ...'
  t0 = omp_get_wtime()

  !$omp parallel do schedule(dynamic) &
  !$omp   default(shared) &
  !$omp   private(ntc, near_idx, tcj, IalphaAsv, Kmat, qpoint, qradii, i, delta_k)
  do k = 1_8, ntri
    qpoint = sum(x(:,:,k), dim=2) / real(nvr, r128)
    qradii = 7.0_r128 * sqrt(sum(w(:,k)))

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

    if (ntc == 0_8) then
      deallocate(near_idx)
      cycle
    end if

    allocate(tcj(3,ntc), IalphaAsv(ntc), Kmat(ntc,nvr))

    do i = 1_8, ntc
      tcj(:,i) = tx(:, near_idx(i))
    end do

    IalphaAsv = 0.0_r128
    call evaluate_solid_angle_integral_r128(ntc, tcj, nvr,             &
        x(:,:,k), nx_s(:,:,k), w(:,k),                                 &
        tri_vert(:,:,k), 3_8*nquad_bdry, xbd(:,:,k), use_nearroot, IalphaAsv)

    call lap3ddlpmat_r128(ntc, tcj, nvr, x(:,:,k), nx_s(:,:,k), w(:,k), Kmat)

    do i = 1_8, ntc
      delta_k = - IalphaAsv(i) / (4.0_r128*pi) - sum(Kmat(i,:))
      !$omp atomic
      K_corr(near_idx(i)) = K_corr(near_idx(i)) + delta_k
    end do

    deallocate(tcj, IalphaAsv, Kmat)
    deallocate(near_idx)
  end do
  !$omp end parallel do

  t1 = omp_get_wtime()
  write(*,'(A,F7.2,A)') '  done in ', t1-t0, ' s'

  do i = 1_8, ntarget
    u(i) = u(i) - K_corr(i)
  end do

  err_max = maxval(abs(u(1:ntarget)))
  err_l2  = sqrt(sum(u(1:ntarget)**2) / real(ntarget, r128))
  write(*,'(/,A)')      '--- Result: DLP(1) = 0 outside ---'
  write(*,'(A,ES14.6)') 'max |u|  = ', real(err_max, 8)
  write(*,'(A,ES14.6)') 'rms |u|  = ', real(err_l2, 8)

  block
    integer(8) :: imax_arr(1), imax
    imax_arr = maxloc(abs(u(1:ntarget)))
    imax = imax_arr(1)
    write(*,'(A,I0,A,3ES23.15)') 'argmax: i=', imax, &
        '  tx=', real(tx(1,imax), 8), real(tx(2,imax), 8), real(tx(3,imax), 8)
    write(*,'(A,ES23.15)') 'u(imax)        = ', real(u(imax), 8)
    write(*,'(A,ES23.15)') 'K_corr(imax)   = ', real(K_corr(imax), 8)
  end block

  ! Cleanup
  deallocate(x, nx_s, w, xbd, tri_vert, tri2face, tri2cell, ptr)
  deallocate(sx_all, snx_all, sw_all, sigma_all)
  deallocate(tx, u, K_corr)

end program test_solid_angle_r128