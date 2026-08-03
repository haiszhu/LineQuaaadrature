
program test_solid_angle_fast
  use linequaaadrature_mod, only: r64
  use ellipsoid_mesh_mod,   only: create_ellipsoid_tri_mesh_r64
  use solidangle_mod,       only: evaluate_solid_angle_integral_r64, &
                                  evaluate_solid_angle_integral_fast_driver_r64
  use lap3d_mod,            only: lap3ddlpmat_r64, lap3ddlp_direct_r64
  use omp_lib,              only: omp_get_wtime
  implicit none

  integer(8), parameter :: order      = 14_8
  integer(8), parameter :: mp         = 8_8
  integer(8), parameter :: np         = 8_8
  real(r64),  parameter :: ratio      = 1.0_r64
  integer(8), parameter :: nquad_bdry = order + 8_8
  integer(8), parameter :: nq         = 0_8
  integer(8), parameter :: nvr        = order*(order+1_8)/2_8
  integer(8), parameter :: ntri       = 12_8*mp*np
  integer(8), parameter :: nbd        = 3_8*nquad_bdry
  integer(8), parameter :: nplotpts   = 20_8

  real(r64), parameter :: pi = 3.14159265358979323846264338327950288_r64

  real(r64),  allocatable :: x(:,:,:), nx_s(:,:,:), w(:,:)
  real(r64),  allocatable :: xbd(:,:,:), tri_vert(:,:,:)
  integer(8), allocatable :: tri2face(:), tri2cell(:,:), ptr(:)

  integer(8)              :: N_src
  real(r64),  allocatable :: sx_all(:,:), snx_all(:,:), sw_all(:), sigma_all(:)

  integer(8)              :: ntarget
  real(r64),  allocatable :: tx(:,:)
  real(r64),  allocatable :: u_base(:), u_fast(:), K_base(:), K_fast(:), u_dir(:)

  integer(8)              :: ntc
  integer(8), allocatable :: near_idx(:)
  real(r64),  allocatable :: tcj(:,:), Ib(:), Ia(:), Kmat(:,:)
  real(r64)               :: db, da

  real(r64) :: tx2(3,2), omega_base(2), omega_fast(2), ialpha2(2)
  real(r64) :: domain(6), px, py, pz, qpoint(3), qradii
  real(r64) :: worst_abs
  real(8)   :: t0, t1
  integer(8) :: k, i, ix, iy, iz

  logical, parameter :: use_nearroot = .true.

  write(*,'(/,A)') '=== test_solid_angle_fast ==='
  write(*,'(A,I0,A,I0,A,I0,A,F4.1)') &
      'order=', order, '  mp=np=', mp, '  ntri=', ntri, '  ratio=', ratio

  allocate(x(3,nvr,ntri), nx_s(3,nvr,ntri), w(nvr,ntri))
  allocate(xbd(3,nbd,ntri), tri_vert(3,3,ntri))
  allocate(tri2face(ntri), tri2cell(2,ntri), ptr(ntri+1))
  call create_ellipsoid_tri_mesh_r64(mp, np, order, nq, ratio, nquad_bdry, nvr, ntri, &
      x, nx_s, w, xbd, tri2face, tri2cell, tri_vert, ptr)

  write(*,'(/,A)') '--- Part A: solid angle sum ---'
  tx2(:,1) = [1.5_r64, 0.0_r64, 0.0_r64]
  tx2(:,2) = [1.1_r64, 0.0_r64, 0.0_r64]

  omega_base = 0.0_r64;  omega_fast = 0.0_r64
  do k = 1_8, ntri
    ialpha2 = 0.0_r64
    call evaluate_solid_angle_integral_r64(2_8, tx2, nvr, x(:,:,k), nx_s(:,:,k), &
        w(:,k), tri_vert(:,:,k), nbd, xbd(:,:,k), use_nearroot, ialpha2)
    omega_base = omega_base + ialpha2

    ialpha2 = 0.0_r64
    call evaluate_solid_angle_integral_fast_driver_r64(2_8, tx2, nvr, x(:,:,k), &
        nx_s(:,:,k), w(:,k), tri_vert(:,:,k), nbd, xbd(:,:,k), ialpha2)
    omega_fast = omega_fast + ialpha2
  end do

  write(*,'(A)')                   '  target        base        fast'
  write(*,'(A,ES11.3,1X,ES11.3)') '  (1.5,0,0)  ', abs(omega_base(1)), abs(omega_fast(1))
  write(*,'(A,ES11.3,1X,ES11.3)') '  (1.1,0,0)  ', abs(omega_base(2)), abs(omega_fast(2))

  write(*,'(/,A)') '--- Part B: DLP accuracy (direct sum + near correction) ---'

  N_src = ntri * nvr
  allocate(sx_all(3,N_src), snx_all(3,N_src), sw_all(N_src), sigma_all(N_src))
  sx_all    = reshape(x,    [3_8, N_src])
  snx_all   = reshape(nx_s, [3_8, N_src])
  sw_all    = reshape(w,    [N_src])
  sigma_all = -1.0_r64

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
  write(*,'(A,I0,A)') '  ntarget = ', ntarget, ' exterior points'

  allocate(u_dir(ntarget), u_base(ntarget), u_fast(ntarget))
  allocate(K_base(ntarget), K_fast(ntarget))
  K_base = 0.0_r64;  K_fast = 0.0_r64

  call lap3ddlp_direct_r64(ntarget, tx(:,1:ntarget), N_src, sx_all, snx_all, &
                           sw_all, sigma_all, u_dir)

  worst_abs = 0.0_r64
  t0 = omp_get_wtime()
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
      allocate(tcj(3,ntc), Ib(ntc), Ia(ntc), Kmat(ntc,nvr))
      do i = 1_8, ntc
        tcj(:,i) = tx(:, near_idx(i))
      end do

      Ib = 0.0_r64
      call evaluate_solid_angle_integral_r64(ntc, tcj, nvr, x(:,:,k), nx_s(:,:,k), &
          w(:,k), tri_vert(:,:,k), nbd, xbd(:,:,k), use_nearroot, Ib)

      Ia = 0.0_r64
      call evaluate_solid_angle_integral_fast_driver_r64(ntc, tcj, nvr, x(:,:,k), &
          nx_s(:,:,k), w(:,k), tri_vert(:,:,k), nbd, xbd(:,:,k), Ia)

      worst_abs = max(worst_abs, maxval(abs(Ia - Ib)))

      call lap3ddlpmat_r64(ntc, tcj, nvr, x(:,:,k), nx_s(:,:,k), w(:,k), Kmat)

      do i = 1_8, ntc
        db = - Ib(i) / (4.0_r64*pi) - sum(Kmat(i,:))
        da = - Ia(i) / (4.0_r64*pi) - sum(Kmat(i,:))
        K_base(near_idx(i)) = K_base(near_idx(i)) + db
        K_fast(near_idx(i)) = K_fast(near_idx(i)) + da
      end do

      deallocate(tcj, Ib, Ia, Kmat)
    end if

    deallocate(near_idx)
  end do
  t1 = omp_get_wtime()
  write(*,'(A,F7.2,A)') '  near-field correction done in ', t1-t0, ' s'

  u_base = u_dir - K_base
  u_fast = u_dir - K_fast

  write(*,'(A,ES10.3)') '  max |fast - base| on IalphaAsvestas = ', worst_abs
  write(*,'(/,A)')      '--- Result: DLP(1) = 0 outside ---'
  write(*,'(A,ES10.3,A,ES10.3)') '  base   max |u| = ', maxval(abs(u_base)), &
       '   rms |u| = ', sqrt(sum(u_base**2)/real(ntarget, r64))
  write(*,'(A,ES10.3,A,ES10.3)') '  fast   max |u| = ', maxval(abs(u_fast)), &
       '   rms |u| = ', sqrt(sum(u_fast**2)/real(ntarget, r64))
  write(*,'(A,ES10.3)') '  max |u_fast - u_base| = ', maxval(abs(u_fast - u_base))

  deallocate(x, nx_s, w, xbd, tri_vert, tri2face, tri2cell, ptr)
  deallocate(sx_all, snx_all, sw_all, sigma_all)
  deallocate(tx, u_dir, u_base, u_fast, K_base, K_fast)

end program test_solid_angle_fast
