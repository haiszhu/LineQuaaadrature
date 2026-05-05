program debug_moments_brf
  use linequaaadrature_mod, only: r64
  use lq_kernel_mod, only: KERNEL_MOMENTS_MN, estimate_nearroot_lengths_r64, &
                           build_nearroot_nodes_r64, bary_row_r64, line_quad_BrF_r64
  use qotential_legacy_mod, only: moments_r64
  implicit none

  integer, parameter :: r128 = 16

  character(len=1024) :: data_dir
  integer(8) :: dims(5)
  integer(8) :: m, nbd, sbdnp, nquad, ncol
  integer(8) :: order, flag, len, lenl, lenr, n_up
  integer(8) :: i, j
  real(r64) :: root_info(4)
  real(r64) :: root_re, root_im
  logical :: root_ok

  real(r64), allocatable :: r0(:), sxbd(:,:), sxpbd(:,:), stangbd(:,:), sspbd(:)
  real(r64), allocatable :: tgl(:), wgl(:), dgl(:,:), w_bclag(:)
  real(r64), allocatable :: legmat(:,:), bclagmatlr(:,:), kdata(:)
  real(r64), allocatable :: t_up(:), w_up(:), row(:), r_up(:,:)
  real(r64), allocatable :: Nk_up(:,:), Mk_up(:,:), f_up(:,:), f_up2(:,:), Br(:,:)
  real(r64) :: max_abs_f, max_rel_f, max_abs_Br, max_abs_t, max_abs_Br_minus_I
  real(r64) :: Br_minus_I
  integer(8) :: imax, qmax

  data_dir = 'test/moments/debug_moments_brf'
  call get_command_argument(1, data_dir)
  data_dir = trim(strip_slash(data_dir))

  call read_int_vec(path(data_dir, 'dims.txt'), dims, 5_8)
  m = dims(1)
  nbd = dims(2)
  sbdnp = dims(3)
  nquad = dims(4)
  ncol = dims(5)
  order = ncol/2_8 - 1_8

  allocate(r0(3), sxbd(3,nbd), sxpbd(3,nbd), stangbd(3,nbd), sspbd(nbd))
  allocate(tgl(nquad), wgl(nquad), dgl(nquad,nquad), w_bclag(nquad))
  allocate(legmat(nquad,nquad), bclagmatlr(nquad,2), kdata(3))

  call read_vec(path(data_dir, 'r0.txt'), r0, 3_8)
  call read_mat(path(data_dir, 'sxbd.txt'), sxbd, 3_8, nbd)
  call read_mat(path(data_dir, 'sxpbd.txt'), sxpbd, 3_8, nbd)
  call read_mat(path(data_dir, 'stangbd.txt'), stangbd, 3_8, nbd)
  call read_vec(path(data_dir, 'sspbd.txt'), sspbd, nbd)
  call read_vec(path(data_dir, 'tgl.txt'), tgl, nquad)
  call read_vec(path(data_dir, 'wgl.txt'), wgl, nquad)
  call read_mat(path(data_dir, 'dgl.txt'), dgl, nquad, nquad)
  call read_vec(path(data_dir, 'w_bclag.txt'), w_bclag, nquad)
  call read_mat(path(data_dir, 'legmat.txt'), legmat, nquad, nquad)
  call read_mat(path(data_dir, 'bclagmatlr.txt'), bclagmatlr, nquad, 2_8)
  call read_vec(path(data_dir, 'kdata.txt'), kdata, 3_8)
  call read_vec(path(data_dir, 'root_kernel.txt'), root_info, 4_8)

  root_re = root_info(1)
  root_im = root_info(2)
  root_ok = (nint(root_info(3)) /= 0)
  flag = nint(kdata(2))

  if (.not. root_ok) error stop 'debug_moments_brf: saved root_ok is false'

  call estimate_nearroot_lengths_r64(root_re, abs(root_im), 12_8, len, lenl, lenr)
  n_up = nquad*(len - 1_8)

  allocate(t_up(n_up), w_up(n_up), row(nquad), r_up(n_up,3))
  allocate(Nk_up(n_up,order+1_8), Mk_up(n_up,order+1_8))
  allocate(f_up(n_up,ncol), f_up2(n_up,ncol), Br(nquad,n_up))

  call build_nearroot_nodes_r64(root_re, nquad, tgl, wgl, len, lenl, lenr, t_up, w_up)

  do i = 1, n_up
    call bary_row_r64(nquad, tgl, w_bclag, t_up(i), row)
    r_up(i,1) = sum(row*sxbd(1,:))
    r_up(i,2) = sum(row*sxbd(2,:))
    r_up(i,3) = sum(row*sxbd(3,:))
  end do

  Nk_up = 0.0_r64
  Mk_up = 0.0_r64
  call moments_r64(r0, n_up, r_up, order, flag, Nk_up, Mk_up)
  f_up(:,1:order+1_8) = Nk_up
  f_up(:,order+2_8:ncol) = Mk_up

  Br = 0.0_r64
  f_up2 = 0.0_r64
  call line_quad_BrF_r64(nquad, n_up, ncol, &
                         sxbd, sxpbd, &
                         tgl, wgl, w_bclag, &
                         t_up, w_up, &
                         r0, kdata, KERNEL_MOMENTS_MN, &
                         Br, f_up2)

  max_abs_f = maxval(abs(f_up - f_up2))
  max_rel_f = maxval(abs(f_up - f_up2) / max(max(abs(f_up), abs(f_up2)), 1.0e-300_r64))
  max_abs_Br = maxval(abs(Br))
  max_abs_t = maxval(abs(t_up - tgl))
  max_abs_Br_minus_I = 0.0_r64
  do i = 1, nquad
    do j = 1, n_up
      Br_minus_I = Br(i,j)
      if (n_up == nquad .and. i == j) Br_minus_I = Br_minus_I - 1.0_r64
      max_abs_Br_minus_I = max(max_abs_Br_minus_I, abs(Br_minus_I))
    end do
  end do
  imax = 1_8
  qmax = 1_8
  do i = 1, n_up
    do j = 1, ncol
      if (abs(f_up(i,j) - f_up2(i,j)) > abs(f_up(imax,qmax) - f_up2(imax,qmax))) then
        imax = i
        qmax = j
      end if
    end do
  end do

  write(*,'(a)') '=== debug_moments_brf ==='
  write(*,'(a,5(i0,1x))') 'dims m/nbd/sbdnp/nquad/ncol = ', m, nbd, sbdnp, nquad, ncol
  write(*,'(a,3(i0,1x))') 'len/lenl/lenr = ', len, lenl, lenr
  write(*,'(a,i0)') 'n_up = ', n_up
  write(*,'(a,i0)') 'order = ', order
  write(*,'(a,i0)') 'flag = ', flag
  write(*,'(a,es24.16,1x,es24.16)') 'root = ', root_re, root_im
  write(*,'(a,es24.16)') 'max abs f_up-f_up2 = ', max_abs_f
  write(*,'(a,es24.16)') 'max rel f_up-f_up2 = ', max_rel_f
  write(*,'(a,es24.16)') 'max abs Br = ', max_abs_Br
  write(*,'(a,es24.16)') 'max abs t_up-tgl = ', max_abs_t
  write(*,'(a,es24.16)') 'max abs Br-I = ', max_abs_Br_minus_I
  write(*,'(a,2(i0,1x))') 'max f diff at i/q = ', imax, qmax
  write(*,'(a,es24.16)') 'f_up  = ', f_up(imax,qmax)
  write(*,'(a,es24.16)') 'f_up2 = ', f_up2(imax,qmax)
  write(*,'(a,es24.16)') 'diff  = ', f_up(imax,qmax) - f_up2(imax,qmax)
  call print_moment_recurrence_terms(imax, order)

contains

  function strip_slash(s) result(out)
    character(len=*), intent(in) :: s
    character(len=1024) :: out
    integer :: n
    out = trim(s)
    n = len_trim(out)
    if (n > 0 .and. out(n:n) == '/') out(n:n) = ' '
  end function strip_slash

  function path(dir, file) result(out)
    character(len=*), intent(in) :: dir, file
    character(len=1024) :: out
    out = trim(dir)//'/'//trim(file)
  end function path

  subroutine read_int_vec(fname, a, n)
    character(len=*), intent(in) :: fname
    integer(8), intent(out) :: a(n)
    integer(8), intent(in) :: n
    integer :: u
    open(newunit=u, file=trim(fname), status='old', action='read')
    read(u,*) a
    close(u)
  end subroutine read_int_vec

  subroutine read_vec(fname, a, n)
    character(len=*), intent(in) :: fname
    real(r64), intent(out) :: a(n)
    integer(8), intent(in) :: n
    integer :: u
    open(newunit=u, file=trim(fname), status='old', action='read')
    read(u,*) a
    close(u)
  end subroutine read_vec

  subroutine read_mat(fname, a, nrow, ncol)
    character(len=*), intent(in) :: fname
    real(r64), intent(out) :: a(nrow,ncol)
    integer(8), intent(in) :: nrow, ncol
    integer(8) :: i
    integer :: u
    open(newunit=u, file=trim(fname), status='old', action='read')
    do i = 1, nrow
      read(u,*) a(i,1:ncol)
    end do
    close(u)
  end subroutine read_mat

  subroutine print_moment_recurrence_terms(ii, order_in)
    integer(8), intent(in) :: ii, order_in
    integer(8) :: kk
    real(r64) :: rt(3), r0dotr, rnorm, rnorm_inv, rnorm2_inv, r0norm
    real(r64) :: r0mr_vec(3), r0mr, r0mr_inv
    real(r64) :: term1, term2, term3, numer, cond
    real(r128) :: rtq(3), r0q(3), r0dotrq, rnormq, rnorm_invq, rnorm2_invq, r0normq
    real(r128) :: r0mr_vecq(3), r0mrq, r0mr_invq
    real(r128) :: Nq(order_in+1_8), Mq(order_in+1_8)
    real(r128) :: r0dotr_over_rnorm2q, r0norm2_over_rnorm2q, r0mr_over_rnorm2q
    real(r128) :: r0normplusrnormq, denominator1q, denominator2q, LMNcommonq, LMcommonq
    real(r128) :: term1q, term2q, term3q, numerq, condq

    rt = r_up(ii,:)
    r0dotr = sum(r0*rt)
    rnorm = sqrt(sum(rt**2))
    rnorm_inv = 1.0_r64/rnorm
    rnorm2_inv = rnorm_inv*rnorm_inv
    r0norm = sqrt(sum(r0**2))
    r0mr_vec = r0 - rt
    r0mr = sqrt(sum(r0mr_vec**2))
    r0mr_inv = 1.0_r64/r0mr

    r0q = real(r0, r128)
    rtq = real(rt, r128)
    r0dotrq = sum(r0q*rtq)
    rnormq = sqrt(sum(rtq**2))
    rnorm_invq = 1.0_r128/rnormq
    rnorm2_invq = rnorm_invq*rnorm_invq
    r0normq = sqrt(sum(r0q**2))
    r0mr_vecq = r0q - rtq
    r0mrq = sqrt(sum(r0mr_vecq**2))
    r0mr_invq = 1.0_r128/r0mrq
    r0dotr_over_rnorm2q = r0dotrq*rnorm2_invq
    r0norm2_over_rnorm2q = r0normq**2*rnorm2_invq
    r0mr_over_rnorm2q = r0mrq*rnorm2_invq
    r0normplusrnormq = r0normq + rnormq
    denominator1q = r0normq + sum(r0q*r0mr_vecq)*r0mr_invq
    denominator2q = rnormq + sum(rtq*r0mr_vecq)*r0mr_invq
    LMNcommonq = r0normplusrnormq/(denominator1q + denominator2q)
    Nq = 0.0_r128
    Mq = 0.0_r128
    Nq(1) = log((r0normplusrnormq+r0mrq)*LMNcommonq*r0mr_invq)*rnorm_invq
    Nq(2) = Nq(1)*r0dotr_over_rnorm2q + (r0mrq-r0normq)*rnorm2_invq
    do kk = 2_8, order_in
      Nq(kk+1_8) = real(2_8*kk-1_8,r128)/real(kk,r128)*r0dotr_over_rnorm2q*Nq(kk) - &
                   real(kk-1_8,r128)/real(kk,r128)*r0norm2_over_rnorm2q*Nq(kk-1_8) + &
                   1.0_r128/real(kk,r128)*r0mr_over_rnorm2q
    end do
    LMcommonq = 1.0_r128/(r0normq*rnormq+r0dotrq)
    Mq(1) = LMNcommonq*LMcommonq*(((r0normplusrnormq)*r0mr_invq+rnormq/r0normq)*r0mr_invq-1.0_r128/r0normq)
    Mq(2) = r0normq*Mq(1)/(r0normq+r0mrq)
    do kk = 2_8, order_in
      Mq(kk+1_8) = (r0dotrq*Mq(kk) + real(kk-1_8,r128)*Nq(kk-1_8) - r0mr_invq)*rnorm2_invq
    end do

    write(*,'(a)') ''
    write(*,'(a,i0)') 'recurrence diagnostics for node i = ', ii
    write(*,'(a,es24.16)') 'r0mr = ', r0mr
    write(*,'(a,es24.16)') 'rnorm2_inv = ', rnorm2_inv
    do kk = 2_8, order_in
      term1 = r0dotr*Mk_up(ii,kk)
      term2 = real(kk-1_8,r64)*Nk_up(ii,kk-1_8)
      term3 = -r0mr_inv
      numer = term1 + term2 + term3
      cond = (abs(term1) + abs(term2) + abs(term3))/max(abs(numer), 1.0e-300_r64)
      term1q = r0dotrq*Mq(kk)
      term2q = real(kk-1_8,r128)*Nq(kk-1_8)
      term3q = -r0mr_invq
      numerq = term1q + term2q + term3q
      condq = (abs(term1q) + abs(term2q) + abs(term3q))/max(abs(numerq), 1.0e-300_r128)
      write(*,'(a,i0)') 'M recurrence k = ', kk
      write(*,'(a,3(es24.16,1x),a,es24.16,a,es12.4)') &
        '  r64  terms/numer/cond = ', term1, term2, term3, ' | ', numer, ' | ', cond
      write(*,'(a,3(es24.16,1x),a,es24.16,a,es12.4)') &
        '  r128 terms/numer/cond = ', real(term1q,r64), real(term2q,r64), real(term3q,r64), &
        ' | ', real(numerq,r64), ' | ', real(condq,r64)
      write(*,'(a,es24.16,1x,es24.16,1x,es24.16)') &
        '  M_r64/M_r128/diff = ', Mk_up(ii,kk+1_8), real(Mq(kk+1_8),r64), &
        Mk_up(ii,kk+1_8)-real(Mq(kk+1_8),r64)
    end do
  end subroutine print_moment_recurrence_terms

end program debug_moments_brf
