! linequaaadrature_mex.f90
! Standalone (non-module) MEX-facing wrappers.
!
! Why needed: gfortran module procedures export as
!   __linequaaadrature_mod_MOD_foo  (module-mangled)
! mwrap expects the plain Fortran symbol _foo_ (trailing underscore, macOS).
!
! Each wrapper imports the module procedure under a local alias and
! delegates.  The fun argument is now type(c_funptr) throughout the
! call chain, so it appears as a regular derived-type dummy in the
! .mod explicit interface — no shim or dispatch module needed.
!
! Calling convention from MATLAB:
!   fptr = kernel_lookup_mex('asvestas_kernel_r64')   % int64 handle
!   [funvals] = line_kernel_eval_mex(..., fptr, kdata, funvals)

! ------------------------------------------------------------------
! gauss
! ------------------------------------------------------------------
subroutine gauss_r64(n, tgl, wgl, Dgl)
  use linequaaadrature_mod, only: gauss_r64_mod => gauss_r64
  implicit none
  integer(8), intent(in)    :: n
  real(8),    intent(inout) :: tgl(n), wgl(n), Dgl(n,n)
  call gauss_r64_mod(n, tgl, wgl, Dgl)
end subroutine gauss_r64

! ------------------------------------------------------------------
! kernel_lookup
! Looks up a Fortran standalone subroutine by name via dlsym.
! Appends '_' to sym_name before the null terminator.
! Returns the function pointer as integer(8); 0 = not found.
! ------------------------------------------------------------------
subroutine kernel_lookup(sym_name, sym_len, fptr_int)
  use iso_c_binding, only: c_char, c_funptr, c_funloc, c_null_char, c_null_funptr, c_associated
  implicit none
  integer(8),        intent(in)  :: sym_len
  character(c_char), intent(in)  :: sym_name(sym_len)
  integer(8),        intent(out) :: fptr_int

  interface
    subroutine asvestas_kernel_r64(r_s, tau_s, r0j, kdata3, val) bind(C)
      use iso_c_binding, only: c_double
      implicit none
      real(c_double), intent(in)    :: r_s(3), tau_s(3), r0j(3), kdata3(3)
      real(c_double), intent(inout) :: val
    end subroutine asvestas_kernel_r64
  end interface

  integer(8) :: i, nlen
  character(len=128) :: name
  type(c_funptr) :: fptr

  name = ''
  nlen = 0_8
  do i = 1_8, sym_len
    if (sym_name(i) == c_null_char) exit
    if (nlen < len(name)) then
      nlen = nlen + 1_8
      name(nlen:nlen) = sym_name(i)
    end if
  end do

  fptr = c_null_funptr
  select case (trim(name))
  case ('asvestas_kernel_r64')
    fptr = c_funloc(asvestas_kernel_r64)
  case default
    fptr = c_null_funptr
  end select

  if (.not. c_associated(fptr)) then
    fptr_int = 0_8
  else
    fptr_int = transfer(fptr, 0_8)
  end if

end subroutine kernel_lookup

! asvestas_kernel_r64 is now exported directly from solidangle_mod with bind(C).
! No standalone wrapper needed.

! ------------------------------------------------------------------
! evaluate_solid_angle_integral  — standalone wrapper for mwrap
! Delegates to solidangle_mod::evaluate_solid_angle_integral_r64.
! ------------------------------------------------------------------
subroutine evaluate_solid_angle_integral(m, tx, n, sx, snx, sw, r_vert, nbd, sxbd_in, IalphaAsvestas)
  use solidangle_mod, only: eval_sa => evaluate_solid_angle_integral_r64
  implicit none
  integer(8), intent(in)    :: m, n, nbd
  real(8),    intent(in)    :: tx(3,m), sx(3,n), snx(3,n), sw(n)
  real(8),    intent(in)    :: r_vert(3,3), sxbd_in(3,nbd)
  real(8),    intent(inout) :: IalphaAsvestas(m)
  call eval_sa(m, tx, n, sx, snx, sw, r_vert, nbd, sxbd_in, IalphaAsvestas)
end subroutine evaluate_solid_angle_integral

! ------------------------------------------------------------------
! create_ellipsoid_tri_mesh  — standalone wrapper for mwrap
! Delegates to ellipsoid_mesh_mod::create_ellipsoid_tri_mesh_r64.
!
! mwrap only supports 1D/2D arrays, so the 3D arrays x, nx, xbd,
! tri_vert are passed flattened (first dims merged):
!   x_flat(3*nvr, ntri)           ← x(3,nvr,ntri)
!   nx_flat(3*nvr, ntri)          ← nx(3,nvr,ntri)
!   xbd_flat(3*3*nquad_bdry,ntri) ← xbd(3,3*nquad_bdry,ntri)
!   tv_flat(9, ntri)              ← tri_vert(3,3,ntri)
!
! tri2face, tri2cell, ptr are integer(8) in the module but MATLAB
! passes/receives double; this wrapper converts in/out.
!
! c_f_pointer is used to remap flat 2D views to 3D without copying.
! ------------------------------------------------------------------
subroutine create_ellipsoid_tri_mesh(mp, np, p, nq, ratio, nquad_bdry, nvr, ntri, &
    x_flat, nx_flat, w, xbd_flat, tri2face_d, tri2cell_d, tv_flat, ptr_d)
  use ellipsoid_mesh_mod, only: create_r64 => create_ellipsoid_tri_mesh_r64
  use iso_c_binding,      only: c_f_pointer, c_loc
  implicit none
  integer(8), intent(in)           :: mp, np, p, nq, nquad_bdry, nvr, ntri
  real(8),    intent(in)           :: ratio
  real(8),    intent(inout), target :: x_flat(3*nvr, ntri)
  real(8),    intent(inout), target :: nx_flat(3*nvr, ntri)
  real(8),    intent(inout)        :: w(nvr, ntri)
  real(8),    intent(inout), target :: xbd_flat(3*3*nquad_bdry, ntri)
  real(8),    intent(inout)        :: tri2face_d(ntri)
  real(8),    intent(inout)        :: tri2cell_d(2, ntri)
  real(8),    intent(inout), target :: tv_flat(9, ntri)
  real(8),    intent(inout)        :: ptr_d(ntri+1)

  real(8), pointer :: x3(:,:,:), nx3(:,:,:), xbd3(:,:,:), tv3(:,:,:)
  integer(8), allocatable :: tri2face(:), tri2cell(:,:), ptr_i(:)
  integer(8) :: k

  ! Remap flat 2D -> 3D without copying
  call c_f_pointer(c_loc(x_flat),   x3,   [3_8, nvr, ntri])
  call c_f_pointer(c_loc(nx_flat),  nx3,  [3_8, nvr, ntri])
  call c_f_pointer(c_loc(xbd_flat), xbd3, [3_8, 3_8*nquad_bdry, ntri])
  call c_f_pointer(c_loc(tv_flat),  tv3,  [3_8, 3_8, ntri])

  allocate(tri2face(ntri), tri2cell(2,ntri), ptr_i(ntri+1))

  call create_r64(mp, np, p, nq, ratio, nquad_bdry, nvr, ntri, &
      x3, nx3, w, xbd3, tri2face, tri2cell, tv3, ptr_i)

  do k = 1_8, ntri
    tri2face_d(k)    = real(tri2face(k),   8)
    tri2cell_d(1,k)  = real(tri2cell(1,k), 8)
    tri2cell_d(2,k)  = real(tri2cell(2,k), 8)
  end do
  do k = 1_8, ntri + 1_8
    ptr_d(k) = real(ptr_i(k), 8)
  end do

  deallocate(tri2face, tri2cell, ptr_i)

end subroutine create_ellipsoid_tri_mesh
