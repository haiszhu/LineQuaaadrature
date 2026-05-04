! ------------------------------------------------------------------
! lqa_root_refine
! Thin MEX-facing wrapper around line_quad_adaptive_mod::lqa_root_refine_r64.
! Complex input/output are split into real/imag parts.
! ifconv is returned as double to keep MATLAB-side handling simple.
! ------------------------------------------------------------------
subroutine lqa_root_refine_mex(xhat, yhat, zhat, n_expa, tx, ty, tz, &
                           tinit, troot, ifconv)
  use lq_adaptive_mod, only: line_quad_root_refine_r64
  implicit none
  integer(8), intent(in)    :: n_expa
  real(8),    intent(in)    :: xhat(n_expa), yhat(n_expa), zhat(n_expa)
  real(8),    intent(in)    :: tx, ty, tz
  complex(8), intent(in)    :: tinit
  complex(8), intent(inout) :: troot
  integer(8), intent(inout) :: ifconv

  call line_quad_root_refine_r64(xhat, yhat, zhat, n_expa, tx, ty, tz, &
                                 tinit, troot, ifconv)
end subroutine lqa_root_refine_mex

subroutine lqa_root_refine_r128_mex(xhat, yhat, zhat, n_expa, tx, ty, tz, &
                           tinit, troot, ifconv)
  use lq_adaptive_mod, only: lqrr_r128 => line_quad_root_refine_r128
  use iso_c_binding, only: c_char, c_int, c_null_char
  implicit none
  interface
    function hdf5_write_string_pair(file, name1, value1, name2, value2) bind(C)
      import c_char, c_int
      character(kind=c_char), intent(in) :: file(*), name1(*), value1(*), name2(*), value2(*)
      integer(c_int) :: hdf5_write_string_pair
    end function hdf5_write_string_pair
  end interface
  integer(8), intent(in)    :: n_expa
  real(8),    intent(in)    :: xhat(n_expa), yhat(n_expa), zhat(n_expa)
  real(8),    intent(in)    :: tx, ty, tz
  complex(8), intent(in)    :: tinit
  complex(8), intent(inout) :: troot
  integer(8), intent(inout) :: ifconv

  integer(c_int) :: h5_ok
  real(16) :: xhat_128(n_expa), yhat_128(n_expa), zhat_128(n_expa)
  real(16) :: tx_128, ty_128, tz_128
  complex(16) :: tinit_128, troot_128
  character(len=64) :: troot_re_str, troot_im_str

  xhat_128 = real(xhat, 16)
  yhat_128 = real(yhat, 16)
  zhat_128 = real(zhat, 16)
  tx_128 = real(tx, 16)
  ty_128 = real(ty, 16)
  tz_128 = real(tz, 16)
  tinit_128 = cmplx(real(tinit, 16), real(aimag(tinit), 16), kind=16)
  troot_128 = cmplx(real(troot, 16), real(aimag(troot), 16), kind=16)

  call lqrr_r128(xhat_128, yhat_128, zhat_128, n_expa, tx_128, ty_128, tz_128, &
                                 tinit_128, troot_128, ifconv)
  write (troot_re_str, '(ES46.36)') real(troot_128, 16)
  write (troot_im_str, '(ES46.36)') aimag(troot_128)
  h5_ok = hdf5_write_string_pair('lqrr_r128.h5'//c_null_char, &
                                 '/troot_re'//c_null_char, trim(adjustl(troot_re_str))//c_null_char, &
                                 '/troot_im'//c_null_char, trim(adjustl(troot_im_str))//c_null_char)
  troot = cmplx(real(troot_128, 8), real(aimag(troot_128), 8), kind=8)

end subroutine lqa_root_refine_r128_mex
