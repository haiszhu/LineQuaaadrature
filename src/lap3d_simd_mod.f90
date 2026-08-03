
module lap3d_simd_mod
  use linequaaadrature_mod, only: r64
  implicit none
  private
  public :: csimd128lap3dsdlpmat_r64, csimd256lap3dsdlpmat_r64, &
            csimd512lap3dsdlpmat_r64

contains

  subroutine csimd128lap3dsdlpmat_r64(m, r0, n, r, rn, w, As, Ad)
    integer(8), intent(in)    :: m, n
    real(r64),  intent(in)    :: r0(*), r(*), rn(*), w(*)
    real(r64),  intent(inout) :: As(*), Ad(*)

    call lq_csimd128lap3dsdlpmat_c(m, r0, n, r, rn, w, As, Ad)

  end subroutine csimd128lap3dsdlpmat_r64

  subroutine csimd256lap3dsdlpmat_r64(m, r0, n, r, rn, w, As, Ad)
    integer(8), intent(in)    :: m, n
    real(r64),  intent(in)    :: r0(*), r(*), rn(*), w(*)
    real(r64),  intent(inout) :: As(*), Ad(*)

    call lq_csimd256lap3dsdlpmat_c(m, r0, n, r, rn, w, As, Ad)

  end subroutine csimd256lap3dsdlpmat_r64

  subroutine csimd512lap3dsdlpmat_r64(m, r0, n, r, rn, w, As, Ad)
    integer(8), intent(in)    :: m, n
    real(r64),  intent(in)    :: r0(*), r(*), rn(*), w(*)
    real(r64),  intent(inout) :: As(*), Ad(*)

    call lq_csimd512lap3dsdlpmat_c(m, r0, n, r, rn, w, As, Ad)

  end subroutine csimd512lap3dsdlpmat_r64

end module lap3d_simd_mod
