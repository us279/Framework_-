! file: roofline_counter.f90
module roofline_counter
   use iso_fortran_env, only : int64
   use precision_mod, only: wp
   implicit none
   integer(int64) :: flop_count = 0_int64
   integer(int64) :: byte_count = 0_int64
contains
   subroutine add_flops(n)   ! count floating-point operations
      integer, intent(in) :: n
      flop_count = flop_count + int(n,int64)
   end subroutine
   subroutine add_bytes(n)   ! count bytes read/written
      integer, intent(in) :: n
      byte_count = byte_count + int(n,int64)
   end subroutine
end module roofline_counter