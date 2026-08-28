!> Precision kinds. The bottom of the dependency chain: uses nothing.
!> Every ESM has a module exactly like this, and every ESM has had a bug
!> caused by someone writing 1.0 instead of 1.0_wp somewhere in it.
module kinds
  use, intrinsic :: iso_fortran_env, only: real64, real32, int32, int64
  implicit none
  public
  integer, parameter :: sp = real32
  integer, parameter :: dp = real64
  integer, parameter :: wp = dp        !< working precision, switchable
end module kinds
