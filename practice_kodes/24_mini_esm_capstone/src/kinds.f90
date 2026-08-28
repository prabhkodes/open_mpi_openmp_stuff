!> Precision kinds. Bottom of the dependency chain.
module kinds
  use, intrinsic :: iso_fortran_env, only: real64, real32
  implicit none
  public
  integer, parameter :: sp = real32
  integer, parameter :: dp = real64
  integer, parameter :: wp = dp        !< working precision
end module kinds
