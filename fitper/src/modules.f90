module common

  real, dimension(:,:), allocatable :: obs_per1, obs_per2, obs_per_dum, obs_per_seq
  real, dimension(:,:), allocatable :: calc_per1, calc_per2, calc_per
  real, dimension(:,:), allocatable :: calc_per1_copy, calc_per2_copy
  real, dimension(:,:), allocatable :: calc_per_copy, calc_per_ell
  real, dimension(:,:), allocatable :: fit_per, fit_per_seq,  fit_per_ell
  real, dimension(15) :: params
  real :: rdum, sigma, xisqr, rnumod, omc
  integer, dimension(4) :: pegseq_ell1, seqsize_ell1, pegseq_ell2, seqsize_ell2
  integer, dimension(2) :: npreseq
  integer :: pegseq, pegseqobs, seqsize, pegseq2, seqsize2, npreseq_ell
  integer :: ncalc, ncalc1, ncalc2, ncalctot, ncalc_ell,numfailed, peg0, peg
  integer :: nper_tot, nper1, nper2, nper_ell
  logical :: success, success1, ifitrest, ifit2ndseq, secondellseq
  logical, parameter :: verbose = .false.

  real, parameter :: huge = 1.e23
  integer, parameter :: npar = 4, imax = 10

end module common

!Column 1 = index number
!Column 2 = observed periods
!Column 3 = ell of calculated period
!Column 4 = k of calculated period
!Column 5 = calculated period
!Column 6 = calculated - observed
!Column 7 = ell specified in obsperiods
