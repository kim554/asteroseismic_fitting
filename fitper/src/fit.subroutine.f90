!This routine fits rich pulsation spectra, where one has consecutive ell=1
!and/or ell=2 sequences.

subroutine calcfit

  use common
  use match_periods_subroutine

  implicit none

  real, dimension(:), allocatable :: matched_per, weights
  real :: sumofweights, denominator
  integer :: i, j
  logical :: stuck

90  format(F12.6,2F4.0)
100 format(F5.0, F12.2, 2F5.0, F12.2, F8.2, F5.0)
!110 format(F9.1,4F7.1,F6.1,2F6.2,3F7.2,5(F10.3,' (',I1,') '),F16.6)
                                       !nper_tot
110 format(15F9.1,9(F10.3,' (',I1,') '),F16.6)
                 !nper tot
  allocate(matched_per(nper_tot),weights(nper_tot))

  success = .true.
  stuck = .false.
  secondellseq = .false.

  sumofweights = 0.

  if (verbose) then
     print *, ''
     print *, 'obs_per1 '
     do i=1,nper1
        print 90, obs_per1(i,:)
     end do
     print *, ""
     print *, "obs_per2"
     do i=1,nper2
        print 90, obs_per2(i,:)
     end do
  end if

do i = 1, nper1
     fit_per(i,2) = obs_per1(i,1)
     fit_per(i,7) = obs_per1(i,2)
     weights(i) = 1./obs_per1(i,3)**2
     sumofweights = sumofweights + weights(i)
  end do
  do i = 1, nper2
     j = i + nper1
     fit_per(j,2) = obs_per2(i,1)
     fit_per(j,7) = obs_per2(i,2)
     weights(j) = 1./obs_per2(i,3)**2
     sumofweights = sumofweights + weights(i)
  end do

  if (verbose) then
     print *, ''
     print *, "fit_per, with just observed periods filled in"
     do i=1,nper_tot
        print 100, fit_per(i,:)
     end do
  end if

!First take care of the ell=1 observed periods list

  if (verbose) then
     print *, ''
     print *, 'Matching the ell=1 modes'
  end if

  call match_periods(pegseq_ell1,seqsize_ell1,1)

  do i=1,nper1
     fit_per(i,:) = fit_per_ell(i,:)
  end do

  if (verbose) then
     print *, ''
     print *, 'fit_per after matching the ell=1 periods'
     do i=1,nper_tot
        print 100, fit_per(i,:)
     end do
     print *, ''
  end if

  deallocate(fit_per_ell)

!Now take care of the ell=2 sequence

  secondellseq = .false.

  if ( (seqsize_ell2(1) .ne. 0) .or. (seqsize_ell2(2) .ne. 0) &
        .or. (seqsize_ell2(3) .ne. 0) .or. (seqsize_ell2(4) .ne. 0) ) then

     if (verbose) then
        print *, ''
        print *, 'Matching the ell=2 modes'
     end if

     call match_periods(pegseq_ell2,seqsize_ell2,2)

     do i=1,nper2
        j = i + nper1
        fit_per(j,:) = fit_per_ell(i,:)
     end do

     if (verbose) then
        print *, ''
        print *, 'fit_per after matching the ell=2 periods'
        do i=1,nper_tot
           print 100, fit_per(i,:)
        end do
        print *, ''
     end if

     deallocate(fit_per_ell)

  end if

  if (verbose) then
     print *, ''
     print *, 'final matching'
  end if

120 format(I5, 2F12.5, F8.2)

  sigma = 0.
  do i=1,nper_tot
     matched_per(i) = fit_per(i,5)
     omc = matched_per(i) - fit_per(i,2)
     if (verbose) print 120, i, fit_per(i,2), matched_per(i), omc
     sigma = sigma + weights(i)*omc**2
  end do

  denominator = real(nper_tot-1)*sumofweights/real(nper_tot)
  sigma = sqrt(sigma/denominator)
  xisqr = real(nper_tot)*log(sigma**2)+real(npar)*log(real(nper_tot))

!  if (params(4) .ge. params(3)+200.) then
      write(3,110) params(:), matched_per(1), int(fit_per(1,3)),&
       matched_per(2), int(fit_per(2,3)), &
       matched_per(3), int(fit_per(3,3)), &
       matched_per(4), int(fit_per(4,3)), &
       matched_per(5), int(fit_per(5,3)), &
       matched_per(6), int(fit_per(6,3)), &
       matched_per(7), int(fit_per(7,3)), &
       matched_per(8), int(fit_per(8,3)), &
       matched_per(9), int(fit_per(9,3)), &
!       matched_per(10), int(fit_per(10,3)), &
!       matched_per(11), int(fit_per(11,3)), &
!       matched_per(12), int(fit_per(12,3)), &
!       matched_per(13), int(fit_per(13,3)), &
!       matched_per(14), int(fit_per(14,3)), &
!       matched_per(15), int(fit_per(15,3)), &
!       matched_per(16), int(fit_per(16,3)), &
!       matched_per(17), int(fit_per(17,3)), &
!       matched_per(18), int(fit_per(18,3)), &
!       matched_per(19), int(fit_per(13,3)), &
!       matched_per(20), int(fit_per(20,3)), &
!matched_per(21), int(fit_per(21,3)), &
!       matched_per(22), int(fit_per(22,3)), matched_per(23), int(fit_per(23,3)), &
!       matched_per(24), int(fit_per(24,3)), matched_per(25), int(fit_per(25,3)), &
!       matched_per(26), int(fit_per(26,3)),  &
       sigma
!   end if

  deallocate(matched_per,weights)

end subroutine calcfit

