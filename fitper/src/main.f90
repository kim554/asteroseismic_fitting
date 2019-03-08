!To use this fitting program, make sure that the obsperiod file has sorted
!lists of periods, separated into different ell. Put the unknown ell periods
!in the ell = 1 list. npreseq, pegseq(:) and seqsize(:) are specified in the input file.

! XXX  0  npreseq = 2
! XXX  0
! XX1  1  pegseq(1) = 3 (third period in list), seqsize(1) = 4
! XX2  1
! XX3  1
! XX4  1
! XXX  0  pegseq(2) = 7, seqsize(2) = 3
! XXX  0
! XXX  0
! XX1  1  pegseq(3) = 10, seqsize(3) = 2
! XX2  1
! XXX  0  pegseq(4) = 12, seqsize(4) = 1

!If part i is missing, set seqsize(i) to zero. pegseq(i) will be ignored by the program so set
!it to zero as well, or any other number that makes sense to you.

!For instance, if there is only one consecutive ell sequence, set seqsize(2) and seqsize(3) to
!zero. There is no middle bit between two sequences, and there is no second sequence.

!The program still works, even if there are no ell=1 or ell=2 seqences. Set all seqsize(:) to zero,
!and set npreseq = nobs.

!Similarly for the ell=2 sequence

program fitper

  use common

  implicit none

  real*8 :: avg
  integer :: i, j, k1max, dum
  integer, parameter :: ncalc1min = 10, ncalc2min = 15

10 format(F8.0,4F6.0,F5.0,F6.1,F6.2,3F7.2,2F16.6)
20 format(3F5.0,F16.6)
!80 format(F8.0,4F6.0,F5.0,F6.1,F6.2,3F7.2)
80 format(15F9.1)
90 format(F12.6,2F4.0)

  open(unit=1,file="obsperiods")
  open(unit=2,file="calcperiods")
  open(unit=3,file="fitnesspars.dat")
  open(unit=5,file="failed_fullinfo.dat")
  open(unit=7,file="failed.dat")

  success = .true.
  success1 = .true.

!Read the list of observed periods

  read(1,*)
  read(1,*) nper1, nper2
  read(1,*)
  read(1,*)
  read(1,*) npreseq(:)
  read(1,*)
  read(1,*) pegseq_ell1(:)
  read(1,*)
  read(1,*) seqsize_ell1(:)
  read(1,*)
  read(1,*) pegseq_ell2(:)
  read(1,*) seqsize_ell2(:)
  read(1,*)
  nper_tot = nper1 + nper2

! Dumb user check
  do i=1,4
     if (seqsize_ell1(i) .lt. 0) then
        print *, 'Check the size of each sequence. At least one of them is negative'
        stop
     end if
     if (seqsize_ell2(i) .lt. 0) then
        print *, 'Check the size of each sequence. At least one of them is negative'
        stop
     end if
  end do

  allocate(obs_per1(nper1,3),obs_per2(nper2,3))

  if (verbose) print *, 'Observed period list read in'

  do i=1,nper1
     read(1,*) obs_per1(i,:)
     if (verbose) print 90, obs_per1(i,:)
  end do

  do i=1,nper2
     read(1,*) obs_per2(i,:)
     if (verbose) print 90, obs_per2(i,:)
  end do

  do

!Fill the calculated period array

!First figure out how many calculated periods there are and allocate array

     ncalc1 = 0
     ncalc2 = 0
     ncalctot = 0

     read(2,*,end=2) params
     if (verbose) then
        print *, ''
        print 80, params(:)

     end if
!     write(4,80) params
     read(2,*) rdum, rdum
     if (rdum .gt. 10) backspace(2)
     do
        read(2,*,end=2) rdum
        dum = int(rdum)
        if (dum .eq. 100000) exit
        if (dum .eq. 1) ncalc1 = ncalc1 + 1
        if (dum .eq. 2) ncalc2 = ncalc2 + 1
     end do

     ncalc1 = ncalc1
     ncalctot = ncalc1 + ncalc2
     if (verbose) then
        print *, ''
        print *, "ncalc1, ncalc2, ncalctot"
        print *, ncalc1, ncalc2, ncalctot
     end if
     k1max = ncalc1
     allocate(calc_per1(ncalc1,4), calc_per2(ncalc2,4), calc_per(ncalctot,4))

!Rewind to the correct spot and read again, this time filling the arrays
!Column 1 = index number
!Column 2 = l
!Column 3 = k (first period of a given l has k=1, the next one k=2, etc...)
!This is not the formal indentification, just there to give a notion of
!low k vs high k)
!Column 4 = calculated period

     do i = 1,ncalctot+2
        backspace(2)
     end do

     do i=1,ncalc1
        calc_per1(i,1) = real(i)
        calc_per1(i,3) = real(i)
        read(2,*,end=2) calc_per1(i,2), calc_per1(i,4)
     end do

     do i=1,ncalc2
        calc_per2(i,1) = real(ncalc1 + i)
        calc_per2(i,3) = real(i)
        read(2,*,end=2) calc_per2(i,2), calc_per2(i,4)
     end do

     ncalc = ncalctot
     do i=1,ncalc1
        do j=1,4
           calc_per(i,j) = calc_per1(i,j)
        enddo
     end do
     do i=ncalc1+1, ncalc
        do j=1,4
           calc_per(i,j) = calc_per2(i-ncalc1,j)
        end do
     end do

     if (verbose) then
        print *, ""
        print *, "calc_per1"
        do i=1,ncalc1
           print 20, calc_per1(i,:)
        end do
        print *, ""
        print *, "calc_per2"
        do i=1,ncalc2
           print 20, calc_per2(i,:)
        end do
     end if

!Position ourselves at beginning of next period listing then do fitting

     success = .true.
     success1 = .true.

     do
        read(2,*) rdum
        dum = int(rdum)
        if (dum .eq. 100000) exit
     end do

!Calculate fit
     if ( (ncalc1 .ge. ncalc1min) .and. (ncalc2 .ge. ncalc2min) &
          .and. (ncalc1 .ge. nper1) .and. (ncalc2 .ge. nper2) .and. &
          (calc_per1(ncalc1,4) .gt. obs_per1(nper1,1)) ) then
!.and. (calc_per2(ncalc2,4) .gt. obs_per2(nper2,1)) ) then
        allocate(fit_per(nper_tot,7))
        do i=1,nper_tot
           do j=1,7
              fit_per(i,j) = 0.
           end do
        end do
        call calcfit
        deallocate(fit_per)
     else
        if (verbose) then
           if ( (ncalc1 .lt. ncalc1min) .or. (ncalc1 .lt. nper1) .or. &
              (calc_per(ncalc1,4) .le. obs_per1(nper1,1)) ) then
              print *, 'model only has ', ncalc1, ' ell=1 periods. Not proceeding'
           else 
              print *, 'model only has ', ncalc2, ' ell=2 periods. Not proceeding'
           end if
        end if
        success = .false.
     end if
     
     if (success) then
        avg = sigma
     else
        write(7,80) params
        write(5,80) params
        write(5,*) rnumod
     end if

1    continue

     rewind(1)

     deallocate(calc_per,calc_per1,calc_per2)
  end do

2 continue

end program fitper
