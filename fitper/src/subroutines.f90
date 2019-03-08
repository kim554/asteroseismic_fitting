module subroutines

contains

!**************************************************************************

  subroutine fillseq(ell_fit)

    use common

    implicit none

    integer, dimension(4) :: seqsize_local
    integer :: i, nperlocal, ell_fit

50 format(F6.0, F12.4, F6.0, 4F12.4)


!The columns of fit_per and fit_per_seq are the following
!Column 1 = index number
!Column 2 = observed period
!Column 3 = ell of calculated period
!Column 4 = k of calculated period
!Column 5 = calculated period
!Column 6 = calculated - observed
!Column 7 = ell specified in obsperiods

    if (verbose) then
       print *, ''
       print *, 'pegseqobs ', pegseqobs
    end if

!Extract the consecutive modes matrice

    select case(ell_fit)
       case(1)
          seqsize_local(:) = seqsize_ell1(:)
          nperlocal = nper1
       case(2)
          seqsize_local(:) = seqsize_ell2(:)
          nperlocal = nper2
       case default
          print *, 'ell assignment not recognized'
          stop
    end select

    if (verbose) then
        print *, ''
        print *, 'peg ', peg
    end if

    if (success1 .eqv. .false.) then
       pegseq = pegseq + seqsize_local(2)
    end if

    if (verbose) print *, ""


    if (ell_fit .eq. 1) then
!Make sure we don't try to go above the top edge of the array
       if (peg .lt. 1) then
          if (verbose) then
             print *, "flirting with the top of the array", peg, "is less than 1"
          end if
          do i = 1, seqsize
             fit_per_seq(i,1) = calc_per_ell(i,1)
             fit_per_seq(i,2) = obs_per_dum(pegseqobs-1+i,1)
             fit_per_seq(i,3) = calc_per_ell(i,2)
             fit_per_seq(i,4) = calc_per_ell(i,3)
             fit_per_seq(i,5) = calc_per_ell(i,4)
             fit_per_seq(i,6) = obs_per_dum(pegseqobs-1+i,1) - &
                  calc_per_ell(i,4)
             fit_per_seq(i,7) = obs_per_dum(pegseqobs-1+i,2)
          end do
       else
          do i = 1, seqsize
!Make sure we don't go past the bottom edge of calc_per_ell
             if (peg-1+i .gt. ncalc_ell) exit
             if (verbose) then
                print *, "pegseq-1+i", pegseq-1+1
                print *, "pegseqobs", pegseqobs
             end if
             fit_per_seq(i,1) = calc_per_ell(peg-1+i,1)
             fit_per_seq(i,2) = obs_per_dum(pegseqobs-1+i+npreseq(ell_fit),1)
             fit_per_seq(i,3) = calc_per_ell(pegseq-1+i,2)
             fit_per_seq(i,4) = calc_per_ell(pegseq-1+i,3)
             fit_per_seq(i,5) = calc_per_ell(peg-1+i,4)
             fit_per_seq(i,6) = obs_per_dum(pegseqobs-1+i,1) - calc_per_ell(peg-1+i,4)
             fit_per_seq(i,7) = obs_per_dum(pegseqobs-1+i,2)
          end do
       end if
    end if


    if (ell_fit .eq. 2) then
!Make sure we don't try to go above the top edge of the array
       if (peg .lt. 1) then
          if (verbose) then
             print *, "flirting with the top of the array", peg, "is less than 1"
          end if
          do i = 1, seqsize
             fit_per_seq(i,1) = calc_per_ell(i,1)
             fit_per_seq(i,2) = obs_per_dum(pegseqobs-1+i,1)
             fit_per_seq(i,3) = calc_per_ell(i,2)
             fit_per_seq(i,4) = calc_per_ell(i,3)
             fit_per_seq(i,5) = calc_per_ell(i,4)
             fit_per_seq(i,6) = obs_per_dum(pegseqobs-1+i,1) - &
                  calc_per_ell(i,4)
             fit_per_seq(i,7) = obs_per_dum(pegseqobs-1+i,2)
          end do
       else
          do i = 1, seqsize
!Make sure we don't go past the bottom edge of calc_per_ell
             if (peg-1+i .gt. ncalc_ell) exit
             if (verbose) then
                print *, "pegseq-1+i", pegseq-1+1
                print *, "pegseqobs", pegseqobs
             end if
             fit_per_seq(i,1) = calc_per_ell(peg-1+i,1)
             fit_per_seq(i,2) = obs_per_dum(pegseqobs-1+i+npreseq(ell_fit),1)
             fit_per_seq(i,3) = calc_per_ell(pegseq-1+i,2)
             fit_per_seq(i,4) = calc_per_ell(pegseq-1+i,3)
             fit_per_seq(i,5) = calc_per_ell(peg-1+i,4)
             fit_per_seq(i,6) = obs_per_dum(pegseqobs-1+i,1) - calc_per_ell(peg-1+i,4)
             fit_per_seq(i,7) = obs_per_dum(pegseqobs-1+i,2)
          end do
       end if
    end if


    if (verbose) then
       print *, ''
       print *, 'seqsize', seqsize
       print *, "fit_per_seq"
       do i=1,seqsize
          print 50, fit_per_seq(i,:)
       end do
    end if

    do i = 1, ncalc_ell
       calc_per_copy(i,:) = calc_per_ell(i,:)
    end do

  end subroutine fillseq

!**************************************************************************
!Fits the periods before the ell=1 sequence following instructions given
!in obsperiods
!0 Whatever fits best
!1 ell=1
!2 ell=2

subroutine fitrest1

    use common

    implicit none

    real :: omcmin1, omcmin2
    real, parameter :: big = 1.0E12
    integer :: i, j, jmin1, jmin2, typeoffit, n

50  format(2I5, 3F5.0, F22.12)
60  format(I5, 3F5.0, F22.12)

!Perform the matchup

    do i=1,npreseq_ell
       omcmin1 = big
       omcmin2 = big
       typeoffit = int(obs_per_dum(i,2))

       select case(typeoffit)
       case(0)
          do j=1,ncalc2
             omc = abs(obs_per_dum(i,1) - calc_per2(j,4))
             if (omc.lt.omcmin2) then
                omcmin2 = omc
                jmin2 = j
             end if
          end do
          do j=1,ncalc1
             omc = abs(obs_per_dum(i,1) - calc_per1(j,4))
             if (omc.lt.omcmin1) then
                omcmin1 = omc
                jmin1 = j
             end if
          end do
       case(1)
          do j=1,ncalc1
             omc = abs(obs_per_dum(i,1) - calc_per1(j,4))
             if (omc.lt.omcmin1) then
                omcmin1 = omc
                jmin1 = j
             end if
          end do
       case(2)
          do j=1,ncalc2
             omc = abs(obs_per_dum(i,1) - calc_per2(j,4))
             if (omc.lt.omcmin2) then
                omcmin2 = omc
                jmin2 = j
             end if
          end do
       case default
          print *, 'ell assignment not recognized'
          stop
       end select

       if (omcmin1 .le. omcmin2) then
          fit_per_ell(i,1) = calc_per1(jmin1,1)
          fit_per_ell(i,2) = obs_per_dum(i,1)
          fit_per_ell(i,3) = calc_per1(jmin1,2)
          fit_per_ell(i,4) = calc_per1(jmin1,3)
          fit_per_ell(i,5) = calc_per1(jmin1,4)
          fit_per_ell(i,6) = obs_per_dum(i,1) - calc_per1(jmin1,4)
          fit_per_ell(i,7) = obs_per_dum(i,2)
       else
          fit_per_ell(i,1) = calc_per2(jmin2,1)
          fit_per_ell(i,2) = obs_per_dum(i,1)
          fit_per_ell(i,3) = calc_per2(jmin2,2)
          fit_per_ell(i,4) = calc_per2(jmin2,3)
          fit_per_ell(i,5) = calc_per2(jmin2,4)
          fit_per_ell(i,6) = obs_per_dum(i,1) - calc_per2(jmin2,4)
          fit_per_ell(i,7) = obs_per_dum(i,2)
       end if

    end do

  end subroutine fitrest1

!**************************************************************************
!Fits the periods after the ell=1 sequence following instructions given
!in obsperiods
!0 Whatever fits best
!1 ell=1
!2 ell=2

  subroutine fitrest2(ell_fit)

    use common

    implicit none

    real, dimension(:,:), allocatable :: calc_per_leftover, obs_per_copy
    real, dimension(4) :: pegseq_local
    real :: omcmin, omcmin1, omcmin2
    real, parameter :: big = 1.0E12
    integer :: i, j, imin, jmin1, jmin2, typeoffit, nperlocal, ncalclocal
    integer :: seqsizelocal, ell_fit, ibest2ndseq, ilastused, marker

!Some prep work to get values of constants we need for this routine
!and to prepare some arrays we want to use.

    select case(ell_fit)
       case(1)
          allocate(obs_per_copy(nper1,3))
          obs_per_copy(:,:) = obs_per1(:,:)
          pegseq_local(:) = pegseq_ell1(:)
          nperlocal = nper1
          ncalclocal = ncalc1
          marker = seqsize_ell1(1)
       case(2)
          allocate(obs_per_copy(nper2,3))
          obs_per_copy(:,:) = obs_per2(:,:)
          pegseq_local(:) = pegseq_ell2(:)
          nperlocal = nper2
          ncalclocal = ncalc2
          marker = seqsize_ell2(1)
       case default
          print *, 'ell assignment not recognized'
          stop
    end select

290 format(F12.6,2F4.0)
280 format(3F5.0,F12.6)

!Messy stuff to take care of only when fitting the central bit
    if (ifitrest .eqv. .false.) then
       if (verbose) then
          print *, ''
          print *, 'Observed period list using now'
          do i=1,nperlocal
             print 290, obs_per_copy(i,:)
          end do
          print *, ''
          print *, 'calc_per_ell'
          do i=1,ncalclocal
             print 280, calc_per_ell(i,:)
          end do
       end if

       if (verbose) then
          print *, ''
          print *, 'pegseq_local(2) ', pegseq_local(2)
          print *, ''
       end if

!We need the index of the last calculated period used to fit the
!First sequence. We call that index ilastused
       ilastused = fit_per_ell(marker,1)

       if (verbose) print *, 'nper_ell, ilastused ', nper_ell, ilastused

!Next figure out index of the period in the relevant calculated
!period list (calc_per_ell) that best fits the first period of the
!second sequence in observed period list. We'll need that. We call
!that index ibest2ndseq.
       omcmin = big
       do i=1,ncalc_ell
          omc = abs(obs_per_copy(pegseq_local(3),1) - calc_per_ell(i,4))
!       if (verbose) then
!          print *, 'omc, omcmin ', omc, omcmin
!       end if
          if (omc .lt. omcmin) then
             omcmin = omc
             ibest2ndseq = i
          end if
       end do

!Safety provision in case the best fit was used in fitting the first sequence
          if (ibest2ndseq .le. ilastused) ibest2ndseq = ilastused + 1
          if (verbose) print *, 'i, ibest2ndseq ', i, ibest2ndseq

300 format(2F12.6, I4)
       if (verbose) then
          print *, ''
          print *, 'First period of second sequence and corresponding fit'
          print 300, obs_per_copy(pegseq_local(3),1), &
                 calc_per_ell(ibest2ndseq,4), ibest2ndseq
       end if
    end if

    deallocate(obs_per_copy)

!If fitting the periods after the 2nd sequence
    if (ifitrest) then
       imin = fit_per_ell(pegseq2-1,1)
       seqsizelocal = ncalc1-imin
    else
!If fitting the periods between the two sequences
       imin = ilastused
       seqsizelocal = ibest2ndseq - ilastused - 1
!Failsafe in case there are not enough periods in the calculated list to fit
!the middle sequence. We will also have to adjust where we start in the
!calculated period list when we fit the second sequence so set flag.
       if (seqsizelocal .lt. seqsize2) then
          seqsizelocal = seqsize2
          success1 = .false.
       end if
    end if

    if (verbose) then
       print *, ''
       print *, 'seqsizelocal, seqsize2', seqsizelocal, seqsize2
    end if

    allocate(calc_per_leftover(seqsizelocal,4))

50  format(I5, 3F5.0, F18.12)

!Load calculated periods in small array
    do i=1,seqsizelocal
       if (ell_fit .eq. 1) then
          calc_per_leftover(i,:) = calc_per1_copy(imin+i,:)
       else
          calc_per_leftover(i,:) = calc_per2_copy(imin+i,:)
       end if
    end do

    if (verbose) then
       print *, ''
       print *, 'calc_per_leftover'
       do i=1,seqsizelocal
          print 280, calc_per_leftover(i,:)
       end do
    end if

!Perform of matchup.

    do i=1,seqsize2
       omcmin1 = big
       omcmin2 = big
       typeoffit = int(obs_per_seq(i,2))

       select case(typeoffit)
       case(0)
          do j=1,ncalc2
             omc = abs(obs_per_seq(i,1) - calc_per2(j,4))
             if (omc.lt.omcmin2) then
                omcmin2 = omc
                jmin2 = j
             end if
          end do
          do j=1,seqsizelocal
             omc = abs(obs_per_seq(i,1) - calc_per_leftover(j,4))
             if (omc.lt.omcmin1) then
                omcmin1 = omc
                jmin1 = j
             end if
          end do
       case(1)
          do j=1,seqsizelocal
             omc = abs(obs_per_seq(i,1) - calc_per_leftover(j,4))
             if (omc.lt.omcmin1) then
                omcmin1 = omc
                jmin1 = j
             end if
          end do
       case(2)
          do j=1,ncalc2
             omc = abs(obs_per_seq(i,1) - calc_per2(j,4))
             if (omc.lt.omcmin2) then
                omcmin2 = omc
                jmin2 = j
             end if
          end do
       case default
          print *, 'ell assignment not recognized'
          stop
       end select

       if (verbose) then
          print *, ''
          print *, 'omcmin1, omcmin2 ', omcmin1, omcmin2
          print *, 'jmin1 ', jmin1
       end if

       if (omcmin1 .le. omcmin2) then
          fit_per_ell(i+pegseq2-1,1) = calc_per_leftover(jmin1,1)
          fit_per_ell(i+pegseq2-1,2) = obs_per_seq(i,1)
          fit_per_ell(i+pegseq2-1,3) = calc_per_leftover(jmin1,2)
          fit_per_ell(i+pegseq2-1,4) = calc_per_leftover(jmin1,3)
          fit_per_ell(i+pegseq2-1,5) = calc_per_leftover(jmin1,4)
          fit_per_ell(i+pegseq2-1,6) = obs_per_seq(i,1) - calc_per_leftover(jmin1,4)
          fit_per_ell(i+pegseq2-1,7) = obs_per_seq(i,2)
       else
          fit_per_ell(i+pegseq2-1,1) = calc_per2(jmin2,1)
          fit_per_ell(i+pegseq2-1,2) = obs_per_seq(i,1)
          fit_per_ell(i+pegseq2-1,3) = calc_per2(jmin2,2)
          fit_per_ell(i+pegseq2-1,4) = calc_per2(jmin2,3)
          fit_per_ell(i+pegseq2-1,5) = calc_per2(jmin2,4)
          fit_per_ell(i+pegseq2-1,6) = obs_per_seq(i,1) - calc_per2(jmin2,4)
          fit_per_ell(i+pegseq2-1,7) = obs_per_seq(i,2)
       end if

    end do

    deallocate(calc_per_leftover)

270 format(F5.0,F12.2,2F5.0,F12.2,F8.2,F5.0)
    if (verbose) then
      print *, ''
      print *, 'fit_per_ell at the end of fitrest2'
      do i=1,nper_ell
         print 270, fit_per_ell(i,:)
      end do
    end if

  end subroutine fitrest2

!**************************************************************************

end module subroutines
