module match_periods_subroutine

contains

!This is a key routine. It performs the period matching. It needs as input
!a period list to match (observed periods), and a number of integers that
!mark where the consecutive ell sequences start and end.

  subroutine match_periods(pegseq_dum,seqsize_dum,ell_fit)

    use common
    use subroutines
    use fitseq_subroutine

    implicit none

    real :: omcmin, period
    integer, dimension(4) :: pegseq_dum, seqsize_dum
    integer :: i, j, ell_fit, imin, jmax, imax2

40  format(3F4.0,F16.8)
60  format(4F16.2)
90  format(F12.6,2F4.0)

    ifit2ndseq = .false.
    allocate(calc_per_copy(ncalc,4),calc_per1_copy(ncalc1,4), calc_per2_copy(ncalc2,4))

    if ( (ell_fit .ne. 1) .and. (ell_fit .ne. 2) ) then
       print *, 'Error, invalid ell_fit ', ell_fit
       stop
    end if

!!$ Fill calculated period arrays we are going to use. One is all of them, the others
!!$ ell=1 or ell=2.

    do i = 1, ncalctot
       calc_per_copy(i,:) = calc_per(i,:)
    end do

    do i = 1, ncalc1
       calc_per1_copy(i,:) = calc_per1(i,:)
    end do

    do i = 1, ncalc2
       calc_per2_copy(i,:) = calc_per2(i,:)
    end do

!A little check. Sometimes, the first ell=1 period is a repeat of a period
!further down the line. We want to eliminate it (or replace it by zero)

    if (calc_per_copy(1,4) .ge. calc_per_copy(2,4)) calc_per_copy(1,4) = 0.

    if (ell_fit .eq. 1) then
       ncalc_ell = ncalc1
       nper_ell = nper1
       npreseq_ell = npreseq(1)
       allocate(calc_per_ell(ncalc1,4),obs_per_dum(nper1,3))
       if (verbose) then
          print *, ""
          print *, "calc_per_ell"
       end if
       do i = 1, ncalc1
          calc_per_ell(i,:) = calc_per1(i,:)
          if (verbose) then
             print 40, calc_per_ell(i,:)
          end if
       end do
       do i = 1, nper1
          obs_per_dum(i,:) = obs_per1(i,:)
       end do
    else
       ncalc_ell = ncalc2
       nper_ell = nper2
       npreseq_ell = npreseq(2)
       allocate(calc_per_ell(ncalc2,4),obs_per_dum(nper2,3))
       do i = 1, ncalc2
          calc_per_ell(i,:) = calc_per2(i,:)
          if (verbose) then
             print 40, calc_per_ell(i,:)
          end if
       end do
       do i = 1, nper2
          obs_per_dum(i,:) = obs_per2(i,:)
       end do
    end if

    allocate(fit_per_ell(nper_ell,7))
    do i=1,nper_ell
       do j=1,7
          fit_per_ell(i,j) = 0.
       end do
    end do

!!$ Fill the consecutive ell sequence of observed periods (if there is one)
!!$ Location in the observed period list of the first period in the sequence and number
!!$ of observed periods in the sequence

    if (seqsize_dum(1) .gt. 0) then

       pegseq = pegseq_dum(1)
       seqsize = seqsize_dum(1)
       allocate(obs_per_seq(seqsize,3))

       if (verbose) then
          print *, ""
          print *, 'obs_per_seq, first ell sequence'
       end if
       do i=1,seqsize
          j=i+pegseq-1
          obs_per_seq(i,:) = obs_per_dum(j,:)
          if (verbose) print 90, obs_per_seq(i,:)
       end do

!!$ Fit the first consecutive ell sequence.

       call fitseq(ell_fit)
       if (success1 .eqv. .false.) pegseq = pegseq+1
       do i=pegseq,pegseq+seqsize-1
          fit_per_ell(i,1) = fit_per_seq(i-pegseq+1,1)
          fit_per_ell(i,2) = obs_per_dum(i,1)
          fit_per_ell(i,3) = fit_per_seq(i-pegseq+1,3)
          fit_per_ell(i,4) = fit_per_seq(i-pegseq+1,4)
          fit_per_ell(i,5) = fit_per_seq(i-pegseq+1,5)
          fit_per_ell(i,6) = obs_per_dum(i,1) - fit_per_seq(i-pegseq+1,5)
          fit_per_ell(i,7) = obs_per_dum(i,2)

!          fit_per_ell(i,1) = fit_per_seq(i-pegseq+1,1)
!          fit_per_ell(i,2) = obs_per_dum(i,1)
!          fit_per_ell(i,3) = fit_per_seq(i-pegseq+1,3)
!          fit_per_ell(i,4) = fit_per_seq(i-pegseq+1,4)
!          fit_per_ell(i,5) = fit_per_seq(i-pegseq+1,5)
!          fit_per_ell(i,6) = obs_per_dum(i,1) - fit_per_seq(i-pegseq+1,5)
!          fit_per_ell(i,7) = obs_per_dum(i,2)
       end do


       if (verbose) then
          print *, ''
          print *, "fit_per_ell, after fitting first ell sequence"
          do i=1,nper_ell
             print 100, fit_per_ell(i,:)
          end do
       end if

       deallocate(fit_per_seq,obs_per_seq)

    end if

!!$ Fit the bit in between consecutive ell sequences (if there is one)

    if (seqsize_dum(2) .gt. 0) then

       pegseq2 = pegseq_dum(2)
       seqsize2 = seqsize_dum(2)
       allocate(obs_per_seq(seqsize2,3))

!!$ Fill the observed periods (it's not really a sequence at this point)

       if (verbose) then
          print *, ''
          print *, 'obs_per_seq for middle bit'
       end if
       do i=1,seqsize2
          j=i+pegseq2-1
          obs_per_seq(i,:) = obs_per_dum(j,:)
          if (verbose) print 90, obs_per_seq(i,:)
       end do

       ifitrest = .false.
       call fitrest2(ell_fit)

       if (verbose) then
          print *, ''
          print *, "fit_per_ell, after fitting the middle bit"
          do i=1,nper_ell
             print 100, fit_per_ell(i,:)
          end do
       end if

       deallocate(obs_per_seq)

    end if
    success1 = .true.

!!$ Work on the next consecutive ell sequence (if there is one)

    if (seqsize_dum(3) .gt. 0) then

       secondellseq = .true.
       ifit2ndseq = .true.
       if (success1) then
          pegseq = pegseq_dum(3)
       else
          pegseq = seqsize_dum(1) + seqsize_dum(2)
       end if
       seqsize = seqsize_dum(3)
       allocate(obs_per_seq(seqsize,3))

!!$ Fill the ell sequence of observed periods

       if (verbose) then
          print *, ''
          print *, 'obs_per_seq for second ell sequence'
       end if

       jmax = pegseq-1
       imin = 1
       imax2 = seqsize

!    do i=imin,imax2
!       j=i+jmax
!       obs_per_seq(i,:) = obs_per_dum(j,:)
!       if (verbose) print 90, obs_per_seq(i,:)
!    end do

       call fitseq(ell_fit)

       if (ell_fit .eq. 1) then

          do i=pegseqobs,pegseqobs+seqsize-1
             fit_per_ell(i,1) = fit_per_seq(i-pegseqobs+1,1)
             fit_per_ell(i,2) = obs_per_dum(i,1)
             fit_per_ell(i,3) = fit_per_seq(i-pegseqobs+1,3)
             fit_per_ell(i,4) = fit_per_seq(i-pegseqobs+1,4)
             fit_per_ell(i,5) = fit_per_seq(i-pegseqobs+1,5)
             fit_per_ell(i,6) = obs_per_dum(i,1) - fit_per_seq(i-pegseqobs+1,5)
             fit_per_ell(i,7) = obs_per_dum(i,2)
          end do

       end if

       if (ell_fit .eq. 2) then
    
          do i=pegseqobs,pegseqobs+seqsize-1
             fit_per_ell(i+npreseq(ell_fit),1) = fit_per_seq(i-pegseqobs+1,1)
             fit_per_ell(i+npreseq(ell_fit),2) = obs_per_dum(i+npreseq(ell_fit),1)
             fit_per_ell(i+npreseq(ell_fit),3) = fit_per_seq(i-pegseqobs+1,3)
             fit_per_ell(i+npreseq(ell_fit),4) = fit_per_seq(i-pegseqobs+1,4)
             fit_per_ell(i+npreseq(ell_fit),5) = fit_per_seq(i-pegseqobs+1,5)
             fit_per_ell(i+npreseq(ell_fit),6) = obs_per_dum(i,1) - fit_per_seq(i-pegseqobs+1,5)
             fit_per_ell(i+npreseq(ell_fit),7) = obs_per_dum(i,2)
          end do
   
       end if

       if (verbose) then
          print *, ''
          print *, "fit_per_ell, after fitting second ell sequence for ell =", ell_fit
          do i=1,nper_ell
             print 100, fit_per_ell(i,:)
          end do
       end if
       deallocate(fit_per_seq, obs_per_seq)

       success1 = .true.

    end if

!!$Work on the rest of the observed period list (if there is one)

    if (seqsize_dum(4) .gt. 0) then

       pegseq2 = pegseq_dum(4)
       seqsize2 = seqsize_dum(4)
       allocate(obs_per_seq(seqsize2,3))

!!$ Fill the observed periods (it's not really a sequence at this point)

       if (verbose) then
          print *, ''
          print *, 'obs_per_seq after second ell sequence'
       end if
       do i=1,seqsize2
          j=i+pegseq2-1
          obs_per_seq(i,:) = obs_per_dum(j,:)
          if (verbose) print 90, obs_per_seq(i,:)
       end do

       ifitrest = .true.
       call fitrest2(ell_fit)
       deallocate(obs_per_seq)

30     format(3F8.1,2I4)
100    format(F5.0, F12.2, 2F5.0, F12.2, F8.2, F5.0)

       if (verbose) then
          print *, ''
          print *, 'fit_per_ell, after fitting periods at the end of the list'
          do i=1,nper_ell
             print 100, fit_per_ell(i,:)
          end do
          print *, ''
       end if
    end if

!!$Last, we need to take care of the period(s) before the first ell=1 sequence
    if (success) then
       if (npreseq_ell .gt. 0) then
          call fitrest1
          if (verbose) then
             print *, ''
             print *, 'fit_per_ell, after fitting periods at the beginning of the list'
             do i=1,nper_ell
                print 100, fit_per_ell(i,:)
             end do
             print *, ''
          end if
       end if
    end if

    deallocate(calc_per_copy,calc_per1_copy, calc_per2_copy)
    deallocate(calc_per_ell,obs_per_dum)

  end subroutine match_periods

end module match_periods_subroutine
