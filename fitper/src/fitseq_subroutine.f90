module fitseq_subroutine

contains

!This routine fits an observed ell sequence with a sequence of same ell
!calculated periods. The observed ell sequence starts at pegseq and contains
!seqsize periods.

  subroutine fitseq(ell_fit)

    use common
    use subroutines

    implicit none


    real :: omcmin, sigma1, sigma2, sigma3
    integer, dimension(4) :: seqsizelocal
    integer :: i, j, ell_fit, nperlocal

    select case(ell_fit)
       case(1)
          nperlocal = nper1
          seqsizelocal = seqsize_ell1
       case(2)
          nperlocal = nper2
          seqsizelocal = seqsize_ell2
       case default
          print *, 'ell assignment not recognized'
          stop
    end select

    allocate(fit_per_seq(seqsize,7))

310 format(F12.5,2F4.0)
    if (verbose) then
       print *, ''
       print *, 'obs_per_dum'
       do i=1,nperlocal
          print 310, obs_per_dum(i,:)
       end do
    end if

!Save the index number of the calculated period closest to observed
!period pegseq. Period pegseq is the period in the observed period list that
!begins the sequence.

    omcmin= 100000000.
    do j = 1,ncalc_ell
       omc = abs(obs_per_dum(pegseq,1)-calc_per_ell(j,4))
1002 format(I3,4F16.6)
    if ( (ell_fit .eq. 2) .and. (verbose) ) then
        print 1002, j, omc, omcmin, obs_per_dum(pegseq,1), calc_per_ell(j,4)
    end if
       if (omc .LT. omcmin) then
          omcmin=omc
          fit_per_ell(pegseq,1) = calc_per_ell(j,3)
       end if
    end do

300 format(F5.0,F12.2,2F5.0,F12.2,F8.2,F5.0)
    if (verbose) then
       print *, ''
       print *, 'fit_per_ell for fit', ell_fit
       do i=1,nper_ell
          print 300, fit_per_ell(i,:)
       end do
    end if

    if (secondellseq) then
       pegseqobs = seqsizelocal(1)+seqsizelocal(2)+1
    else
       pegseqobs = 1
    end if

    if (verbose) then
       print *, ''
       print *, 'success1 ', success1
       print *, 'pegseq ', pegseq
       print *, 'fit_per_ell(pegseq, 4) ', fit_per_ell(pegseq,4)
       print *, 'seqsizelocal(2) ', seqsizelocal(2)
    end if

    if (success1) then
       peg0 = fit_per_ell(pegseq,1)
    else
       peg0 = fit_per_ell(pegseq,4) + seqsizelocal(2) + 1
    end if

    if (verbose) then
       print *, ''
       print *, 'peg0 ', peg0
    end if

    peg = peg0
    call fillseq(ell_fit)

    sigma1 = 0.
    do i=1,seqsize
       sigma1 = sigma1 + (fit_per_seq(i,5)-fit_per_seq(i,2))**2
    end do

!Try to "slide" the sequence up or down by one and see if either works better
    peg = peg0 - 1
    call fillseq(ell_fit)
    sigma2 = 0.
    do i=1,seqsize
       sigma2 = sigma2 + (fit_per_seq(i,5)-fit_per_seq(i,2))**2
    end do
!If we are fitting the second sequence, we have to make sure that we are not
!trying to reuse a period already used to fit the middle bit.
    if (ifit2ndseq) then
       if (verbose) then
          print *, ''
          print *, 'fit_per_seq(1,1), fit_per_ell(pegseq-1,4)'
          print *, fit_per_seq(1,1), fit_per_ell(pegseq-1,4)
       end if
       if (fit_per_seq(1,1) .eq. fit_per_ell(pegseq-1,4)) sigma2=huge
    end if

    peg = peg0 + 1
    call fillseq(ell_fit)
    sigma3 = 0.
    do i=1,seqsize
       sigma3 = sigma3 + (fit_per_seq(i,5)-fit_per_seq(i,2))**2
    end do

    if ( (sigma1 .le. sigma2) .and. (sigma1 .lt. sigma3) ) peg = peg0
    if ( (sigma2 .lt. sigma1) .and. (sigma2 .lt. sigma3) ) peg = peg0 - 1
    if ( (sigma3 .lt. sigma1) .and. (sigma3 .lt. sigma2) ) peg = peg0 + 1

    sigma = min(sigma1, sigma2, sigma3)
    call fillseq(ell_fit)

  end subroutine fitseq

end module fitseq_subroutine
