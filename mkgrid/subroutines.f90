module subroutines

contains

!*****************************************************************************

  subroutine assign_data_arrays

    use common_values

    implicit none

! For subroutine calc_menv    

!                   n=4 n=3 n=2        n=1                  n=0
    coeffs_menv = (/ 0., 0., 0., 0.399661453893467, -77.5354961002561 /)

! For subroutine calc_mhe
!                 n=4       n=3                  n=2                n=1                  n=0
    coeffs_mhe = (/0., -3.3311415683546e-6, 0.007656935511203, -5.40570588350895, 1420.699975197 /)

! For subroutine calc_mh
!                 n=4 n=3 n=2                n=1                  n=0
    coeffs_mh =  (/0., 0., 0., 0.409728844615393, 145.723487268456 /)
    
! For subroutine calc_xhebar
!                    n=4 n=3          n=2                n=1                  n=0
    coeffs_xhebar = (/0., 0., 0.00166075407339736, -1.99185174214316, 666.135254017898 /)

! For subroutine calc_alph1
!                    n=4            n=3              n=2                n=1                  n=0
    coeffs_alph1 = (/ 0., -4.951628484337e-7, 0.0011739125053575, -0.940072408489723, 261.65313829671 /)
    
! For subroutine calc_alph2
!                    n=4            n=3              n=2                n=1                  n=0
    coeffs_alph2 = (/ 0.,  -9.835868293707e-7, 0.002279872026178,  -1.7667474777928429, 463.175286878752 /)

! For subroutine calc_h1
    
    mass_points_for_h1 = (/ 550., 570., 600., 610., 660., 700., 770., 830. /)
!           n=4       n=3                      n=2                n=1                  n=0
    cy1 = (/ 0., 0.,                  0.,                -0.0434782608695652, 93.8260869565217 /)
    cy2 = (/ 0., 0.,                  0.00202190543193,  -2.3067405143413,    726.911824468814 /)
    cy3 = (/ 0., 0.,                 -0.003958104981123,  5.04244306418272,  -1530.86189258329 /)
    cy4 = (/ 0., -9.267138030101e-06, 0.021199481108407, -16.138083337802,    4153.90614722919 /)
    cy5 = (/ 0., 0.,                 -0.000212984245514,  0.313510528029639, -49.1564251078307 /)

    cfa = (/ 0., 0.000096822212763,  -0.161422897647,     89.647724215,      -16514.7125992    /)
    cfb = (/ 0., 0.00025264008691,   -0.4514403308,       268.97756988,      -53367.63825      /)
    cfc = (/ 0., 0.00003162606767,   -0.06344945469,      42.24207945,       -9260.55792       /)
    cfd = (/ 0., 0.00000915347912,   -0.02257703047,      18.51299278,       -4981.948353      /)

! For subroutine calc_h2    
!                 n=4          n=3                 n=2                n=1                  n=0
    coeffs_h2 = (/ 0., 7.104859260031e-7, -0.0014507181916943, 0.941204705302367, -137.842876253142/)

! For subroutine calc_h3    
!                 n=4         n=3                  n=2                n=1               n=0
    coeffs_h3 = (/ 0., -3.398627130502e-7, 0.000927762414475, -0.779601007843499, 244.44712006846 /)

! For subroutine calc_w1
!                          n=4               n=3             n=2           n=1            n=0
    coeffs_w1 = (/ 0.0000000059452392, -0.000016437679, 0.017367411, -8.2702334447, 1532.2685599 /)
 
! For subroutine calc_w2
!                 n=4 n=3 n=2           n=1                n=0
    coeffs_w2 = (/ 0., 0., 0., -0.0325440793775398, 29.1388488682968 /)

! For subroutine calc_w3
!                          n=4               n=3              n=2           n=1             n=0
    coeffs_w3 = (/ -0.0000000061747678, 0.000017131709, -0.0182384149, 8.8248047556, -1586.2564762 /)
    
! For subroutine calc_w1
!                 n=4 n=3              n=2              n=1                 n=0
    coeffs_w4 = (/ 0., 0., -0.000153229977608543, 0.123946811311292, -17.8603370485237 /)

! For subroutine calc_oxygen_profile_crude     
    mass_points = (/ 525. ,536.5, 559., 581.5, 601., 620.5, 646., 682.5, 737.5, 803.5, 857. /)
    
    hs_fid(1,:) = (/68., 70., 66., 67., 72., 76., 73., 66., 65., 62., 62. /) 
    hs_fid(2,:) = (/56., 53., 56., 56., 58., 57., 55., 54., 51., 50., 51. /)
    hs_fid(3,:) = (/42., 40., 38., 38., 36., 37., 37., 37., 39., 42., 45. /)
    
    ws_fid(1,:) = (/58., 45., 52., 52., 40., 39., 41., 44., 49., 53., 66. /)
    ws_fid(2,:) = (/ 2., 17.,  3.,  2.,  9., 10., 10.,  4.,  2.,  3.,  1. /)
    ws_fid(3,:) = (/31., 30., 39., 42., 48., 49., 47., 50., 46., 42., 31. /)
    ws_fid(4,:) = (/ 5.,  4.,  3.,  2.,  1.,  1.,  1.,  1.,  1.,  1.,  1. /)
    
  end subroutine assign_data_arrays
    
!*****************************************************************************

  subroutine calc_menv(input,mass, menv)

    use common_values

    implicit none

    real, intent(in) :: mass,input
    real, intent(out) :: menv
    real :: minmenv
    real, dimension(5) :: coeffs

    coeffs = coeffs_menv

    if (input .lt. 0.d0) then
       menv = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    else
       menv = input
    end if

    if (mass .le. 600) minmenv = 120. ! Can go down to 120.
    if ((mass .gt. 600) .and. (mass .le. 700)) minmenv = 140. ! Can go down to 140.
    if (mass .gt. 700) minmenv = 180. ! can go down to 180.
       
    if (menv .lt. minmenv) menv = minmenv

  end subroutine calc_menv

    
!*****************************************************************************
    
  
  subroutine calc_mhe(input, mass, menv, mhe)

    use common_values

    implicit none

    real, intent(in) :: input, menv, mass
    real, intent(out) :: mhe
    real, parameter :: minmhe = 210.
    real, dimension(5) :: coeffs

    coeffs = coeffs_mhe
    
    if (input .lt. 0.d0) then
       if (mass .lt. 600) then
          mhe = minmhe
       else
          mhe = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
       end if
    else
       mhe = input
    end if
    
    if (mhe .lt. menv) mhe = menv

  end subroutine calc_mhe
  
!*****************************************************************************

  subroutine calc_mh(input, mass, mhe, mh)

    use common_values

    real, intent(in) :: input, mass, mhe
    real, intent(out) :: mh
    real :: minmh
    real, dimension(5) :: coeffs

    coeffs = coeffs_mh
    minmh = mhe + 150.
    
    if (input .lt. 0.d0) then
       mh = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    else
       mh = input
    end if
    
    if (mh .lt. minmh) mh = minmh
     
  end subroutine calc_mh
  
!*****************************************************************************
  subroutine calc_xhebar(input, mass, xhebar)

    use common_values

    implicit none

    real, intent(in) :: input, mass
    real, intent(out) :: xhebar
    real, dimension(5) :: coeffs

    coeffs = coeffs_xhebar

    if (input .lt. 0.d0) then
       if (mass .lt. 600) then
          xhebar = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
       else
          xhebar = 69.
       end if
    else
       xhebar = input
    end if
    
    if (xhebar .gt. 100.) xhebar = 100.

  end subroutine calc_xhebar

!*****************************************************************************
  subroutine calc_alph1(input, mass, alph1)

    use common_values

    implicit none

    real, intent(in) :: input, mass
    real, intent(out) :: alph1
    real, dimension(5) :: coeffs
    real, parameter :: minalph1 = 4.d0, maxalph1 = 20.d0
    
    coeffs = coeffs_alph1
    
    if (input .lt. 0.d0) then
       alph1 = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    else
       alph1 = input
    end if

    if (alph1 .lt. minalph1) alph1 = minalph1
    if (alph1 .gt. maxalph1) alph1 = maxalph1

  end subroutine calc_alph1

!*****************************************************************************

  subroutine calc_alph2(input, mass, alph2)

    use common_values

    implicit none

    real, intent(in) :: input, mass
    real, intent(out) :: alph2
    real, dimension(5) :: coeffs
    real, parameter :: minalph2 = 4.d0, maxalph2 = 20.d0, avg_alph2 = 10.

    coeffs = coeffs_alph2

    if (input .lt. 0.d0) then
       if (mass .gt. 600) then
          alph2 = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
       else
          alph2 = avg_alph2
       end if
    else
       alph2 = input
    end if

    if (alph2 .lt. minalph2) alph2 = minalph2
    if (alph2 .gt. maxalph2) alph2 = maxalph2

  end subroutine calc_alph2
  
!*****************************************************************************
  subroutine calc_oxygen_profile_crude(mass,hs,ws)

    use common_values

    implicit none

    real, intent(in) :: mass
    real, dimension(3), intent(out) :: hs
    real, dimension(4), intent(out) :: ws
    integer :: i
    
    if (mass .lt. mass_points(1)) then
       do i=1,3
          hs(i) = hs_fid(i,1)
          ws(i) = ws_fid(i,1)
       end do
       ws(4) = ws_fid(4,1)
    end if
    if ( (mass .ge. mass_points(1)) .and. (mass .lt. mass_points(2)) ) then
       do i=1,3
          hs(i) = hs_fid(i,1)
          ws(i) = ws_fid(i,1)
       end do
       ws(4) = ws_fid(4,1)
    end if
    if ( (mass .ge. mass_points(2)) .and. (mass .lt. mass_points(3)) ) then
       do i=1,3
          hs(i) = hs_fid(i,2)
          ws(i) = ws_fid(i,2)
       end do
       ws(4) = ws_fid(4,2)
    end if
    if ( (mass .ge. mass_points(3)) .and. (mass .lt. mass_points(4)) ) then
       do i=1,3
          hs(i) = hs_fid(i,3)
          ws(i) = ws_fid(i,3)          
       end do
       ws(4) = ws_fid(4,3)
    end if
    if ( (mass .ge. mass_points(4)) .and. (mass .lt. mass_points(5)) ) then
       do i=1,3
          hs(i) = hs_fid(i,4)
          ws(i) = ws_fid(i,4)
       end do
       ws(4) = ws_fid(4,4)
    end if
    if ( (mass .ge. mass_points(5)) .and. (mass .lt. mass_points(6)) ) then
       do i=1,3
          hs(i) = hs_fid(i,5)
          ws(i) = ws_fid(i,5)
       end do
       ws(4) = ws_fid(4,5)
    end if
    if ( (mass .ge. mass_points(6)) .and. (mass .lt. mass_points(7)) ) then
       do i=1,3
          hs(i) = hs_fid(i,6)
          ws(i) = ws_fid(i,6)
       end do
       ws(4) = ws_fid(4,6)
    end if
    if ( (mass .ge. mass_points(7)) .and. (mass .lt. mass_points(8)) ) then
       do i=1,3
          hs(i) = hs_fid(i,7)
          ws(i) = ws_fid(i,7)
       end do
       ws(4) = ws_fid(4,7)
    end if
    if ( (mass .ge. mass_points(8)) .and. (mass .lt. mass_points(9)) ) then
       do i=1,3
          hs(i) = hs_fid(i,8)
          ws(i) = ws_fid(i,8)
       end do
       ws(4) = ws_fid(4,8)
    end if
    if ( (mass .ge. mass_points(9)) .and. (mass .lt. mass_points(10)) ) then
       do i=1,3
          hs(i) = hs_fid(i,9)
          ws(i) = ws_fid(i,9)
       end do
       ws(4) = ws_fid(4,9)
    end if
    if ( (mass .ge. mass_points(10)) .and. (mass .lt. mass_points(11)) ) then
       do i=1,3
          hs(i) = hs_fid(i,10)
          ws(i) = ws_fid(i,10)
       end do
       ws(4) = ws_fid(4,10)
    end if
    if ( mass .gt. mass_points(11) ) then
       do i=1,3
          hs(i) = hs_fid(i,11)
          ws(i) = ws_fid(i,11)
       end do
       ws(4) = ws_fid(4,11)
    end if
    
  end subroutine calc_oxygen_profile_crude

!*****************************************************************************
  subroutine calc_oxygen_profile(mass,hs,ws)

    use common_values

    implicit none

    real, intent(in) :: mass
    real, dimension(3), intent(out) :: hs
    real, dimension(4), intent(out) :: ws
    integer :: i

    call calc_h1(mass, hs(1))
    call calc_h2(mass, hs(2))
    call calc_h3(mass, hs(3))

    call calc_w1(mass, ws(1))
    call calc_w2(mass, ws(2))
    call calc_w3(mass, ws(3))
    call calc_w4(mass, ws(4))    
   
  end subroutine calc_oxygen_profile

!*****************************************************************************
! See AK research notebook 04/01/21 and associated spreadsheet
  
  subroutine calc_h1(mass, h)

    use common_values

    implicit none

    real, intent(in) :: mass
    real :: h
    real, dimension(5) :: coeffs
    real :: minh1 = 70.

    if (mass .le. mass_points_for_h1(1)) then
       !h = minh1
       coeffs = cy1
       h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    end if
    if ( (mass .gt. mass_points_for_h1(1)) .and. (mass .le. mass_points_for_h1(2)) ) then
       coeffs = cfa
       h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    end if
    if ( (mass .gt. mass_points_for_h1(2)) .and. (mass .le. mass_points_for_h1(3)) ) then
       coeffs = cy2
       h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    end if
    if ( (mass .gt. mass_points_for_h1(3)) .and. (mass .le. mass_points_for_h1(4)) ) then
       coeffs = cfb
       h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    end if
    if ( (mass .gt. mass_points_for_h1(4)) .and. (mass .le. mass_points_for_h1(5)) ) then
       coeffs = cy3
       h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    end if   
    if ( (mass .gt. mass_points_for_h1(5)) .and. (mass .le. mass_points_for_h1(6)) ) then
       coeffs = cfc
       h =  coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    end if
    if ( (mass .gt. mass_points_for_h1(6)) .and. (mass .le. mass_points_for_h1(7)) ) then
       coeffs = cy4
       h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    end if
    if ( (mass .gt. mass_points_for_h1(7)) .and. (mass .le. mass_points_for_h1(8)) ) then
       coeffs = cfd
       h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    end if
    if (mass .gt. mass_points_for_h1(8)) then
       coeffs = cy5
       h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    end if

  end subroutine calc_h1

!*****************************************************************************
  
  subroutine calc_h2(mass, h)

    use common_values

    implicit none

    real, intent(in) :: mass
    real :: h
    real, dimension(5) :: coeffs

    coeffs = coeffs_h2

    if (mass .lt. 600) then
       h = 56.
    else
       h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    end if

  end subroutine calc_h2

!*****************************************************************************
  
  subroutine calc_h3(mass, h)

    use common_values

    implicit none

    real, intent(in) :: mass
    real :: h
    real, dimension(5) :: coeffs

    coeffs = coeffs_h3

    h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    
  end subroutine calc_h3

!*****************************************************************************
  
  subroutine calc_w1(mass, h)

    use common_values
    
    implicit none

    real, intent(in) :: mass
    real :: h
    real, dimension(5) :: coeffs
    real, parameter :: w1_lowmass = 55., w1max = 70.

    coeffs = coeffs_w1

    if (mass .lt. 500) then
       h = w1_lowmass
    else
       h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    end if

    if (h .gt. w1max) h = w1max
    
  end subroutine calc_w1

!*****************************************************************************
  
  subroutine calc_w2(mass, h)

    use common_values

    implicit none

    real, intent(in) :: mass
    real :: h
    real, dimension(5) :: coeffs

    coeffs = coeffs_w2

    h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    
  end subroutine calc_w2

!*****************************************************************************
! The value for w3 calculated by this routine does not actually gets used,
! but including just in case.

  subroutine calc_w3(mass, h)

    use common_values

    implicit none

    real, intent(in) :: mass
    real :: h
    real, dimension(5) :: coeffs

    coeffs = coeffs_w3

    h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    
  end subroutine calc_w3

!*****************************************************************************
  
  subroutine calc_w4(mass, h)

    use common_values

    implicit none

    real, intent(in) :: mass
    real :: h
    real, dimension(5) :: coeffs

    coeffs = coeffs_w4

    h = coeffs(1)*mass**4 + coeffs(2)*mass**3 + coeffs(3)*mass**2 + coeffs(4)*mass + coeffs(5)
    
  end subroutine calc_w4

!*****************************************************************************

end module subroutines
