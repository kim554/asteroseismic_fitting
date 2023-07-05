!Build a 15 parameter gridparameter file
! 1     2     3     4    5   6      7       8      9     10  11  12  13  14  15  16
!teff, mass, menv, mhe, mh, Xhe, alphenv, alphhe, alpha, h1, h2, h3, w1, w2, w3  w4
!add to that tiny1, in 16th place

program getpar

  use subroutines
  use common_values

  implicit none

  real, dimension(36) :: datain
  real, dimension(16) :: minval, maxval, step
  real, dimension(4) :: ws
  real, dimension(3) :: hs
  real :: fitness, ff, x, tiny1, mass, sum
  real :: menv
  real, parameter :: eps=1.d-4, rmax = 96.5, minw = 1.
  logical :: parameterize_w3
  
  open(unit=101,file='inputparameters')
  open(unit=102,file='gridparameters')

  call assign_data_arrays
  
  !bounds and step sizes for problem at hand
  read(101,*)
  read(101,*) minval(:)
  read(101,*)
  read(101,*) maxval(:)
  read(101,*)
  read(101,*) step(:)
  datain(1:16) = minval(:)

  !  print *, "it's alive!"

  !  print 10, datain(1:15)
  !  print 10, minval
  !  print 10, maxval
  !  print 10, step

10 format(16F9.1)

  !Run grid

  loop01: do
     !print *, "par1", datain(1), maxval(1)+eps
     if (datain(1) .gt. (maxval(1)+eps)) exit
     loop02: do
        !print *, "par2", datain(2), maxval(2)+eps
        if (datain(2) .gt. (maxval(2)+eps)) exit
        mass = datain(2)
        call  calc_oxygen_profile(mass,hs,ws)
        loop03: do
           !print *, "par3", datain(3), maxval(3)+eps
           if (datain(3) .gt. (maxval(3)+eps)) exit
           call calc_menv(datain(3),mass,menv)
           datain(3) = menv
           loop04: do
              !print *, "par4", datain(4), maxval(4)+eps
              if (datain(4) .gt. (maxval(4)+eps)) exit
              call calc_mhe(datain(4),mass,menv,x)
              datain(4) = x
              loop05: do
                 !print *, "par5", datain(5), maxval(5)+eps
                 if (datain(5) .gt. (maxval(5)+eps)) exit
                 call calc_mh(datain(5),mass,datain(4),x)
                 datain(5) = x
                 loop06: do
                    !print *, "par6", datain(6), maxval(6)+eps
                    if (datain(6) .gt. (maxval(6)+eps)) exit
                    call calc_xhebar(datain(6),mass,x)
                    datain(6) = x
                    loop07: do
                       !print *, "par7", datain(7), maxval(7)+eps
                       if (datain(7) .gt. (maxval(7)+eps)) exit
                       call calc_alph1(datain(7),mass,x)
                       datain(7) = x
                       loop08: do
                          !print *, "par8", datain(8), maxval(8)+eps
                          if (datain(8) .gt. (maxval(8)+eps)) exit
                          call calc_alph2(datain(8),mass,x)
                          datain(8) = x
                          loop09: do
                             !print *, "par9", datain(9), maxval(9)+eps
                             if (datain(9) .gt. (maxval(9)+eps)) exit
                             loop10: do
                                !print*, "par10", datain(10), maxval(10)+eps
                                if (datain(10) .gt. (maxval(10)+eps)) exit
                                if (datain(10) .lt. 0) then
                                   datain(10) = hs(1)
                                end if
                                loop11: do
                                   !print *, "par11", datain(11), maxval(11)+eps
                                   if (datain(11) .gt. (maxval(11)+eps)) exit
                                   if (datain(11) .lt. 0) then
                                      datain(11) = hs(2)
                                   end if
                                   if (datain(11) .gt. datain(10)) then
                                      datain(11) = datain(10)
                                   end if
                                   loop12: do
                                      !print *, "par12", datain(12), maxval(12)+eps
                                      if (datain(12) .gt. (maxval(12)+eps)) exit
                                      if (datain(12) .lt. 0) then
                                         datain(12) = hs(3)
                                      end if
                                      if (datain(12) .gt. datain(11)) then
                                         datain(12) = datain(11)
                                      end if
                                      loop13: do
                                         !print *, "par13", datain(13), maxval(13)+eps
                                         if (datain(13) .gt. (maxval(13)+eps))  exit
                                         if (datain(13) .lt. 0) then
                                            datain(13) = ws(1)
                                         end if
                                         loop14: do
                                            !print *, "par14", datain(14), maxval(14)+eps
                                            if (datain(14) .gt. (maxval(14)+eps)) exit
                                            if (datain(14) .lt. 0.d0) then
                                               datain(14) = ws(2)
                                            end if
                                            if (datain(14) .lt. minw) datain(14) = minw
                                            loop15: do
                                               !print *, "par15", datain(15), maxval(15)+eps
                                               if (datain(15) .gt. (maxval(15)+eps)) exit
                                               if (datain(15) .lt. 0.d0) then
                                                  parameterize_w3 = .true.
                                               end if
                                               !print *, "par16", datain(16), maxval(16)+eps
                                                  loop16: do
                                                     if (datain(16) .gt. (maxval(16)+eps)) exit
                                                     if (datain(16) .lt. 0.d0) then
                                                        datain(16) = ws(4)
                                                     end if
                                                     if (datain(16) .lt. minw) then
                                                        datain(16) = minw
                                                     end if
                                                     if (parameterize_w3) then
                                                        datain(15) = rmax - (datain(13) + datain(14) + datain(16))
                                                        !datain(15) = ws(3)
                                                     end if
                                                     sum = datain(13) + datain(14) + datain(15) + datain(16)
                                                     if (sum .gt. rmax + eps) then
                                                        print *, "Sum(w1:w4) greater than ", rmax
                                                        stop
                                                     end if
                                                     write(102,10) datain(1:16)
                                                     if (datain(16) .gt. 0.) datain(16) = datain(16) + step(16)
                                                  end do loop16
                                               if (datain(16) .gt. 0.) datain(16) = minval(16)
                                               if (datain(15) .gt. 0.) datain(15) = datain(15) + step(15) 
                                               end do loop15
                                            if (datain(15) .gt. 0.) datain(15) = minval(15)
                                            if (datain(14) .gt. 0.) datain(14) = datain(14) + step(14)
                                         end do loop14
                                         if (datain(14) .gt. 0.) datain(14) = minval(14)
                                         if (datain(13) .gt. 0.) datain(13) = datain(13) + step(13)  
                                      end do loop13
                                      if (datain(13) .gt. 0. ) datain(13) = minval(13)
                                      if (datain(12) .gt. 0.) datain(12) = datain(12) + step(12)             
                                   end do loop12
                                   if (datain(12) .gt. 0.) datain(12) = minval(12)
                                   if (datain(11) .gt. 0.) datain(11) = datain(11) + step(11)
                                end do loop11
                                if (datain(11) .gt. 0.) datain(11) = minval(11)
                                if (datain(10) .gt. 0.) datain(10) = datain(10) + step(10)   
                             end do loop10
                             if (datain(10) .gt. 0.)  datain(10) = minval(10)
                             if (datain(9) .gt. 0.) datain(9) = datain(9) + step(9)
                          end do loop09
                          if (datain(9) .gt. 0.) datain(9) = minval(9)
                          if (datain(8) .gt. 0.) datain(8) = datain(8) + step(8)     
                       end do loop08
                       if (datain(8) .gt. 0.) datain(8) = minval(8)
                       if (datain(7) .gt. 0.) datain(7) = datain(7) + step(7) 
                    end do loop07
                    if (datain(7) .gt. 0.) datain(7) = minval(7)
                    if (datain(6) .gt. 0.) datain(6) = datain(6) + step(6)  
                 end do loop06
                 if (datain(6) .gt. 0.) datain(6) = minval(6)
                 if (datain(5) .gt. 0.) datain(5) = datain(5) + step(5) 
              end do loop05
              if (datain(5) .gt. 0.) datain(5) = minval(5)
              if (datain(4) .gt. 0.) datain(4) = datain(4) + step(4)
           end do loop04
           if (datain(4) .gt. 0.) datain(4) = minval(4)
           if (datain(3) .gt. 0.) datain(3) = datain(3) + step(3)   
        end do loop03
        if (datain(3) .gt. 0.) datain(3) = minval(3)
        datain(2) = datain(2) + step(2)                                          
     end do loop02
     datain(2) = minval(2)
     datain(1) = datain(1) + step(1)                                          
  end do loop01

end program getpar
