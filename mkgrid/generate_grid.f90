!Build a 15 parameter gridparameter file
! 1     2     3     4    5   6      7       8      9     10  11  12  13  14  15
!teff, mass, menv, mhe, mh, Xhe, alphenv, alphhe, alpha, h1, h2, h3, w1, w2, w3

program getpar

  implicit none

  real*8, dimension(36) :: datain
  real*8, dimension(15) :: minval, maxval, step
  real*8 :: fitness, ff
  real*8, parameter :: eps=1.d-4, rmax = 95.

  open(unit=101,file='inputparameters')
  open(unit=102,file='gridparameters')

!bounds and step sizes for problem at hand
  read(101,*)
  read(101,*) minval(:)
  read(101,*)
  read(101,*) maxval(:)
  read(101,*)
  read(101,*) step(:)
  datain(1:15) = minval(:)

  print *, "it's alive!"
 
  print 10, datain(1:15)
  print 10, minval
  print 10, maxval
  print 10, step

10 format(15F9.1)
 
!Run grid
  do
     if (datain(1) .gt. (maxval(1)+eps)) exit
     do
        if (datain(2) .gt. (maxval(2)+eps)) exit
        do
           if (datain(3) .gt. (maxval(3)+eps)) exit
           do
              if (datain(4) .gt. (maxval(4)+eps)) exit
              do
                 if (datain(5) .gt. (maxval(5)+eps)) exit
                 do
                    if (datain(6) .gt. (maxval(6)+eps)) exit
                    do
                       if (datain(7) .gt. (maxval(7)+eps)) exit
                       do
                          if (datain(8) .gt. (maxval(8)+eps)) exit
                          do
                             if (datain(9) .gt. (maxval(9)+eps)) exit
                             do
                                if (datain(10) .gt. (maxval(10)+eps)) exit
                                do
                                   if (datain(11) .gt. (maxval(11)+eps)) exit
                                   do
                                      if (datain(12) .gt. (maxval(12)+eps)) exit
                                      do
                                         if (datain(13) .gt. (maxval(13)+eps)) exit
                                         do
                                            if (datain(14) .gt. (maxval(14)+eps)) exit
                                            do
                                               if (datain(15) .ge. (maxval(15)+eps)) exit 
                                               if (((datain(4)-datain(3)) .ge. 0.0d0) .and. &
                                                    ((datain(5)-datain(4)) .ge. 200.) .and. &
                                                    (datain(13) + datain(14) + datain(15)) .lt. rmax) &
                                                    then
                                               write(102,10) datain(1:15)
                                               endif
                                               datain(15) = datain(15) + step(15)
                                            end do
                                            datain(15) = minval(15)
                                            datain(14) = datain(14) + step(14)
                                         end do
                                         datain(14) = minval(14)
                                         datain(13) = datain(13) + step(13)
                                      end do
                                      datain(13) = minval(13)
                                      datain(12) = datain(12) + step(12)
                                   end do
                                   datain(12) = minval(12)
                                   datain(11) = datain(11) + step(11)
                                end do
                                datain(11) = minval(11)
                                datain(10) = datain(10) + step(10)
                             end do
                             datain(10) = minval(10)
                             datain(9) = datain(9) + step(9)
                          end do
                          datain(9) = minval(9)
                          datain(8) = datain(8) + step(8)
                       end do
                       datain(8) = minval(8)
                       datain(7) = datain(7) + step(7)
                    end do
                    datain(7) = minval(7)
                    datain(6) = datain(6) + step(6)
                 end do
                 datain(6) = minval(6)
                 datain(5) = datain(5) + step(5)
              end do
              datain(5) = minval(5)
              datain(4) = datain(4) + step(4)
           end do
           datain(4) = minval(4)
           datain(3) = datain(3) + step(3)
        end do
        datain(3) = minval(3)
        datain(2) = datain(2) + step(2)
     end do
     datain(2) = minval(2)
     datain(1) = datain(1) + step(1)
  end do
  datain(1) = minval(1)
  
end program getpar

