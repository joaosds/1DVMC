program energies

 !B.D.: Small program to calculate the exact energies for different chain lengths for further comparison

implicit none
 integer :: N(6), i
 real*8 :: gs
 real*8, parameter ::  pi = 4.0d0*DATAN(1.0d0) 
 
 N(1) = 32
 N(2) = 64
 N(3) = 96
 N(4) = 128
 N(5) = 256
 N(6) = 512
 
 do i = 1, 6
  gs = (pi**2/(real(24,8)*real(N(i),8)))*(real(N(i),8)+5/real(N(i),8))
  write(*,*) gs
 enddo
 
   gs = (pi**2/(real(24,8)*real(1000000,8)))*(real(1000000,8)+5/real(1000000,8))
  write(*,*) gs
 
end program energies
