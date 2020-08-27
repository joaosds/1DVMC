! ---------------------- Variational Monte Carlo 1D --------------------
!
! Author: João Augusto Sobral da Silva / Github: @joaosds
!
! The code is commented, and each subroutine has a brief description (B.D.)
! of its functionality. For any discussions you can contact me at Github!
!
! History:
!
!  v1.0 (08/27/20) - First implementation 
!
! ---------------------- Informations ---------------------------------
!
! File 'random.f90' - All codes responsible for random generating numbers
! and related topics/ we use the implementation 'rkiss05.f90' by Thomas Vojta,
! vojta@mst.edu
!
! File 'vmcvariables.f90' - Global variables
!
! File 'modelvmc.f90' - retains the lattice definition, energy and spin-spin
! calculation and initial wavefunction.
! 
! File 'inputpar.dat': Input parameters as following:
! N, seed -> Chain Length and random seed;
! nz, sweep -> Number of next neighbors considered in the model / values
! used after the thermalization for the mean calculation;
! thermalization -> steps for themalization of the MC algorithm;
! Nbin - > Bin number for calculation of the means.
!
! File 'meanene.dat': saves the energy obtained from nonlocal+local method
! and only local methods (best suitable for large N).
!
! File 'spincor.dat': saves the z component of the spin-spin correlations;
!
! Note: One can change line 116 and 117 to inversesvd if the input matrices
! are ill-conditioned. This may improve a little the accuracy of the algorithm.
!
! ----------------------------------------------------------------------

 include 'vmcvariables.f90'
 include 'random.f90'
 include 'modelvmc.f90'

program VMC

 ! B.D.: Main structure of the algorithm. Here the pseudo-code of my
 ! master thesis is fully implemented, with some adjustments for mean
 ! values and the standard deviation;
 
 use vmcvariables
 implicit none
 integer :: lt, Nbin, pq
 real*8, allocatable :: meancor(:), meancors(:), meancorsigma(:)
 real*8 :: enfinalsigma,enlfinalsigma
 real*8 :: meanen, meanenl, meanens, meanenls
 
! ------------------------ Initial Conditions -------------------------!

 open(unit = 3,File = "inputpar.dat", Status = 'OLD')
 
 read(3,*) N, seed
 read(3,*) nz, sweep
 read(3,*) thermalization
 read(3,*) Nbin
 close(3)
 
 L = N 
 N = N/2 !!!!!!! Attention !!!!!!!!!!!
 ncor = N
  
 call kissinit(seed) ! Start ran subroutine
 
 allocate(AupinvBU(N,N),AdowninvBU(N,N))
 allocate(Aupinv(N,N), Adowninv(N,N), updif(N), downdif(N))
 allocate(psi(N,L), Aup(N,N), Adown(N,N), Aupnew(N,N), Adownnew(N,N))
 allocate(chain(L), pup(N), pdown(N),nneighbors(L,nz))
 allocate(zcor(ncor),zcormean(ncor))
 allocate(meancor(ncor))
 allocate(meancors(ncor), meancorsigma(ncor))
 
 meanen = 0.0d0
 meanenl = 0.0d0
 meancor = 0.0d0
 meanens = 0.0d0
 meanenls = 0.0d0
 meancors = 0.0d0
 
 zcormean = 0.0d0
 nneighbors = 0
 zcor = 0.0d0
 
 psi = zero
 Aup = zero
 Adown = zero
 Aupinv = zero
 Adowninv = zero
 Aupnew = zero
 Adownnew = zero 
 updif = zero
 downdif = zero
 Aupinv = zero
 Adowninv = zero
 AupinvBU = zero
 AdowninvBU = zero
 psuc = 0.0d0
 pfail = 0.0d0
 
! ---------------------------------------------------------------------! 
 
 call wavefunction
 call initialconfig
 
 deallocate(psi)
! -------------- Backup of initial configurations -------------------- !

  Aupinv = Aup
  Adowninv = Adown
  
  call inverse(Aupinv, N) ! Aupinv and Adowninv are the inverse matrices
  call inverse(Adowninv, N)

 do lt = 1, thermalization
  call mainvmc
 enddo
 
sweepbin = sweep/Nbin

meanen = 0.0d0
meanenl = 0.0d0

do pq = 1, Nbin

 eav1 = 0.0d0
 eav2 = 0.0d0
 eav3 = 0.0d0
 eav4 = 0.0d0
 zcormean = 0.0d0
 
 do lt = 1, sweepbin
  call mainvmc
  
  ! ---------------------- Spin-Spin correlations ---------------------!
  
  call spincorrelations(zcor,ncor)
  zcormean = zcormean + zcor
  
  ! -------------------------- Energy Configuration -------------------!
  call energyconfiguration(elocal,etotal)
  eav1 = eav1 + real(etotal,8)
  eav3 = eav3 + elocal
  
 enddo
 
 zcormean = zcormean/real(sweepbin,8)
 eav3 = eav3/real(sweepbin,8)
 eav1 = eav1/real(sweepbin,8)

 meancor = meancor + zcormean
 meancors = meancors + zcormean**2
 
 meanen = meanen + eav1
 meanenl = meanenl + eav3
 meanens = meanens + eav1**2
 meanenls = meanenls + eav3**2
 
enddo
! -------------------------- Local +  Nonlocal energy -----------------!

 meanen = meanen/Nbin   !final energy
 meanenl = meanenl/Nbin !final energy
 
 meanens = meanens/Nbin
 meanenls = meanenls/Nbin
 
 enfinalsigma = (meanens-meanen**2)/(Nbin-1)
 enfinalsigma = dsqrt(enfinalsigma)
 
 enlfinalsigma = (meanenls-meanenl**2)/(Nbin-1)
 enlfinalsigma = dsqrt(enlfinalsigma)
 
 ! ------------------- Mean and Standard Deviation --------------------!
 
 meancor = meancor/Nbin
 meancors = meancors/Nbin
 meancorsigma = (meancors-meancor**2)/(Nbin-1)
 meancorsigma = dsqrt(meancorsigma)
 
 open(unit = 10, file = "spincor.dat", status = 'old')
 open(unit = 40, file = "meanene.dat", status = 'old')
 
 do lt = 1, ncor
  write(10,*) lt, meancor(lt), meancorsigma(lt)
 enddo
 
 write(40,*) meanen, enfinalsigma
 write(40,*) 3.0d0*meanenl, enlfinalsigma
 
 close(10)
 close(40)

 ! -------------- Safety measure: Deallocation and Files Closure ------!
 
 deallocate(meancor, meancors, meancorsigma)
 deallocate(AupinvBU,AdowninvBU)
 deallocate(Aupinv, Adowninv, updif, downdif)
 deallocate(Aup, Adown, Aupnew, Adownnew)
 deallocate(chain, pup, pdown, nneighbors)
 deallocate(zcor, zcormean)

end program vmc

subroutine mainvmc

 ! B.D.: Set a backup of the initial configuration and swap one column
 ! of the old one. Use the Sherman routine to get both the new inverse
 ! and determinat. At the end, call the metropolis algorithm.
 
 use vmcvariables
 implicit none
 integer :: pq, ipos
 
 pq = 1
 do while (pq .le. N) 
  pq = pq+1
  AupinvBU = Aupinv      ! backup
  AdowninvBU = Adowninv  ! backup
  
  Aupnew = Aup     ! new configuration
  Adownnew = Adown ! new configuration 
  
  kthup = ipos(N) ! sort a column of ups, the position will be pup(kthup)
  kthdown = ipos(N) ! position pdown(kthdown)
 
  call columnswap(kthup, kthdown) ! Aupnew an Adownnew are the newconfigurations
 
  call ShermanLemma(N, kthup, lambdaup, AupinvBU , updif) 
  call ShermanLemma(N, kthdown, lambdadown, AdowninvBU , downdif)
  ! ---------------------------- Sherman Lemma ------------------------ !
  Lambdatemp = lambdaup*lambdadown
  Lambda = real(lambdatemp*conjg(lambdatemp),8)

  call metropolis
 enddo
  
 return
endsubroutine mainvmc

subroutine wavefunction

 ! B.F.: Calculates the initial configuration for the Fermi Sea of
 ! the spinons ansatz (with the gutzwiller projection to half-filling).
 
 use vmcvariables
 implicit none
 integer :: i, j
 real*8 :: k

 do i = 1,N !line N = N/2
  k = (2.0d0*real(i,8)-1.0d0)*pi/real(L,8)-pi/2.0d0
  do j = 1,L ! column, runs through chain
   Psi(i,j) = cos(k*real(j,8))+Im*sin(k*real(j,8))
  enddo
 enddo
 
 return
end subroutine wavefunction

Subroutine sort(nn,a)

! B.D.: Sort a vector a(n) in increasing order

 Implicit None
 Integer*4 :: nn,i,j
 Integer*4 :: a(nn), temp

 Do i=1,nn-1
  Do j=i+1,nn
   If (a(i) .gt. a(j)) Then
    temp = a(i)
    a(i) = a(j)
    a(j) = temp
   End If
  End Do
 End Do

 Return
 
End

subroutine inverse(matrix, na)

 ! B.F.: Calculates the inverse matrix from the LU decomposition, using
 ! the routines ZGETRF and ZGETRI of LAPACK
 
 implicit none
 integer :: na
 complex*16 :: matrix(na,na)
 complex*16, allocatable :: work(:)
 integer :: ipiv(na)
 integer ::  info, lwork
 complex*16 :: zero = (0.0d0,0.0d0)
 
 external ZGETRF, ZGETRI
 
 ipiv = 0
 info = 0
 lwork = 33*2*na
 
 allocate(work(lwork))
 work = zero
 
  ! --------- LU factorization of matrix Aup and Adown ------------------!

  call ZGETRF(na, na, matrix, na, ipiv, info)

! ------------------ Inverse of both matrices -------------------------!
! The inverses are stored at the input matrices

 call ZGETRI(na, matrix, na, ipiv, WORK, lwork, info)

 return
end subroutine inverse

subroutine inversesvd (matrix, na)

 ! B.F.: Calculates the pseudo-inverse matrix from the SVD decomposition, using
 ! the routines ZGESVD, ZGETRF, ZGETRI of LAPACK
 
 implicit none
 integer :: na, i
 complex*16, allocatable :: work(:)
 Real*8 :: work2(2*na), work3(2*na), conditional
 integer :: ipiv(na), info
 complex*16 :: Aupinv2(na,na), Aupinv3(na,na)
 complex*16 :: matrix(na,na)
 complex*16 :: zero = (0.0d0,0.0d0)
 
 !ZGESVD variables
  Real*8 ::  S(na), rworksvd(5*na) ! s is sorted
  Complex*16 :: U(na,na), VT(na,na), UT(na,na), V(na,na)
  integer :: lworksvd
  Real*8 :: Spseudo(na,na)

 rworksvd = 0.0d0 
 Spseudo = 0.0d0
 VT = zero
 U = zero
 UT = zero
 V = zero
 S = 0.0d0
 ipiv = 0
 info = 0
 conditional = 0.d0
 
 lworksvd = 66*na
 
 allocate(work(lworksvd))
 
 work = zero
 work2 = 0.d0
 work3 = 0.d0

 Aupinv2 = matrix
! ------------------------------------ SVD decomposition --------------!

 Call ZGESVD ('A','A',na,na,Aupinv2, na, S, U, na, VT,na, work, lworksvd,rworksvd,info)

 Spseudo = 0.0d0
 do i = 1, na
  if (S(i) .LT. 1e-14) THEN
  Spseudo(i,i) = 0.0d0
  else
  Spseudo(i,i) = 1.0d0/S(i)
  endif
 enddo
 
 Aupinv2 = zero
 Aupinv3 = zero
 
 call ZGETRF(na, na, U, na, ipiv, info)
 call ZGETRI(na, U, na, ipiv, WORK, lworksvd, info)
 UT = U
 ipiv = 0
 call ZGETRF(na, na, VT, na, ipiv, info)
 call ZGETRI(na, VT, na, ipiv, WORK, lworksvd, info)
 V = VT
 Aupinv2 = matmul(V,Spseudo)
 Aupinv3 = matmul(Aupinv2, UT)
 
 matrix = Aupinv3
 return
end subroutine inversesvd

subroutine columnswap(kthup1,kthdown1)

 ! B.D.: Randomly changes the column of the up and down matrices of
 ! configuration.
 
 use vmcvariables
 implicit none
 integer :: j, kthup1, kthdown1
 complex*16 :: temp1, temp2
 
 temp1 = zero
 temp2 = zero


 do j = 1, N
  temp1 = Aupnew(j, kthup1)
  temp2 = Adownnew(j, kthdown1)
  Aupnew(j,kthup1) = temp2
  Adownnew(j,kthdown1) = temp1
  updif(j) = Aupnew(j,kthup1) - temp1
  downdif(j) = Adownnew(j,kthdown1) - temp2
 end do

 return
end subroutine columnswap

subroutine ShermanLemma(na, fkth, lambdaq, matrixinv, u)

 ! B.D.: Implements the Sherman-Morrison formula and the Matrix 
 ! determinant Lemma to calculate both the inverse matrix and the
 ! determinant ratio with the new configuration.
 
 use vmcvariables
 implicit none
 integer :: na, fkth, i ,j
 complex*16 :: w(na), z(na), lambdaq, matrixinv(na,na), u(na), temp1

 w = zero
 z = zero
 lambdaq = zero
 
 ! 2nd index = column of configuration matrix
 
 ! z = Ainv*u
 
 do i = 1, na !column
 temp1 = zero
  do j = 1, na
  temp1 = temp1 + matrixinv(i,j)*u(j)
  enddo
 z(i) = temp1 ! column vector
 enddo

 do j = 1, na
  w(j) = matrixinv(fkth,j) ! row vector
 enddo
 
 Lambdaq = z(fkth)
 
 ! ----------------------- Update of Inverse and det ----------------- !
 Lambdaq = one + Lambdaq + smallzero
 
 do i = 1, na
  do j = 1, na ! column
   matrixinv(i,j) =  matrixinv(i,j) -z(i)*w(j)/Lambdaq
  enddo
 enddo
 
 return
end subroutine ShermanLemma

subroutine metropolis

 ! B.D.: Metropolis-Hasting algorithm implementation. We can also save
 ! the chances of sucess when swapping the columns (psuc) and failure
 ! (pfail) to check at the end.
 
 use vmcvariables
 implicit none
 integer :: tempint
 real*8 :: w, rkiss05
 
 w = rkiss05()

  if (lambda .ge. w) then ! Accept the new configuration
  Aupinv = AupinvBU
  Adowninv = AdowninvBU
  Aup = Aupnew
  Adown = Adownnew
  
    ! Swap position vectors and update chain positions
    
  tempint = pup(kthup)
  pup(kthup) = pdown(kthdown)
  pdown(kthdown) = tempint

  chain(pup(kthup)) = 1 !up
  chain(pdown(kthdown)) = 0 !down
  psuc = psuc + 1.0d0
    
  else  !Mantain old configuration
   AupinvBU = zero
   AdowninvBU = zero
   Aupnew = zero
   Adownnew = zero
   pfail = pfail +1.0d0
  endif

 return
end subroutine metropolis

