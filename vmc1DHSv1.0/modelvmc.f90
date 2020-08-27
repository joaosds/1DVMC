
subroutine initialconfig
 use vmcvariables
 implicit none
 integer :: cond1, cond2, ir , ipos, ind
 integer :: t, u, i, j, tempt
 
 chain = 0
 pup = 0
 pdown = 0

 
 ! ----------------- Fix chain in half filling -----------------------!

  do i = 1,N ! Loop for spin up in half-filling
  cond1 = 0
   do while (cond1 .eq. 0)
    ir = ipos(L) !sorteia um número  aleatório entre 1 e L
    cond2 = 0   !reseta o valor em cada loop 
    if (chain(ir) .eq. 1) then  ! occupied state
     cond2 = 1  ! breaks while loop for half-filling condition
    endif          
    if (cond2 .eq. 0) then 
     pup(i) = ir           ! Position of spin up
     chain(ir) = 1         ! Spin up at ir
     cond1 = 1             ! Necessary for half-filling condition
    endif            
   enddo
 enddo
 
  call sort(n,pup) ! from pup construction it's out of order
  
 j = 1
 do i = 1,L
   !it's automatically organized from line 49
  if (chain(i) .eq. 0) then !all other positions are up
   pdown(j) = i  !position of spin
   j = j+1
  endif
 enddo

 ! --------- Assign values of psi to matrix configuration ------------ !

 t = 0 ! at each column of A values will be assigned based on the W.F.
 u = 0 ! so that it's ordered by the position of spins up (or down)
     
 ! Run through the chain and check the spin. Get the position and assign
 ! the column of matrices A() with values of Psi     


 do i = 1, N
  do j = 1, N
  
  ind = pup(j)
  Aup(i,j) = psi(i,ind)
  
  ind = pdown(j)
  Adown(i,j) = psi(i,ind)
  enddo
 enddo 
 
 
 ! ---------------------------- Definition of neighbors vector --------!
 
 do j = 1, nz ! 2nd index is for the 2nd, 3rd... neighbor;
  do i = 1, L
   tempt = i+j
   if (tempt .gt. L) then
   tempt = tempt - L
   nneighbors(i,j) = tempt
   else
   nneighbors(i,j) = tempt
   endif 
  enddo
 enddo
 
 return  
end subroutine initialconfig



subroutine spincorrelations(zcor2,ncor2)
 use vmcvariables
 implicit none
 integer :: ncor2, i, j, t
 real*8 :: zcor2(ncor)

 zcor2 = 0.0d0
 do i = 1, N !half-filling and periodic boundary conditions! N is the maximum |i-j|
  t = i
  do j = 1, ncor2 ! Percorre até distância máxima de L/2 -1
   t = t+1 ! right neighboor
   zcor2(j)= zcor2(j) +0.250d0*(2.0d0*real(chain(i),8)-1.0d0)*(2.0d0*real(chain(t),8)-1.0d0)
  enddo
 enddo
 
 zcor2 = zcor2/real(N,8)
 
 return
end subroutine

subroutine energyconfiguration(elocal3,elocal2)
 use vmcvariables
 implicit none
 integer ::  i, j, t, kthupnew, kthdownnew, chain2(L)
 real*8 :: elocal1, elocal3, dist
 complex*16 :: elocal2, Lambdaq
 
 ! It's necessary to build a new chain since exchanges will change it
 
 do i = 1, N ! Not L, we assign two values to chain2 at each iteration!
  t = pup(i)
  chain2(t) = i
  t = pdown(i)
  chain2(t) = i
 enddo
 !chain2 is necessary to determine the given position of each spin
 
 ! --------------------------- Local Contribution ---------------------!

 elocal2 = zero
 elocal3 = 0.0d0
 
 do i = 1, L
  do j = 1, nz
   
   ! Local energy !
    
   t = nneighbors(i,j) ! right neighboor
   elocal1 = 0.250d0*(2.0d0*real(chain(i),8)-1.0d0)*(2.0d0*real(chain(t),8)-1.0d0)
   dist = i-t
   dist = real(dist,8)
   dist = dist*Pi/real(L,8)
   !dist = 2*dist**2
   dist = 2*sin(dist)**2
   elocal1 = elocal1/dist
   elocal3 = elocal3+elocal1
   elocal2 = elocal2 + cmplx(elocal1, 0.0d0, 8)
   ! Non Local energy !
   
   ! Only hopping between neighbors with different spin is permitted;
  
   if (chain(i) .ne. chain(t)) then ! different spin
    if (chain(i) .eq. 1) then
     kthupnew = chain2(i)
     kthdownnew = chain2(t)
    else
     kthupnew = chain2(t)
     kthdownnew = chain2(i)
   endif
   
   AupinvBU = Aupinv      ! backup
   AdowninvBU = Adowninv  ! backup
   
   Aupnew = Aup     ! new configuration
   Adownnew = Adown ! new configuration 
   
   call columnswap(kthupnew,kthdownnew)
   
   call ShermanLemma(N, kthupnew, lambdaup, AupinvBU , updif) 
   call ShermanLemma(N, kthdownnew, lambdadown, AdowninvBU , downdif)
  
  Lambdaq = lambdaup*lambdadown
  elocal2 = elocal2 -  cmplx(0.50d0,0.0d0,8)*Lambdaq/cmplx(dist, 0.0d0, 8)
  ! After doing one hopping we must return to the initial configuration;(line 560-564)
  endif
  enddo
 enddo
 
 elocal3 = elocal3/real(L,8)
 elocal2 = elocal2/real(L,8)
 elocal3 = elocal3*(Pi/L)**2
 elocal2 = elocal2*(Pi/L)**2
 return
end subroutine energyconfiguration


