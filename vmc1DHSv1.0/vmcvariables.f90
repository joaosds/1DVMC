module vmcvariables

 complex*16 :: zero = (0.0d0,0.0d0), one = (1.0d0,0.0d0), onec = (1.0d0,1.0d0), Im = (0.0d0, 1.0d0) ! complex constants
 complex*16 :: smallzero = (1e-15,1e-15) ! small zero offset
 integer :: N, L, seed, sweep, thermalization, ncor, nz, sweepbin ! loop and chain variables
 complex*16, allocatable :: psi(:,:) ! wavefunctions
 real*8, parameter ::  pi = 4.0d0*DATAN(1.0d0) ! pi in double precision * error at last digit (1, not 2)
 integer, allocatable :: chain(:), pup(:), pdown(:), nneighbors(:,:)! position/chain vectors
 complex*16, allocatable :: Aup(:,:), Adown(:,:), Aupinv(:,:), Adowninv(:,:), Aupnew(:,:), Adownnew(:,:)
 integer :: kthup, kthdown ! Column swap index
 real*8 :: psuc, pfail, ptot
 Real*8, allocatable :: zcor(:), zcormean(:)! z spin correlations
 complex*16, allocatable :: AupinvBU(:,:), AdowninvBU(:,:), updif(:), downdif(:) 
 complex*16 :: lambdaup, lambdadown, lambdatemp
 real*8 :: lambda
 real*8 :: elocal, eav1, eav2, eav3, eav4
 complex*16 :: etotal ! total energy: local + nonlocal
 
end module vmcvariables


