integer function ipos(cm)

 ! B.D.: Generates number between 1 and cm
 
 implicit none
 integer :: intl,intl2,cm
 real*8 :: r,rkiss05

 r = rkiss05()*cm
 intl = int(r) + 1
 
 if (intl .gt. cm) then 
  intl2 = cm
  ipos = intl2
 else
 ipos = intl
 endif
!
end function ipos

function rkiss05()
 implicit none
 integer,parameter      :: r8b = SELECTED_REAL_KIND(P = 14,R = 99)   ! 8-byte reals
 integer,parameter      :: i4b = SELECTED_INT_KIND(8)            ! 4-byte integers 
 real(r8b),parameter    :: am = 4.656612873077392578d-10       ! multiplier 1/2^31
 real(r8b)              :: rkiss05  
 integer(i4b)           :: kiss
 integer(i4b)           :: x,y,z,w              ! working variables for the four generators
 common /kisscom/x,y,z,w 

 x = 69069 * x + 1327217885
 y = ieor (y, ishft (y, 13)); y = ieor (y, ishft (y, -17)); y = ieor (y, ishft (y, 5))
 z = 18000 * iand (z, 65535) + ishft (z, - 16)
 w = 30903 * iand (w, 65535) + ishft (w, - 16)
 kiss = ishft(x + y + ishft (z, 16) + w , -1)
 rkiss05 = kiss*am
end function rkiss05

SUBROUTINE kissinit(iinit)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 

      integer(i4b) idum,ia,im,iq,ir,iinit
      integer(i4b) k,x,y,z,w,c1
      real(r8b)    rkiss05,rdum
      parameter (ia=16807,im=2147483647,iq=127773,ir=2836)
      common /kisscom/x,y,z,w

      !!! Test integer representation !!!
      c1=-8
      c1=ishftc(c1,-3)
!     print *,c1
      if (c1.ne.536870911) then
         print *,'Nonstandard integer representation. Stoped.'
         stop
      endif

      idum=iinit
      idum= abs(1099087573 * idum)               ! 32-bit LCG to shuffle seeds
      if (idum.eq.0) idum=1
      if (idum.ge.IM) idum=IM-1

      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         x=idum+1 
      else 
         x=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then 
         y=idum+1 
      else 
         y=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         z=idum+1 
      else 
         z=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         w=idum+1 
      else 
         w=idum
      endif

      rdum=rkiss05()
      
      return
end subroutine kissinit

