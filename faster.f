 c   DESCRIPTION. This is an implementation of the Ising method following
c   the algorithm developed by Metropolis et al. (1953). The actual program
c   starts at line 100. The preceding lines of code are a module used in the
c   program to generate random numbers.

      module randomnumber
      integer*8::ip,iq,is,np,nbit,ic
      parameter(ip=1279)
      parameter(iq=418)
      parameter(is=ip-iq)
      parameter(np=14)
      parameter(nbit=31)
      integer*8,dimension(1:ip)::ix
      save ic,ix
      contains

      real*8 function gaus()
      real*8::v1,v2,fac,r,gset
      integer*8::iset
      save iset,gset

      if(iset.eq.0)then
      r=2.d0
      do while(r>1.or.r.eq.0.d0)
        v1=2.d0*dran_u()-1.d0
        v2=2.d0*dran_u()-1.d0
        r=v1*v1+v2*v2
      enddo
      fac=dsqrt(-2.d0*dlog(r)/r)
      gset=v1*fac
      gaus=v2*fac
      iset=1
      else
        gaus=gset
        iset=0
      endif
      end function gaus

      subroutine dran_ini(iseed0)
      integer*8::m,np1,nn,nn1,i,j
      real*8::dseed,p,t,x
      dseed=iseed0
      do i=1,ip
        ix(i)=0
        do j=0,nbit-1
            if(rand_xx(dseed).lt.0.5d0) ix(i)=ibset(ix(i),j)
      enddo
      enddo
      ic=0
      end subroutine dran_ini

        subroutine dran_read(iunit)
      integer*8::i
      read(iunit,*)ic
        read(iunit,*)(ix(i),i=1,ip)
        end subroutine dran_read

        subroutine dran_write(iunit)
        integer*8::i
      write(iunit,*) ic
        write(iunit,*) (ix(i),i=1,ip)
        end subroutine dran_write

        integer*8 function i_dran(n)
        integer*8::i_ran,n
      ic=ic+1
        if(ic.gt.ip) ic=1
      if(ic.gt.iq)then
      ix(ic)=ieor(ix(ic),ix(ic-iq))
      else
          ix(ic)=ieor(ix(ic),ix(ic+is))
        endif
        i_ran=ix(ic)
        if(n.gt.0)i_dran=mod(i_ran,n)+1
      end function i_dran


      real*8 function dran_u()
      real*8::rmax
        parameter (rmax=2147483647.0)
      ic=ic+1
        if(ic.gt.ip) ic=1
      if(ic.gt.iq)then
      ix(ic)=ieor(ix(ic),ix(ic-iq))
      else
          ix(ic)=ieor(ix(ic),ix(ic+is))
        endif
      dran_u=dble(ix(ic))/rmax
      end function dran_u

        real*8 function rand_xx(dseed)
        real*8:: a,c,xm,rm,dseed
        parameter (xm=2.d0**32,rm=1.d0/xm,a=69069.d0,c=1.d0)
        dseed=mod(dseed*a+c,xm)
        rand_xx=dseed*rm
        end function rand_xx
      end module randomnumber

c THE ACTUAL PROGRAM STARTS HERE. (Explain variables here).

      program metropolis
      use randomnumber
      implicit none
c     Declaration of variables
      integer i,j,k
      integer lattice, steps, m, n
      integer randseed(1:8)
      real*8 temp, energy, p, ksi
      real*8 start, finish
      integer, dimension(:), allocatable :: spin
      integer, dimension(:), allocatable :: nl,nr,nu,nd

c     Set size of lattice (nodes per side), temperature & MC steps
      lattice = 500
      temp = 0.d0
      steps = 30

c     Set seed by reading current time (array of dim. 8) and multiplying every
c     number (yr*mo*day*(...)*miliseconds) to get a different seed every time
      call date_and_time(values=randseed)
      call dran_ini(abs(product(randseed)))

c     Save results onto a file for every MC step.
      open(70,file="results.txt",status="unknown")

c     Allocate memory
      allocate(spin(1:lattice**2))
      allocate(nr(1:lattice**2),nl(1:lattice**2))
      allocate(nu(1:lattice**2),nd(1:lattice**2))

c     Fill spin array with uniform +1 and -1 values
      do m=1,lattice
        do n=1,lattice
          i=(n-1)*lattice+m
          spin(i) = 1
          if(dran_u() .lt. 0.5) spin(i) = -1
          write(70,*) spin(i)
        end do
      end do

c     Save neighbours of each node just once (improves performance)
      do m=1,lattice
        do n=1,lattice
          i = (n-1)*lattice+m

          !Right neighbours
          nr(i) = m+1
          if(nr(i).eq.lattice+1) nr(i) = 1
          nr(i)=(n-1)*lattice+nr(i)

          !Left neighbours
          nl(i) = m-1
          if(nl(i).eq.0) nl(i) = lattice
          nl(i)=(n-1)*lattice+nl(i)

          !Up neighbours
          nu(i) = n+1
          if(nu(i).eq.lattice+1) nu(i) = 1
          nu(i) = (nu(i)-1)*lattice+m

          !Down neighbours
          nd(i) = n-1
          if(nd(i).eq.0) nd(i) = lattice
          nd(i) = (nd(i)-1)*lattice+m
        end do
      end do

c     Metropoli's algorithm
      call cpu_time(start)
      do i=1,steps*lattice**2
        energy = 0.d0
        !Go to a random node
        j = nint((lattice**2-1)*dran_u()+1)
        !Compute associated difference in energy
        energy=2.d0*float(spin(j)*(spin(nu(j))+spin(nd(j))+
     &   spin(nr(j))+spin(nl(j))))
        !Compute probability
        p = min(1.d0,exp(-(energy/temp)))
        ksi = dran_u()
        !Compare with uniform random value and take decision
        if(ksi .lt. p) spin(j) = -spin(j)
        !Save results every MC iteration
        if(mod(i,lattice**2) .eq. 0) write(70,*) spin
      end do
      call cpu_time(finish)

      close(70)
      write(*,*) "Results saved in file ./results.txt"
      write(*,*) "Runtime: ",finish-start,"seconds."
      end program metropolis
