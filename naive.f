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
c
c
c
c
      program metropolis
      use randomnumber
      implicit none
c     Declaration of variables
      integer i,j,k
      integer m, n, lattice, pos, steps
      integer randseed(1:8)
      real*8 temp, energy, p, ksi
      real*8 start, finish
      integer, dimension(:,:), allocatable :: spin

c     Set size of lattice (nodes per side), temperature & MC steps
      lattice = 50
      temp = 0.d0
      steps = 1

c     Set seed by reading current time (array of dim. 8) and multiplying every
c     number (yr*mo*day*(...)*miliseconds) to get a different seed every time
      call date_and_time(values=randseed)
      call dran_ini(abs(product(randseed)))

c     Save results onto a file for every MC step.
      open(70,file="results.txt",status="unknown")
      call cpu_time(start)
c     Initialize the spin lattice with random values
      allocate(spin(1:lattice,1:lattice))
      do i=1,lattice
        do j=1,lattice
          spin(i,j) = 1
          if(dran_u() .lt. 0.5) spin(i,j) = -1
          write(70,"(I3)",advance='no') i
          write(70,"(I3)",advance='no') j
          write(70,"(I3)") spin(i,j)
        end do
      end do

c     Metropolis et. al's algorithm:
      do i=1,(steps*lattice**2)
        n = nint((lattice-1)*dran_u()+1)
        m = nint((lattice-1)*dran_u()+1)

c       Compute the difference in energy (normalized units)
        energy = 0.d0
        energy = energy + float(spin(pos(n+1,lattice),m) +
     &  spin(pos(n-1,lattice),m) + spin(n,pos(m+1,lattice)) +
     &  spin(n,pos(m-1,lattice)))
        energy = 2.d0*spin(n,m)*energy

c       Compute the value of p and decide whether to change sign or not
        p = min(1.d0,exp(-(energy/temp)))
        ksi = dran_u()
        if(ksi .lt. p) spin(n,m) = -spin(n,m)

c       Save the results every MC iteration
        if(mod(i,steps) .eq. 0) then
          do j=1,lattice
            do k=1,lattice
              write(70,"(I3)",advance='no') j
              write(70,"(I3)",advance='no') k
              write(70,"(I3)") spin(j,k)
            end do
          end do
        end if
      end do
      call cpu_time(finish)
      close(70)

      write(*,*) "Results saved in file ./results.txt"
      write(*,*) "Runtime: ",finish-start," seconds."
      end program metropolis

c     Implementation of periodic boundary conditions
      function pos(node,lattice) result(new)
        integer node, new, lattice
        new = mod(node,lattice)
        if(new .eq. 0) new = lattice
      end function pos
