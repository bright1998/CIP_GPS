c    -----------------
      subroutine init
c    -----------------
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine chck(da,damin,is,ie,js,je)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: damin
            integer :: is,ie,js,je
         end subroutine chck

         subroutine convertx(da,dam,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: da
            double precision,dimension(:,:),intent(inout) :: dam
            integer :: is,ie,js,je
         end subroutine convertx

         subroutine converty(da,dam,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: da
            double precision,dimension(:,:),intent(inout) :: dam
            integer :: is,ie,js,je
         end subroutine converty
      end interface

C OT-Vortex
      do j=1,j0
      do i=1,i0
         ro(i,j) = 25.d0/36.d0/pi
         pr(i,j) = 5.d0/12.d0/pi
         vz(i,j) = 0.d0
         bz(i,j) = 0.d0
         vx(i,j) =-Sin(2.d0*pi*y(j))
         bx(i,j) =-Sin(2.d0*pi*y(j))
         vy(i,j) = Sin(2.d0*pi*x(i))
         by(i,j) = Sin(4.d0*pi*x(i))
      enddo
      enddo

      call chck(ro,dmin,1,i0,1,j0)
      call chck(pr,pmin,1,i0,1,j0)
      call convertx(vx,vxm,1,i0-1,1,j0)
      call convertx(bx,bxm,1,i0-1,1,j0)
      call converty(vy,vym,1,i0,1,j0-1)
      call converty(by,bym,1,i0,1,j0-1)

C Calculate the space-derivative of the physical quantities
      do j=1,j0
      do i=2,i0-1
         rodx(i,j) = (-ro(i-1,j) + ro(i+1,j))/(dx(i-1) + dx(i))
         prdx(i,j) = (-pr(i-1,j) + pr(i+1,j))/(dx(i-1) + dx(i))
         vzdx(i,j) = (-vz(i-1,j) + vz(i+1,j))/(dx(i-1) + dx(i))
      enddo
      enddo

      do j=1,j0
      do i=2,i0-2
         vxdxm(i,j) = (-vxm(i-1,j) + vxm(i+1,j))/(dxm(i-1) + dxm(i))
      enddo
      enddo

      do j=1,j0-1
      do i=2,i0-1
         vydxm(i,j) = (-vym(i-1,j) + vym(i+1,j))/(dx(i-1) + dx(i))
      enddo
      enddo

      do j=2,j0-1
      do i=1,i0
         rody(i,j) = (-ro(i,j-1) + ro(i,j+1))/(dy(j-1) + dy(j))
         prdy(i,j) = (-pr(i,j-1) + pr(i,j+1))/(dy(j-1) + dy(j))
         vzdy(i,j) = (-vz(i,j-1) + vz(i,j+1))/(dy(j-1) + dy(j))
      enddo
      enddo

      do j=2,j0-1
      do i=1,i0-1
         vxdym(i,j) = (-vxm(i,j-1) + vxm(i,j+1))/(dy(j-1) + dy(j))
      enddo
      enddo

      do j=2,j0-2
      do i=1,i0
         vydym(i,j) = (-vym(i,j-1) + vym(i,j+1))/(dym(j-1) + dym(j))
      enddo
      enddo

C Calculate temperature [(gm - 1)*C_p*T]
      do j=1,j0
      do i=1,i0
         te(i,j) = gm*pr(i,j)/ro(i,j)
      enddo
      enddo

C Preserve the initial values
      do j=1,j0
      do i=1,i0
         ro0(i,j) = ro(i,j)
         vx0(i,j) = vx(i,j)
         vy0(i,j) = vy(i,j)
         vz0(i,j) = vz(i,j)
         bx0(i,j) = bx(i,j)
         by0(i,j) = by(i,j)
         bz0(i,j) = bz(i,j)
         pr0(i,j) = pr(i,j)
         te0(i,j) = te(i,j)
      enddo
      enddo

      do j=1,j0
      do i=1,i0-1
         vxm0(i,j) = vxm(i,j)
         bxm0(i,j) = bxm(i,j)
      enddo
      enddo

      do j=1,j0-1
      do i=1,i0
         vym0(i,j) = vym(i,j)
         bym0(i,j) = bym(i,j)
      enddo
      enddo

      do j=1,j0
      do i=1,i0
         rodx0(i,j) = rodx(i,j)
         prdx0(i,j) = prdx(i,j)
         vzdx0(i,j) = vzdx(i,j)
         rody0(i,j) = rody(i,j)
         prdy0(i,j) = prdy(i,j)
         vzdy0(i,j) = vzdy(i,j)
      enddo
      enddo

      do j=1,j0
      do i=1,i0-1
         vxdxm0(i,j) = vxdxm(i,j)
         vxdym0(i,j) = vxdym(i,j)
      enddo
      enddo

      do j=1,j0-1
      do i=1,i0
         vydxm0(i,j) = vydxm(i,j)
         vydym0(i,j) = vydym(i,j)
      enddo
      enddo

      call bndc
      call bndc_mag

C Preserve the initial values
      do j=1,j0
      do i=1,i0
         ro0(i,j) = ro(i,j)
         vx0(i,j) = vx(i,j)
         vy0(i,j) = vy(i,j)
         vz0(i,j) = vz(i,j)
         bx0(i,j) = bx(i,j)
         by0(i,j) = by(i,j)
         bz0(i,j) = bz(i,j)
         pr0(i,j) = pr(i,j)
         te0(i,j) = te(i,j)
      enddo
      enddo

      do j=1,j0
      do i=1,i0-1
         vxm0(i,j) = vxm(i,j)
         bxm0(i,j) = bxm(i,j)
      enddo
      enddo

      do j=1,j0-1
      do i=1,i0
         vym0(i,j) = vym(i,j)
         bym0(i,j) = bym(i,j)
      enddo
      enddo

      do j=1,j0
      do i=1,i0
         rodx0(i,j) = rodx(i,j)
         prdx0(i,j) = prdx(i,j)
         vzdx0(i,j) = vzdx(i,j)
         rody0(i,j) = rody(i,j)
         prdy0(i,j) = prdy(i,j)
         vzdy0(i,j) = vzdy(i,j)
      enddo
      enddo

      do j=1,j0
      do i=1,i0-1
         vxdxm0(i,j) = vxdxm(i,j)
         vxdym0(i,j) = vxdym(i,j)
      enddo
      enddo

      do j=1,j0-1
      do i=1,i0
         vydxm0(i,j) = vydxm(i,j)
         vydym0(i,j) = vydym(i,j)
      enddo
      enddo

      return
      end subroutine init
