c    -----------------
      subroutine mocz
c    -----------------
c MOC-Interpolation for updating Vz
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine vanLeer(dl,qp,qm,q,cp,cm,is,ie)
            double precision,dimension(:),intent(out) :: qp,qm
            double precision,dimension(:),intent(in) :: dl,q,cp,cm
            integer :: is,ie
         end subroutine vanLeer
      end interface

      dimension :: cp(1:ijmax),cm(1:ijmax),vanv(1:ijmax),vanb(1:ijmax),
     &             dl(1:ijmax)
      dimension :: vp(1:ijmax),vm(1:ijmax),bp(1:ijmax),bm(1:ijmax)

      do j=2,j0-1
c Interpolation to x-direction
         do i=2,i0-2
c Define the characteristic velocity at (i+1/2,j)
            avero = 0.5d0*(roh(i,j) + roh(i+1,j))
            cp(i) = bxm(i,j)/dsqrt(pi4*avero)
            cm(i) =-cp(i)
         enddo

         do i=2,i0-1
            vanv(i) = vzh(i,j)
            vanb(i) = bz(i,j)
         enddo

         do i=1,i0-1
            dl(i) = dx(i)
         enddo

         call vanLeer(dl,vp,vm,vanv,cp,cm,3,i0-3)
         call vanLeer(dl,bp,bm,vanb,cp,cm,3,i0-3)

         do i=3,i0-3
            averop = roh(i,j)
            averom = roh(i+1,j)
            sqropi = 1.d0/dsqrt(pi4*averop)
            sqromi = 1.d0/dsqrt(pi4*averom)
            bzmx(i,j) = (-vp(i) + vm(i) + bp(i)*sqropi + bm(i)*sqromi)
     &                /(sqropi + sqromi)
         enddo
      enddo

      do i=2,i0-1
c Interpolation to y-direction
         do j=2,j0-2
c Define the characteristic velocity at (i,j+1/2)
            avero = 0.5d0*(roh(i,j) + roh(i,j+1))
            cp(j) = bym(i,j)/dsqrt(pi4*avero)
            cm(j) =-cp(j)
         enddo

         do j=2,j0-1
            vanv(j) = vzh(i,j)
            vanb(j) = bz(i,j)
         enddo

         do j=1,j0-1
            dl(j) = dy(j)
         enddo

         call vanLeer(dl,vp,vm,vanv,cp,cm,3,j0-3)
         call vanLeer(dl,bp,bm,vanb,cp,cm,3,j0-3)

         do j=3,j0-3
            averop = roh(i,j)
            averom = roh(i,j+1)
            sqropi = 1.d0/dsqrt(pi4*averop)
            sqromi = 1.d0/dsqrt(pi4*averom)
            bzmy(i,j) = (-vp(j) + vm(j) + bp(j)*sqropi + bm(j)*sqromi)
     &                /(sqropi + sqromi)
         enddo
      enddo

      return
      end subroutine mocz
