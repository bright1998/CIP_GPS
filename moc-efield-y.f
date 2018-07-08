c    --------------------
      subroutine moceley
c    --------------------
c MOC-Interpolation for calculating Ey
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

      do j=1,j0
c Interpolation to x-direction
         do i=1,i0-1
C Define the characteristic velocity at (i+1/2,j)
            avero = 0.5d0*(ro(i,j) + ro(i+1,j))
            ca = bxm(i,j)/dsqrt(pi4*avero)
            cp(i) = vxm(i,j) + ca
            cm(i) = vxm(i,j) - ca
         enddo

         do i=1,i0
            vanv(i) = vz(i,j)
            vanb(i) = bz(i,j)
         enddo

         do i=1,i0-1
            dl(i) = dx(i)
         enddo

         call vanLeer(dl,vp,vm,vanv,cp,cm,2,i0-2)
         call vanLeer(dl,bp,bm,vanb,cp,cm,2,i0-2)

         do i=2,i0-2
            sqrop = dsqrt(pi4*ro(i,j))
            sqrom = dsqrt(pi4*ro(i+1,j))
            sqropi = 1.d0/sqrop
            sqromi = 1.d0/sqrom
            vzmx(i,j) = (vp(i)*sqrop + vm(i)*sqrom - bp(i) + bm(i))
     &                /(sqrop + sqrom)
            bzmx(i,j) = (-vp(i) + vm(i) + bp(i)*sqropi + bm(i)*sqromi)
     &                /(sqropi + sqromi)
         enddo
      enddo

      do j=1,j0
      do i=2,i0-2
         ey(i,j) = -(vzmx(i,j)*bxm(i,j) - vxm(i,j)*bzmx(i,j))
      enddo
      enddo

      return
      end subroutine moceley
