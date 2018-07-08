c    --------------------
      subroutine mocelex
c    --------------------
c MOC-Interpolation for calculating Ex
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

      do i=1,i0
c Interpolation to y-direction
         do j=1,j0-1
C Define the characteristic velocity at (i,j+1/2)
            avero = 0.5d0*(ro(i,j) + ro(i,j+1))
            ca = bym(i,j)/dsqrt(pi4*avero)
            cp(j) = vym(i,j) + ca
            cm(j) = vym(i,j) - ca
         enddo

         do j=1,j0
            vanv(j) = vz(i,j)
            vanb(j) = bz(i,j)
         enddo

         do j=1,j0-1
            dl(j) = dy(j)
         enddo

         call vanLeer(dl,vp,vm,vanv,cp,cm,2,j0-2)
         call vanLeer(dl,bp,bm,vanb,cp,cm,2,j0-2)

         do j=2,j0-2
            sqrop = dsqrt(pi4*ro(i,j))
            sqrom = dsqrt(pi4*ro(i,j+1))
            sqropi = 1.d0/sqrop
            sqromi = 1.d0/sqrom
            vzmy(i,j) = (vp(j)*sqrop + vm(j)*sqrom - bp(j) + bm(j))
     &                /(sqrop + sqrom)
            bzmy(i,j) = (-vp(j) + vm(j) + bp(j)*sqropi + bm(j)*sqromi)
     &                /(sqropi + sqromi)
         enddo
      enddo

      do j=2,j0-2
      do i=1,i0
         ex(i,j) = -(vym(i,j)*bzmy(i,j) - vzmy(i,j)*bym(i,j))
      enddo
      enddo

      return
      end subroutine mocelex
