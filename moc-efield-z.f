c    --------------------
      subroutine mocelez
c    --------------------
c MOC-Interpolation for calculating Ez
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

      do i=1,i0-1
c Interpolation to y-direction
         do j=1,j0-1
c Define the characteristic velocity at (i+1/2,j+1/2)
            avero = 0.25d0*(ro(i,j)   + ro(i+1,j)
     &                    + ro(i,j+1) + ro(i+1,j+1))
            avevy = 0.5d0*(vym(i,j) + vym(i+1,j))
            aveby = 0.5d0*(bym(i,j) + bym(i+1,j))
            ca = aveby/dsqrt(pi4*avero)
            cp(j) = avevy + ca
            cm(j) = avevy - ca
         enddo

         do j=1,j0
            vanv(j) = vxm(i,j)
            vanb(j) = bxm(i,j)
         enddo

         do j=1,j0-1
            dl(j) = dy(j)
         enddo

         call vanLeer(dl,vp,vm,vanv,cp,cm,2,j0-2)
         call vanLeer(dl,bp,bm,vanb,cp,cm,2,j0-2)

         do j=2,j0-2
            averop = 0.5d0*(ro(i,j) +   ro(i+1,j))
            averom = 0.5d0*(ro(i,j+1) + ro(i+1,j+1))
            sqrop = dsqrt(pi4*averop)
            sqrom = dsqrt(pi4*averom)
            sqropi = 1.d0/sqrop
            sqromi = 1.d0/sqrom
            vxmy(i,j) = (vp(j)*sqrop + vm(j)*sqrom - bp(j) + bm(j))
     &                /(sqrop + sqrom)
            bxmy(i,j) = (-vp(j) + vm(j) + bp(j)*sqropi + bm(j)*sqromi)
     &                /(sqropi + sqromi)
         enddo
      enddo

      do j=1,j0-1
c Interpolation to x-direction
         do i=1,i0-1
c Define the characteristic velocity at (i+1/2,j+1/2)
            avero = 0.25d0*(ro(i,j)   + ro(i+1,j)
     &                    + ro(i,j+1) + ro(i+1,j+1))
            avevx = 0.5d0*(vxm(i,j) + vxm(i,j+1))
            avebx = 0.5d0*(bxm(i,j) + bxm(i,j+1))
            ca = avebx/dsqrt(pi4*avero)
            cp(i) = avevx + ca
            cm(i) = avevx - ca
         enddo

         do i=1,i0
            vanv(i) = vym(i,j)
            vanb(i) = bym(i,j)
         enddo

         do i=1,i0-1
            dl(i) = dx(i)
         enddo

         call vanLeer(dl,vp,vm,vanv,cp,cm,2,i0-2)
         call vanLeer(dl,bp,bm,vanb,cp,cm,2,i0-2)

         do i=2,i0-2
            averop = 0.5d0*(ro(i,j) +   ro(i,j+1))
            averom = 0.5d0*(ro(i+1,j) + ro(i+1,j+1))
            sqrop = dsqrt(pi4*averop)
            sqrom = dsqrt(pi4*averom)
            sqropi = 1.d0/sqrop
            sqromi = 1.d0/sqrom
            vymx(i,j) = (vp(i)*sqrop + vm(i)*sqrom - bp(i) + bm(i))
     &                /(sqrop + sqrom)
            bymx(i,j) = (-vp(i) + vm(i) + bp(i)*sqropi + bm(i)*sqromi)
     &                /(sqropi + sqromi)
         enddo
      enddo

      do j=2,j0-2
      do i=2,i0-2
         ez(i,j) = -(vxmy(i,j)*bymx(i,j) - vymx(i,j)*bxmy(i,j))
      enddo
      enddo

      return
      end subroutine mocelez
