c    -----------------
      subroutine mocx
c    -----------------
c MOC-Interpolation for updating Vx
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

      do i=2,i0-2
c Interpolate to y-direction
         do j=2,j0-2
c Define the characteristic velocity at (i+1/2,j+1/2)
            avero = 0.25d0*(roh(i,j)   + roh(i+1,j)
     &                    + roh(i,j+1) + roh(i+1,j+1))
            aveby = 0.5d0*(bym(i,j) + bym(i+1,j))
            cp(j) = aveby/dsqrt(pi4*avero)
            cm(j) =-cp(j)
         enddo

         do j=2,j0-1
            vanv(j) = vxmh(i,j)
            vanb(j) = bxm(i,j)
         enddo

         do j=1,j0-1
            dl(j) = dy(j)
         enddo

         call vanLeer(dl,vp,vm,vanv,cp,cm,3,j0-3)
         call vanLeer(dl,bp,bm,vanb,cp,cm,3,j0-3)

         do j=3,j0-3
            averop = 0.5d0*(roh(i,j) +   roh(i+1,j))
            averom = 0.5d0*(roh(i,j+1) + roh(i+1,j+1))
            sqropi = 1.d0/dsqrt(pi4*averop)
            sqromi = 1.d0/dsqrt(pi4*averom)
            bxmy(i,j) = (-vp(j) + vm(j) + bp(j)*sqropi + bm(j)*sqromi)
     &                /(sqropi + sqromi)
         enddo
      enddo

      return
      end subroutine mocx
