c    -----------------------------------------------------------
      subroutine cipadv(cadx,cady,da,dan,dadx,dady,dadxn,dadyn,
     &u,v,is,ie,js,je)
c    -----------------------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:),intent(in) :: cadx,cady
      double precision,dimension(:,:),intent(in) :: da,dadx,dady,u,v
      double precision,dimension(:,:),intent(inout) :: dan,dadxn,dadyn
      integer :: is,ie,js,je

      do j=js,je
      do i=is,ie
         if(u(i,j) .le. 0.d0) then
            iup = i+1
            dx1 = cadx(i)
         else
            iup = i-1
            dx1 =-cadx(i-1)
         endif

         if(v(i,j) .le. 0.d0) then
            jup = j+1
            dy1 = cady(j)
         else
            jup = j-1
            dy1 =-cady(j-1)
         endif

         dx2 = dx1**2
         dx3 = dx1**3
         rdx1 = 1.d0/dx1
         rdx2 = rdx1*rdx1
         rdx3 = rdx2*rdx1

         dy2 = dy1**2
         dy3 = dy1**3
         rdy1 = 1.d0/dy1
         rdy2 = rdy1*rdy1
         rdy3 = rdy2*rdy1

         a8 = da(i,j) - da(iup,j) - da(i,jup) + da(iup,jup)
          g = (-a8 + (dadx(i,jup) - dadx(i,j))*dx1
     &      + (dady(iup,j) - dady(i,j))*dy1)*rdx1*rdy1

          a = (dadx(i,j) + dadx(iup,j))*rdx2
     &      + 2.d0*(da(i,j) - da(iup,j))*rdx3
          c = (a8 - (dadx(i,jup) - dadx(i,j))*dx1)*rdx2*rdy1
          e = 3.d0*(da(iup,j) - da(i,j))*rdx2
     &      - (2.d0*dadx(i,j) + dadx(iup,j))*rdx1

          b = (dady(i,j) + dady(i,jup))*rdy2
     &      + 2.d0*(da(i,j) - da(i,jup))*rdy3
          d = (a8 - (dady(iup,j) - dady(i,j))*dy1)*rdx1*rdy2
          f = 3.d0*(da(i,jup) - da(i,j))*rdy2
     &      - (2.d0*dady(i,j) + dady(i,jup))*rdy1

         xx =-u(i,j)*dt
         yy =-v(i,j)*dt

         dan(i,j) = da(i,j) + g*xx*yy
     &            + ((a*xx + c*yy + e)*xx + dadx(i,j))*xx
     &            + ((b*yy + d*xx + f)*yy + dady(i,j))*yy
         dadxn(i,j) = dadx(i,j)
     &              + (3.d0*a*xx + 2.d0*c*yy + 2.d0*e)*xx
     &              + (g + d*yy)*yy
         dadyn(i,j) = dady(i,j)
     &              + (3.d0*b*yy + 2.d0*d*xx + 2.d0*f)*yy
     &              + (g + c*xx)*xx

         damin = dmin1(da(i,j),da(iup,j),da(i,jup))
         damax = dmax1(da(i,j),da(iup,j),da(i,jup))
         if(dan(i,j) .lt. damin) dan(i,j) = damin
         if(dan(i,j) .gt. damax) dan(i,j) = damax
      enddo
      enddo

      return
      end subroutine cipadv
