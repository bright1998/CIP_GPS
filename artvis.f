c    --------------------------------------------------------------
      subroutine artvis(DIV,q,qx,qy,qdivv,xvisl,yvisl,is,ie,js,je)
c    --------------------------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(out) :: DIV,q,qx,qy,qdivv,
     &xvisl,yvisl
      integer :: is,ie,js,je

      r3 = 1.d0/3.d0

      do j=js,je
      do i=is,ie
         q(i,j) = 0.d0
         qx(i,j) = 0.d0
         qy(i,j) = 0.d0
         qdivv(i,j) = 0.d0
         xvisl(i,j) = 0.d0
         yvisl(i,j) = 0.d0
         Cs = dsqrt(gm*prh(i,j)/roh(i,j))
         delvx =-vxmh(i-1,j) + vxmh(i,j)
         delvy =-vymh(i,j-1) + vymh(i,j)
         DIV(i,j) = delvx*dxmi(i-1) + delvy*dymi(j-1)
C Ogata & Yabe 1999 Computer Physics Communications 119 179
         if(DIV(i,j) .lt. 0.d0) then
            delv = DIV(i,j)*dmin1(dxm(i-1),dym(j-1))
            q(i,j) = fv*(-Cs + 0.5d0*(gm + 1.d0)*delv)*roh(i,j)*delv

C Stone & Norman
C Tensor Term
            dd = dmax1(dxm(i-1),dym(j-1))
            vis = fvt*dd**2*roh(i,j)*DIV(i,j)
            qx(i,j) = vis*(delvx*dxmi(i-1) - DIV(i,j)*r3)
            qy(i,j) = vis*(delvy*dymi(j-1) - DIV(i,j)*r3)
            qdivv(i,j) = vis*r3*((delvx*dxmi(i-1) - delvy*dymi(j-1))**2
     &                         + (delvy*dymi(j-1) - 0.d0           )**2
     &                         + (0.d0            - delvx*dxmi(i-1))**2)
         endif

C Linear Term
         if(delvx .lt. 0.d0) then
            xvisl(i,j) =-fvl*roh(i,j)*delvx*Cs
         endif
         if(delvy .lt. 0.d0) then
            yvisl(i,j) =-fvl*roh(i,j)*delvy*Cs
         endif
      enddo
      enddo

      return
      end subroutine artvis
