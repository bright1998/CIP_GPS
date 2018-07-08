c    -----------------------
      subroutine normal_cip
c    -----------------------
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine artvis(DIV,q,qx,qy,qdivv,xvisl,yvisl,is,ie,js,je)
            double precision,dimension(:,:),intent(out) :: DIV,q,qx,qy,
     &qdivv,xvisl,yvisl
            integer :: is,ie,js,je
         end subroutine artvis
      end interface

      dimension :: DIV(1:i0,1:j0),DIVh(1:i0,1:j0)
      dimension :: q(1:i0,1:j0),qx(1:i0,1:j0),qy(1:i0,1:j0),
     &             qdivv(1:i0,1:j0),xvisl(1:i0,1:j0),yvisl(1:i0,1:j0)

C Artificial Viscosity
      call artvis(DIV,q,qx,qy,qdivv,xvisl,yvisl,3,i0-2,3,j0-2)

C MOC-Interpolation
      call mocx
      call mocy
      call mocz

C Vx
      do j=4,j0-3
      do i=3,i0-3
         aveby = 0.5d0*(by(i,j) + by(i+1,j))
         avebz = 0.5d0*(bz(i,j) + bz(i+1,j))
         averoi = 2.d0/(roh(i,j) + roh(i+1,j))
         vxmn(i,j) = vxmh(i,j) + dt*(
c gas pressure gradient term
     &-(-(prh(i,j)   + q(i,j)   + qx(i,j)   + xvisl(i,j))
     & + (prh(i+1,j) + q(i+1,j) + qx(i+1,j) + xvisl(i+1,j))
c magnetic pressure term
     & + pi4i*(aveby*(-by(i,j) + by(i+1,j))
     &       + avebz*(-bz(i,j) + bz(i+1,j))))*averoi*dxi(i)
c magnetic tension term
     & + pi4i*averoi*aveby*(-bxmy(i,j-1) + bxmy(i,j))*dymi(j-1)
c external force
     & + Fexxm(i,j))
      enddo
      enddo

C Vy
      do j=3,j0-3
      do i=4,i0-3
         avebx = 0.5d0*(bx(i,j) + bx(i,j+1))
         avebz = 0.5d0*(bz(i,j) + bz(i,j+1))
         averoi = 2.d0/(roh(i,j) + roh(i,j+1))
         vymn(i,j) = vymh(i,j) + dt*(
c gas pressure gradient term
     &-(-(prh(i,j)   + q(i,j)   + qy(i,j)   + yvisl(i,j))
     & + (prh(i,j+1) + q(i,j+1) + qy(i,j+1) + yvisl(i,j+1))
c magnetic pressure term
     & + pi4i*(avebx*(-bx(i,j) + bx(i,j+1))
     &       + avebz*(-bz(i,j) + bz(i,j+1))))*averoi*dyi(j)
c magnetic tension force
     & + pi4i*averoi*avebx*(-bymx(i-1,j) + bymx(i,j))*dxmi(i-1)
c external force
     & + Fexym(i,j))
      enddo
      enddo

C Vz
      do j=4,j0-3
      do i=4,i0-3
         vzn(i,j) = vzh(i,j) + dt*(
c magnetic tension force
     & pi4i/roh(i,j)*(bx(i,j)*(-bzmx(i-1,j) + bzmx(i,j))*dxmi(i-1)
     &              + by(i,j)*(-bzmy(i,j-1) + bzmy(i,j))*dymi(j-1))
     & + Fexz(i,j))
      enddo
      enddo

C Density & Pressure
      do j=4,j0-3
      do i=4,i0-3
         delvx =-vxmn(i-1,j) + vxmn(i,j)
         delvy =-vymn(i,j-1) + vymn(i,j)
         DIVh(i,j) = delvx*dxmi(i-1) + delvy*dymi(j-1)
         ron(i,j) = roh(i,j)*(1.0d0 - dt*DIVh(i,j))

         prn(i,j) = prh(i,j) - dt*(gm*(prh(i,j) + q(i,j))*DIVh(i,j)
     & + (gm - 1.0d0)*qdivv(i,j) + xvisl(i,j)*delvx*dxmi(i-1)
     &                           + yvisl(i,j)*delvy*dymi(j-1))
      enddo
      enddo

      return
      end subroutine normal_cip
