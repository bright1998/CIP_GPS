c    ---------------------
      subroutine allocate
c    ---------------------
      use common
      implicit double precision(a-h,o-z)

      allocate(x(1:i0),dx(1:i0),dxi(1:i0),
     &         y(1:j0),dy(1:j0),dyi(1:j0))

      allocate( ro(1:i0,1:j0),ro0(1:i0,1:j0),
     &         roh(1:i0,1:j0),ron(1:i0,1:j0),
     &          vx(1:i0,1:j0),vx0(1:i0,1:j0),vxh(1:i0,1:j0),
     &          vy(1:i0,1:j0),vy0(1:i0,1:j0),vyh(1:i0,1:j0),
     &          vz(1:i0,1:j0),vz0(1:i0,1:j0),
     &         vzh(1:i0,1:j0),vzn(1:i0,1:j0),
     &          bx(1:i0,1:j0),bx0(1:i0,1:j0),
     &          by(1:i0,1:j0),by0(1:i0,1:j0),
     &          bz(1:i0,1:j0),bz0(1:i0,1:j0),
     &          pr(1:i0,1:j0),pr0(1:i0,1:j0),
     &         prh(1:i0,1:j0),prn(1:i0,1:j0),
     &          te(1:i0,1:j0),te0(1:i0,1:j0))

      allocate(xm(1:i0),dxm(1:i0),dxmi(1:i0),
     &         ym(1:j0),dym(1:j0),dymi(1:j0))

      allocate( vxm(1:i0,1:j0),vxm0(1:i0,1:j0),
     &         vxmh(1:i0,1:j0),vxmn(1:i0,1:j0),
     &          vym(1:i0,1:j0),vym0(1:i0,1:j0),
     &         vymh(1:i0,1:j0),vymn(1:i0,1:j0),
     &          bxm(1:i0,1:j0),bxm0(1:i0,1:j0),
     &          bym(1:i0,1:j0),bym0(1:i0,1:j0))

C For CIP Method
      allocate(  rodx(1:i0,1:j0), rodx0(1:i0,1:j0),
     &          rodxh(1:i0,1:j0), rodxn(1:i0,1:j0),
     &          vxdxm(1:i0,1:j0),vxdxm0(1:i0,1:j0),
     &         vxdxmh(1:i0,1:j0),vxdxmn(1:i0,1:j0),
     &          vydxm(1:i0,1:j0),vydxm0(1:i0,1:j0),
     &         vydxmh(1:i0,1:j0),vydxmn(1:i0,1:j0),
     &           vzdx(1:i0,1:j0), vzdx0(1:i0,1:j0),
     &          vzdxh(1:i0,1:j0), vzdxn(1:i0,1:j0),
     &           prdx(1:i0,1:j0), prdx0(1:i0,1:j0),
     &          prdxh(1:i0,1:j0), prdxn(1:i0,1:j0))

      allocate(  rody(1:i0,1:j0), rody0(1:i0,1:j0),
     &          rodyh(1:i0,1:j0), rodyn(1:i0,1:j0),
     &          vxdym(1:i0,1:j0),vxdym0(1:i0,1:j0),
     &         vxdymh(1:i0,1:j0),vxdymn(1:i0,1:j0),
     &          vydym(1:i0,1:j0),vydym0(1:i0,1:j0),
     &         vydymh(1:i0,1:j0),vydymn(1:i0,1:j0),
     &           vzdy(1:i0,1:j0), vzdy0(1:i0,1:j0),
     &          vzdyh(1:i0,1:j0), vzdyn(1:i0,1:j0),
     &           prdy(1:i0,1:j0), prdy0(1:i0,1:j0),
     &          prdyh(1:i0,1:j0), prdyn(1:i0,1:j0))

      allocate(vxmy(1:i0,1:j0),vymx(1:i0,1:j0),
     &         vzmx(1:i0,1:j0),vzmy(1:i0,1:j0),
     &         bxmy(1:i0,1:j0),bymx(1:i0,1:j0),
     &         bzmx(1:i0,1:j0),bzmy(1:i0,1:j0),
     &         ex(1:i0,1:j0),ey(1:i0,1:j0),ez(1:i0,1:j0))

      allocate(Fexx(1:i0,1:j0),Fexy(1:i0,1:j0),Fexz(1:i0,1:j0))
      allocate(Fexxm(1:i0,1:j0),Fexym(1:i0,1:j0))

      do j=1,j0
      do i=1,i0
         Fexx(i,j) = 0.d0
         Fexy(i,j) = 0.d0
         Fexz(i,j) = 0.d0
      enddo
      enddo

      do j=1,j0
      do i=1,i0-1
         Fexxm(i,j) = 0.d0
      enddo
      enddo

      do j=1,j0-1 
      do i=1,i0
         Fexym(i,j) = 0.d0
      enddo
      enddo

      return
      end subroutine allocate
