c    ----------------
      subroutine cip
c    ----------------
      use common
      implicit double precision(a-h,o-z)

      interface
         subroutine cipadv(cadx,cady,da,dan,dadx,dady,dadxn,dadyn,
     &u,v,is,ie,js,je)
            double precision,dimension(:),intent(in) :: cadx,cady
            double precision,dimension(:,:),intent(in) :: da,dadx,dady,
     &u,v
            double precision,dimension(:,:),intent(inout) :: dan,dadxn,
     &dadyn
            integer :: is,ie,js,je
         end subroutine cipadv

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

         subroutine convertxm(dam,da,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: dam
            double precision,dimension(:,:),intent(inout) :: da
            integer :: is,ie,js,je
         end subroutine convertxm

         subroutine converty(da,dam,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: da
            double precision,dimension(:,:),intent(inout) :: dam
            integer :: is,ie,js,je
         end subroutine converty

         subroutine convertym(dam,da,is,ie,js,je)
            double precision,dimension(:,:),intent(in) :: dam
            double precision,dimension(:,:),intent(inout) :: da
            integer :: is,ie,js,je
         end subroutine convertym

         subroutine bdfrdx(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdfrdx

         subroutine bdfrex(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdfrex

         subroutine bdperx(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdperx

         subroutine bdsymx(da,margin,mbnd,men,coff)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: coff
            integer :: margin,mbnd,men
         end subroutine bdsymx

         subroutine bdfrdy(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdfrdy

         subroutine bdfrey(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdfrey

         subroutine bdpery(da,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            integer :: margin,mbnd,men
         end subroutine bdpery

         subroutine bdsymy(da,margin,mbnd,men,coff)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: coff
            integer :: margin,mbnd,men
         end subroutine bdsymy

         subroutine tconduct(is,ie,js,je)
            integer :: is,ie,js,je
         end subroutine tconduct
      end interface

      dimension :: DIV(1:i0,1:j0),DIVh(1:i0,1:j0)
      dimension :: s(1:i0,1:j0),c(1:i0,1:j0),d(1:i0,1:j0),e(1:i0,1:j0)

C Advection Phase

C Density
C Log convert
      do j=1,j0
      do i=1,i0
         rodx(i,j) = rodx(i,j)/ro(i,j)
         rody(i,j) = rody(i,j)/ro(i,j)
         ro(i,j) = dlog(ro(i,j))
      enddo
      enddo

      call cipadv(dx,dy,ro,roh,rodx,rody,rodxh,rodyh,vx,vy,
     &2,i0-1,2,j0-1)

C Exp convert
      do j=2,j0-1
      do i=2,i0-1
         roh(i,j) = dexp(roh(i,j))
         rodxh(i,j) = roh(i,j)*rodxh(i,j)
         rodyh(i,j) = roh(i,j)*rodyh(i,j)
      enddo
      enddo

C Pressure
C Log convert
      do j=1,j0
      do i=1,i0
         prdx(i,j) = prdx(i,j)/pr(i,j)
         prdy(i,j) = prdy(i,j)/pr(i,j)
         pr(i,j) = dlog(pr(i,j))
      enddo
      enddo

      call cipadv(dx,dy,pr,prh,prdx,prdy,prdxh,prdyh,vx,vy,
     &2,i0-1,2,j0-1)

C Exp convert
      do j=2,j0-1
      do i=2,i0-1
         prh(i,j) = dexp(prh(i,j))
         prdxh(i,j) = prh(i,j)*prdxh(i,j)
         prdyh(i,j) = prh(i,j)*prdyh(i,j)
      enddo
      enddo

      call convertx(vy,vymx,1,i0-1,1,j0)
      call converty(vx,vxmy,1,i0,1,j0-1)

C Vx
      call cipadv(dxm,dy,vxm,vxmh,vxdxm,vxdym,vxdxmh,vxdymh,vxm,vymx,
     &2,i0-2,2,j0-1)
C Vy
      call cipadv(dx,dym,vym,vymh,vydxm,vydym,vydxmh,vydymh,vxmy,vym,
     &2,i0-1,2,j0-2)
C Vz
      call cipadv(dx,dy,vz,vzh,vzdx,vzdy,vzdxh,vzdyh,vx,vy,
     &2,i0-1,2,j0-1)

      call chck(roh,dmin,2,i0-1,2,j0-1)
      call chck(prh,pmin,2,i0-1,2,j0-1)

      call convertxm(vxmh,vxh,3,i0-2,2,j0-1)
      call convertym(vymh,vyh,2,i0-1,3,j0-2)

C Non-advection Phase
      call normal_cip

      call chck(ron,dmin,4,i0-3,4,j0-3)
      call chck(prn,pmin,4,i0-3,4,j0-3)

      call convertx(vyh,vymx,2,i0-2,3,j0-2)
      call converty(vxh,vxmy,3,i0-2,2,j0-2)

C update space-derivatives of physical variables
      call cipdsrc

C update
      do j=5,j0-4
      do i=5,i0-4
         ro(i,j) = ron(i,j)
         pr(i,j) = prn(i,j)
         vz(i,j) = vzn(i,j)
         rodx(i,j) = rodxn(i,j)
         prdx(i,j) = prdxn(i,j)
         vzdx(i,j) = vzdxn(i,j)
         rody(i,j) = rodyn(i,j)
         prdy(i,j) = prdyn(i,j)
         vzdy(i,j) = vzdyn(i,j)
      enddo
      enddo

      do j=5,j0-4
      do i=5,i0-5
         vxm(i,j) = vxmn(i,j)
         vxdxm(i,j) = vxdxmn(i,j)
         vxdym(i,j) = vxdymn(i,j)
      enddo
      enddo

      do j=5,j0-5
      do i=5,i0-4
         vym(i,j) = vymn(i,j)
         vydxm(i,j) = vydxmn(i,j)
         vydym(i,j) = vydymn(i,j)
      enddo
      enddo

      call bndc

      call convertxm(vxm,vx,2,i0-1,1,j0)
      call convertym(vym,vy,1,i0,2,j0-1)

      if(nbnd1 .eq. 1) then
         call bdfrdx(vx,vx0,1,0,0)
      elseif(nbnd1 .eq. 2) then
         call bdfrex(vx,1,0,0)
      elseif(nbnd1 .eq. 3) then
         call bdperx(vx,1,0,0)
      elseif(nbnd1 .eq. 4) then
         call bdsymx(vx,1,0,0,-1.d0)
      endif

      if(nbnd2 .eq. 1) then
         call bdfrdx(vx,vx0,1,1,0)
      elseif(nbnd2 .eq. 2) then
         call bdfrex(vx,1,1,0)
      elseif(nbnd2 .eq. 3) then
         call bdperx(vx,1,1,0)
      elseif(nbnd2 .eq. 4) then
         call bdsymx(vx,1,1,0,-1.d0)
      endif

      if(nbnd3 .eq. 1) then
         call bdfrdy(vy,vy0,1,0,0)
      elseif(nbnd3 .eq. 2) then
         call bdfrey(vy,1,0,0)
      elseif(nbnd3 .eq. 3) then
         call bdpery(vy,1,0,0)
      elseif(nbnd3 .eq. 4) then
         call bdsymy(vy,1,0,0,-1.d0)
      endif

      if(nbnd4 .eq. 1) then
         call bdfrdy(vy,vy0,1,1,0)
      elseif(nbnd4 .eq. 2) then
         call bdfrey(vy,1,1,0)
      elseif(nbnd4 .eq. 3) then
         call bdpery(vy,1,1,0)
      elseif(nbnd4 .eq. 4) then
         call bdsymy(vy,1,1,0,-1.d0)
      endif

C Electric Field
      call mocelex
      call moceley
      call mocelez

C Bx
      do j=3,j0-2
      do i=2,i0-2
C bxm(i,j) = bxm(i,j) - dt*(dEz/dy - dEy/dz)
C dEy/dz = 0
C Ez(i+1/2,j+1/2) = -(VxB)_z = -(Vx By - Vy Bx)
         bxm(i,j) = bxm(i,j) - dt*(ez(i,j) - ez(i,j-1))*dymi(j-1)
      enddo
      enddo

C By
      do j=2,j0-2
      do i=3,i0-2
C bym(i,j) = bym(i,j) - dt*(dEx/dz - dEz/dx)
C dEx/dz = 0
C Ez(i+1/2,j+1/2) = -(VxB)_z = -(Vx By - Vy Bx)
         bym(i,j) = bym(i,j) + dt*(ez(i,j) - ez(i-1,j))*dxmi(i-1)
      enddo
      enddo

C Bz
      do j=3,j0-2
      do i=3,i0-2
C bz(i,j) = bz(i,j) - dt*(dEy/dx - dEx/dy)
C Ex(i,j+1/2) = -(VxB)_x = -(Vy Bz - Vz By)
C Ey(i+1/2,j) = -(VxB)_y = -(Vz Bx - Vx Bz)
         bz(i,j) = bz(i,j) - dt*((ey(i,j) - ey(i-1,j))*dxmi(i-1)
     &                         - (ex(i,j) - ex(i,j-1))*dymi(j-1))
      enddo
      enddo

      call bndc_mag

      call convertxm(bxm,bx,2,i0-1,1,j0)
      call convertym(bym,by,1,i0,2,j0-1)

      if(nbnd1 .eq. 1) then
         call bdfrdx(bx,bx0,1,0,0)
      elseif(nbnd1 .eq. 2) then
         call bdfrex(bx,1,0,0)
      elseif(nbnd1 .eq. 3) then
         call bdperx(bx,1,0,0)
      elseif(nbnd1 .eq. 4) then
         call bdsymx(bx,1,0,0,-1.d0)
      endif

      if(nbnd2 .eq. 1) then
         call bdfrdx(bx,bx0,1,1,0)
      elseif(nbnd2 .eq. 2) then
         call bdfrex(bx,1,1,0)
      elseif(nbnd2 .eq. 3) then
         call bdperx(bx,1,1,0)
      elseif(nbnd2 .eq. 4) then
         call bdsymx(bx,1,1,0,-1.d0)
      endif

      if(nbnd3 .eq. 1) then
         call bdfrdy(by,by0,1,0,0)
      elseif(nbnd3 .eq. 2) then
         call bdfrey(by,1,0,0)
      elseif(nbnd3 .eq. 3) then
         call bdpery(by,1,0,0)
      elseif(nbnd3 .eq. 4) then
         call bdsymy(by,1,0,0,-1.d0)
      endif

      if(nbnd4 .eq. 1) then
         call bdfrdy(by,by0,1,1,0)
      elseif(nbnd4 .eq. 2) then
         call bdfrey(by,1,1,0)
      elseif(nbnd4 .eq. 3) then
         call bdpery(by,1,1,0)
      elseif(nbnd4 .eq. 4) then
         call bdsymy(by,1,1,0,-1.d0)
      endif

C Calculate temperature [(gm - 1)*C_p*T]
      do j=1,j0
      do i=1,i0
         te(i,j) = gm*pr(i,j)/ro(i,j)
      enddo
      enddo
      if(ntcon .eq. 1) call tconduct(2,i0-1,2,j0-1)

      return
      end subroutine cip
