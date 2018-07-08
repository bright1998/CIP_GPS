c    -----------------
      subroutine bndc
c    -----------------
      use common
      implicit double precision(a-h,o-z)

      interface
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

         subroutine bdconx(da,margin,mbnd,men,const)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: const
            integer :: margin,mbnd,men
         end subroutine bdconx

         subroutine bdinix(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdinix

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

         subroutine bdcony(da,margin,mbnd,men,const)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: const
            integer :: margin,mbnd,men
         end subroutine bdcony

         subroutine bdiniy(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdiniy
      end interface

      if(nbnd1 .eq. 1) then
         call bdfrdx(ro,ro0,4,0,0)
         call bdinix(rodx,rodx0,4,0,0)
         call bdfrdx(rody,rody0,4,0,0)

         call bdfrdx(vxm,vxm0,4,0,1)
         call bdinix(vxdxm,vxdxm0,4,0,1)
         call bdfrdx(vxdym,vxdym0,4,0,1)

         call bdfrdx(vym,vym0,4,0,0)
         call bdinix(vydxm,vydxm0,4,0,0)
         call bdfrdx(vydym,vydym0,4,0,0)

         call bdfrdx(vz,vz0,4,0,0)
         call bdinix(vzdx,vzdx0,4,0,0)
         call bdfrdx(vzdy,vzdy0,4,0,0)

         call bdfrdx(pr,pr0,4,0,0)
         call bdinix(prdx,prdx0,4,0,0)
         call bdfrdx(prdy,prdy0,4,0,0)
      elseif(nbnd1 .eq. 2) then
         call bdfrex(ro,4,0,0)
         call bdconx(rodx,4,0,0,0.d0)
         call bdfrex(rody,4,0,0)

         call bdfrex(vxm,4,0,1)
         call bdconx(vxdxm,4,0,1,0.d0)
         call bdfrex(vxdym,4,0,1)

         call bdfrex(vym,4,0,0)
         call bdconx(vydxm,4,0,0,0.d0)
         call bdfrex(vydym,4,0,0)

         call bdfrex(vz,4,0,0)
         call bdconx(vzdx,4,0,0,0.d0)
         call bdfrex(vzdy,4,0,0)

         call bdfrex(pr,4,0,0)
         call bdconx(prdx,4,0,0,0.d0)
         call bdfrex(prdy,4,0,0)
      elseif(nbnd1 .eq. 3) then
         call bdperx(ro,4,0,0)
         call bdperx(rodx,4,0,0)
         call bdperx(rody,4,0,0)

         call bdperx(vxm,4,0,1)
         call bdperx(vxdxm,4,0,1)
         call bdperx(vxdym,4,0,1)

         call bdperx(vym,4,0,0)
         call bdperx(vydxm,4,0,0)
         call bdperx(vydym,4,0,0)

         call bdperx(vz,4,0,0)
         call bdperx(vzdx,4,0,0)
         call bdperx(vzdy,4,0,0)

         call bdperx(pr,4,0,0)
         call bdperx(prdx,4,0,0)
         call bdperx(prdy,4,0,0)
      elseif(nbnd1 .eq. 4) then
         call bdsymx(ro,4,0,0, 1.d0)
         call bdsymx(rodx,4,0,0,-1.d0)
         call bdsymx(rody,4,0,0, 1.d0)

         call bdsymx(vxm,4,0,1,-1.d0)
         call bdsymx(vxdxm,4,0,1, 1.d0)
         call bdsymx(vxdym,4,0,1,-1.d0)

         call bdsymx(vym,4,0,0, 1.d0)
         call bdsymx(vydxm,4,0,0,-1.d0)
         call bdsymx(vydym,4,0,0, 1.d0)

         call bdsymx(vz,4,0,0, 1.d0)
         call bdsymx(vzdx,4,0,0,-1.d0)
         call bdsymx(vzdy,4,0,0, 1.d0)

         call bdsymx(pr,4,0,0, 1.d0)
         call bdsymx(prdx,4,0,0,-1.d0)
         call bdsymx(prdy,4,0,0, 1.d0)
      endif

      if(nbnd2 .eq. 1) then
         call bdfrdx(ro,ro0,4,1,0)
         call bdinix(rodx,rodx0,4,1,0)
         call bdfrdx(rody,rody0,4,1,0)

         call bdfrdx(vxm,vxm0,4,1,1)
         call bdinix(vxdxm,vxdxm0,4,1,1)
         call bdfrdx(vxdym,vxdym0,4,1,1)

         call bdfrdx(vym,vym0,4,1,0)
         call bdinix(vydxm,vydxm0,4,1,0)
         call bdfrdx(vydym,vydym0,4,1,0)

         call bdfrdx(vz,vz0,4,1,0)
         call bdinix(vzdx,vzdx0,4,1,0)
         call bdfrdx(vzdy,vzdy0,4,1,0)

         call bdfrdx(pr,pr0,4,1,0)
         call bdinix(prdx,prdx0,4,1,0)
         call bdfrdx(prdy,prdy0,4,1,0)
      elseif(nbnd2 .eq. 2) then
         call bdfrex(ro,4,1,0)
         call bdconx(rodx,4,1,0,0.d0)
         call bdfrex(rody,4,1,0)

         call bdfrex(vxm,4,1,1)
         call bdconx(vxdxm,4,1,1,0.d0)
         call bdfrex(vxdym,4,1,1)

         call bdfrex(vym,4,1,0)
         call bdconx(vydxm,4,1,0,0.d0)
         call bdfrex(vydym,4,1,0)

         call bdfrex(vz,4,1,0)
         call bdconx(vzdx,4,1,0,0.d0)
         call bdfrex(vzdy,4,1,0)

         call bdfrex(pr,4,1,0)
         call bdconx(prdx,4,1,0,0.d0)
         call bdfrex(prdy,4,1,0)
      elseif(nbnd2 .eq. 3) then
         call bdperx(ro,4,1,0)
         call bdperx(rodx,4,1,0)
         call bdperx(rody,4,1,0)

         call bdperx(vxm,4,1,1)
         call bdperx(vxdxm,4,1,1)
         call bdperx(vxdym,4,1,1)

         call bdperx(vym,4,1,0)
         call bdperx(vydxm,4,1,0)
         call bdperx(vydym,4,1,0)

         call bdperx(vz,4,1,0)
         call bdperx(vzdx,4,1,0)
         call bdperx(vzdy,4,1,0)

         call bdperx(pr,4,1,0)
         call bdperx(prdx,4,1,0)
         call bdperx(prdy,4,1,0)
      elseif(nbnd2 .eq. 4) then
         call bdsymx(ro,4,1,0, 1.d0)
         call bdsymx(rodx,4,1,0,-1.d0)
         call bdsymx(rody,4,1,0, 1.d0)

         call bdsymx(vxm,4,1,1,-1.d0)
         call bdsymx(vxdxm,4,1,1, 1.d0)
         call bdsymx(vxdym,4,1,1,-1.d0)

         call bdsymx(vym,4,1,0, 1.d0)
         call bdsymx(vydxm,4,1,0,-1.d0)
         call bdsymx(vydym,4,1,0, 1.d0)

         call bdsymx(vz,4,1,0, 1.d0)
         call bdsymx(vzdx,4,1,0,-1.d0)
         call bdsymx(vzdy,4,1,0, 1.d0)

         call bdsymx(pr,4,1,0, 1.d0)
         call bdsymx(prdx,4,1,0,-1.d0)
         call bdsymx(prdy,4,1,0, 1.d0)
      endif

      if(nbnd3 .eq. 1) then
         call bdfrdy(ro,ro0,4,0,0)
         call bdfrdy(rodx,rodx0,4,0,0)
         call bdiniy(rody,rody0,4,0,0)

         call bdfrdy(vxm,vxm0,4,0,0)
         call bdfrdy(vxdxm,vxdxm0,4,0,0)
         call bdiniy(vxdym,vxdym0,4,0,0)

         call bdfrdy(vym,vym0,4,0,1)
         call bdfrdy(vydxm,vydxm0,4,0,1)
         call bdiniy(vydym,vydym0,4,0,1)

         call bdfrdy(vz,vz0,4,0,0)
         call bdfrdy(vzdx,vzdx0,4,0,0)
         call bdiniy(vzdy,vzdy0,4,0,0)

         call bdfrdy(pr,pr0,4,0,0)
         call bdfrdy(prdx,prdx0,4,0,0)
         call bdiniy(prdy,prdy0,4,0,0)
      elseif(nbnd3 .eq. 2) then
         call bdfrey(ro,4,0,0)
         call bdfrey(rodx,4,0,0)
         call bdcony(rody,4,0,0,0.d0)

         call bdfrey(vxm,4,0,0)
         call bdfrey(vxdxm,4,0,0)
         call bdcony(vxdym,4,0,0,0.d0)

         call bdfrey(vym,4,0,1)
         call bdfrey(vydxm,4,0,1)
         call bdcony(vydym,4,0,1,0.d0)

         call bdfrey(vz,4,0,0)
         call bdfrey(vzdx,4,0,0)
         call bdcony(vzdy,4,0,0,0.d0)

         call bdfrey(pr,4,0,0)
         call bdfrey(prdx,4,0,0)
         call bdcony(prdy,4,0,0,0.d0)
      elseif(nbnd3 .eq. 3) then
         call bdpery(ro,4,0,0)
         call bdpery(rodx,4,0,0)
         call bdpery(rody,4,0,0)

         call bdpery(vxm,4,0,0)
         call bdpery(vxdxm,4,0,0)
         call bdpery(vxdym,4,0,0)

         call bdpery(vym,4,0,1)
         call bdpery(vydxm,4,0,1)
         call bdpery(vydym,4,0,1)

         call bdpery(vz,4,0,0)
         call bdpery(vzdx,4,0,0)
         call bdpery(vzdy,4,0,0)

         call bdpery(pr,4,0,0)
         call bdpery(prdx,4,0,0)
         call bdpery(prdy,4,0,0)
      elseif(nbnd3 .eq. 4) then
         call bdsymy(ro,4,0,0, 1.d0)
         call bdsymy(rodx,4,0,0, 1.d0)
         call bdsymy(rody,4,0,0,-1.d0)

         call bdsymy(vxm,4,0,0, 1.d0)
         call bdsymy(vxdxm,4,0,0, 1.d0)
         call bdsymy(vxdym,4,0,0,-1.d0)

         call bdsymy(vym,4,0,1,-1.d0)
         call bdsymy(vydxm,4,0,1,-1.d0)
         call bdsymy(vydym,4,0,1, 1.d0)

         call bdsymy(vz,4,0,0, 1.d0)
         call bdsymy(vzdx,4,0,0, 1.d0)
         call bdsymy(vzdy,4,0,0,-1.d0)

         call bdsymy(pr,4,0,0, 1.d0)
         call bdsymy(prdx,4,0,0, 1.d0)
         call bdsymy(prdy,4,0,0,-1.d0)
      endif

      if(nbnd4 .eq. 1) then
         call bdfrdy(ro,ro0,4,1,0)
         call bdfrdy(rodx,rodx0,4,1,0)
         call bdiniy(rody,rody0,4,1,0)

         call bdfrdy(vxm,vxm0,4,1,0)
         call bdfrdy(vxdxm,vxdxm0,4,1,0)
         call bdiniy(vxdym,vxdym0,4,1,0)

         call bdfrdy(vym,vym0,4,1,1)
         call bdfrdy(vydxm,vydxm0,4,1,1)
         call bdiniy(vydym,vydym0,4,1,1)

         call bdfrdy(vz,vz0,4,1,0)
         call bdfrdy(vzdx,vzdx0,4,1,0)
         call bdiniy(vzdy,vzdy0,4,1,0)

         call bdfrdy(pr,pr0,4,1,0)
         call bdfrdy(prdx,prdx0,4,1,0)
         call bdiniy(prdy,prdy0,4,1,0)
      elseif(nbnd4 .eq. 2) then
         call bdfrey(ro,4,1,0)
         call bdfrey(rodx,4,1,0)
         call bdcony(rody,4,1,0,0.d0)

         call bdfrey(vxm,4,1,0)
         call bdfrey(vxdxm,4,1,0)
         call bdcony(vxdym,4,1,0,0.d0)

         call bdfrey(vym,4,1,1)
         call bdfrey(vydxm,4,1,1)
         call bdcony(vydym,4,1,1,0.d0)

         call bdfrey(vz,4,1,0)
         call bdfrey(vzdx,4,1,0)
         call bdcony(vzdy,4,1,0,0.d0)

         call bdfrey(pr,4,1,0)
         call bdfrey(prdx,4,1,0)
         call bdcony(prdy,4,1,0,0.d0)
      elseif(nbnd4 .eq. 3) then
         call bdpery(ro,4,1,0)
         call bdpery(rodx,4,1,0)
         call bdpery(rody,4,1,0)

         call bdpery(vxm,4,1,0)
         call bdpery(vxdxm,4,1,0)
         call bdpery(vxdym,4,1,0)

         call bdpery(vym,4,1,1)
         call bdpery(vydxm,4,1,1)
         call bdpery(vydym,4,1,1)

         call bdpery(vz,4,1,0)
         call bdpery(vzdx,4,1,0)
         call bdpery(vzdy,4,1,0)

         call bdpery(pr,4,1,0)
         call bdpery(prdx,4,1,0)
         call bdpery(prdy,4,1,0)
      elseif(nbnd4 .eq. 4) then
         call bdsymy(ro,4,1,0, 1.d0)
         call bdsymy(rodx,4,1,0, 1.d0)
         call bdsymy(rody,4,1,0,-1.d0)

         call bdsymy(vxm,4,1,0, 1.d0)
         call bdsymy(vxdxm,4,1,0, 1.d0)
         call bdsymy(vxdym,4,1,0,-1.d0)

         call bdsymy(vym,4,1,1,-1.d0)
         call bdsymy(vydxm,4,1,1,-1.d0)
         call bdsymy(vydym,4,1,1, 1.d0)

         call bdsymy(vz,4,1,0, 1.d0)
         call bdsymy(vzdx,4,1,0, 1.d0)
         call bdsymy(vzdy,4,1,0,-1.d0)

         call bdsymy(pr,4,1,0, 1.d0)
         call bdsymy(prdx,4,1,0, 1.d0)
         call bdsymy(prdy,4,1,0,-1.d0)
      endif

      return
      end subroutine bndc

c    ---------------------
      subroutine bndc_mag
c    ---------------------
      use common
      implicit double precision(a-h,o-z)

      interface
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

         subroutine bdconx(da,margin,mbnd,men,const)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: const
            integer :: margin,mbnd,men
         end subroutine bdconx

         subroutine bdinix(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdinix

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

         subroutine bdcony(da,margin,mbnd,men,const)
            double precision,dimension(:,:),intent(inout) :: da
            double precision :: const
            integer :: margin,mbnd,men
         end subroutine bdcony

         subroutine bdiniy(da,da0,margin,mbnd,men)
            double precision,dimension(:,:),intent(inout) :: da
            double precision,dimension(:,:),intent(in) :: da0
            integer :: margin,mbnd,men
         end subroutine bdiniy
      end interface

      if(nbnd1 .eq. 1) then
         call bdfrdx(bxm,bxm0,4,0,1)
         call bdfrdx(bym,bym0,4,0,0)
         call bdfrdx(bz,bz0,4,0,0)
      elseif(nbnd1 .eq. 2) then
         call bdfrex(bxm,4,0,1)
         call bdfrex(bym,4,0,0)
         call bdfrex(bz,4,0,0)
      elseif(nbnd1 .eq. 3) then
         call bdperx(bxm,4,0,1)
         call bdperx(bym,4,0,0)
         call bdperx(bz,4,0,0)
      elseif(nbnd1 .eq. 4) then
         call bdsymx(bxm,4,0,1,-1.d0)
         call bdsymx(bym,4,0,0, 1.d0)
         call bdsymx(bz,4,0,0, 1.d0)
      endif

      if(nbnd2 .eq. 1) then
         call bdfrdx(bxm,bxm0,4,1,1)
         call bdfrdx(bym,bym0,4,1,0)
         call bdfrdx(bz,bz0,4,1,0)
      elseif(nbnd2 .eq. 2) then
         call bdfrex(bxm,4,1,1)
         call bdfrex(bym,4,1,0)
         call bdfrex(bz,4,1,0)
      elseif(nbnd2 .eq. 3) then
         call bdperx(bxm,4,1,1)
         call bdperx(bym,4,1,0)
         call bdperx(bz,4,1,0)
      elseif(nbnd2 .eq. 4) then
         call bdsymx(bxm,4,1,1,-1.d0)
         call bdsymx(bym,4,1,0, 1.d0)
         call bdsymx(bz,4,1,0, 1.d0)
      endif

      if(nbnd3 .eq. 1) then
         call bdfrdy(bxm,bxm0,4,0,0)
         call bdfrdy(bym,bym0,4,0,1)
         call bdfrdy(bz,bz0,4,0,0)
      elseif(nbnd3 .eq. 2) then
         call bdfrey(bxm,4,0,0)
         call bdfrey(bym,4,0,1)
         call bdfrey(bz,4,0,0)
      elseif(nbnd3 .eq. 3) then
         call bdpery(bxm,4,0,0)
         call bdpery(bym,4,0,1)
         call bdpery(bz,4,0,0)
      elseif(nbnd3 .eq. 4) then
         call bdsymy(bxm,4,0,0, 1.d0)
         call bdsymy(bym,4,0,1,-1.d0)
         call bdsymy(bz,4,0,0, 1.d0)
      endif

      if(nbnd4 .eq. 1) then
         call bdfrdy(bxm,bxm0,4,1,0)
         call bdfrdy(bym,bym0,4,1,1)
         call bdfrdy(bz,bz0,4,1,0)
      elseif(nbnd4 .eq. 2) then
         call bdfrey(bxm,4,1,0)
         call bdfrey(bym,4,1,1)
         call bdfrey(bz,4,1,0)
      elseif(nbnd4 .eq. 3) then
         call bdpery(bxm,4,1,0)
         call bdpery(bym,4,1,1)
         call bdpery(bz,4,1,0)
      elseif(nbnd4 .eq. 4) then
         call bdsymy(bxm,4,1,0, 1.d0)
         call bdsymy(bym,4,1,1,-1.d0)
         call bdsymy(bz,4,1,0, 1.d0)
      endif

      return
      end subroutine bndc_mag

c    -------------------------------------------
      subroutine bdfrdx(da,da0,margin,mbnd,men)
c    -------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision,dimension(:,:),intent(in) :: da0
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         ibnd = 1+margin
         do j=1,j0
         do i=1,margin
            da(ibnd-i,j) = (da0(ibnd-i,j) - da0(ibnd,j)) + da(ibnd,j)
         enddo
         enddo
      else
         ibnd = i0-margin-men
         do j=1,j0
         do i=1,margin
            da(ibnd+i,j) = (da0(ibnd+i,j) - da0(ibnd,j)) + da(ibnd,j)
         enddo
         enddo
      endif

      return
      end subroutine bdfrdx

c    ---------------------------------------
      subroutine bdfrex(da,margin,mbnd,men)
c    ---------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         do j=1,j0
         do i=1,margin
            da(i,j) = da(margin+1,j)
         enddo
         enddo
      else
         ibnd1 = i0+1-men
         do j=1,j0
         do i=1,margin
            da(ibnd1-i,j) = da(ibnd1-margin-1,j)
         enddo
         enddo
      endif

      return
      end subroutine bdfrex

c    ---------------------------------------
      subroutine bdperx(da,margin,mbnd,men)
c    ---------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         ibnd = 1+margin
         do i=1,margin
            do j=1,j0
               da(ibnd-i,j) = da(i0-margin-i,j)
            enddo
         enddo
      else
         ibnd = i0-margin-men
         do i=1,margin
            do j=1,j0
               if(men .eq. 0) then
                  da(ibnd+i,j) = da(margin+1+i,j)
               else
                  da(ibnd+i,j) = da(margin+i,j)
               endif
            enddo
         enddo
      endif

      return
      end subroutine bdperx

c    --------------------------------------------
      subroutine bdsymx(da,margin,mbnd,men,coff)
c    --------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision :: coff
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         ibnd1 = 1+margin !(4 for men=0,1)
         ibnd2 = 1+margin-men !(4 for men=0, 3 for men=1)
         do i=1,margin
            do j=1,j0
               da(ibnd1-i,j) = coff*da(ibnd2+i,j)
            enddo
         enddo
      else
         ibnd1 = i0-margin !(i0-3 for men=0,1)
         ibnd2 = i0-margin+men !(i0-3 for men=0, i0-2 for men=1)
         do i=1,margin
            do j=1,j0
               da(ibnd1+i,j) = coff*da(ibnd2-i,j)
            enddo
         enddo
      endif

      end subroutine bdsymx

c    ---------------------------------------------
      subroutine bdconx(da,margin,mbnd,men,const)
c    ---------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision :: const
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         do j=1,j0
         do i=1,margin
            da(i,j) = const
         enddo
         enddo
      else
         ibnd1 = i0+1-men
         do j=1,j0
         do i=1,margin
            da(ibnd1-i,j) = const
         enddo
         enddo
      endif

      return
      end subroutine bdconx

c    -------------------------------------------
      subroutine bdinix(da,da0,margin,mbnd,men)
c    -------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision,dimension(:,:),intent(in) :: da0
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         ibnd = 1+margin
         do j=1,j0
         do i=1,margin
            da(ibnd-i,j) = da0(ibnd-i,j)
         enddo
         enddo
      else
         ibnd = i0-margin-men
         do j=1,j0
         do i=1,margin
            da(ibnd+i,j) = da0(ibnd+i,j)
         enddo
         enddo
      endif

      return
      end subroutine bdinix

c    -------------------------------------------
      subroutine bdfrdy(da,da0,margin,mbnd,men)
c    -------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision,dimension(:,:),intent(in) :: da0
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         jbnd = 1+margin
         do j=1,margin
         do i=1,i0
            da(i,jbnd-j) = (da0(i,jbnd-j) - da0(i,jbnd)) + da(i,jbnd)
         enddo
         enddo
      else
         jbnd = j0-margin-men
         do j=1,margin
         do i=1,i0
            da(i,jbnd+j) = (da0(i,jbnd+j) - da0(i,jbnd)) + da(i,jbnd)
         enddo
         enddo
      endif

      return
      end subroutine bdfrdy

c    ---------------------------------------
      subroutine bdfrey(da,margin,mbnd,men)
c    ---------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         do j=1,margin
         do i=1,i0
            da(i,j) = da(i,margin+1)
         enddo
         enddo
      else
         jbnd1 = j0+1-men
         do j=1,margin
         do i=1,i0
            da(i,jbnd1-j) = da(i,jbnd1-margin-1)
         enddo
         enddo
      endif

      return
      end subroutine bdfrey

c    ---------------------------------------
      subroutine bdpery(da,margin,mbnd,men)
c    ---------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         jbnd = 1+margin
         do j=1,margin
            do i=1,i0
               da(i,jbnd-j) = da(i,j0-margin-j)
            enddo
         enddo
      else
         jbnd = j0-margin-men
         do j=1,margin
            do i=1,i0
               if(men .eq. 0) then
                  da(i,jbnd+j) = da(i,margin+1+j)
               else
                  da(i,jbnd+j) = da(i,margin+j)
               endif
            enddo
         enddo
      endif

      return
      end subroutine bdpery

c    --------------------------------------------
      subroutine bdsymy(da,margin,mbnd,men,coff)
c    --------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision :: coff
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         jbnd1 = 1+margin
         jbnd2 = 1+margin-men
         do j=1,margin
            do i=1,i0
               da(i,jbnd1-j) = coff*da(i,jbnd2+j)
            enddo
         enddo
      else
         jbnd1 = j0-margin
         jbnd2 = j0-margin+men
         do j=1,margin
            do i=1,i0
               da(i,jbnd1+j) = coff*da(i,jbnd2-j)
            enddo
         enddo
      endif

      end subroutine bdsymy

c    ---------------------------------------------
      subroutine bdcony(da,margin,mbnd,men,const)
c    ---------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision :: const
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         do j=1,margin
         do i=1,i0
            da(i,j) = const
         enddo
         enddo
      else
         jbnd1 = j0+1-men
         do j=1,margin
         do i=1,i0
            da(i,jbnd1-j) = const
         enddo
         enddo
      endif

      return
      end subroutine bdcony

c    -------------------------------------------
      subroutine bdiniy(da,da0,margin,mbnd,men)
c    -------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(inout) :: da
      double precision,dimension(:,:),intent(in) :: da0
      integer :: margin,mbnd,men

      if(mbnd .eq. 0) then
         jbnd = 1+margin
         do j=1,margin
         do i=1,i0
            da(i,jbnd-j) = da0(i,jbnd-j)
         enddo
         enddo
      else
         jbnd = j0-margin-men
         do j=1,margin
         do i=1,i0
            da(i,jbnd+j) = da0(i,jbnd+j)
         enddo
         enddo
      endif

      return
      end subroutine bdiniy
