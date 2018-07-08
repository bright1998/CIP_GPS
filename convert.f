c    -----------------------------------------
      subroutine convertx(da,dam,is,ie,js,je)
c    -----------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(in) :: da
      double precision,dimension(:,:),intent(inout) :: dam
      integer :: is,ie,js,je

      do j=js,je
      do i=is,ie
         tx = (xm(i) - x(i))/(x(i+1) - x(i))
         dam(i,j) = (1.d0 - tx)*da(i,j) + tx*da(i+1,j)
      enddo
      enddo

      return
      end subroutine convertx

c    ------------------------------------------
      subroutine convertxm(dam,da,is,ie,js,je)
c    ------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(in) :: dam
      double precision,dimension(:,:),intent(inout) :: da
      integer :: is,ie,js,je

      do j=js,je
      do i=is,ie
         tx = (xm(i) - x(i))/(-xm(i-1) + xm(i))
         da(i,j) = tx*dam(i-1,j) + (1.d0 - tx)*dam(i,j)
      enddo
      enddo

      return
      end subroutine convertxm

c    -----------------------------------------
      subroutine converty(da,dam,is,ie,js,je)
c    -----------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(in) :: da
      double precision,dimension(:,:),intent(inout) :: dam
      integer :: is,ie,js,je

      do j=js,je
      do i=is,ie
         ty = (ym(j) - y(j))/(y(j+1) - y(j))
         dam(i,j) = (1.d0 - ty)*da(i,j) + ty*da(i,j+1)
      enddo
      enddo

      return
      end subroutine converty

c    ------------------------------------------
      subroutine convertym(dam,da,is,ie,js,je)
c    ------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:,:),intent(in) :: dam
      double precision,dimension(:,:),intent(inout) :: da
      integer :: is,ie,js,je

      do j=js,je
      do i=is,ie
         ty = (ym(j) - y(j))/(-ym(j-1) + ym(j))
         da(i,j) = ty*dam(i,j-1) + (1.d0 - ty)*dam(i,j)
      enddo
      enddo

      return
      end subroutine convertym
