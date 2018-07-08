c    --------------------------------------------
      subroutine vanLeer(dl,qp,qm,q,cp,cm,is,ie)
c    --------------------------------------------
      use common
      implicit double precision(a-h,o-z)

      double precision,dimension(:),intent(out) :: qp,qm
      double precision,dimension(:),intent(in) :: dl,q,cp,cm
      integer :: is,ie

      dimension :: delq(1:ijmax),delqds(1:ijmax)

      do i=is-1,ie+1
         delq(i) = (-q(i) + q(i+1))/dl(i)
      enddo

      do i=is,ie+1
         if(delq(i-1)*delq(i) .gt. 0.d0) then
            delqds(i) = 2.d0*delq(i-1)*delq(i)/(delq(i-1) + delq(i))
         else
            delqds(i) = 0.d0
         endif
      enddo

      do i=is,ie
         if(cp(i) .gt. 0.d0) then
            qp(i) = q(i) + 0.5d0*(dl(i) - cp(i)*dt)*delqds(i)
         else
            qp(i) = q(i+1) - 0.5d0*(dl(i) + cp(i)*dt)*delqds(i+1)
         endif

         if(cm(i) .gt. 0.d0) then
            qm(i) = q(i) + 0.5d0*(dl(i) - cm(i)*dt)*delqds(i)
         else
            qm(i) = q(i+1) - 0.5d0*(dl(i) + cm(i)*dt)*delqds(i+1)
         endif
      enddo

      return
      end subroutine vanLeer
