c    --------------------
      subroutine cipdsrc
c    --------------------
      use common
      implicit double precision(a-h,o-z)

c Density, Pressure, Vz
      do j=5,j0-4
      do i=5,i0-4
         rodxn(i,j) = rodxh(i,j)
     & + (-ron(i-1,j) + ron(i+1,j) + roh(i-1,j) - roh(i+1,j))
     &  /(dx(i-1) + dx(i))
     & - (rodxh(i,j)*(-vxh(i-1,j) + vxh(i+1,j))
     &  + rodyh(i,j)*(-vyh(i-1,j) + vyh(i+1,j)))*dt/(dx(i-1) + dx(i))

         prdxn(i,j) = prdxh(i,j)
     & + (-prn(i-1,j) + prn(i+1,j) + prh(i-1,j) - prh(i+1,j))
     &  /(dx(i-1) + dx(i))
     & - (prdxh(i,j)*(-vxh(i-1,j) + vxh(i+1,j))
     &  + prdyh(i,j)*(-vyh(i-1,j) + vyh(i+1,j)))*dt/(dx(i-1) + dx(i))

         vzdxn(i,j) = vzdxh(i,j)
     & + (-vzn(i-1,j) + vzn(i+1,j) + vzh(i-1,j) - vzh(i+1,j))
     &  /(dx(i-1) + dx(i))
     & - (vzdxh(i,j)*(-vxh(i-1,j) + vxh(i+1,j))
     &  + vzdyh(i,j)*(-vyh(i-1,j) + vyh(i+1,j)))*dt/(dx(i-1) + dx(i))
      enddo
      enddo

      do j=5,j0-4
      do i=5,i0-4
         rodyn(i,j) = rodyh(i,j)
     & + (-ron(i,j-1) + ron(i,j+1) + roh(i,j-1) - roh(i,j+1))
     &  /(dy(j-1) + dy(j))
     & - (rodxh(i,j)*(-vxh(i,j-1) + vxh(i,j+1))
     &  + rodyh(i,j)*(-vyh(i,j-1) + vyh(i,j+1)))*dt/(dy(j-1) + dy(j))

         prdyn(i,j) = prdyh(i,j)
     & + (-prn(i,j-1) + prn(i,j+1) + prh(i,j-1) - prh(i,j+1))
     &  /(dy(j-1) + dy(j))
     & - (prdxh(i,j)*(-vxh(i,j-1) + vxh(i,j+1))
     &  + prdyh(i,j)*(-vyh(i,j-1) + vyh(i,j+1)))*dt/(dy(j-1) + dy(j))

         vzdyn(i,j) = vzdyh(i,j)
     & + (-vzn(i,j-1) + vzn(i,j+1) + vzh(i,j-1) - vzh(i,j+1))
     &  /(dy(j-1) + dy(j))
     & - (vzdxh(i,j)*(-vxh(i,j-1) + vxh(i,j+1))
     &  + vzdyh(i,j)*(-vyh(i,j-1) + vyh(i,j+1)))*dt/(dy(j-1) + dy(j))
      enddo
      enddo

      do j=4,j0-3
      do i=4,i0-4
         vxdxmn(i,j) = vxdxmh(i,j)
     & + (-vxmn(i-1,j) + vxmn(i+1,j) + vxmh(i-1,j) - vxmh(i+1,j))
     &  /(dxm(i-1) + dxm(i))
     & - (vxdxmh(i,j)*(-vxmh(i-1,j) + vxmh(i+1,j))
     &  + vxdymh(i,j)*(-vymx(i-1,j) + vymx(i+1,j)))*dt
     &  /(dxm(i-1) + dxm(i))
      enddo
      enddo

      do j=3,j0-3
      do i=5,i0-4
         vydxmn(i,j) = vydxmh(i,j)
     & + (-vymn(i-1,j) + vymn(i+1,j) + vymh(i-1,j) - vymh(i+1,j))
     &  /(dx(i-1) + dx(i))
     & - (vydxmh(i,j)*(-vxmy(i-1,j) + vxmy(i+1,j))
     &  + vydymh(i,j)*(-vymh(i-1,j) + vymh(i+1,j)))*dt/(dx(i-1) + dx(i))
      enddo
      enddo

      do j=5,j0-4
      do i=3,i0-3
         vxdymn(i,j) = vxdymh(i,j)
     & + (-vxmn(i,j-1) + vxmn(i,j+1) + vxmh(i,j-1) - vxmh(i,j+1))
     &  /(dy(j-1) + dy(j))
     & - (vxdxmh(i,j)*(-vxmh(i,j-1) + vxmh(i,j+1))
     &  + vxdymh(i,j)*(-vymx(i,j-1) + vymx(i,j+1)))*dt/(dy(j-1) + dy(j))
      enddo
      enddo

      do j=4,j0-4
      do i=4,i0-3
         vydymn(i,j) = vydymh(i,j)
     & + (-vymn(i,j-1) + vymn(i,j+1) + vymh(i,j-1) - vymh(i,j+1))
     &  /(dym(j-1) + dym(j))
     & - (vydxmh(i,j)*(-vxmy(i,j-1) + vxmy(i,j+1))
     &  + vydymh(i,j)*(-vymh(i,j-1) + vymh(i,j+1)))*dt
     &  /(dym(j-1) + dym(j))
      enddo
      enddo

      return
      end subroutine cipdsrc
