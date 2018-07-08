c    -----------------
      subroutine read
c    -----------------
      use common
      implicit double precision(a-h,o-z)

C Read the data for Restart Calculation
      open(22,file='out/restart.data',form='unformatted')
      read(22)tm,ns,nf,nbcast
      read(22)
     &((  ro(i,j),i=1,i0),j=1,j0),
     &((rodx(i,j),i=1,i0),j=1,j0),
     &((rody(i,j),i=1,i0),j=1,j0),
     &((  vx(i,j),i=1,i0),j=1,j0),
     &((  vy(i,j),i=1,i0),j=1,j0),
     &((  vz(i,j),i=1,i0),j=1,j0),
     &((vzdx(i,j),i=1,i0),j=1,j0),
     &((vzdy(i,j),i=1,i0),j=1,j0),
     &((  bx(i,j),i=1,i0),j=1,j0),
     &((  by(i,j),i=1,i0),j=1,j0),
     &((  bz(i,j),i=1,i0),j=1,j0),
     &((  pr(i,j),i=1,i0),j=1,j0),
     &((prdx(i,j),i=1,i0),j=1,j0),
     &((prdy(i,j),i=1,i0),j=1,j0),
     &((  te(i,j),i=1,i0),j=1,j0)
      read(22)
     &((  vxm(i,j),i=1,i0-1),j=1,j0),
     &((vxdxm(i,j),i=1,i0-1),j=1,j0),
     &((vxdym(i,j),i=1,i0-1),j=1,j0),
     &((  bxm(i,j),i=1,i0-1),j=1,j0)
      read(22)
     &((  vym(i,j),i=1,i0),j=1,j0-1),
     &((vydxm(i,j),i=1,i0),j=1,j0-1),
     &((vydym(i,j),i=1,i0),j=1,j0-1),
     &((  bym(i,j),i=1,i0),j=1,j0-1)
      close(22)

      write(6,*)' '
      write(6,*)
     &'read data for restart calculation'

      return
      end subroutine read
