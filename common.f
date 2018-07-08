c    ---------------
      module common
c    ---------------
      implicit double precision(a-h,o-z)

c main.f
      double precision :: tm,etm,pi,pi4,pi4i
      integer :: ns,nf,nbcast

c prms.f
      integer :: nsmax,ic,i0,jc,j0,ncflconst,ntcon,itemax,ijmax,
     &           nrestart
      integer :: nbnd1,nbnd2,nbnd3,nbnd4
      double precision :: xcmin,xcmax,rax
      double precision :: ycmin,ycmax,ray
      double precision :: tmmax,ftm,dtlim,cn,dtcon,gm
      double precision :: pec,peci,epsit
      double precision :: dmin,pmin
      double precision :: fv,fvt,fvl

c cflc.f
      double precision :: dt

C Arraies on Normal Grids
      double precision,allocatable :: x(:),dx(:),dxi(:)
      double precision,allocatable :: y(:),dy(:),dyi(:)

      double precision,allocatable ::  ro(:,:),ro0(:,:),
     &                                roh(:,:),ron(:,:)
      double precision,allocatable :: vx(:,:),vx0(:,:),vxh(:,:)
      double precision,allocatable :: vy(:,:),vy0(:,:),vyh(:,:)
      double precision,allocatable ::  vz(:,:),vz0(:,:),
     &                                vzh(:,:),vzn(:,:)
      double precision,allocatable :: bx(:,:),bx0(:,:)
      double precision,allocatable :: by(:,:),by0(:,:)
      double precision,allocatable :: bz(:,:),bz0(:,:)
      double precision,allocatable ::  pr(:,:),pr0(:,:),
     &                                prh(:,:),prn(:,:)
      double precision,allocatable :: te(:,:),te0(:,:)

C Arraies on Staggered Grids
      double precision,allocatable :: xm(:),dxm(:),dxmi(:)
      double precision,allocatable :: ym(:),dym(:),dymi(:)

      double precision,allocatable ::  vxm(:,:),vxm0(:,:),
     &                                vxmh(:,:),vxmn(:,:)
      double precision,allocatable ::  vym(:,:),vym0(:,:),
     &                                vymh(:,:),vymn(:,:)
      double precision,allocatable :: bxm(:,:),bxm0(:,:)
      double precision,allocatable :: bym(:,:),bym0(:,:)

C For CIP Method
      double precision,allocatable ::  rodx(:,:),rodx0(:,:),
     &                                rodxh(:,:),rodxn(:,:)
      double precision,allocatable ::  vxdxm(:,:),vxdxm0(:,:),
     &                                vxdxmh(:,:),vxdxmn(:,:)
      double precision,allocatable ::  vydxm(:,:),vydxm0(:,:),
     &                                vydxmh(:,:),vydxmn(:,:)
      double precision,allocatable ::  vzdx(:,:),vzdx0(:,:),
     &                                vzdxh(:,:),vzdxn(:,:)
      double precision,allocatable ::  prdx(:,:),prdx0(:,:),
     &                                prdxh(:,:),prdxn(:,:)

      double precision,allocatable ::  rody(:,:),rody0(:,:),
     &                                rodyh(:,:),rodyn(:,:)
      double precision,allocatable ::  vxdym(:,:),vxdym0(:,:),
     &                                vxdymh(:,:),vxdymn(:,:)
      double precision,allocatable ::  vydym(:,:),vydym0(:,:),
     &                                vydymh(:,:),vydymn(:,:)
      double precision,allocatable ::  vzdy(:,:),vzdy0(:,:),
     &                                vzdyh(:,:),vzdyn(:,:)
      double precision,allocatable ::  prdy(:,:),prdy0(:,:),
     &                                prdyh(:,:),prdyn(:,:)

      double precision,allocatable :: vxmy(:,:),vymx(:,:),
     &                                vzmx(:,:),vzmy(:,:)
      double precision,allocatable :: bxmy(:,:),bymx(:,:),
     &                                bzmx(:,:),bzmy(:,:)
      double precision,allocatable :: ex(:,:),ey(:,:),ez(:,:)

C Source Term
      double precision,allocatable :: Fexx(:,:),Fexy(:,:),Fexz(:,:)
      double precision,allocatable :: Fexxm(:,:),Fexym(:,:)

      end module common
