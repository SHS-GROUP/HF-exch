!=======================================================================
      program xcneo_hybrid

!=======================================================================
      implicit none
!     include 'mpif.h'
      include 'omp_lib.h'

!     integer nproc,myid,ierr


!     call mpi_set(nproc,myid,ierr)

!     call xcneo_hybrid_driver(nproc,myid)
      call xcneo_driver

!     call mpi_final(ierr)

      end
!=======================================================================
      subroutine xcneo_driver

!=======================================================================
      implicit none
!     include 'mpif.h'
      include 'omp_lib.h'

!     integer nproc,myid
!-------Basis Set Info-------(
      integer,allocatable :: ELCAM(:,:)  ! Angular mom for electrons
      double precision,allocatable :: ELCEX(:) ! Exponents: elec basis
      double precision,allocatable :: ELCBFC(:,:) ! Basis centers: elec basis
      integer,allocatable :: AMPEB2C(:) ! Map primitive index to contracted
      double precision,allocatable :: AGEBFCC(:) ! Map prim index to contract coef
      integer,allocatable :: KPESTR(:)  ! Map contracted index to primitive start
      integer,allocatable :: KPEEND(:)  ! Map contracted index to primitive end
      double precision,allocatable :: zan(:) ! Classical nuclear charges
      double precision,allocatable :: cat(:,:) ! XYZ Coordinates of classical atoms
      integer nat
      integer nelec
      integer NAE               ! Number of regular electrons
      integer NBE               ! Number of special electrons
      integer NAalpE,NAbetE
      double precision pmass    ! Mass of nonelectron quantum particle 
!-------Basis Set Info-------)
      integer i,j,idum,istat
      integer NUCST
      logical read_CE
      logical LDBG
      logical LSOSCF
      logical LOCBSE

      double precision a2bohr,bohr2a
      parameter(bohr2a=0.529177249d+00)
      parameter(a2bohr=1.0d+00/0.529177249d+00)

      integer npebf,nebf,npbf
      integer nebf2,npbf2,NPR,NEBFLT
      integer NPRA,NPRB

      integer ngee

      integer junk

      double precision wtime,wtime1,wtime2


      wtime = omp_get_wtime()


!-------READ-INPUT-FILE-AND-ALLOCATE-MEMORY-FOR-BASIS-SET--------------(
      open(unit=9,file='basis_definition.inp')

      read(9,*)nat
      if(allocated(zan)) deallocate(zan)
      allocate( zan(nat),stat=istat )
      if(allocated(cat)) deallocate(cat)
      allocate( cat(3,nat),stat=istat )
      do i=1,nat
         read(9,*)zan(i),cat(1,i),cat(2,i),cat(3,i)
         cat(1,i)=a2bohr*cat(1,i)
         cat(2,i)=a2bohr*cat(2,i)
         cat(3,i)=a2bohr*cat(3,i)
      end do

      read(9,*)nebf
      read(9,*)npebf
      if(allocated(AMPEB2C)) deallocate(AMPEB2C)
      allocate( AMPEB2C(npebf),stat=istat )
      if(allocated(ELCEX)) deallocate(ELCEX)
      allocate( ELCEX(npebf),stat=istat )
      if(allocated(AGEBFCC)) deallocate(AGEBFCC)
      allocate( AGEBFCC(npebf),stat=istat )
      if(allocated(ELCAM)) deallocate(ELCAM)
      allocate( ELCAM(npebf,3),stat=istat )
      if(allocated(ELCBFC)) deallocate(ELCBFC)
      allocate( ELCBFC(npebf,3),stat=istat )
      do i=1,npebf
         read(9,*)idum,AMPEB2C(i),ELCAM(i,1),ELCAM(i,2),ELCAM(i,3),
     x ELCEX(i),AGEBFCC(i),ELCBFC(i,1),ELCBFC(i,2),ELCBFC(i,3)
         do j=1,3
            ELCBFC(i,j)=a2bohr*ELCBFC(i,j)
         end do
      end do

      read(9,*)nelec
      read(9,*)NAE
      read(9,*)NBE
      read(9,*) read_CE
      read(9,*) LDBG
      read(9,*) LSOSCF
      read(9,*) LOCBSE

      close(9)

      if(allocated(KPESTR)) deallocate(KPESTR)
      allocate( KPESTR(nebf),stat=istat )
      if(allocated(KPEEND)) deallocate(KPEEND)
      allocate( KPEEND(nebf),stat=istat )
      call make_KPE(nebf,npebf,AMPEB2C,KPESTR,KPEEND)

      ngee=nebf*nebf*nebf*nebf

      write(*,*)
      write(*,*)'nat   =',nat
      write(*,*)'npebf =',npebf
      write(*,*)'nebf  =',nebf
      write(*,*)'npbf  =',npbf
      write(*,*)'ngee  =',ngee
      write(*,*)'nelec   =',nelec
      write(*,*)'NAE     =',NAE
      write(*,*)'NBE     =',NBE
      write(*,*)'read_CE =',read_CE
      write(*,*)'LDBG    =',LDBG
      write(*,*)'LSOSCF  =',LSOSCF
      write(*,*)'LOCBSE  =',LOCBSE

      write(*,*)
      write(*,*)' CHECK CONTRACTED ELECTRONIC BASIS FUNCTIONS '
      write(*,*)'CONT INDEX    KPESTR     KPEEND'
      do i=1,nebf
         write(*,8000) i,KPESTR(i),KPEEND(i)
      end do

      WRITE(*,*)
      WRITE(*,*)'ELECTRONIC BASIS FUNCTIONS:'
      WRITE(*,*)
      WRITE(*,*)'PRIM  CONT    ANG       EXPONENT CONTRACT  -X- -Y- -Z-'
      WRITE(*,*)'INDEX INDEX   MOM                  COEF'
      DO i=1,npebf
        WRITE(*,9000) i,AMPEB2C(i),ELCAM(i,1),ELCAM(i,2),ELCAM(i,3),
     x ELCEX(i),AGEBFCC(i),ELCBFC(i,1),ELCBFC(i,2),ELCBFC(i,3)
      END DO
      call ELCNORM3(npebf,nebf,
     x              AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

      WRITE(*,*)'ELECTRONIC BASIS FUNCTIONS:'
      WRITE(*,*)'CONTRACT COEFF HAVE BEEN NORMALIZED'
      WRITE(*,*)
      WRITE(*,*)'PRIM  CONT    ANG       EXPONENT CONTRACT  -X- -Y- -Z-'
      WRITE(*,*)'INDEX INDEX   MOM                  COEF'
      DO i=1,npebf
        WRITE(*,9000) i,AMPEB2C(i),ELCAM(i,1),ELCAM(i,2),ELCAM(i,3),
     x ELCEX(i),AGEBFCC(i),ELCBFC(i,1),ELCBFC(i,2),ELCBFC(i,3)
      END DO

!-------READ-INPUT-FILE-AND-ALLOCATE-MEMORY-FOR-BASIS-SET--------------)

      call class_nuc_rep(nat,zan,cat)

      call elec_ovlap(npebf,nebf,nebf*nebf,
     x                AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

      call check_elec_ovlap(nebf)

      call calc_GAM_ecore(nebf,npebf,nebf*nebf,
     x                    nat,zan,cat,
     x                    AMPEB2C,AGEBFCC,
     x                    ELCEX,ELCAM,ELCBFC)

      call calc_GAM_ee(nebf,npebf,ngee,
     x                 AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

      nebf2=nebf*nebf
      if(nelec.gt.1) then
         NPR=(nebf-(nelec/2))*(nelec/2) !Number OCC-VIR PAIRS
      else
         NPR=nelec*(nebf-nelec)
      end if
      if(nae.gt.1) then
         NPRA=(nebf-(nae/2))*(nae/2)
      else
         NPRA=nae*(nebf-nae)
      end if
      if(nbe.gt.1) then
         NPRB=(nebf-(nbe/2))*(nbe/2)
      else
         NPRB=nbe*(nebf-nbe)
      end if
      NEBFLT=nebf*(nebf+1)/2

      wtime1 = omp_get_wtime() - wtime
      wtime  = omp_get_wtime()

      call scf(nelec,nae,nbe,npra,nprb,nebflt,
     x         npebf,nebf,nebf2,ngee,
     x         read_CE,
     x         LSOSCF,LOCBSE,LDBG,
     x         nat,pmass,cat,zan,
     x         KPESTR,KPEEND,AMPEB2C,AGEBFCC,
     x         ELCEX,ELCAM,ELCBFC)

         wtime2 = omp_get_wtime() - wtime

         write(*,3000) wtime1,wtime2

 3000 FORMAT(/8X,'  +--------------------------------------+',/,
     X        8X,'  |    TIMING SUMMARY FOR CALCULATION    |',/,
     x        8X,'  +--------------------------------------+',/,
     x        8X,'    TIME TO EVALUATE INTEGRALS:',1X,F12.4/
     x        8X,'                  TIME FOR SCF:',1X,F12.4/)

 8000 format(/1X,I3,I6,I5)
 9000 format(/1X,I3,I6,I5,I3,I3,F12.6,F10.6,F10.6,F10.6,F10.6)


      return
      end
!=======================================================================
      subroutine mpi_set(nproc,myid,ierr)
! Setup for MPI
!=======================================================================
      implicit none
!     include 'mpif.h'
! Variables Returned
      integer nproc,myid,ierr

!     call MPI_INIT(ierr)
!     call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
!     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

      return
      end
!=======================================================================
      subroutine mpi_final(ierr)

!=======================================================================
      implicit none
!     include 'mpif.h'
! Variables Returned
      integer ierr

!     call MPI_FINALIZE(ierr)

      return
      end

!======================================================================
      subroutine make_KPE(nebf,npebf,AMPEB2C,KPESTR,KPEEND)
 
!     Create a map between contracted EBF index and:
!     1) Beginning primitive index of the contracted shell:  KPESTR()
!     2) Ending primitive index of the contracted shell: KPEEND()
!======================================================================
      implicit none
! Input Variables
      integer nebf,npebf
      integer AMPEB2C(npebf)
! Variables Returned
      integer KPESTR(nebf)
      integer KPEEND(nebf)
! Local Variables
      integer ichk,ist,iend
      integer iep,iec


      ichk=1
      ist=1
      do iep=1,npebf
         iec=AMPEB2C(iep)
         if(iec.ne.ichk) then
            ichk=ichk+1
            iend=iep-1
            KPESTR(iec-1)=ist
            KPEEND(iec-1)=iend
            ist=iep
         end if
         if(iep.eq.npebf) then
            iend=iep
            KPESTR(iec)=ist
            KPEEND(iec)=iend
         end if
      end do


      return
      end

