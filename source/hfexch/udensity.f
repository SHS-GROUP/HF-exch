!=======================================================================
      subroutine construct_uDE(nelec,nebf,CE,DE)

! Form electronic density matrix for single-spin
!=======================================================================
      implicit none

! Input Variables
      integer nelec,nebf
      double precision CE(nebf,nebf)

! Variables Returned
      double precision DE(nebf,nebf)

! Local Variables
      integer i,j,k

      DO i=1,nebf
         DO j=1,nebf

          DE(j,i)=0.0d+00
          do k=1,nelec
            DE(j,i) = DE(j,i) + CE(j,k)*CE(i,k)
          end do

         END DO
      END DO


      return
      end
!=======================================================================
      subroutine read_uCE(nebf,NAalpE,NAbetE,NBalpE,NBbetE,
     x                    DAalpE,DAbetE,DBalpE,DBbetE,
     x                    CAalpE,CAbetE,CBalpE,CBbetE)
!=======================================================================
      implicit none
! Input Variables
      integer nebf
      integer NAalpE,NAbetE
      integer NBalpE,NBbetE
! Variables Returned
      double precision DAalpE(nebf,nebf)
      double precision DAbetE(nebf,nebf)
      double precision DBalpE(nebf,nebf)
      double precision DBbetE(nebf,nebf)
      double precision CAalpE(nebf,nebf)
      double precision CAbetE(nebf,nebf)
      double precision CBalpE(nebf,nebf)
      double precision CBbetE(nebf,nebf)
! Local Variables
      integer ia
      integer ie1,je1
      integer I,J,NMOS,IC,IMIN,IMAX,NUM1,JJ,ICC,NSTM,MODJ,MODIC
      double precision VEC(nebf,nebf)
      integer k
      double precision ans


 9040 FORMAT(I2,I3,5E15.8)
 9060 FORMAT(1X,'*** ERROR IN READ_CAE:   PROBLEM READING ORBITALS!'/
     *       1X,'POSSIBLY A DAMAGED OR MANGLED ORBITAL INPUT GROUP?'/
     *       1X,'ERROR OCCURED AT ORBITAL=',I6,' (MODULUS 100=',I4,'),'/
     *       1X,'         ITS LINE NUMBER=',I6/
     *       1X,'DATA READ FROM INPUT WAS ORBITAL=',I6,' LINE=',I6)

      if (NAalpE.gt.0) then

       open(862,file='guessCAalpE.inp',status='unknown')
       NSTM=0
       NMOS=nebf
       NUM1=nebf

       DO 280 J = 1,NMOS
        IMAX = 0
        IC = 0
  240   CONTINUE
        IMIN = IMAX+1
        IMAX = IMAX+5
        IC = IC+1
        IF(IMAX .GT. NUM1) IMAX = NUM1
        READ(862,9040) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
     *                                    IMAX+NSTM)
!       READ(IR,9040,END=300,ERR=300) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
!    *                                                   IMAX+NSTM)
        MODJ  = MOD(J ,100 )
        MODIC = MOD(IC,1000)
        IF(JJ.EQ.MODJ . AND.  ICC.EQ.MODIC) GO TO 260
           WRITE(*,9060) J,MODJ,IC,JJ,ICC
           STOP
  260   CONTINUE
        IF(IMAX .LT. NUM1) GO TO 240
  280  CONTINUE

       close(862)

       CAalpE=VEC

       DO ie1=1,nebf
       DO je1=1,nebf

        DAalpE(je1,ie1)=0.0d+00

        do k=1,NAalpE
         DAalpE(je1,ie1) = DAalpE(je1,ie1) + CAalpE(je1,k)*CAalpE(ie1,k)
        end do

       END DO
       END DO

      end if

      if (NAbetE.gt.0) then

       open(862,file='guessCAbetE.inp',status='unknown')
       NSTM=0
       NMOS=nebf
       NUM1=nebf

       DO 380 J = 1,NMOS
        IMAX = 0
        IC = 0
  340   CONTINUE
        IMIN = IMAX+1
        IMAX = IMAX+5
        IC = IC+1
        IF(IMAX .GT. NUM1) IMAX = NUM1
        READ(862,9040) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
     *                                    IMAX+NSTM)
!       READ(IR,9040,END=300,ERR=300) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
!    *                                                   IMAX+NSTM)
        MODJ  = MOD(J ,100 )
        MODIC = MOD(IC,1000)
        IF(JJ.EQ.MODJ . AND.  ICC.EQ.MODIC) GO TO 360
           WRITE(*,9060) J,MODJ,IC,JJ,ICC
           STOP
  360   CONTINUE
        IF(IMAX .LT. NUM1) GO TO 340
  380  CONTINUE

       close(862)

       CAbetE=VEC

       DO ie1=1,nebf
       DO je1=1,nebf

        DAbetE(je1,ie1)=0.0d+00

        do k=1,NAalpE
         DAbetE(je1,ie1) = DAbetE(je1,ie1) + CAbetE(je1,k)*CAbetE(ie1,k)
        end do

       END DO
       END DO

      end if

      if (NBalpE.gt.0) then

       open(862,file='guessCBalpE.inp',status='unknown')
       NSTM=0
       NMOS=nebf
       NUM1=nebf

       DO 480 J = 1,NMOS
        IMAX = 0
        IC = 0
  440   CONTINUE
        IMIN = IMAX+1
        IMAX = IMAX+5
        IC = IC+1
        IF(IMAX .GT. NUM1) IMAX = NUM1
        READ(862,9040) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
     *                                    IMAX+NSTM)
!       READ(IR,9040,END=300,ERR=300) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
!    *                                                   IMAX+NSTM)
        MODJ  = MOD(J ,100 )
        MODIC = MOD(IC,1000)
        IF(JJ.EQ.MODJ . AND.  ICC.EQ.MODIC) GO TO 460
           WRITE(*,9060) J,MODJ,IC,JJ,ICC
           STOP
  460   CONTINUE
        IF(IMAX .LT. NUM1) GO TO 440
  480  CONTINUE

       close(862)

       CBalpE=VEC

       DO ie1=1,nebf
       DO je1=1,nebf

        DBalpE(je1,ie1)=0.0d+00

        do k=1,NAalpE
         DBalpE(je1,ie1) = DBalpE(je1,ie1) + CBalpE(je1,k)*CBalpE(ie1,k)
        end do

       END DO
       END DO

      end if

      if (NBbetE.gt.0) then

       open(862,file='guessCBbetE.inp',status='unknown')
       NSTM=0
       NMOS=nebf
       NUM1=nebf

       DO 580 J = 1,NMOS
        IMAX = 0
        IC = 0
  540   CONTINUE
        IMIN = IMAX+1
        IMAX = IMAX+5
        IC = IC+1
        IF(IMAX .GT. NUM1) IMAX = NUM1
        READ(862,9040) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
     *                                    IMAX+NSTM)
!       READ(IR,9040,END=300,ERR=300) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
!    *                                                   IMAX+NSTM)
        MODJ  = MOD(J ,100 )
        MODIC = MOD(IC,1000)
        IF(JJ.EQ.MODJ . AND.  ICC.EQ.MODIC) GO TO 560
           WRITE(*,9060) J,MODJ,IC,JJ,ICC
           STOP
  560   CONTINUE
        IF(IMAX .LT. NUM1) GO TO 540
  580  CONTINUE

       close(862)

       CBbetE=VEC

       DO ie1=1,nebf
       DO je1=1,nebf

        DBbetE(je1,ie1)=0.0d+00

        do k=1,NAalpE
         DBbetE(je1,ie1) = DBbetE(je1,ie1) + CBbetE(je1,k)*CBbetE(ie1,k)
        end do

       END DO
       END DO

      end if

      return
      end
!======================================================================
      subroutine guess_uA_elec(NAE,NBE,nebf,xxse,GAM_ecore,
     x                         DAE,DBE,CAE,CBE)
 
!     Diagonalize the core electron Hamiltonian
!     to construct initial alpha/beta electron guess density
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer NAE,NBE
      double precision xxse(nebf,nebf)
      double precision GAM_ecore(nebf,nebf)
! Variables Returned
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
! Local variables
      double precision CAE(nebf,nebf)
      double precision CBE(nebf,nebf)
      double precision EVF(nebf)


      if (NAE.gt.1) then
       call UROOTHAN(CAE,EVF,xxse,GAM_ecore,nebf)
       call construct_uDE(NAE,nebf,CAE,DAE)
      end if
      if (NBE.gt.1) then
       call UROOTHAN(CBE,EVF,xxse,GAM_ecore,nebf)
       call construct_uDE(NBE,nebf,CBE,DBE)
      end if

      return
      end
!======================================================================
      subroutine guess_uelec(nae,nbe,nebf,xxse,GAM_ecore,
     x                       DAE,DBE,CAE,CBE)
 
!     Diagonalize the core electron Hamiltonian
!     to construct initial regular and special electronic guess density
!     for one particular spin (e.g. alpha B constr to be orth to alpha A)
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nae
      integer nbe
      double precision xxse(nebf,nebf)
      double precision GAM_ecore(nebf,nebf)
! Variables Returned
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision CAE(nebf,nebf)
      double precision CBE(nebf,nebf)
! Local variables
      integer i,j
      double precision C(nebf,nebf)
      double precision EVF(nebf)
      double precision zero
      parameter(zero=0.0d+00)

      DAE=zero
      DBE=zero
      CAE=zero
      CBE=zero

      call UROOTHAN(C,EVF,xxse,GAM_ecore,nebf)

! If nae>0 store first nae evectors in CAE
      if (nae.gt.0) then

       do i=1,nae
         do j=1,nebf
           CAE(j,i)=C(j,i)
         end do
       end do

! nae>0 so if nbe>0 store remaining nebf-nae evectors in CBE
       if (nbe.gt.0) then
        do i=nae+1,nebf
          do j=1,nebf
            CBE(j,i-nae)=C(j,i)
          end do
        end do
       end if

      else

! nae=0 so if nbe>0 store first nbe evectors in CBE
       if (nbe.gt.0) then
        do i=1,nbe
          do j=1,nebf
            CBE(j,i)=C(j,i)
          end do
        end do
       end if

      end if


      call construct_uDE(NAE,nebf,CAE,DAE)
      call construct_uDE(NBE,nebf,CBE,DBE)

      return
      end

