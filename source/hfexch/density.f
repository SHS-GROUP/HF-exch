!=======================================================================
      subroutine construct_DE(nelec,nebf,CE,DE)

! Form electronic density matrix 
!=======================================================================
      implicit none

! Input Variables
      integer nelec,nebf
      double precision CE(nebf,nebf)

! Variables Returned
      double precision DE(nebf,nebf)

! Local Variables
      integer i,j,k
      integer nocc
      double precision coeff
      double precision one,two
      parameter(one=1.0d+00,two=2.0d+00)

      if(nelec.gt.1) then
       coeff=two
       nocc=nelec/2
      else
       coeff=one
       nocc=1
      end if

      DO i=1,nebf
         DO j=1,nebf

          DE(j,i)=0.0d+00
          do k=1,nocc
            DE(j,i) = DE(j,i) + coeff*CE(j,k)*CE(i,k)
          end do

         END DO
      END DO


      return
      end
!=======================================================================
      subroutine read_CAE(nebf,NAE,DAE,C)
!=======================================================================
      implicit none
! Input Variables
      integer nebf,NAE
! Variables Returned
      double precision DAE(nebf,nebf)
! Local Variables
      integer ia
      integer nocc
      integer ie1,je1
      integer I,J,NMOS,IC,IMIN,IMAX,NUM1,JJ,ICC,NSTM,MODJ,MODIC
      integer k
      double precision ans
      double precision coeff
      double precision C(nebf,nebf)
      double precision VEC(nebf,nebf)
      double precision one,two
      parameter(one=1.0d+00,two=2.0d+00)


!--------FORMATTING----------------------------------------------------(
 9040 FORMAT(I2,I3,5E15.8)
 9060 FORMAT(1X,'*** ERROR IN READ_CAE:   PROBLEM READING ORBITALS!'/
     *       1X,'POSSIBLY A DAMAGED OR MANGLED ORBITAL INPUT GROUP?'/
     *       1X,'ERROR OCCURED AT ORBITAL=',I6,' (MODULUS 100=',I4,'),'/
     *       1X,'         ITS LINE NUMBER=',I6/
     *       1X,'DATA READ FROM INPUT WAS ORBITAL=',I6,' LINE=',I6)
!--------FORMATTING----------------------------------------------------)

      open(862,file='guessCAE.inp',status='unknown')

      NSTM=0

      NMOS=nebf
      NUM1=nebf

!-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------(
      DO 280 J = 1,NMOS
         IMAX = 0
         IC = 0
  240    CONTINUE
            IMIN = IMAX+1
            IMAX = IMAX+5
            IC = IC+1
            IF(IMAX .GT. NUM1) IMAX = NUM1
            READ(862,9040) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
     *                                        IMAX+NSTM)
!           READ(IR,9040,END=300,ERR=300) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
!    *                                                       IMAX+NSTM)
            MODJ  = MOD(J ,100 )
            MODIC = MOD(IC,1000)
            IF(JJ.EQ.MODJ . AND.  ICC.EQ.MODIC) GO TO 260
               WRITE(*,9060) J,MODJ,IC,JJ,ICC
               STOP
  260       CONTINUE
         IF(IMAX .LT. NUM1) GO TO 240
  280 CONTINUE
!-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------)

      close(862)

      C=VEC

!-----FORM-DENSITY-MATRIX----------------------------------------------(
      if(NAE.gt.1) then
       coeff=two
       nocc=nae/2
      else
       coeff=one
       nocc=1
      end if

      DO ie1=1,nebf
         DO je1=1,nebf

            DAE(je1,ie1)=0.0d+00

            do k=1,nocc
               DAE(je1,ie1) = DAE(je1,ie1) + coeff*C(je1,k)*C(ie1,k)
            end do

         END DO
      END DO
!-----FORM-DENSITY-MATRIX----------------------------------------------)


      return
      end
!=======================================================================
      subroutine read_CBE(nebf,NBE,DBE,C)
!=======================================================================
      implicit none
! Input Variables
      integer nebf,NBE
! Variables Returned
      double precision DBE(nebf,nebf)
! Local Variables
      integer ia
      integer nocc
      integer ie1,je1
      integer I,J,NMOS,IC,IMIN,IMAX,NUM1,JJ,ICC,NSTM,MODJ,MODIC
      integer K
      double precision ans
      double precision coeff
      double precision C(nebf,nebf)
      double precision VEC(nebf,nebf)
      double precision one,two
      parameter(one=1.0d+00,two=2.0d+00)


!--------FORMATTING----------------------------------------------------(
 9040 FORMAT(I2,I3,5E15.8)
 9060 FORMAT(1X,'*** ERROR IN READ_CBE:   PROBLEM READING ORBITALS!'/
     *       1X,'POSSIBLY A DAMAGED OR MANGLED ORBITAL INPUT GROUP?'/
     *       1X,'ERROR OCCURED AT ORBITAL=',I6,' (MODULUS 100=',I4,'),'/
     *       1X,'         ITS LINE NUMBER=',I6/
     *       1X,'DATA READ FROM INPUT WAS ORBITAL=',I6,' LINE=',I6)
!--------FORMATTING----------------------------------------------------)

      open(863,file='guessCBE.inp',status='unknown')

      NSTM=0

      NMOS=nebf
      NUM1=nebf

!-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------(
      DO 280 J = 1,NMOS
         IMAX = 0
         IC = 0
  240    CONTINUE
            IMIN = IMAX+1
            IMAX = IMAX+5
            IC = IC+1
            IF(IMAX .GT. NUM1) IMAX = NUM1
            READ(863,9040) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
     *                                        IMAX+NSTM)
!           READ(IR,9040,END=300,ERR=300) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
!    *                                                       IMAX+NSTM)
            MODJ  = MOD(J ,100 )
            MODIC = MOD(IC,1000)
            IF(JJ.EQ.MODJ . AND.  ICC.EQ.MODIC) GO TO 260
               WRITE(*,9060) J,MODJ,IC,JJ,ICC
               STOP
  260       CONTINUE
         IF(IMAX .LT. NUM1) GO TO 240
  280 CONTINUE
!-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------)

      close(863)

      C=VEC

!-----FORM-DENSITY-MATRIX----------------------------------------------(
      if(NBE.gt.1) then
       coeff=two
       nocc=nbe/2
      else
       coeff=one
       nocc=1
      end if

      DO ie1=1,nebf
         DO je1=1,nebf

            DBE(je1,ie1)=0.0d+00

            do k=1,nocc
               DBE(je1,ie1) = DBE(je1,ie1) + coeff*C(je1,k)*C(ie1,k)
            end do

         END DO
      END DO
!-----FORM-DENSITY-MATRIX----------------------------------------------)


      return
      end
!======================================================================
      SUBROUTINE guess_A_elec(NAE,nebf,xxse,GAM_ecore,DAE,CAE)
 
!     Diagonalize the core electron Hamiltonian
!     to construct initial regular electron guess density
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer NAE
      double precision xxse(nebf,nebf)
      double precision GAM_ecore(nebf,nebf)
! Variables Returned
      double precision DAE(nebf,nebf)
! Local variables
      double precision CAE(nebf,nebf)
      double precision EVF(nebf)


      call UROOTHAN(CAE,EVF,xxse,GAM_ecore,nebf)

      call construct_DE(NAE,nebf,CAE,DAE)


      RETURN
      END
!======================================================================
      subroutine guess_elec(nae,nbe,nebf,xxse,GAM_ecore,
     x                      DAE,DBE,CAE,CBE)
 
!     Diagonalize the core electron Hamiltonian
!     to construct initial regular and special electronic guess density
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
      integer nocca,noccb
      double precision C(nebf,nebf)
      double precision EVF(nebf)
      double precision zero
      parameter(zero=0.0d+00)

      if (nae.gt.1) then
       nocca=nae/2
      else
       nocca=nae
      end if
      if (nbe.gt.1) then
       noccb=nbe/2
      else
       noccb=nbe
      end if

      DAE=zero
      DBE=zero
      CAE=zero
      CBE=zero

      call UROOTHAN(C,EVF,xxse,GAM_ecore,nebf)

! Store first nocca evectors as occ reg elec vectors
      do i=1,nocca
        do j=1,nebf
          CAE(j,i)=C(j,i)
        end do
      end do

! Store nebf-nocca remaining evectors as occ and virt spec elec vectors
      do i=nocca+1,nebf
        do j=1,nebf
          CBE(j,i-nocca)=C(j,i)
        end do
      end do

      call construct_DE(NAE,nebf,CAE,DAE)
      call construct_DE(NBE,nebf,CBE,DBE)

      return
      end
C======================================================================
      SUBROUTINE DENDIF(D0,D1,NB,DIFF)
C
C     ----- CALCULATE THE LARGEST ABSOLUTE CHANGE IN THE
C           DENSITY MATRIX -----
C======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION D0(NB,NB),D1(NB,NB)
C
      DIFF = 0.0D+00
      DO I=1,NB
         DO J=1,NB
            DIFIJ=ABS(D0(I,J)-D1(I,J))
            IF(DIFIJ.GT.DIFF) DIFF=DIFIJ
         END DO
      END DO
C
      RETURN
      END
C
C======================================================================
      SUBROUTINE COPYDEN(D0,D,N)
C
C     ---> COPY CURRENT DENSITY TO OLD DENSITY
C
C     D0 :: OLD DENSITY
C     D  :: CURRENT DENSITY
C     N  :: NUMBER OF BASIS FUNCTIONS
C
C======================================================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D0(N,N),D(N,N)

      DO I=1,N
         DO J=1,N
            D0(I,J)=D(I,J)
         END DO
      END DO
C
      RETURN
      END

