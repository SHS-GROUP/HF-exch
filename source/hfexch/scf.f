!======================================================================
      subroutine scf(nelec,NAE,NBE,NPRA,NPRB,NEBFLT,
     x               npebf,nebf,nebf2,ngee,
     x               read_CE,
     x               LHF,LSOSCF,LOCBSE,LDBG,
     x               nat,cat,zan,
     x               KPESTR,KPEEND,AMPEB2C,AGEBFCC,
     x               ELCEX,ELCAM,ELCBFC,GAM_ee)

!
! PERFORM A NUCLEAR-ELECTRONIC RESTRICTED XC HARTREE-FOCK CALCULATION
! FOR AN NAE-REGULAR ELECTRON NBE-SPECIAL ELECTRON ONE-PROTON SYSTEM
!
!     **DEFINITIONS:
!
!     *FOR REGULAR ELECTRONS:
!     NAE    ::  NUMBER OF REGULAR ELECTRONS
!     DAE    ::  NEW REGULAR ELECTRON DENSITY MATRIX
!     DAE0   ::  OLD REGULAR ELECTRON DENSITY MATRIX
!     VECAE  ::  REGULAR ELECTRON MOS
!     VECAE0 ::  OLD REGULAR ELECTRON MOS
!     AEE    ::  REGULAR ELECTRON ORBITAL EIGENVALUES
!
!     * FOR SPECIAL ELECTRONS:
!     NBE    ::  NUMBER OF SPECIAL ELECTRONS
!     DBE    ::  NEW SPECIAL ELECTRON DENSITY MATRIX
!     DBE0   ::  OLD SPECIAL ELECTRON DENSITY MATRIX
!     VECBE  ::  SPECIAL ELECTRON MOS
!     VECBE0 ::  OLD SPECIAL ELECTRON MOS
!     BEE    ::  SPECIAL ELECTRON ORBITAL EIGENVALUES
!
!     * FOR HF CALCULATION:
!     USE AE VARIABLES
!
!======================================================================
      implicit none
! Input Variables
      logical LOCBSE   ! Use OCBSE scheme as is (restricted variational freedom for reg/sp elecs)
      logical LOCBSE2  ! Use modified OCBSE scheme (complete variational freedom for reg elecs)
      logical read_CE
      logical LHF
      logical LDBG
      integer nelec
      integer NAE,NBE
      integer NPRA,NPRB
      integer NEBFLT
      integer nebf
      integer nebf2
      integer ngee
!-----DIRECT-SCF-RELATED-----------------------------------------------(
      integer npebf
      integer nat
!-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      integer KPESTR(nebf)  ! Map contracted index to primitive start
      integer KPEEND(nebf)  ! Map contracted index to primitive end
!-------Basis Set Info-------)
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
!-----DIRECT-SCF-RELATED-----------------------------------------------)

! Local variables
      double precision zero,one
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00) 

      double precision xxse(nebf,nebf)  ! Elec overlap matrix
      double precision GAM_ecore(nebf2)
      double precision GAM_ee(ngee)

      integer maxit

      integer i,j,k,l

      double precision TOLE
      double precision diffE
      double precision diffAE
      double precision diffBE

      double precision E_total

      double precision E_A
      double precision E_A_ecore
      double precision E_A_ee

      double precision E_B
      double precision E_B_ecore
      double precision E_B_ee

      double precision E_AB

      double precision E_nuc

      double precision DAE(NEBF,NEBF)
      double precision DAE0(NEBF,NEBF)
      double precision DBE(NEBF,NEBF)
      double precision DBE0(NEBF,NEBF)
      double precision VECAE(NEBF,NEBF)
      double precision VECAE0(NEBF,NEBF)
      double precision VECBE(NEBF,NEBF)
      double precision VECBE0(NEBF,NEBF)
      double precision AEE(NEBF)
      double precision BEE(NEBF)
      double precision FAE(nebf,nebf)
      double precision XFAE(nebf,nebf)
      double precision FBE(nebf,nebf)
      double precision XFBE(nebf,nebf)

      double precision FAEint(nebf,nebf)
      double precision FBEint(nebf,nebf)

      double precision E_total_old
      double precision Delta_E_tot

      logical LDIFFE

!--------SOSCF-RELATED-VARIABLES------------(
      logical LSOSCF
      logical EIGAVL
      integer NA
      integer ITER
      integer ITSOA ! SOSCF iteration counter
      integer ITSOB ! SOSCF iteration counter
      integer L0,L1
      integer NFT15
      integer NFT16
      double precision FLT(NEBFLT) !FLT: Lower triangle focke
      double precision HSTARTA(NPRA)
C      double precision HSTARTB(NPRB)
      double precision GRADA(NPRA)
C      double precision GRADB(NPRB)
      double precision PGRADA(NPRA)
C      double precision PGRADB(NPRB)
      double precision DISPLIA(NPRA)
C      double precision DISPLIB(NPRB)
      double precision DGRADA(NPRA)  ! WRK1
C      double precision DGRADB(NPRB)  ! WRK1
      double precision DISPLA(NPRA)  ! WRK2
C      double precision DISPLB(NPRB)  ! WRK2
      double precision UPDTA(NPRA)   ! WRK3
C      double precision UPDTB(NPRB)   ! WRK3
      double precision DISPLNA(NPRA) ! WRK1+NPR
C      double precision DISPLNB(NPRB) ! WRK1+NPR
      double precision DGRADIA(NPRA) ! WRK2+NPR
C      double precision DGRADIB(NPRB) ! WRK2+NPR
      double precision UPDTIA(NPRA)  ! WRK3+NPR
C      double precision UPDTIB(NPRB)  ! WRK3+NPR
      double precision ORBGRDA
C      double precision ORBGRDB
      double precision SMALL
      double precision SOGTOL ! ORBGRAD TOL to activate soscf
      double precision XA(NPRA)
C      double precision XB(NPRB)
      double precision GA(nebf,nebf) !G(L0,L0)
      double precision GB(nebf,nebf) !G(L0,L0)
      double precision WRK(nebf) !WRK(L0)
!cc   double precision CCC(nebf,nebf) !WRK(L0)
!cc   NPR=(L0-NA)*NA ! Line 2134 RHFCL ?NA is NUM ALPHA E?
!--------SOSCF-RELATED-VARIABLES------------)

!--------OUTPUT-FORMATTING---------------------------------------------(
 9000 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS ')

 9050 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        ORBGRAD_A ')

 9100 FORMAT(1X,I3,F20.10,F17.10,2F17.10)

 9150 FORMAT(1X,I3,F20.10,F17.10,3F17.10)

 9200 FORMAT(/1X,'FINAL HF ENERGY IS',F20.10,' AFTER',I4,
     *           ' ITERATIONS')

 9300 FORMAT(/6X,'-----------------------------------------------',/,
     x        6X,'         HF(exch) ENERGETIC COMPONENTS         ',/,
     x        6X,'-----------------------------------------------',/,
     x       12X,'    E_NUC=',1X,F20.10/
     x       12X,'---------------------------------',/,
     x       12X,' E_A_CORE=',1X,F20.10/
     x       12X,'   E_A_EE=',1X,F20.10/
     x       12X,'      E_A=',1X,F20.10/
     x       12X,'---------------------------------',/,
     x       12X,' E_B_CORE=',1X,F20.10/
     x       12X,'   E_B_EE=',1X,F20.10/
     x       12X,'      E_B=',1X,F20.10/
     x       12X,'---------------------------------',/,
     x       12X,'     E_AB=',1X,F20.10/
     x       12X,'---------------------------------',/,
     x       12X,'  E_TOTAL=',1X,F20.10/
     x        6X,'-----------------------------------------------',/)

 9400 FORMAT(/1X,'          INITIAL GUESS ENERGETICS:            ')

 9500 FORMAT(/6X,' ** BEGIN SELF-CONSISTENT-FIELD CALCULATION **')

 9610 FORMAT(/1X,' A ELECTRONIC ORBITALS AND EIGENVALUES:       ')

 9620 FORMAT(/1X,' B ELECTRONIC ORBITALS AND EIGENVALUES:       ')

 9800 FORMAT(10X,15(1H-),'START SECOND ORDER SCF',15(1H-))
                                           
!--------OUTPUT-FORMATTING---------------------------------------------)

C      LOCBSE2=LOCBSE
      LOCBSE2=.false.
      if(LOCBSE2) then
       LOCBSE=.false.
      end if
      if (LHF) then
       write(*,*) "Performing HF calculation"
       LOCBSE=.false.
       LOCBSE2=.false.
       nae=nelec
      end if
      if(LOCBSE) write(*,*) "Using LOCBSE"
      if(LOCBSE2) write(*,*) "Using LOCBSE2"

      DAE=0.0d+00
      DBE=0.0d+00
      DAE0=0.0d+00
      DBE0=0.0d+00
      vecAE=0.0d+00
      vecBE=0.0d+00
      vecAE0=0.0d+00
      vecBE0=0.0d+00

!----------CALCULATE-CLASSICAL-NUCLEAR-REPULSION-ENERGY----------------(
!      call class_nuc_rep(nat,zan,cat,E_nuc)
      open(800,file='ENUCRP.dat',status='unknown')
      read(800,*) E_nuc
      close(800)
      write(*,*)'READ IN NUC REPULSION'
!----------CALCULATE-CLASSICAL-NUCLEAR-REPULSION-ENERGY----------------)

!--------------READ-INTEGRALS-NEEDED-FOR-NEO-HF------------------------(
      nebf2=nebf*nebf
      write(*,*)
      call read_elec_ovlap(nebf,xxse)
      write(*,*)'READ IN ELEC OVLAP'
      call read_GAM_ecore(nebf,nebf2,GAM_ecore)
      write(*,*)'READ IN GAM_ECORE'
C      call read_GAM_ee(nebf,ngee,GAM_ee)
C      write(*,*)'READ IN GAM_EE'
      write(*,*)
!--------------READ-INTEGRALS-NEEDED-FOR-NEO-HF------------------------)

!-------------INITIAL-GUESSES------------------------------------------(
      if(read_CE) then
!        READ IN GUESS FOR E:
         call read_CAE(nebf,NAE,DAE,VECAE0)
         if(.not.LHF) call read_CBE(nebf,NBE,DBE,VECBE0)
      else
!       STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
        write(*,*)'ABOUT TO CALL guess_A_elec'
!       call guess_elec(nelec,nebf,xxse,GAM_ecore,DE)
        if ((LOCBSE).or.(LOCBSE2)) then
          call guess_elec(nae,nbe,nebf,xxse,GAM_ecore,
     x                          DAE,DBE,VECAE0,VECBE0)
          write(*,*)'BACK FROM guess_elec for OCBSE'
        else
         call guess_A_elec(NAE,nebf,xxse,GAM_ecore,DAE,VECAE0)
         if(.not.LHF) then
          call guess_A_elec(NBE,nebf,xxse,GAM_ecore,DBE,VECBE0)
         end if
         write(*,*)'BACK FROM guess_elec'
         write(*,*)
        end if
      end if

C ARS( debug: print out initial guess MOs here
      if (LDBG) then
       write(*,*)
       write(*,*) "------------------"
       write(*,*) "INITIAL GUESS MOs:"
       write(*,*) "------------------"
       write(*,*)
       WRITE(*,9610)
       call PREVNU(vecAE0,AEE,nebf,nebf,nebf)
       WRITE(*,9620)
       call PREVNU(vecBE0,BEE,nebf,nebf,nebf)
       write(*,*)
      end if
C )
!-------------INITIAL-GUESSES------------------------------------------)

!-------------SETUP-FOR-POSSIBLE-SOSCF---------------------------------(
!     LSOSCF=.FALSE.
      if((nae.eq.1).or.LOCBSE) then
         LSOSCF=.FALSE.
      end if
      if(LSOSCF) THEN
         NFT15=15
         OPEN(NFT15, FILE='WORK15', STATUS='UNKNOWN',
     *        ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         SOGTOL=1.0d+00
         SMALL=1.0D-06
         ITSOA=0
         L0=nebf
         L1=nebf
         NA=nae/2
      end if
!-------------SETUP-FOR-POSSIBLE-SOSCF---------------------------------)

!
!     SET CONVERGENCE CRITERIA AND MAXIMUM ITERATIONS 
!
      TOLE = 1.0D-06
      maxit=100
      if(LOCBSE) maxit=400
!
!     BEGIN XCSCF ITERATIONS
      WRITE(*,9500)
!     WRITE(*,9000)
!
      E_total=0.0d+00
      E_A_ecore=0.0d+00
      E_A_ee=0.0d+00
      E_A=0.0d+00
      E_B_ecore=0.0d+00
      E_B_ee=0.0d+00
      E_B=0.0d+00
      E_AB=0.0d+00

      E_total_old=0.0d+00
      ORBGRDA=0.0d+00
C      ORBGRDB=0.0d+00

      DO I=1,MAXIT

!--------------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------(
!     write(*,*)
!     write(*,*)'IN: xcuscf '
!     write(*,*)'Before call to UHF_FOCK, GAM_EE='
!     write(*,*)GAM_ee
!     write(*,*)
       if((.not.((locbse).or.(locbse2))).or.(i.eq.1)) then

C Call HF Fock build for NAE electrons
         call fock_hf(LDBG,nebf,nebf2,NAE,ngee,
     x                DAE,GAM_ecore,GAM_ee,
     x                FAE,E_A,E_A_ecore,E_A_ee)

         if (.not.LHF) then
C Call HF Fock build for NBE electrons
          call fock_hf(LDBG,nebf,nebf2,NBE,ngee,
     x                 DBE,GAM_ecore,GAM_ee,
     x                 FBE,E_B,E_B_ecore,E_B_ee)

C Call interaction Fock build for all particles
          call fock_int(LDBG,nelec,NAE,NBE,
     x                  nebf,nebf2,ngee,
     x                  DAE,DBE,GAM_ee,
     x                  FAEint,FBEint, 
     x                  E_AB)

          call add2fock(nebf,FAEint,FAE)
          call add2fock(nebf,FBEint,FBE)
         end if

          IF (LDBG) then
           write(*,*)
           write(*,*) "FAE:"
           call prt_lower_triangle(nebf,nebflt,FAE)
           write(*,*)
           if (.not.LHF) then
            write(*,*) "FBE:"
            call prt_lower_triangle(nebf,nebflt,FBE)
           end if
           write(*,*)
          END IF  

          E_total=E_A+E_B+E_AB+E_nuc

       end if

!--------------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------)
         if(I.eq.1) then
            WRITE(*,9400)

      WRITE(*,9300) E_nuc,
     x              E_A_ecore,E_A_ee,E_A,
     x              E_B_ecore,E_B_ee,E_B,
     x              E_AB,
     x              E_total

            if(LSOSCF) then 
               WRITE(*,9050)
            else
               WRITE(*,9000)
            end if
         end if

         if (LOCBSE) then
! Do OCBSE procedure (restricted solutions for regular and special electrons)

           call OCBSE(nebf,nae,nbe,vecAE0,vecBE0,FAE,FBE,xxse,
     x                vecAE,vecBE,AEe,BEe)

! Form regular electronic density matrix and store stuff for next it
           call construct_DE(NAE,nebf,vecAE,DAE)
           CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
           CALL COPYDEN(DAE0,DAE,NEBF)
           CALL COPYDEN(vecAE0,vecAE,NEBF)

! Form special electronic density matrix and store stuff for next it
           call construct_DE(NBE,nebf,vecBE,DBE)
           CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
           CALL COPYDEN(DBE0,DBE,NEBF)
           CALL COPYDEN(vecBE0,vecBE,NEBF)

! Calculate energy for this it and Fock matrices for next it

C Call HF Fock build for NAE electrons
           call fock_hf(LDBG,nebf,nebf2,NAE,ngee,
     x                  DAE,GAM_ecore,GAM_ee,
     x                  FAE,E_A,E_A_ecore,E_A_ee)

C Call HF Fock build for NBE electrons
           call fock_hf(LDBG,nebf,nebf2,NBE,ngee,
     x                  DBE,GAM_ecore,GAM_ee,
     x                  FBE,E_B,E_B_ecore,E_B_ee)

C Call interaction Fock build for all particles
           call fock_int(LDBG,nelec,NAE,NBE,
     x                   nebf,nebf2,ngee,
     x                   DAE,DBE,GAM_ee,
     x                   FAEint,FBEint, 
     x                   E_AB)

            call add2fock(nebf,FAEint,FAE)
            call add2fock(nebf,FBEint,FBE)

            IF (LDBG) then
             write(*,*)
             write(*,*) "FAE:"
             call prt_lower_triangle(nebf,nebflt,FAE)
             write(*,*)
             write(*,*) "FBE:"
             call prt_lower_triangle(nebf,nebflt,FBE)
             write(*,*)
            END IF  

            E_total=E_A+E_B+E_AB+E_nuc

         else if (LOCBSE2) then
! Do OCBSE2 procedure (restricted solutions for special electrons)

! Regular electrons
         if(LSOSCF) THEN
          ITER=I
          EIGAVL = ITER.GT.1
         end if
         IF(LSOSCF .AND.  EIGAVL) THEN
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAE,FLT)
          call SOGRAD(GRADA,FLT,vecAE,WRK,NPRA,NA,L0,L1,NEBFLT,ORBGRDA)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 900  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDA.LT.SOGTOL  .OR.  ITSOA.GT.0) THEN
              IF(ITSOA.EQ.0) THEN
              WRITE(*,9800)
                 call SOHESS(HSTARTA,AEE,NPRA,L0,NA,NA)
              END IF
              ITSOA = ITSOA+1
           call SONEWT(HSTARTA,GRADA,PGRADA,DISPLIA,DGRADA,DISPLA,UPDTA,
     *                 DISPLNA,DGRADIA,UPDTIA,ORBGRDA,NPRA,ITSOA,NFT15)
            call SOTRAN(DISPLIA,vecAE,GA,WRK,NPRA,L0,L1,NA,NA,ORBGRDA)
             CALL DCOPY(NPRA,GRADA,1,PGRADA,1)
              call construct_DE(NAE,nebf,vecAE,DAE)
              GO TO 950  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------)

  900 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAE,AEE,xxse,FAE,nebf)
         call construct_DE(NAE,nebf,vecAE,DAE)

  950 CONTINUE
!        --> FIND LARGEST CHANGE IN Alpha E DENSITY
         CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
         CALL COPYDEN(DAE0,DAE,NEBF)

! Do OCBSE procedure for special electrons
           call OCBSE2(nebf,npebf,nae,nbe,vecAE,vecBE0,FBE,
     x                 xxse,elcam,elcbfc,elcex,
     x                 ampeb2c,agebfcc,kpestr,kpeend,
     x                 vecBE,BEe)

! Form special electronic density matrix and store stuff for next it
           call construct_DE(NBE,nebf,vecBE,DBE)
           CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
           CALL COPYDEN(DBE0,DBE,NEBF)
           CALL COPYDEN(vecAE0,vecAE,NEBF)
           CALL COPYDEN(vecBE0,vecBE,NEBF)

! Calculate energy for this it and Fock matrices for next it

C Call HF Fock build for NAE electrons
           call fock_hf(LDBG,nebf,nebf2,NAE,ngee,
     x                  DAE,GAM_ecore,GAM_ee,
     x                  FAE,E_A,E_A_ecore,E_A_ee)

C Call HF Fock build for NBE electrons
           call fock_hf(LDBG,nebf,nebf2,NBE,ngee,
     x                  DBE,GAM_ecore,GAM_ee,
     x                  FBE,E_B,E_B_ecore,E_B_ee)

C Call interaction Fock build for all particles
           call fock_int(LDBG,nelec,NAE,NBE,
     x                   nebf,nebf2,ngee,
     x                   DAE,DBE,GAM_ee,
     x                   FAEint,FBEint, 
     x                   E_AB)

            call add2fock(nebf,FAEint,FAE)
            call add2fock(nebf,FBEint,FBE)

            IF (LDBG) then
             write(*,*)
             write(*,*) "FAE:"
             call prt_lower_triangle(nebf,nebflt,FAE)
             write(*,*)
             write(*,*) "FBE:"
             call prt_lower_triangle(nebf,nebflt,FBE)
             write(*,*)
            END IF  

            E_total=E_A+E_B+E_AB+E_nuc

         else

!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------(
         if(LSOSCF) THEN
          ITER=I
          EIGAVL = ITER.GT.1
         end if
         IF(LSOSCF .AND.  EIGAVL) THEN
          write(*,*) "nae:",nae
          write(*,*) "na:",na
          write(*,*) "npra:",npra
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAE,FLT)
          call SOGRAD(GRADA,FLT,vecAE,WRK,NPRA,NA,L0,L1,NEBFLT,ORBGRDA)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 700  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDA.LT.SOGTOL  .OR.  ITSOA.GT.0) THEN
              IF(ITSOA.EQ.0) THEN
              WRITE(*,9800)
                 call SOHESS(HSTARTA,AEE,NPRA,L0,NA,NA)
              END IF
              ITSOA = ITSOA+1
           call SONEWT(HSTARTA,GRADA,PGRADA,DISPLIA,DGRADA,DISPLA,UPDTA,
     *                 DISPLNA,DGRADIA,UPDTIA,ORBGRDA,NPRA,ITSOA,NFT15)
            call SOTRAN(DISPLIA,vecAE,GA,WRK,NPRA,L0,L1,NA,NA,ORBGRDA)
             CALL DCOPY(NPRA,GRADA,1,PGRADA,1)
              call construct_DE(NAE,nebf,vecAE,DAE)
              GO TO 750  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------)

  700 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAE,AEE,xxse,FAE,nebf)
         call construct_DE(NAE,nebf,vecAE,DAE)

  750 CONTINUE
!        --> FIND LARGEST CHANGE IN Alpha E DENSITY
         CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
         CALL COPYDEN(DAE0,DAE,NEBF)

!-----------------------POSSIBLE-SOSCF-BETA----------------------------(
!        if(LSOSCF) THEN
!        ITER=I
!        EIGAVL = ITER.GT.1
C         IF(LSOSCF .AND.  EIGAVL) THEN
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
C           call pack_LT(nebf,nebfLT,FBE,FLT)
C          call SOGRAD(GRADB,FLT,vecBE,WRK,NPRB,NBE,L0,L1,NEBFLT,ORBGRDB)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 700  ! Check on convergence behavior
!!!!!!      END IF
C            IF(ORBGRDB.LT.SOGTOL  .OR.  ITSOB.GT.0) THEN
C              IF(ITSOB.EQ.0) THEN
!             WRITE(*,9800)
C                 call SOHESS(HSTARTB,BEE,NPRB,L0,NBE,NBE)
C              END IF
C              ITSOB = ITSOB+1
C           call SONEWT(HSTARTB,GRADB,PGRADB,DISPLIB,DGRADB,DISPLB,UPDTB,
C     *                 DISPLNB,DGRADIB,UPDTIB,ORBGRDB,NPRB,ITSOB,NFT16)
C            call SOTRAN(DISPLIB,vecBE,GB,WRK,NPRB,L0,L1,NBE,NBE,ORBGRDB)
C             CALL DCOPY(NPRB,GRADB,1,PGRADB,1)
C              call construct_DAE(NBE,nebf,vecBE,DBE)
C              GO TO 850  ! Use the new C's to form new density (change)
C            END IF
C         END IF
!-----------------------POSSIBLE-SOSCF-BETA----------------------------)

  800 CONTINUE
!        call ROOTHAN(DBE,vecBE,BEE,xxse,FBE,nebf,nelec,1,NUCST)
        if (.not.LHF) then
         call UROOTHAN(vecBE,BEE,xxse,FBE,nebf)
         call construct_DE(NBE,nebf,vecBE,DBE)
        end if

  850 CONTINUE
!        --> FIND LARGEST CHANGE IN Beta E DENSITY
         CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
         CALL COPYDEN(DBE0,DBE,NEBF)

         end if ! end if for not ocbse or ocbse2

!        --> CALCULATE CHANGE IN TOTAL ENERGY
         Delta_E_tot=E_total-E_total_old
         E_total_old=E_total

!        --> PRINT SUMMARY OF THIS ITERATION
         if(LSOSCF) then
            WRITE(*,9150) I,E_total,Delta_E_tot,
     x                    DIFFAE,DIFFBE,ORBGRDA
         else
            WRITE(*,9150) I,E_total,Delta_E_tot,
     x                    DIFFAE,DIFFBE
         end if
C ARS( debug: print out MOs here
      if (LDBG) then
       WRITE(*,9610)
       call PREVNU(vecAE,AEE,nebf,nebf,nebf)
       WRITE(*,9620)
       call PREVNU(vecBE,BEE,nebf,nebf,nebf)
      end if
C )
! Output the vectors for this iteration for restart if necessary:
         call write_MOs(860,nebf,VECAE)
         call write_MOs(861,nebf,VECBE)

         LDIFFE=( (DIFFAE.LT.TOLE).and.(DIFFBE.LT.TOLE) )
         IF(LDIFFE) GOTO 100
         IF(I.EQ.MAXIT) GOTO 10

      END DO
 
  10  CONTINUE
!     IF WE GET HERE SOMETHING WENT WRONG

      if(LSOSCF) THEN
         close(NFT15)
      end if

      WRITE(*,*)
      WRITE(*,*)'WARNING:  ITERATION LIMIT EXCEEDED'
      E_total=zero
      WRITE(*,*)
!     STOP
!
  100 CONTINUE
!     IF WE GET HERE WE ARE DONE - CONVERGENCE ACHIEVED

      if(LSOSCF) THEN
         close(NFT15)
      end if

!     PRINT FINAL ENERGY AND PUNCH THE ORBITALS
      WRITE(*,9200) E_total,I

      WRITE(*,9300) E_nuc,
     x              E_A_ecore,E_A_ee,E_A,
     x              E_B_ecore,E_B_ee,E_B,
     x              E_AB,
     x              E_total

!  OUTPUT ELEC AND NUC EIGENVALUES AND EIGENVECTORS
      WRITE(*,9610)
      call PREVNU(vecAE,AEE,nebf,nebf,nebf)
      WRITE(*,9620)
      call PREVNU(vecBE,BEE,nebf,nebf,nebf)

! PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------(
!     subroutine write_MOs(IFIL,nbf,VEC)
!     IFIL=852 :: FinalCE.dat
!     IFIL=853 :: FinalCP.dat
!     IFIL=860 :: FinalCAE.dat
!     IFIL=861 :: FinalCBE.dat
      call write_MOs(860,nebf,VECAE)
      call write_MOs(861,nebf,VECBE)
! PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------)
!

      RETURN
      END
C======================================================================
      SUBROUTINE UROOTHAN(C,EVF,S,F,NB)
C
C     SOLVES THE ROOTHAN EQUATIONS:
C     ==============================
C     ON INPUT:
C     S         ::  OVERLAP MATRIX
C     F         ::  FOCK MATRIX
C     NB        ::  NUMBER OF BASIS FUNCTIONS
C     ==============================
C     WORKING:
C     EV        ::  EIGENVALUES OF S MATRIX
C     EVECS     ::  EIGENVECTORS OF S MATRIX
C     EVECST    ::  TRANSPOSE OF EIGENVECTORS OF S MATRIX
C     X         ::  TRANSFORMATION MATRIX
C     XP        ::  TRANSPOSE OF TRANFORMATION MATRIX
C     FP        ::  TRANSFORMED FOCK MATRIX
C     CP        ::  EIGENVECTORS FROM FOCK MATRIX
C     FV1       ::  WORK SPACE
C     FV2       ::  WORK SPACE
C     FV3       ::  WORK SPACE
C     FV4       ::  WORK SPACE
C     ==============================
C     OUTPUT:
C     C         ::  TRANSFORMED EIGENVECTORS
C     EVF       ::  EIGENVALUES FROM FOCK MATRIX
C     ============================== 
C
C======================================================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00)
      DIMENSION S(NB,NB),F(NB,NB),FP(NB,NB),EVF(NB),X(NB,NB),XP(NB,NB),
     *          CP(NB,NB),C(NB,NB),D(NB,NB),EVS(NB),EVECS(NB,NB),
     *          FV1(NB),FV2(NB),FV3(NB,NB),FV4(NB,NB),DEVS(NB,NB),
     *          EVECST(NB,NB),xmat(NB,NB)
      LOGICAL DEBUG
C
C     CONSTRUCT ORTHONORMALIZING TRANSFORMATION MATRIX
C
C     ---> DIAGONALIZE OVERLAP MATRIX
C
C     DEBUG=.TRUE.
      DEBUG=.FALSE.
      CALL RS(NB,NB,S,EVS,2,EVECS,FV1,FV2,IERR)
      IF(DEBUG) THEN
         WRITE(*,*)
         WRITE(*,*)'---- OVERLAP MATRIX --> EIGENVALUES ----'
         WRITE(*,*)
         DO I=1,NB
               WRITE(*,*)'EVS(',I,')=',EVS(I)
         END DO
         WRITE(*,*)
         WRITE(*,*)'---- OVERLAP MATRIX --> EIGENVECTORS ----'
         WRITE(*,*)
         DO I=1,NB
            DO J=1,NB
               WRITE(*,*)'EVECS(',I,J,')=',EVECS(I,J)
            END DO
         END DO
      END IF
C           
C     ---- SYMMETRIC ORTHOGONALIZATION SCHEME ----
C
C     ---> FORM DIAGONAL EIGENVALUE MATRIX
C     ---> GET TRANSPOSE OF VECTORS FROM DIAGONALIZATION OF S
C
      DO I=1,NB
         DO J=1,NB
            DEVS(I,J)=ZERO
            EVECST(I,J)=EVECS(J,I)
         END DO
         DEVS(I,I)=1/SQRT(EVS(I))
      END DO
      IF(DEBUG) THEN
         WRITE(*,*)
         WRITE(*,*)'---- S^-(1/2) EIGENVALUE MATRIX ----'
         WRITE(*,*)
         DO I=1,NB
            DO J=1,NB
               WRITE(*,*)'DEVS(',I,J,')=',DEVS(I,J)
            END DO
         END DO
      END IF
C
C     ---> FORM TRANSFORMATION MATRIX 
C
C     MULTIPLY EVECS.DEVS.EVECST = X
C     
      CALL MATMULT(NB,NB,NB,NB,DEVS,EVECST,FV3)
      CALL MATMULT(NB,NB,NB,NB,EVECS,FV3,X)
      IF(DEBUG) THEN
         WRITE(*,*)
         WRITE(*,*)'---- X MATRIX ----'
         WRITE(*,*)
         DO I=1,NB
            DO J=1,NB
               WRITE(*,*)'X(',I,J,')=',X(I,J)
            END DO
         END DO
      END IF
C
C     FORM TRANSPOSE OF TRANSFORMATION MATRIX
C
      DO I=1,NB
         DO J=1,NB
            XP(I,J)=X(J,I)
         END DO
      END DO
C
C     TRANSFORM FOCK MATRIX 
C 
C     MULTIPLY XP.F.X = FP
C
      CALL MATMULT(NB,NB,NB,NB,F,X,FV4)
      CALL MATMULT(NB,NB,NB,NB,XP,FV4,FP)
C
C     DIAGONALIZE FOCK MATRIX 
C
      CALL RS(NB,NB,FP,EVF,2,CP,FV1,FV2,IERR)
C
C     TRANSFORM MO COEFFICIENTS
C      C = X.CP  
C   
      CALL MATMULT(NB,NB,NB,NB,X,CP,C)

CDDD
c     write(*,*) 
c     write(*,*)'--------ORTHO-CHEC----------' 
c     do i=1,NB
c     do j=1,NB
c     xmat(i,j)=0.0d0
c     do k1=1,NB
c     do k2=1,NB
c     xmat(i,j)=xmat(i,j)+(C(k1,i)*S(k1,k2)*C(k2,j))
c     end do
c     end do
c     write(*,*) i,j,xmat(i,j)
c     end do
c     end do
c     write(*,*)'--------END OF ORTHO-CHEC----------' 
C
      RETURN
      END
!======================================================================
      subroutine OCBSE(nebf,nae,nbe,vecAE0,vecBE0,FAE,FBE,
     x                 Selec,vecAE,vecBE,AEen,BEen)
! 
!     Perform OCBSE procedure
!       vecAE0: Regular electron coefficients from previous iteration
!       vecBE0: Special electron coefficients from previous iteration
!       vecAE:  Regular electron coefficients from current iteration
!       vecBE:  Special electron coefficients from current iteration
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nae,nbe
      double precision vecAE0(nebf,nebf),vecBE0(nebf,nebf)
      double precision FAE(nebf,nebf),FBE(nebf,nebf)
      double precision Selec(nebf,nebf)
! Variables Returned
      double precision vecAE(nebf,nebf),vecBE(nebf,nebf)
      double precision AEen(nebf),BEen(nebf)
! Local variables
      integer nocca,noccb,nvirt 
      integer noccvirta,noccvirtb
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
      nvirt=nebf-nocca-noccb
      noccvirta=nebf-noccb
      noccvirtb=nebf-nocca

      vecAE=zero
      vecBE=zero
      AEen=zero
      BEen=zero

      call OCBSE_driver(nebf,nae,nbe,nocca,noccb,
     x                  nvirt,noccvirta,noccvirtb,
     x                  vecAE0,vecBE0,FAE,FBE,Selec,
     x                  vecAE,vecBE,AEen,BEen)

      return
      end
!======================================================================
      subroutine OCBSE_driver(nebf,nae,nbe,nocca,noccb,
     x                        nvirt,noccvirta,noccvirtb,
     x                        vecAE0,vecBE0,FAE,FBE,Selec,
     x                        vecAE,vecBE,AEen,BEen)
!
!     OCBSE Procedure:
!       - Construct transformation matrix for regular electrons
!       - Transform FAE and diagonalize to obtain vecAE
!       - Construct transformation matrix for special electrons
!       - Transform FBE and diagonalize to obtain vecBE
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nae,nbe
      double precision vecAE0(nebf,nebf),vecBE0(nebf,nebf)
      double precision FAE(nebf,nebf),FBE(nebf,nebf)
      double precision Selec(nebf,nebf)
! Variables Returned
      double precision vecAE(nebf,nebf),vecBE(nebf,nebf)
      double precision AEen(nebf),BEen(nebf)
! Local variables
      integer i,j
      integer nocca,noccb,nvirt
      integer noccvirta,noccvirtb
      double precision WA(nebf,noccvirta),WB(nebf,noccvirtb)
      double precision WAtrans(noccvirta,nebf),WBtrans(noccvirtb,nebf)
      double precision wFAEw(noccvirta,noccvirta)
      double precision wFBEw(noccvirtb,noccvirtb)
      double precision wSAw(noccvirta,noccvirta)
      double precision wSBw(noccvirtb,noccvirtb)
      double precision AUXA(nebf,noccvirta)
      double precision AUXB(nebf,noccvirtb)
      double precision xAEen(noccvirta)
      double precision xBEen(noccvirtb)
      double precision xvecAE(noccvirta,noccvirta)
      double precision xvecBE(noccvirtb,noccvirtb)
      double precision blockvecAE(nebf,noccvirta)
      double precision blockvecBE(nebf,noccvirtb)
      double precision zero
      parameter(zero=0.0d+00)

      logical debug
      debug=.false.

! Initialize
      WA=zero
      WB=zero
      WAtrans=zero
      WBtrans=zero
      wFAEw=zero
      wFBEw=zero
      wSAw=zero
      wSBw=zero
      AUXA=zero
      AUXB=zero
      xAEen=zero
      xBEen=zero
      xvecAE=zero
      xvecBE=zero
      blockvecAE=zero
      blockvecBE=zero

C ARS(
      if (debug) then
      write(*,*) "MATRIX Previous vecAE:"
      call PREVNU(vecAE0,AEen,nebf,nebf,nebf)
      write(*,*) "MATRIX Previous vecBE:"
      call PREVNU(vecBE0,BEen,nebf,nebf,nebf)
      write(*,*) "MATRIX Previous S:"
      call PREVNU(Selec,BEen,nebf,nebf,nebf)
      end if
C )

!!!!!!!!!!!!!!!! Solve for regular electronic solution !!!!!!!!!!!!!!!!

! Form regular electronic transformation matrix
      do i=1,nvirt
        do j=1,nebf
          WA(j,i)=vecBE0(j,noccb+i)
        end do
      end do

      do i=1,nocca
        do j=1,nebf
          WA(j,i+nvirt)=vecAE0(j,i)
        end do
      end do

C ARS(
      if (debug) then
      write(*,*) "MATRIX WA:"
      call PREVNU(WA,xAEen,noccvirta,nebf,nebf)
      end if
C )

! Form transpose
      do i=1,noccvirta
        do j=1,nebf
          WAtrans(i,j)=WA(j,i)
        end do
      end do

! Transform FAE as Wtrans * FAE * W
      call matmult(nebf,nebf,nebf,noccvirta,
     x                   FAE,WA,AUXA)
      call matmult(noccvirta,nebf,nebf,noccvirta,
     x                   WAtrans,AUXA,wFAEw)

! Transform AO overlap matrix as Wtrans * S * W
      call matmult(nebf,nebf,nebf,noccvirta,
     x                   Selec,WA,AUXA)
      call matmult(noccvirta,nebf,nebf,noccvirta,
     x                   WAtrans,AUXA,wSAw)

C ARS(
      if (debug) then
      write(*,*) "MATRIX wSAw:"
      call PREVNU(wSAw,xAEen,noccvirta,noccvirta,noccvirta)
      end if
C )

! Diagonalize transformed FAE to obtain solutions in new basis
      call UROOTHAN(xvecAE,xAEen,wSAw,wFAEw,noccvirta)

C ARS(
      if (debug) then
      WRITE(*,*) "MATRIX xvecAE:"
      call PREVNU(xvecAE,xAEen,noccvirta,noccvirta,noccvirta)
      end if
C )

! Transform evectors to AO basis as vecAE = W * xvecAE
      call matmult(nebf,noccvirta,noccvirta,noccvirta,
     x                   WA,xvecAE,blockvecAE)

! Pass evectors to output variables (zeros for unfilled noccb part)
      do i=1,noccvirta
        AEen(i)=xAEen(i)
        do j=1,nebf
          vecAE(j,i)=blockvecAE(j,i)
        end do
      end do

C ARS(
      if (debug) then
      WRITE(*,*) "MATRIX New vecAE:"
      call PREVNU(vecAE,AEen,nebf,nebf,nebf)
      end if
C )

!!!!!!!!!!!!!!!! Solve for special electronic solution !!!!!!!!!!!!!!!!

! Form regular electronic transformation matrix
      do i=1,nvirt
        do j=1,nebf
          WB(j,i)=vecAE(j,nocca+i)
        end do
      end do

      do i=1,noccb
        do j=1,nebf
          WB(j,i+nvirt)=vecBE0(j,i)
        end do
      end do

C ARS(
      if (debug) then
      write(*,*) "MATRIX WB:"
      call PREVNU(WB,xBEen,noccvirtb,nebf,nebf)
      end if
C )

! Form transpose
      do i=1,noccvirtb
        do j=1,nebf
          WBtrans(i,j)=WB(j,i)
        end do
      end do

! Transform FBE as Wtrans * FBE * W
      call matmult(nebf,nebf,nebf,noccvirtb,
     x                   FBE,WB,AUXB)
      call matmult(noccvirtb,nebf,nebf,noccvirtb,
     x                   WBtrans,AUXB,wFBEw)

! Transform AO overlap matrix as Wtrans * S * W
      call matmult(nebf,nebf,nebf,noccvirtb,
     x                   Selec,WB,AUXB)
      call matmult(noccvirtb,nebf,nebf,noccvirtb,
     x                   WBtrans,AUXB,wSBw)

C ARS(
      if (debug) then
      write(*,*) "MATRIX wSBw:"
      call PREVNU(wSBw,xBEen,noccvirtb,noccvirtb,noccvirtb)
      end if
C )

! Diagonalize transformed FBE to obtain solutions in new basis
      call UROOTHAN(xvecBE,xBEen,wSBw,wFBEw,noccvirtb)

C ARS(
      if (debug) then
      WRITE(*,*) "MATRIX xvecBE:"
      call PREVNU(xvecBE,xBEen,noccvirtb,noccvirtb,noccvirtb)
      end if
C )

! Transform evectors to AO basis as vecBE = W * xvecBE
      call matmult(nebf,noccvirtb,noccvirtb,noccvirtb,
     x                   WB,xvecBE,blockvecBE)

! Pass evectors to output variables (zeros for unfilled noccb part)
      do i=1,noccvirtb
        BEen(i)=xBEen(i)
        do j=1,nebf
          vecBE(j,i)=blockvecBE(j,i)
        end do
      end do

C ARS(
      if (debug) then
      WRITE(*,*) "MATRIX New vecBE:"
      call PREVNU(vecBE,BEen,nebf,nebf,nebf)
      end if
C )

      return
      end
!======================================================================
      subroutine OCBSE2(nebf,npebf,nae,nbe,vecAE,vecBE0,FBE,
     x                  Selec,elcam,elcbfc,elcex,
     x                  ampeb2c,agebfcc,kpestr,kpeend,
     x                  vecBE,BEen)
! 
!     Perform OCBSE procedure for special electrons only
!       vecBE0: Special electron coefficients from previous iteration
!       vecAE:  Regular electron coefficients from current iteration
!       vecBE:  Special electron coefficients from current iteration
!
!======================================================================
      implicit none
! Input Variables
      integer nebf,npebf
      integer nae,nbe
      double precision vecAE(nebf,nebf),vecBE0(nebf,nebf)
      double precision FBE(nebf,nebf)
      double precision Selec(nebf,nebf)
      integer ampeb2c(npebf)               ! Map prim index to contr index
      integer kpestr(nebf)                 ! Map contr index to prim start
      integer kpeend(nebf)                 ! Map contr index to prim end
      integer elcam(npebf,3)               ! Angular mom for electrons
      double precision agebfcc(npebf)      ! Map prim index to contr coeff
      double precision elcex(npebf)        ! Exponents: elec basis
      double precision elcbfc(npebf,3)     ! Basis centers: elec basis
! Variables Returned
      double precision vecBE(nebf,nebf)
      double precision BEen(nebf)
! Local variables
      integer nocca,noccb,nvirt 
      integer noccvirta,noccvirtb
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
      nvirt=nebf-nocca-noccb
      noccvirta=nebf-noccb
      noccvirtb=nebf-nocca

      vecBE=zero
      BEen=zero

      call OCBSE2_driver(nebf,npebf,nae,nbe,nocca,noccb,
     x                   nvirt,noccvirta,noccvirtb,
     x                   vecAE,vecBE0,FBE,Selec,
     x                   elcam,elcbfc,elcex,
     x                   ampeb2c,agebfcc,kpestr,kpeend,
     x                   vecBE,BEen)

      return
      end
!======================================================================
      subroutine OCBSE2_driver(nebf,npebf,nae,nbe,nocca,noccb,
     x                               nvirt,noccvirta,noccvirtb,
     x                               vecAE,vecBE0,FBE,Selec,
     x                               elcam,elcbfc,elcex,
     x                               ampeb2c,agebfcc,kpestr,kpeend,
     x                               vecBE,BEen)
!
!     OCBSE Procedure:
!       - Construct transformation matrix for special electrons
!       - Transform FBE and diagonalize to obtain vecBE
!
!  Switching to virtual space (no remnant of occupied in proj basis)
!
!======================================================================
      implicit none
! Input Variables
      integer nebf,npebf
      integer nae,nbe
      integer nocca,noccb,nvirt
      integer noccvirta,noccvirtb
      double precision vecAE(nebf,nebf),vecBE0(nebf,nebf)
      double precision FBE(nebf,nebf)
      double precision Selec(nebf,nebf)
      integer ampeb2c(npebf)               ! Map prim index to contr index
      integer kpestr(nebf)                 ! Map contr index to prim start
      integer kpeend(nebf)                 ! Map contr index to prim end
      integer elcam(npebf,3)               ! Angular mom for electrons
      double precision agebfcc(npebf)      ! Map prim index to contr coeff
      double precision elcex(npebf)        ! Exponents: elec basis
      double precision elcbfc(npebf,3)     ! Basis centers: elec basis
! Variables Returned
      double precision vecBE(nebf,nebf)
      double precision BEen(nebf)
! Local variables
      integer i,j
      integer ie1,je1
      integer iec1,jec1
      integer ie1_start,je1_start
      integer ie1_end,je1_end
      integer mo1,mo2,motodiscard
      integer i1,j1,k1
      integer l1,m1,n1
      double precision a1,b1
      double precision cof_ie1,cof_je1
      double precision ans,ovlap,ovlapcheck,maxovlap
      double precision Amat1(3), Bmat1(3)
      double precision ovlaparr(nebf),projorb(nebf) ! fix2
      double precision WB(nebf,noccvirtb)
      double precision WBtrans(noccvirtb,nebf)
      double precision wFBEw(noccvirtb,noccvirtb)
      double precision wSBw(noccvirtb,noccvirtb)
      double precision AUXB(nebf,noccvirtb)
      double precision xAEen(noccvirta)
      double precision xBEen(noccvirtb)
      double precision xvecBE(noccvirtb,noccvirtb)
      double precision blockvecBE(nebf,noccvirtb)
      double precision zero
      parameter(zero=0.0d+00)

      logical debug
      debug=.false.

! Initialize
      WB=zero
      WBtrans=zero
      wFBEw=zero
      wSBw=zero
      AUXB=zero
      xBEen=zero
      xvecBE=zero
      blockvecBE=zero

C ARS(
      if (debug) then
      write(*,*) "nae,nbe:",nae,nbe
      write(*,*) "nocca,noccb:",nocca,noccb
      write(*,*) "MATRIX vecAE:"
      call PREVNU(vecAE,BEen,nebf,nebf,nebf)
      write(*,*) "MATRIX Previous vecBE:"
      call PREVNU(vecBE0,BEen,nebf,nebf,nebf)
      write(*,*) "MATRIX Previous S:"
      call PREVNU(Selec,BEen,nebf,nebf,nebf)
      end if
C )

!!!!!!!!!!!!!!!! Solve for special electronic solution !!!!!!!!!!!!!!!!

! Form special electronic transformation matrix
      do i=1,noccvirtb
        do j=1,nebf
          WB(j,i)=vecAE(j,nocca+i)
        end do
      end do

C ARS(
      if (debug) then
      write(*,*) "MATRIX WB:"
      call PREVNU(WB,xBEen,noccvirtb,nebf,nebf)
      end if
C )

! Form transpose
      do i=1,noccvirtb
        do j=1,nebf
          WBtrans(i,j)=WB(j,i)
        end do
      end do

! Transform FBE as Wtrans * FBE * W
      call matmult(nebf,nebf,nebf,noccvirtb,
     x                   FBE,WB,AUXB)
      call matmult(noccvirtb,nebf,nebf,noccvirtb,
     x                   WBtrans,AUXB,wFBEw)

! Transform AO overlap matrix as Wtrans * S * W
      call matmult(nebf,nebf,nebf,noccvirtb,
     x                   Selec,WB,AUXB)
      call matmult(noccvirtb,nebf,nebf,noccvirtb,
     x                   WBtrans,AUXB,wSBw)

C ARS(
      if (debug) then
      write(*,*) "MATRIX wSBw:"
      call PREVNU(wSBw,xBEen,noccvirtb,noccvirtb,noccvirtb)
      end if
C )

! Diagonalize transformed FBE to obtain solutions in new basis
      call UROOTHAN(xvecBE,xBEen,wSBw,wFBEw,noccvirtb)

C ARS(
      if (debug) then
      WRITE(*,*) "MATRIX xvecBE:"
      call PREVNU(xvecBE,xBEen,noccvirtb,noccvirtb,noccvirtb)
      end if
C )

! Transform evectors to AO basis as vecBE = W * xvecBE
      call matmult(nebf,noccvirtb,noccvirtb,noccvirtb,
     x                   WB,xvecBE,blockvecBE)

! Pass evectors to output variables (zeros for unfilled noccb part)
      do i=1,noccvirtb
        BEen(i)=xBEen(i)
        do j=1,nebf
          vecBE(j,i)=blockvecBE(j,i)
        end do
      end do

C ARS(
      if (debug) then
      WRITE(*,*) "MATRIX New vecBE:"
      call PREVNU(vecBE,BEen,nebf,nebf,nebf)
      end if
C )

      return
      end
!======================================================================
      subroutine matmult(rowsa,colsa,rowsb,colsb,A,B,C)
!
!     Computes C = A * B
!======================================================================
      implicit none
! Input Variables
      integer rowsa,colsa,rowsb,colsb
      double precision A(rowsa,colsa)
      double precision B(rowsb,colsb)
! Output Variables
      double precision C(rowsa,colsb)
! Local Variables
      integer i,j,k
      double precision zero
      parameter(zero=0.0d+00)

      C=zero

      if (colsa.ne.rowsb) then
       write(*,*) "ERROR IN MATMULT, QUITTING"
       STOP
      end if

      do i=1,colsb
        do j=1,rowsa
          do k=1,colsa
            C(j,i)=C(j,i)+A(j,k)*B(k,i)
          end do
        end do
      end do

      return
      end

