!======================================================================
      subroutine uscf(nelec,NAE,NBE,NAalpE,NAbetE,NBalpE,NBbetE,
     x                NPRA,NPRB,
     x                npebf,nebf,nebf2,nebflt,ngee,
     x                read_CE,
     x                LHF,LSOSCF,LOCBSE,LDBG,
     x                nat,cat,zan,
     x                KPESTR,KPEEND,AMPEB2C,AGEBFCC,
     x                ELCEX,ELCAM,ELCBFC,GAM_ee)

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
      integer NAalpE,NAbetE
      integer NBalpE,NBbetE
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
C      double precision diffE
C      double precision diffAE
C      double precision diffBE
      double precision DIFFAalpE,DIFFAbetE
      double precision DIFFBalpE,DIFFBbetE

      double precision E_total

      double precision E_A
      double precision E_A_ecore
      double precision E_A_ee

      double precision E_B
      double precision E_B_ecore
      double precision E_B_ee

      double precision E_AB

      double precision E_nuc

      double precision DAalpE(NEBF,NEBF)
      double precision DAalpE0(NEBF,NEBF)
      double precision DAbetE(NEBF,NEBF)
      double precision DAbetE0(NEBF,NEBF)
      double precision DBalpE(NEBF,NEBF)
      double precision DBalpE0(NEBF,NEBF)
      double precision DBbetE(NEBF,NEBF)
      double precision DBbetE0(NEBF,NEBF)
      double precision vecAalpE(NEBF,NEBF)
      double precision vecAalpE0(NEBF,NEBF)
      double precision vecAbetE(NEBF,NEBF)
      double precision vecAbetE0(NEBF,NEBF)
      double precision vecBalpE(NEBF,NEBF)
      double precision vecBalpE0(NEBF,NEBF)
      double precision vecBbetE(NEBF,NEBF)
      double precision vecBbetE0(NEBF,NEBF)
      double precision AalpEE(NEBF)
      double precision AbetEE(NEBF)
      double precision BalpEE(NEBF)
      double precision BbetEE(NEBF)
      double precision FAalpE(nebf,nebf)
      double precision FAbetE(nebf,nebf)
      double precision FBalpE(nebf,nebf)
      double precision FBbetE(nebf,nebf)
      double precision XFAalpE(nebf,nebf)
      double precision XFAbetE(nebf,nebf)
      double precision XFBalpE(nebf,nebf)
      double precision XFBbetE(nebf,nebf)

      double precision FAalpEint(nebf,nebf)
      double precision FAbetEint(nebf,nebf)
      double precision FBalpEint(nebf,nebf)
      double precision FBbetEint(nebf,nebf)

      double precision E_total_old
      double precision Delta_E_tot

      logical LDIFFE

!--------SOSCF-RELATED-VARIABLES------------(
      logical LSOSCF,LSOSCFA,LSOSCFB
      logical EIGAVL
      integer ITER
      integer ITSOA ! SOSCF iteration counter
      integer ITSOB ! SOSCF iteration counter
      integer L0,L1
      integer NFT15
      integer NFT16
      double precision FLT(NEBFLT) !FLT: Lower triangle focke
      double precision HSTARTA(NPRA)
      double precision HSTARTB(NPRB)
      double precision GRADA(NPRA)
      double precision GRADB(NPRB)
      double precision PGRADA(NPRA)
      double precision PGRADB(NPRB)
      double precision DISPLIA(NPRA)
      double precision DISPLIB(NPRB)
      double precision DGRADA(NPRA)  ! WRK1
      double precision DGRADB(NPRB)  ! WRK1
      double precision DISPLA(NPRA)  ! WRK2
      double precision DISPLB(NPRB)  ! WRK2
      double precision UPDTA(NPRA)   ! WRK3
      double precision UPDTB(NPRB)   ! WRK3
      double precision DISPLNA(NPRA) ! WRK1+NPR
      double precision DISPLNB(NPRB) ! WRK1+NPR
      double precision DGRADIA(NPRA) ! WRK2+NPR
      double precision DGRADIB(NPRB) ! WRK2+NPR
      double precision UPDTIA(NPRA)  ! WRK3+NPR
      double precision UPDTIB(NPRB)  ! WRK3+NPR
      double precision ORBGRDA
      double precision ORBGRDB
      double precision SMALL
      double precision SOGTOL ! ORBGRAD TOL to activate soscf
      double precision XA(NPRA)
      double precision XB(NPRB)
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

 9150 FORMAT(1X,I3,F20.10,F17.10,6F17.10)

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

 9610 FORMAT(/1X,' A alpha ELECTRONIC ORBITALS AND EIGENVALUES: ')
 9611 FORMAT(/1X,' A  beta ELECTRONIC ORBITALS AND EIGENVALUES: ')

 9620 FORMAT(/1X,' B alpha ELECTRONIC ORBITALS AND EIGENVALUES: ')
 9621 FORMAT(/1X,' B  beta ELECTRONIC ORBITALS AND EIGENVALUES: ')

 9800 FORMAT(10X,15(1H-),'START SECOND ORDER SCF',15(1H-))
                                           
!--------OUTPUT-FORMATTING---------------------------------------------)

C      LOCBSE2=LOCBSE
      LOCBSE2=.false.
      if(LOCBSE2) then
       LOCBSE=.false.
      end if
C      if (LHF) then
C       write(*,*) "Performing HF calculation"
C       LOCBSE=.false.
C       LOCBSE2=.false.
C       nae=nelec
C      end if
      if(LOCBSE) write(*,*) "Using LOCBSE"
      if(LOCBSE2) write(*,*) "Using LOCBSE2"

      if (.not.(LOCBSE.or.LOCBSE2)) then
       write(*,*) "Only LOCBSE/LOCBSE2 supported..."
       return
      end if

      DAalpE=0.0d+00
      DAbetE=0.0d+00
      DBalpE=0.0d+00
      DBbetE=0.0d+00
      DAalpE0=0.0d+00
      DAbetE0=0.0d+00
      DBalpE0=0.0d+00
      DBbetE0=0.0d+00
      vecAalpE=0.0d+00
      vecAbetE=0.0d+00
      vecBalpE=0.0d+00
      vecBbetE=0.0d+00
      vecAalpE0=0.0d+00
      vecAbetE0=0.0d+00
      vecBalpE0=0.0d+00
      vecBbetE0=0.0d+00
      AalpEE=0.0d+00
      AbetEE=0.0d+00
      BalpEE=0.0d+00
      BbetEE=0.0d+00

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
         call read_uCE(nebf,NAalpE,NAbetE,NBalpE,NBbetE,
     x                 DAalpE,DAbetE,DBalpE,DBbetE,
     x                 vecAalpE0,vecAbetE0,vecBalpE0,vecBbetE0)
      else
!       STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
        write(*,*)'ABOUT TO CALL guess_A_elec'
!       call guess_elec(nelec,nebf,xxse,GAM_ecore,DE)
        if ((LOCBSE).or.(LOCBSE2)) then
          call guess_uelec(NAalpE,NBalpE,nebf,xxse,GAM_ecore,
     x                    DAalpE,DBalpE,vecAalpE0,vecBalpE0)
          call guess_uelec(NAbetE,NBbetE,nebf,xxse,GAM_ecore,
     x                    DAbetE,DBbetE,vecAbetE0,vecBbetE0)
          write(*,*)'BACK FROM guess_elec for OCBSE'
        else
         call guess_uA_elec(NAalpE,NAbetE,nebf,xxse,GAM_ecore,
     x                      DAalpE,DAbetE,vecAalpE0,vecAbetE0)
         call guess_uA_elec(NBalpE,NBbetE,nebf,xxse,GAM_ecore,
     x                      DBalpE,DBbetE,vecBalpE0,vecBbetE0)
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
       call PREVNU(vecAalpE0,AalpEE,nebf,nebf,nebf)
       WRITE(*,9611)
       call PREVNU(vecAbetE0,AbetEE,nebf,nebf,nebf)
       WRITE(*,9620)
       call PREVNU(vecBalpE0,BalpEE,nebf,nebf,nebf)
       WRITE(*,9621)
       call PREVNU(vecBbetE0,BbetEE,nebf,nebf,nebf)
       write(*,*)
      end if
C )
!-------------INITIAL-GUESSES------------------------------------------)

!-------------SETUP-FOR-POSSIBLE-SOSCF---------------------------------(
!     LSOSCF=.FALSE.
C B set of electrons never have SOSCF enabled since OCBSE2 should be invoked
C If OCBSE is requested then no SOSCF should be done at all
      if(LOCBSE) LSOSCF=.FALSE.
      if(LSOSCF) then
         SOGTOL=1.0d+00
         SMALL=1.0D-06
         L0=nebf
         L1=nebf
         LSOSCFA=.true.
         LSOSCFB=.true.
         if(NAalpE.le.1) LSOSCFA=.false.
         if(NAbetE.le.1) LSOSCFB=.false.
      else
         LSOSCFA=.false.
         LSOSCFB=.false.
      end if
      if(LSOSCFA) then  ! A alpha SOSCF
         NFT15=15
         OPEN(NFT15, FILE='WORK15', STATUS='UNKNOWN',
     *        ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         ITSOA=0
      end if
      if(LSOSCFB) then  ! A beta SOSCF
         NFT16=16
         OPEN(NFT16, FILE='WORK16', STATUS='UNKNOWN',
     *        ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         ITSOB=0
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
      ORBGRDB=0.0d+00
      PGRADA=0.0d+00
      PGRADB=0.0d+00

      DO I=1,MAXIT

!--------------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------(
!     write(*,*)
!     write(*,*)'IN: xcuscf '
!     write(*,*)'Before call to UHF_FOCK, GAM_EE='
!     write(*,*)GAM_ee
!     write(*,*)
       if((.not.((locbse).or.(locbse2))).or.(i.eq.1)) then

C Call HF Fock build for NAalpE+NAbetE=NAE electrons
         call fock_uhf(LDBG,nebf,nebf2,NAalpE,NAbetE,ngee,
     x                 DAalpE,DAbetE,GAM_ecore,GAM_ee,
     x                 FAalpE,FAbetE,E_A,E_A_ecore,E_A_ee)

         if (.not.LHF) then
C Call HF Fock build for NBalpE+NBbetE=NBE electrons
          call fock_uhf(LDBG,nebf,nebf2,NBalpE,NBbetE,ngee,
     x                  DBalpE,DBbetE,GAM_ecore,GAM_ee,
     x                  FBalpE,FBbetE,E_B,E_B_ecore,E_B_ee)

C Call interaction Fock build for all electrons
          call fock_uint(LDBG,nelec,NAalpE,NAbetE,NBalpE,NBbetE,
     x                   nebf,nebf2,ngee,
     x                   DAalpE,DAbetE,DBalpE,DBbetE,GAM_ee,
     x                   FAalpEint,FAbetEint,FBalpEint,FBbetEint,
     x                   E_AB)

          call add2fock(nebf,FAalpEint,FAalpE)
          call add2fock(nebf,FAbetEint,FAbetE)
          call add2fock(nebf,FBalpEint,FBalpE)
          call add2fock(nebf,FBbetEint,FBbetE)
         end if

C          IF (LDBG) then
C           write(*,*)
C           write(*,*) "FAE:"
C           call prt_lower_triangle(nebf,nebflt,FAE)
C           write(*,*)
C           if (.not.LHF) then
C            write(*,*) "FBE:"
C            call prt_lower_triangle(nebf,nebflt,FBE)
C           end if
C           write(*,*)
C          END IF  

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
! Do OCBSE procedure (restricted solutions for A and B electrons)

! Alpha electrons
           if ((NAalpE.gt.0).and.(NBalpE.gt.0)) then
            call uOCBSE(nebf,NAalpE,NBalpE,vecAalpE0,vecBalpE0,
     x                  FAalpE,FBalpE,xxse,
     x                  vecAalpE,vecBalpE,AalpEE,BalpEE)
           else if ((NAalpE.gt.0).and.(NBalpE.eq.0)) then
            call UROOTHAN(vecAalpE,AalpEE,xxse,FAalpE,nebf)
           else if ((NAalpE.eq.0).and.(NBalpE.gt.0)) then
            call UROOTHAN(vecBalpE,BalpEE,xxse,FBalpE,nebf)
           end if

           call construct_uDE(NAalpE,nebf,vecAalpE,DAalpE)
           CALL DENDIF(DAalpE0,DAalpE,NEBF,DIFFAalpE)
           CALL COPYDEN(DAalpE0,DAalpE,NEBF)
           CALL COPYDEN(vecAalpE0,vecAalpE,NEBF)

           call construct_uDE(NBalpE,nebf,vecBalpE,DBalpE)
           CALL DENDIF(DBalpE0,DBalpE,NEBF,DIFFBalpE)
           CALL COPYDEN(DBalpE0,DBalpE,NEBF)
           CALL COPYDEN(vecBalpE0,vecBalpE,NEBF)

! Beta electrons
           if ((NAbetE.gt.0).and.(NBbetE.gt.0)) then
            call uOCBSE(nebf,NAbetE,NBbetE,vecAbetE0,vecBbetE0,
     x                  FAbetE,FBbetE,xxse,
     x                  vecAbetE,vecBbetE,AbetEE,BbetEE)
           else if ((NAbetE.gt.0).and.(NBbetE.eq.0)) then
            call UROOTHAN(vecAbetE,AbetEE,xxse,FAbetE,nebf)
           else if ((NAbetE.eq.0).and.(NBbetE.gt.0)) then
            call UROOTHAN(vecBbetE,BbetEE,xxse,FBbetE,nebf)
           end if

           call construct_uDE(NAbetE,nebf,vecAbetE,DAbetE)
           CALL DENDIF(DAbetE0,DAbetE,NEBF,DIFFAbetE)
           CALL COPYDEN(DAbetE0,DAbetE,NEBF)
           CALL COPYDEN(vecAbetE0,vecAbetE,NEBF)

           call construct_uDE(NBbetE,nebf,vecBbetE,DBbetE)
           CALL DENDIF(DBbetE0,DBbetE,NEBF,DIFFBbetE)
           CALL COPYDEN(DBbetE0,DBbetE,NEBF)
           CALL COPYDEN(vecBbetE0,vecBbetE,NEBF)

! Calculate energy for this it and Fock matrices for next it

C Call HF Fock build for NAalpE+NAbetE=NAE electrons
         call fock_uhf(LDBG,nebf,nebf2,NAalpE,NAbetE,ngee,
     x                 DAalpE,DAbetE,GAM_ecore,GAM_ee,
     x                 FAalpE,FAbetE,E_A,E_A_ecore,E_A_ee)

C Call HF Fock build for NBalpE+NBbetE=NBE electrons
         call fock_uhf(LDBG,nebf,nebf2,NBalpE,NBbetE,ngee,
     x                 DBalpE,DBbetE,GAM_ecore,GAM_ee,
     x                 FBalpE,FBbetE,E_B,E_B_ecore,E_B_ee)

C Call interaction Fock build for all electrons
          call fock_uint(LDBG,nelec,NAalpE,NAbetE,NBalpE,NBbetE,
     x                   nebf,nebf2,ngee,
     x                   DAalpE,DAbetE,DBalpE,DBbetE,GAM_ee,
     x                   FAalpEint,FAbetEint,FBalpEint,FBbetEint,
     x                   E_AB)

          call add2fock(nebf,FAalpEint,FAalpE)
          call add2fock(nebf,FAbetEint,FAbetE)
          call add2fock(nebf,FBalpEint,FBalpE)
          call add2fock(nebf,FBbetEint,FBbetE)

C            IF (LDBG) then
C             write(*,*)
C             write(*,*) "FAE:"
C             call prt_lower_triangle(nebf,nebflt,FAE)
C             write(*,*)
C             write(*,*) "FBE:"
C             call prt_lower_triangle(nebf,nebflt,FBE)
C             write(*,*)
C            END IF  

            E_total=E_A+E_B+E_AB+E_nuc

         else if (LOCBSE2) then
! Do OCBSE2 procedure (restricted solutions for B electrons)

! A alpha electrons
         if(LSOSCFA) THEN
          ITER=I
          EIGAVL = ITER.GT.1
         end if
         IF(LSOSCFA .AND.  EIGAVL) THEN
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAalpE,FLT)
          call SOGRAD(GRADA,FLT,vecAalpE,WRK,NPRA,NAalpE,
     x                L0,L1,NEBFLT,ORBGRDA)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 900  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDA.LT.SOGTOL  .OR.  ITSOA.GT.0) THEN
              IF(ITSOA.EQ.0) THEN
              WRITE(*,9800)
                 call SOHESS(HSTARTA,AalpEE,NPRA,L0,NAalpE,NAalpE)
              END IF
              ITSOA = ITSOA+1
           call SONEWT(HSTARTA,GRADA,PGRADA,DISPLIA,DGRADA,DISPLA,UPDTA,
     *                 DISPLNA,DGRADIA,UPDTIA,ORBGRDA,NPRA,ITSOA,NFT15)
            call SOTRAN(DISPLIA,vecAalpE,GA,WRK,NPRA,
     x                  L0,L1,NAalpE,NAalpE,ORBGRDA)
             CALL DCOPY(NPRA,GRADA,1,PGRADA,1)
              call construct_uDE(NAalpE,nebf,vecAalpE,DAalpE)
              GO TO 950  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------)

  900 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAalpE,AalpEE,xxse,FAalpE,nebf)
         call construct_uDE(NAalpE,nebf,vecAalpE,DAalpE)

  950 CONTINUE
!        --> FIND LARGEST CHANGE IN Alpha E DENSITY
         CALL DENDIF(DAalpE0,DAalpE,NEBF,DIFFAalpE)
         CALL COPYDEN(DAalpE0,DAalpE,NEBF)

! Do OCBSE procedure for B alpha electrons
           call uOCBSE2(nebf,npebf,NAalpE,NBalpE,vecAalpE,vecBalpE0,
     x                  FBalpE,xxse,elcam,elcbfc,elcex,
     x                  ampeb2c,agebfcc,kpestr,kpeend,
     x                  vecBalpE,BalpEE)

           call construct_uDE(NBalpE,nebf,vecBalpE,DBalpE)
           CALL DENDIF(DBalpE0,DBalpE,NEBF,DIFFBalpE)
           CALL COPYDEN(DBalpE0,DBalpE,NEBF)
           CALL COPYDEN(vecAalpE0,vecAalpE,NEBF)
           CALL COPYDEN(vecBalpE0,vecBalpE,NEBF)

! A beta electrons
         if(LSOSCFB) THEN
          ITER=I
          EIGAVL = ITER.GT.1
         end if
         IF(LSOSCFB .AND.  EIGAVL) THEN
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAbetE,FLT)
          call SOGRAD(GRADB,FLT,vecAbetE,WRK,NPRB,NAbetE,
     x                L0,L1,NEBFLT,ORBGRDB)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 910  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDB.LT.SOGTOL  .OR.  ITSOB.GT.0) THEN
              IF(ITSOB.EQ.0) THEN
              WRITE(*,9800)
                 call SOHESS(HSTARTB,AbetEE,NPRB,L0,NAbetE,NAbetE)
              END IF
              ITSOB = ITSOB+1
           call SONEWT(HSTARTB,GRADB,PGRADB,DISPLIB,DGRADB,DISPLB,UPDTB,
     *                 DISPLNB,DGRADIB,UPDTIB,ORBGRDB,NPRB,ITSOB,NFT16)
            call SOTRAN(DISPLIB,vecAbetE,GB,WRK,NPRB,
     x                  L0,L1,NAbetE,NAbetE,ORBGRDB)
             CALL DCOPY(NPRB,GRADB,1,PGRADB,1)
              call construct_uDE(NAbetE,nebf,vecAbetE,DAbetE)
              GO TO 960  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------)

  910 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAbetE,AbetEE,xxse,FAbetE,nebf)
         call construct_uDE(NAbetE,nebf,vecAbetE,DAbetE)

  960 CONTINUE
!        --> FIND LARGEST CHANGE IN Alpha E DENSITY
         CALL DENDIF(DAbetE0,DAbetE,NEBF,DIFFAbetE)
         CALL COPYDEN(DAbetE0,DAbetE,NEBF)

! Do OCBSE procedure for B beta electrons
           call uOCBSE2(nebf,npebf,NAbetE,NBbetE,vecAbetE,vecBbetE0,
     x                  FBbetE,xxse,elcam,elcbfc,elcex,
     x                  ampeb2c,agebfcc,kpestr,kpeend,
     x                  vecBbetE,BbetEE)

           call construct_uDE(NBbetE,nebf,vecBbetE,DBbetE)
           CALL DENDIF(DBbetE0,DBbetE,NEBF,DIFFBbetE)
           CALL COPYDEN(DBbetE0,DBbetE,NEBF)
           CALL COPYDEN(vecAbetE0,vecAbetE,NEBF)
           CALL COPYDEN(vecBbetE0,vecBbetE,NEBF)

! Calculate energy for this it and Fock matrices for next it

C Call HF Fock build for NAalpE+NAbetE=NAE electrons
         call fock_uhf(LDBG,nebf,nebf2,NAalpE,NAbetE,ngee,
     x                 DAalpE,DAbetE,GAM_ecore,GAM_ee,
     x                 FAalpE,FAbetE,E_A,E_A_ecore,E_A_ee)

C Call HF Fock build for NBalpE+NBbetE=NBE electrons
         call fock_uhf(LDBG,nebf,nebf2,NBalpE,NBbetE,ngee,
     x                 DBalpE,DBbetE,GAM_ecore,GAM_ee,
     x                 FBalpE,FBbetE,E_B,E_B_ecore,E_B_ee)

C Call interaction Fock build for all electrons
          call fock_uint(LDBG,nelec,NAalpE,NAbetE,NBalpE,NBbetE,
     x                   nebf,nebf2,ngee,
     x                   DAalpE,DAbetE,DBalpE,DBbetE,GAM_ee,
     x                   FAalpEint,FAbetEint,FBalpEint,FBbetEint,
     x                   E_AB)

          call add2fock(nebf,FAalpEint,FAalpE)
          call add2fock(nebf,FAbetEint,FAbetE)
          call add2fock(nebf,FBalpEint,FBalpE)
          call add2fock(nebf,FBbetEint,FBbetE)

C            IF (LDBG) then
C             write(*,*)
C             write(*,*) "FAE:"
C             call prt_lower_triangle(nebf,nebflt,FAE)
C             write(*,*)
C             write(*,*) "FBE:"
C             call prt_lower_triangle(nebf,nebflt,FBE)
C             write(*,*)
C            END IF  

            E_total=E_A+E_B+E_AB+E_nuc

         end if ! end if for ocbse/ocbse2

!        --> CALCULATE CHANGE IN TOTAL ENERGY
         Delta_E_tot=E_total-E_total_old
         E_total_old=E_total

!        --> PRINT SUMMARY OF THIS ITERATION
            WRITE(*,9150) I,E_total,Delta_E_tot,
     x                    DIFFAalpE,DIFFAbetE,
     x                    DIFFBalpE,DIFFBbetE,
     x                    ORBGRDA,ORBGRDB

C ARS( debug: print out MOs here
      if (LDBG) then
       WRITE(*,9610)
       call PREVNU(vecAalpE,AalpEE,nebf,nebf,nebf)
       WRITE(*,9611)
       call PREVNU(vecAbetE,AbetEE,nebf,nebf,nebf)
       WRITE(*,9620)
       call PREVNU(vecBalpE,BalpEE,nebf,nebf,nebf)
       WRITE(*,9621)
       call PREVNU(vecBbetE,BbetEE,nebf,nebf,nebf)
      end if
C )

! Output the vectors for this iteration for restart if necessary:
         call write_MOs(870,nebf,vecAalpE)
         call write_MOs(871,nebf,vecAbetE)
         call write_MOs(880,nebf,vecBalpE)
         call write_MOs(881,nebf,vecBbetE)

         LDIFFE=( (DIFFAalpE.LT.TOLE).and.(DIFFAbetE.LT.TOLE)
     x      .and. (DIFFBalpE.LT.TOLE).and.(DIFFBbetE.LT.TOLE) )
         IF(LDIFFE) GOTO 100
         IF(I.EQ.MAXIT) GOTO 10

      END DO
 
  10  CONTINUE
!     IF WE GET HERE SOMETHING WENT WRONG

      if(LSOSCFA) close(NFT15)
      if(LSOSCFB) close(NFT16)

      WRITE(*,*)
      WRITE(*,*)'WARNING:  ITERATION LIMIT EXCEEDED'
      E_total=zero
      WRITE(*,*)
!     STOP
!
  100 CONTINUE
!     IF WE GET HERE WE ARE DONE - CONVERGENCE ACHIEVED

      if(LSOSCFA) close(NFT15)
      if(LSOSCFB) close(NFT16)

!     PRINT FINAL ENERGY AND PUNCH THE ORBITALS
      WRITE(*,9200) E_total,I

      WRITE(*,9300) E_nuc,
     x              E_A_ecore,E_A_ee,E_A,
     x              E_B_ecore,E_B_ee,E_B,
     x              E_AB,
     x              E_total

!  OUTPUT ELEC AND NUC EIGENVALUES AND EIGENVECTORS
      WRITE(*,9610)
      call PREVNU(vecAalpE,AalpEE,nebf,nebf,nebf)
      WRITE(*,9611)
      call PREVNU(vecAbetE,AbetEE,nebf,nebf,nebf)
      WRITE(*,9620)
      call PREVNU(vecBalpE,BalpEE,nebf,nebf,nebf)
      WRITE(*,9621)
      call PREVNU(vecBbetE,BbetEE,nebf,nebf,nebf)

! PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------(
!     subroutine write_MOs(IFIL,nbf,VEC)
!     IFIL=852 :: FinalCE.dat
!     IFIL=853 :: FinalCP.dat
!     IFIL=860 :: FinalCAE.dat
!     IFIL=861 :: FinalCBE.dat
C     IFIL=870 :: FinalCAalpE.dat
C     IFIL=871 :: FinalCAbetE.dat
C     IFIL=880 :: FinalCBalpE.dat
C     IFIL=881 :: FinalCBbetE.dat
      call write_MOs(870,nebf,vecAalpE)
      call write_MOs(871,nebf,vecAbetE)
      call write_MOs(880,nebf,vecBalpE)
      call write_MOs(881,nebf,vecBbetE)
! PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------)
!

      RETURN
      END
!======================================================================
      subroutine uOCBSE(nebf,nae,nbe,vecAE0,vecBE0,FAE,FBE,
     x                  Selec,vecAE,vecBE,AEen,BEen)
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

      nocca=nae
      noccb=nbe
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
      subroutine uOCBSE2(nebf,npebf,nae,nbe,vecAE,vecBE0,FBE,
     x                   Selec,elcam,elcbfc,elcex,
     x                   ampeb2c,agebfcc,kpestr,kpeend,
     x                   vecBE,BEen)
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

      nocca=nae
      noccb=nbe
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
