C======================================================================
      subroutine fock_uhf(LDBG,nebf,nebf2,NAE,NBE,ngee,
     x                    DAE,DBE,GAM_ecore,GAM_ee,
     x                    FAE,FBE,E_total,E_ecore,E_ee)

C HF Fock procedure
C======================================================================
      implicit none

C Input variables
      logical           LDBG
      integer           nebf
      integer           nebf2
      integer           ngee
      integer           nae,nbe
      double precision  DAE(nebf,nebf) ! Alpha electrons
      double precision  DBE(nebf,nebf) ! Beta electrons
      double precision  GAM_ecore(nebf2)
      double precision  GAM_ee(ngee)

C Output variables
      double precision  FAE(nebf,nebf)
      double precision  FBE(nebf,nebf)
      double precision  E_total
      double precision  E_ecore
      double precision  E_ee

C Local variables
      integer           nebflt
      double precision  DEtot(nebf,nebf)
      double precision  E_A_ecore
      double precision  E_B_ecore
      double precision  E_A_ee
      double precision  E_B_ee
      double precision  zero
      parameter(zero=0.0d+00)

C Initialize values
      E_total=zero
      E_ecore=zero
      E_A_ecore=zero
      E_B_ecore=zero
      E_ecore=zero
      E_A_ee=zero
      E_B_ee=zero
      E_ee=zero
      FAE=zero
      FBE=zero
      DEtot=zero

C Form total electronic density
      call add2fock(nebf,DAE,DEtot)
      call add2fock(nebf,DBE,DEtot)

C One-electron energy
      if (NAE.gt.0) then
       call E_from_GAM_ecore(nebf,nebf2,GAM_ecore,DAE,FAE,E_A_ecore)
      end if
      if (NBE.gt.0) then
       call E_from_GAM_ecore(nebf,nebf2,GAM_ecore,DBE,FBE,E_B_ecore)
      end if

C Two-electron energy
      if((NAE+NBE).gt.1) then
         call uE_from_GAM_ee(nebf,ngee,GAM_ee,DAE,DEtot,FAE,E_A_ee)
         call uE_from_GAM_ee(nebf,ngee,GAM_ee,DBE,DEtot,FBE,E_B_ee)
      end if

      E_ecore = E_A_ecore + E_B_ecore
      E_ee    = E_A_ee + E_B_ee
      E_total = E_ecore + E_ee

      if(LDBG) then

       write(*,*) "E_A_ecore:",E_A_ecore
       write(*,*) "E_B_ecore:",E_B_ecore
       write(*,*) "E_ecore:",E_ecore
       write(*,*)
       write(*,*) "E_A_ee:",E_A_ee
       write(*,*) "E_B_ee:",E_B_ee
       write(*,*) "E_ee:",E_ee

       call UFM_sym_check2(nebf,FAE)
       call UFM_sym_check2(nebf,FBE)

       nebflt=nebf*(nebf+1)/2
       write(*,*)
       write(*,*) "FAE HF:"
       call prt_lower_triangle(nebf,nebflt,FAE)
       write(*,*)
       write(*,*) "FBE HF:"
       call prt_lower_triangle(nebf,nebflt,FBE)
       write(*,*)

      end if

      return
      end
C======================================================================
      subroutine fock_uint(LDBG,nelec,NAalpE,NAbetE,NBalpE,NBbetE,
     x                     nebf,nebf2,ngee,
     x                     DAalpE,DAbetE,DBalpE,DBbetE,GAM_ee,
     x                     FAalpEint,FAbetEint,FBalpEint,FBbetEint,
     x                     E_AB)

C HF interaction Fock procedure
C======================================================================
      implicit none

C Input variables
      logical           LDBG
      integer           nebf
      integer           nebf2
      integer           ngee
      integer           nelec
      integer           NAalpE,NAbetE,NBalpE,NBbetE
      double precision  DAalpE(nebf,nebf) ! Alpha A electrons
      double precision  DAbetE(nebf,nebf) ! Beta  A electrons
      double precision  DBalpE(nebf,nebf) ! Alpha B electrons
      double precision  DBbetE(nebf,nebf) ! Beta  B electrons
      double precision  GAM_ee(ngee)

C Output variables
      double precision  FAalpEint(nebf,nebf)
      double precision  FAbetEint(nebf,nebf)
      double precision  FBalpEint(nebf,nebf)
      double precision  FBbetEint(nebf,nebf)
      double precision  E_AB

C Local variables
      integer           nebflt
      double precision  DAEtot(nebf,nebf)
      double precision  DBEtot(nebf,nebf)
      double precision  E_AB_alp
      double precision  E_AB_bet
      double precision  zero
      parameter(zero=0.0d+00)


C Initialize values
      E_AB=zero
      E_AB_alp=zero
      E_AB_bet=zero
      FAalpEint=zero
      FAbetEint=zero
      FBalpEint=zero
      FBbetEint=zero
      DAEtot=zero
      DBEtot=zero

C Form total electronic densities
      call add2fock(nebf,DAalpE,DAEtot)
      call add2fock(nebf,DAbetE,DAEtot)
      call add2fock(nebf,DBalpE,DBEtot)
      call add2fock(nebf,DBbetE,DBEtot)

C Two-electron energy
      if((NAalpE.gt.0).or.(NBalpE.gt.0)) then
         call E_from_interaction(nebf,ngee,GAM_ee,DAEtot,DBEtot,
     x                           FAalpEint,FBalpEint,E_AB_alp)
         if(NAalpE.eq.0) FAalpEint=zero
         if(NBalpE.eq.0) FBalpEint=zero
         E_AB=E_AB_alp ! same as calculated below
      end if
      if((NAbetE.gt.0).or.(NBbetE.gt.0)) then
         call E_from_interaction(nebf,ngee,GAM_ee,DAEtot,DBEtot,
     x                           FAbetEint,FBbetEint,E_AB_bet)
         if(NAbetE.eq.0) FAbetEint=zero
         if(NBbetE.eq.0) FBbetEint=zero
         E_AB=E_AB_bet ! same as calculated above
      end if

      if(LDBG) then

       write(*,*) "E_AB_alp:",E_AB_alp
       write(*,*) "E_AB_bet:",E_AB_bet

       call UFM_sym_check2(nebf,FAalpEint)
       call UFM_sym_check2(nebf,FAbetEint)
       call UFM_sym_check2(nebf,FBalpEint)
       call UFM_sym_check2(nebf,FBbetEint)

       nebflt=nebf*(nebf+1)/2
       write(*,*) "FAalpE interaction:"
       call prt_lower_triangle(nebf,nebflt,FAalpEint)
       write(*,*)
       write(*,*) "FAbetE interaction:"
       call prt_lower_triangle(nebf,nebflt,FAbetEint)
       write(*,*)
       write(*,*) "FBalpE interaction:"
       call prt_lower_triangle(nebf,nebflt,FBalpEint)
       write(*,*)
       write(*,*) "FBbetE interaction:"
       call prt_lower_triangle(nebf,nebflt,FBbetEint)
       write(*,*)

      end if

      return
      end
!======================================================================
      subroutine uE_from_GAM_ee(nebf,ngee,
     * GAM_ee,DAE,DE,focke,E_ee)
!
!======================================================================
      implicit none

      double precision zero,two,half
      parameter(zero=0.0d+00,two=2.0d+00,half=5.0d-01)
      integer nebf
      integer ngee
      double precision GAM_ee(ngee)
      double precision DAE(nebf,nebf) ! Alpha/beta density matrix
      double precision DE(nebf,nebf)  ! Total density matrix
      double precision focke(nebf,nebf)
      double precision E_ee

! Local variables
      integer nebfLT
      integer ia1,ia2
      integer ie1,je1,ie2,je2
      double precision vee1
      double precision vee2
      double precision val
      double precision xfocke(nebf,nebf)

      E_ee=zero

      call zero_out(nebf,xfocke)

      do ie2=1,nebf
         do je2=1,nebf

            do ie1=1,nebf
               do je1=1,nebf

                  call pack_4D(nebf,nebf,nebf,
     x                         je2,ie2,je1,ie1,ia1)

                  call pack_4D(nebf,nebf,nebf,
     x                         je1,ie2,je2,ie1,ia2)
   
                  vee1=GAM_ee(ia1)
                  vee2=GAM_ee(ia2)

                  xfocke(ie1,je1)=xfocke(ie1,je1)+DE(ie2,je2)*vee1
     x                                          -DAE(ie2,je2)*vee2
                  E_ee=E_ee+half*DAE(ie1,je1)*
     x                  (DE(ie2,je2)*vee1-DAE(ie2,je2)*vee2)

               end do
            end do

         end do
      end do

!     nebfLT=(nebf+nebf*nebf)/2
!     write(*,*)'IN E_from_GAM_ee: EE Contribution to Elec FOCK Matrix:'
!     call print_my_fock(nebf,nebfLT,xfocke)

!  Update the full electronic Fock matrix
      call add2fock(nebf,xfocke,focke)

      return
      end

