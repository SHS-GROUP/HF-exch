C======================================================================
      subroutine fock_hf(LDBG,nebf,nebf2,nelec,ngee,
     x                   DE,GAM_ecore,GAM_ee,
     x                   focke,E_total,E_ecore,E_ee)

C HF Fock procedure
C======================================================================
      implicit none

C Input variables
      logical           LDBG
      integer           nebf
      integer           nebf2
      integer           ngee
      integer           nelec
      double precision  DE(nebf,nebf)
      double precision  GAM_ecore(nebf2)
      double precision  GAM_ee(ngee)

C Output variables
      double precision  focke(nebf,nebf)
      double precision  E_total
      double precision  E_ecore
      double precision  E_ee

C Local variables
      integer           nebflt
      double precision  zero
      parameter(zero=0.0d+00)


C Initialize values
      E_total=zero
      E_ecore=zero
      E_ee=zero
      focke=zero

C One-electron energy
      call E_from_GAM_ecore(nebf,nebf2,GAM_ecore,DE,focke,E_ecore)

C Two-electron energy
      if(nelec.gt.1) then
         call E_from_GAM_ee(nebf,ngee,GAM_ee,DE,focke,E_ee)
      end if

      E_total=E_ecore+E_ee

      if(LDBG) then

       call UFM_sym_check2(nebf,focke)

       nebflt=nebf*(nebf+1)/2
       write(*,*) "FE HF:"
       call prt_lower_triangle(nebf,nebflt,focke)
       write(*,*)

      end if

      return
      end
C======================================================================
      subroutine fock_int(LDBG,nelec,nae,nbe,
     x                    nebf,nebf2,ngee,
     x                    DAE,DBE,GAM_ee,
     x                    FAEint,FBEint, 
     x                    E_AB)

C HF interaction Fock procedure
C======================================================================
      implicit none

C Input variables
      logical           LDBG
      integer           nebf
      integer           nebf2
      integer           ngee
      integer           nelec
      integer           nae,nbe
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  GAM_ee(ngee)

C Output variables
      double precision  FAEint(nebf,nebf)
      double precision  FBEint(nebf,nebf)
      double precision  E_AB

C Local variables
      integer           nebflt
      double precision  zero
      parameter(zero=0.0d+00)


C Initialize values
      E_AB=zero
      FAEint=zero
      FBEint=zero

C Two-electron energy
      if((nae.ge.1).and.(nbe.ge.1)) then
         call E_from_interaction(nebf,ngee,GAM_ee,
     x                           DAE,DBE,FAEint,FBEint,E_AB)
      end if

      if(LDBG) then

       call UFM_sym_check2(nebf,FAEint)
       call UFM_sym_check2(nebf,FBEint)

       nebflt=nebf*(nebf+1)/2
       write(*,*) "FAE interaction:"
       call prt_lower_triangle(nebf,nebflt,FAEint)
       write(*,*)
       write(*,*) "FBE interaction:"
       call prt_lower_triangle(nebf,nebflt,FBEint)
       write(*,*)

      end if

      return
      end
!=======================================================================
      subroutine UFM_sym_check2(nbf,F)

! Checks that the Fock matrices are symmetric
!=======================================================================
      implicit none

! Input Variables
      integer nbf
      double precision F(nbf,nbf)

! No Variables Returned

! Local variables
      logical LSYM
      logical LOUTPUT
      integer i,j
      double precision val,tolerance,maxdif
      parameter(tolerance=1.0d-10)

      maxdif=0.0d+00
      LSYM=.true.
      LOUTPUT=.false.

      do i=1,nbf
         do j=1,i
            val=F(i,j)-F(j,i)
            if(abs(val).gt.tolerance) then
              maxdif=val
              if(LOUTPUT) then
              write(*,*)' FOCK MATRIX IS NOT SYMMETRIC FOR IJ=',I,J
              write(*,*)'>>>> FM(IJ)=',F(i,j)
              write(*,*)'>>>> FM(JI)=',F(j,i)
              end if
              LSYM=.false.
            end if
         end do
      end do

      if(LSYM) then
         write(*,*)'FOCK MATRIX IS SYMMETRIC'
      else
         write(*,*)' FOCK MATRIX IS NOT SYMMETRIC, MAXDIFF =',maxdif
      end if


      return
      end
!======================================================================
      subroutine E_from_GAM_ecore(nebf,nebf2,
     * GAM_ecore,DE,focke,E_ecore)

!======================================================================
      implicit none

! Input Variables
      integer nebf,nebf2

      double precision focke(nebf,nebf)
      double precision DE(nebf,nebf)
!     double precision GAM_ecore(nebf,nebf)
      double precision GAM_ecore(nebf2)
      double precision E_ecore

! Local variables
      integer ie1,je1,ia
      double precision val_ecore
      double precision zero
      parameter(zero=0.0d+00)


      E_ecore=zero

      do ie1=1,nebf
         do je1=1,nebf

            call pack_2D(nebf,je1,ie1,ia)
!           val_ecore=GAM_ecore(ie1,je1)
            val_ecore=GAM_ecore(ia)
            focke(ie1,je1)=val_ecore
            E_ecore=E_ecore+DE(ie1,je1)*val_ecore

         end do
      end do


      return
      end
!======================================================================
      subroutine E_from_GAM_ee(nebf,ngee,
     * GAM_ee,DE,focke,E_ee)
!
!======================================================================
      implicit none

      double precision zero,two,half
      parameter(zero=0.0d+00,two=2.0d+00,half=5.0d-01)
      integer nebf
      integer ngee
      double precision GAM_ee(ngee)
      double precision DE(nebf,nebf)
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
!                 vee1=GAM_ee(ie1,je1,ie2,je2)
!                 vee2=GAM_ee(ie1,je2,ie2,je1)
                  val=(vee1-half*vee2)*half
                  xfocke(ie1,je1)=xfocke(ie1,je1)+two*DE(ie2,je2)*val
                  E_ee=E_ee+DE(ie1,je1)*DE(ie2,je2)*val

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
!======================================================================
      subroutine E_from_interaction(nebf,ngee,GAM_ee,
     *                              DAE,DBE,FAEint,FBEint,E_AB)
!
!======================================================================
      implicit none

      integer nebf
      integer ngee
      double precision GAM_ee(ngee)
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision FAEint(nebf,nebf)
      double precision FBEint(nebf,nebf)
      double precision E_AB

! Local variables
      integer nebfLT
      integer ia1
      integer ie1,je1,ie2,je2
      double precision val
      double precision xFAEint(nebf,nebf)
      double precision xFBEint(nebf,nebf)
      double precision zero,half
      parameter(zero=0.0d+00,half=0.50d+00)

      E_AB=zero

      call zero_out(nebf,xFAEint)
      call zero_out(nebf,xFBEint)

      do ie2=1,nebf
         do je2=1,nebf

            do ie1=1,nebf
               do je1=1,nebf

                  call pack_4D(nebf,nebf,nebf,
     x                         je2,ie2,je1,ie1,ia1)

                  val=GAM_ee(ia1)
                  xFAEint(ie1,je1)=xFAEint(ie1,je1)+DBE(ie2,je2)*val
                  xFBEint(ie1,je1)=xFBEint(ie1,je1)+DAE(ie2,je2)*val
                  E_AB=E_AB+DAE(ie1,je1)*DBE(ie2,je2)*val*half

               end do
            end do

         end do
      end do

      call add2fock(nebf,xFAEint,FAEint)
      call add2fock(nebf,xFBEint,FBEint)

      return
      end
!======================================================================
      subroutine zero_out(nbf,amat)
!
!======================================================================

      parameter(zero=0.0d+00)

      integer nbf
      double precision amat(nbf,nbf)

! Local variables
      integer i
      integer j

      do i=1,nbf
         do j=1,nbf

            amat(i,j)=zero

         end do
      end do


      return
      end
!======================================================================
      subroutine add2fock(nbf,f_part,f_full)
!
!======================================================================

      integer nbf
      double precision f_full(nbf,nbf)
      double precision f_part(nbf,nbf)

! Local variables
      integer i
      integer j

      do i=1,nbf
         do j=1,nbf

            f_full(i,j)=f_full(i,j)+f_part(i,j)

         end do
      end do


      return
      end

