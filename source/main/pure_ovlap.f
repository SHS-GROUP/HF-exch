C======================================================================
      subroutine elec_ovlap(npebf,nebf,nebf2,
     x                      AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)
C Testing routine:  Calculates Overlap Matrix (electronic)
C over contracted basis functions for comparison to GAMESS'
C======================================================================
      implicit none

      integer npebf,nebf,nebf2
C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
c     integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
c     double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
c     double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
c     double precision AGNBFCC(npbf)  ! Nuclear contract coef
C-------Basis Set Info-------)
C Local variables
C--------------------------------(
C Basis set-related local variables
      integer I1,J1,K1
      integer L1,M1,N1
      double precision A1,Amat1(3)
      double precision B1,Bmat1(3)
C--------------------------------)
      integer ia,ie1,je1,iec1,jec1
      double precision GS(nebf2)
      double precision ans
      double precision Cof_ie1,Cof_je1
      double precision zero
      parameter(zero=0.0d+00)

C  zero out the 1-D array to hold the contracted integrals
      do ia=1,nebf2
         GS(ia)=zero
      end do

      do ie1=1,npebf
         do je1=1,npebf

C  Map from primitive BF indices to contracted indices
c           call MPEB2C(ie1,iec1)
c           call MPEB2C(je1,jec1)
            iec1=AMPEB2C(ie1)
            jec1=AMPEB2C(je1)

C  Get primitive Electron Basis Function Contraction Coefficients
c           call GEBFCC(ie1,Cof_ie1)
c           call GEBFCC(je1,Cof_je1)
            Cof_ie1=AGEBFCC(ie1)
            Cof_je1=AGEBFCC(je1)

C  Map the 2-index contracted integral to 1-D:
            call pack_2D(nebf,jec1,iec1,ia)

C Calculate the primitive integral
            A1=ELCEX(ie1)
            I1=ELCAM(ie1,1)
            J1=ELCAM(ie1,2)
            K1=ELCAM(ie1,3)
            Amat1(1)=ELCBFC(ie1,1)
            Amat1(2)=ELCBFC(ie1,2)
            Amat1(3)=ELCBFC(ie1,3)

            B1=ELCEX(je1)
            L1=ELCAM(je1,1)
            M1=ELCAM(je1,2)
            N1=ELCAM(je1,3)
            Bmat1(1)=ELCBFC(je1,1)
            Bmat1(2)=ELCBFC(je1,2)
            Bmat1(3)=ELCBFC(je1,3)

            call gfovlap(I1,J1,K1,A1,Amat1,
     2                   L1,M1,N1,B1,Bmat1,
     3                   ans)

            GS(ia)=GS(ia)+ans*Cof_ie1*Cof_je1

         end do
      end do

C Write to file
      open(814,file='eovlap.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do iec1=1,nebf
         do jec1=1,nebf

C  Map the 2-index contracted integral to 1-D:
            call pack_2D(nebf,jec1,iec1,ia)

            write(814,REC=ia) GS(ia)

         end do
      end do

      close(814)


      return
      end
C======================================================================
      subroutine check_elec_ovlap(nebf)
C======================================================================
      implicit none
C Input Variables
      integer nebf
C Local Variables
      integer nbfLT
      double precision xxse(nebf,nebf)

      call read_elec_ovlap(nebf,xxse)
      write(*,*)
      write(*,*)'CONTRACTED ELECTRONIC OVERLAP MATRIX'
      write(*,*)
      nbfLT=(nebf+nebf*nebf)/2
      call prt_lower_triangle(nebf,nbfLT,xxse)


      return
      end
C======================================================================
      subroutine ELCNORM3(npebf,nebf,AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

C Normalize the contraction coefficients of electronic basis set
C======================================================================
      implicit none

      integer npebf,nebf
C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
c     integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
c     double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
c     double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
c     double precision AGNBFCC(npbf)  ! Nuclear contract coef
C-------Basis Set Info-------)
C Local variables
C--------------------------------(
C Basis set-related local variables
      integer I1,J1,K1
      integer L1,M1,N1
      double precision A1,Amat1(3)
      double precision B1,Bmat1(3)
C--------------------------------)
      integer ia,ie1,je1,iec1,jec1
      integer i_pebf
      integer i_cebf
      double precision norm_fac
c     double precision GS(nebf2)
      double precision kappa(nebf)
      double precision ans
      double precision Cof_ie1,Cof_je1
      double precision zero
      parameter(zero=0.0d+00)

C----------------------------------------------------------------------(
      do ie1=1,npebf

           A1=ELCEX(ie1)
           I1=ELCAM(ie1,1)
           J1=ELCAM(ie1,2)
           K1=ELCAM(ie1,3)
           Amat1(1)=ELCBFC(ie1,1)
           Amat1(2)=ELCBFC(ie1,2)
           Amat1(3)=ELCBFC(ie1,3)

           call gfovlap(I1,J1,K1,A1,Amat1,
     x                  I1,J1,K1,A1,Amat1,
     x                  ans)

c       write(*,*)'CHECK NORM PRIM: iep=  ans=',ie1,ans/AGEBFCC(ie1)**2
           AGEBFCC(ie1)=AGEBFCC(ie1)/sqrt(ans)

      end do
C----------------------------------------------------------------------)

C  zero out the 1-D array to hold the contracted integral overlaps
      do ia=1,nebf
         kappa(ia)=zero
      end do

      do ie1=1,npebf
         do je1=1,npebf
c        do je1=1,ie1

C  Map from primitive BF indices to contracted indices
            iec1=AMPEB2C(ie1)
            jec1=AMPEB2C(je1)

C  Check to see if primitives belong to the same contracted BF
            if(iec1.eq.jec1) then
C  Get primitive Electron Basis Function Contraction Coefficients
               Cof_ie1=AGEBFCC(ie1)
               Cof_je1=AGEBFCC(je1)

C  Map the 2-index contracted integral to 1-D:
c              call pack_2D(nebf,iec1,iec1,ia)

C Calculate the primitive integral
               A1=ELCEX(ie1)
               I1=ELCAM(ie1,1)
               J1=ELCAM(ie1,2)
               K1=ELCAM(ie1,3)
               Amat1(1)=ELCBFC(ie1,1)
               Amat1(2)=ELCBFC(ie1,2)
               Amat1(3)=ELCBFC(ie1,3)

               B1=ELCEX(je1)
               L1=ELCAM(je1,1)
               M1=ELCAM(je1,2)
               N1=ELCAM(je1,3)
               Bmat1(1)=ELCBFC(je1,1)
               Bmat1(2)=ELCBFC(je1,2)
               Bmat1(3)=ELCBFC(je1,3)

               call gfovlap(I1,J1,K1,A1,Amat1,
     x                      L1,M1,N1,B1,Bmat1,
     x                      ans)

c              GS(ia)=GS(ia)+ans*Cof_ie1*Cof_ie1
               kappa(iec1)=kappa(iec1)+ans*Cof_ie1*Cof_je1
   
            end if  ! end if for iec1 eq jec1
         end do
      end do

C Normalize the contraction coefficients
      do i_pebf=1,npebf
         i_cebf=AMPEB2C(i_pebf)
         norm_fac=1.0d+00/sqrt(kappa(i_cebf))
         AGEBFCC(i_pebf)=AGEBFCC(i_pebf)*norm_fac
      end do

      return
      end

C======================================================================
      subroutine ELCNORM2(npebf,nebf,AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

C Normalize the contraction coefficients of electronic basis set
C======================================================================
      implicit none

      integer npebf,nebf
C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
c     integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
c     double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
c     double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
c     double precision AGNBFCC(npbf)  ! Nuclear contract coef
C-------Basis Set Info-------)
C Local variables
C--------------------------------(
C Basis set-related local variables
      integer I1,J1,K1
      integer L1,M1,N1
      double precision A1,Amat1(3)
      double precision B1,Bmat1(3)
C--------------------------------)
      integer ia,ie1,je1,iec1,jec1
      integer i_pebf
      integer i_cebf
      double precision norm_fac
c     double precision GS(nebf2)
      double precision kappa(nebf)
      double precision ans
      double precision Cof_ie1,Cof_je1
      double precision zero
      parameter(zero=0.0d+00)

C  zero out the 1-D array to hold the contracted integral overlaps
      do ia=1,nebf
         kappa(ia)=zero
      end do

      do ie1=1,npebf
c        do je1=1,npebf
c        do je1=1,ie1

C  Map from primitive BF indices to contracted indices
c           iec1=AMPEB2C(ie1)
c           jec1=AMPEB2C(je1)

C  Check to see if primitives belong to the same contracted BF
c           if(iec1.eq.jec1) then
C  Get primitive Electron Basis Function Contraction Coefficients
c              Cof_ie1=AGEBFCC(ie1)
c              Cof_je1=AGEBFCC(je1)

C  Map the 2-index contracted integral to 1-D:
c              call pack_2D(nebf,iec1,iec1,ia)

C Calculate the primitive integral
               A1=ELCEX(ie1)
               I1=ELCAM(ie1,1)
               J1=ELCAM(ie1,2)
               K1=ELCAM(ie1,3)
               Amat1(1)=ELCBFC(ie1,1)
               Amat1(2)=ELCBFC(ie1,2)
               Amat1(3)=ELCBFC(ie1,3)

c              B1=ELCEX(je1)
c              L1=ELCAM(je1,1)
c              M1=ELCAM(je1,2)
c              N1=ELCAM(je1,3)
c              Bmat1(1)=ELCBFC(je1,1)
c              Bmat1(2)=ELCBFC(je1,2)
c              Bmat1(3)=ELCBFC(je1,3)

               call gfovlap(I1,J1,K1,A1,Amat1,
     x                      I1,J1,K1,A1,Amat1,
     x                      ans)

c              GS(ia)=GS(ia)+ans*Cof_ie1*Cof_ie1
c              kappa(iec1)=kappa(iec1)+ans*Cof_ie1*Cof_ie1

               AGEBFCC(ie1)=AGEBFCC(ie1)*1.0d+00/sqrt(ans)
   
c           end if  ! end if for ie1 eq je1
c        end do
      end do

C Normalize the contraction coefficients
c     do i_pebf=1,npebf
c        i_cebf=AMPEB2C(i_pebf)
c        norm_fac=1.0d+00/sqrt(kappa(i_cebf))
c        AGEBFCC(i_pebf)=AGEBFCC(i_pebf)*norm_fac
c     end do

      return
      end

