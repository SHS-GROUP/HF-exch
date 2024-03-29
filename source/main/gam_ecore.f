C======================================================================
      subroutine calc_GAM_ecore(nebf,npebf,nebf2,
     x                          nat,zan,cat,
     x                          AMPEB2C,AGEBFCC,
     x                          ELCEX,ELCAM,ELCBFC)

C======================================================================
      implicit none
C Input Variables
      integer npebf,nebf,nebf2,nat
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
C-------Basis Set Info-------)
C Variables Returned
      double precision GAM_ecore(nebf2)
C Local variables
C--------------------------------(
C Basis set-related local variables
      integer I1,J1,K1
      integer L1,M1,N1
      double precision A1,Amat1(3)
      double precision B1,Bmat1(3)
C--------------------------------)
      integer ia,ip,jp,ie1,je1,iec1,jec1
      double precision ans
      double precision xmass
      double precision coul_sign
      double precision Cof_ie1,Cof_je1
      double precision Cof_ip,Cof_jp
      double precision zero
      parameter(zero=0.0d+00)

      write(*,*)
      write(*,*)'**************************************'
      write(*,*)'    CALCULATING ELEC CORE INTEGRALS   '
      write(*,*)'**************************************'
      write(*,*)
c     call flshbf(6)

C  zero out the 1-D array to hold the contracted integrals
      do ia=1,nebf2
         GAM_ecore(ia)=zero
      end do

      do ie1=1,npebf
         do je1=1,npebf

C  Map from primitive BF indices to contracted indices
c            call MPEB2C(ie1,iec1)
c            call MPEB2C(je1,jec1)
            iec1=AMPEB2C(ie1)
            jec1=AMPEB2C(je1)

C  Get primitive Electron Basis Function Contraction Coefficients
c               call GEBFCC(ie1,Cof_ie1)
c               call GEBFCC(je1,Cof_je1)
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

c           write(*,*)'ia=',ia
c           write(*,*)'Cof_ie1=',Cof_ie1
c           write(*,*)'Cof_je1=',Cof_je1

            coul_sign=-1.0d+00
            xmass=1.0d+00
            call xcalc_GAM_epcore(I1,J1,K1,A1,Amat1,
     x                            L1,M1,N1,B1,Bmat1,
     x                            nat,xmass,cat,zan,coul_sign,ans)

            GAM_ecore(ia)=GAM_ecore(ia)+ans*Cof_ie1*Cof_je1
c           write(*,*)'ecore_integral=',ans
c           write(*,*)'GAM_ecore=',GAM_ecore(ia)
c           write(*,*)'>>>>>>>>>>>>>>>>>>>>>>>>'

         end do
      end do

C Write electron core integrals to file
      open(812,file='GAM_ecr.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do iec1=1,nebf
         do jec1=1,nebf

C  Map the 2-index contracted integral to 1-D:
            call pack_2D(nebf,jec1,iec1,ia)

            write(812,REC=ia) GAM_ecore(ia)

         end do
      end do

      close(812)

      return
      end

C======================================================================
      subroutine xcalc_GAM_epcore(I1,J1,K1,A1,Amat1,
     x                            L1,M1,N1,B1,Bmat1,
     x                            nat,xmass,cat,zan,coul_sign,ans)

C======================================================================
      implicit none
C Input Variables
      integer nat
      integer I1,J1,K1
      integer L1,M1,N1
      double precision A1,Amat1(3)
      double precision B1,Bmat1(3)
      double precision xmass
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision zan(nat) ! Classical nuclear charges
      double precision coul_sign
C Variables Returned
      double precision ans
C Local Variables
      integer iii
      double precision ans_ke
      double precision ans_vec
      double precision xans
      double precision Cmat(3)
      double precision ZNUC 

      logical debug


C     debug=.true.
      debug=.false.


      call gfke(I1,J1,K1,A1,Amat1,
     2          L1,M1,N1,B1,Bmat1,
     3          xmass,ans_ke)


      ans_vec=0.d+00

      DO III=1,NAT
         Cmat(1)=cat(1,III)
         Cmat(2)=cat(2,III)
         Cmat(3)=cat(3,III)
         ZNUC=ZAN(III)

         call gfvec(I1,J1,K1,A1,Amat1,
     2              L1,M1,N1,B1,Bmat1,
     3              Cmat,xans)

         ans_vec = ans_vec + coul_sign * znuc * xans

      END DO

      ans=ans_ke+ans_vec

C>>>>>>>>>> debug
      if(debug) then
         write(*,*)'>>> GAM_epcore <<<'
c        write(*,*)'ijkl=',ie1,je1,ip,jp
         write(*,*)'I1 J1 K1 Alpha1 Amat1',I1,J1,K1,A1,Amat1
         write(*,*)'L1 M1 N1 Beta1 Bmat1 ',L1,M1,N1,B1,Bmat1
         write(*,*)'ans=',ans
         write(*,*)
      end if

      return
      end

