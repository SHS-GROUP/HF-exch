C----------DIRECT-ACCESS-FILES-FOR-XCHF--------------------------------(
c     ----------XC-NEO------------
c          801  ::   GAM_1.ufm  
c          802  ::   GAM_1s.ufm 
c          803  ::   GAM_2.ufm  
c          804  ::   GAM_2s.ufm 
c          805  ::   GAM_3.ufm  
c          806  ::   GAM_4.ufm  
c     -------STANDARD-NEO---------
c          810  ::   GAM_ee.ufm 
c          811  ::   GAM_ep.ufm 
c          812  ::   GAM_ecr.ufm
c          813  ::   GAM_pcr.ufm
c          814  ::   eovlap.ufm 
c          815  ::   novlap.ufm 
C----------DIRECT-ACCESS-FILES-FOR-XCHF--------------------------------)
C
C          800  ::   ENUCRP.dat  classical nuclear repulsion energy
C
C---------FORMATTED-FILES-FOR-MO-COEFFICIENTS--------------------------(
c          850  ::   guessCE.inp               
c          851  ::   guessCP.inp               
c          852  ::   FinalCE.dat
c          853  ::   FinalCP.dat
C---------FORMATTED-FILES-FOR-MO-COEFFICIENTS--------------------------)
c
C=======================================================================
      subroutine read_GAM_ee(ne,ngee,GAM_ee)
C=======================================================================
      implicit none
C Input Variables
      integer ne,ngee
C Variables Returned
      double precision GAM_ee(ngee)
C Local Variables
      integer ia
      integer ip,jp
      integer ie1,je1
      integer ie2,je2
      double precision ans

      open(810,file='GAM_ee.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do ie1=1,ne
       do je1=1,ne
        do ie2=1,ne
         do je2=1,ne

            call pack_4D(ne,ne,ne,
     x                   je2,ie2,je1,ie1,ia)

            read(810,REC=ia) ans

            GAM_ee(ia)=ans

         end do
        end do
       end do
      end do

      close(810)


      return
      end
C=======================================================================
      subroutine read_GAM_ecore(ne,ne2,GAM_ecore)
C=======================================================================
      implicit none
C Input Variables
      integer ne,ne2
C Variables Returned
      double precision GAM_ecore(ne2)
C Local Variables
      integer ia
      integer ie,je
      double precision ans

      open(812,file='GAM_ecr.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do ie=1,ne
         do je=1,ne

C  Map the 2-index contracted integral to 1-D:
            call pack_2D(ne,je,ie,ia)

            read(812,REC=ia) GAM_ecore(ia)

         end do
      end do

      close(812)



      return
      end
C=======================================================================
      subroutine read_elec_ovlap(ne,GS)
C=======================================================================
      implicit none
C Input Variables
      integer ne
C Variables Returned
      double precision GS(ne,ne)
C Local Variables
      integer ia
      integer ie1,je1
      double precision ans

      open(814,file='eovlap.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do ie1=1,ne
         do je1=1,ne

C  Map the 2-index contracted integral to 1-D:
            call pack_2D(ne,je1,ie1,ia)

            read(814,REC=ia) ans 
            GS(ie1,je1)=ans

         end do
      end do

      close(814)



      return
      end
C=======================================================================
      subroutine write_MOs(IFIL,nbf,VEC)
C     IFIL=852 :: FinalCE.dat
C     IFIL=853 :: FinalCP.dat
C     IFIL=860 :: FinalCAE.dat
C     IFIL=861 :: FinalCBE.dat
C     IFIL=870 :: FinalCAalpE.dat
C     IFIL=871 :: FinalCAbetE.dat
C     IFIL=880 :: FinalCBalpE.dat
C     IFIL=881 :: FinalCBbetE.dat
C=======================================================================
      implicit none
C Input Variables
      integer IFIL,nbf
      double precision VEC(nbf,nbf)
C Variables Returned
C Local Variables

      call PUSQLF(IFIL,VEC,nbf,nbf,nbf)


      return
      end
C*MODULE MTHLIB  *DECK PUSQLF
      SUBROUTINE PUSQLF(LUFILE,V,M,N,LDV)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
c     LOGICAL GOPARR,DSKWRK,MASWRK
C
      DIMENSION V(LDV,M)
C
c     COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     ----- PUNCH A RECTANGULAR MATRIX WITH ORDERING LABELS -----
C     -V- IS -N- ROWS BY -M- COLUMNS, WITH TRUE LEAD DIMENSION -LDV-
C     NOTE THAT -PUSQLF- IS AN EXACT CLONE OF -PUSQL- EXCEPT A UNIT
C     NUMBER IS TO BE GIVEN AS AN ARGUMENT.
C
c     IF (.NOT.MASWRK) RETURN

      if(LUFILE.eq.852) then
         open(LUFILE,file='FinalCE.dat',status='unknown')
      else if(LUFILE.eq.853) then
         open(LUFILE,file='FinalCP.dat',status='unknown')
      else if(LUFILE.eq.860) then
         open(LUFILE,file='FinalCAE.dat',status='unknown')
      else if(LUFILE.eq.861) then
         open(LUFILE,file='FinalCBE.dat',status='unknown')
      else if(LUFILE.eq.870) then
         open(LUFILE,file='FinalCAalpE.dat',status='unknown')
      else if(LUFILE.eq.871) then
         open(LUFILE,file='FinalCAbetE.dat',status='unknown')
      else if(LUFILE.eq.880) then
         open(LUFILE,file='FinalCBalpE.dat',status='unknown')
      else if(LUFILE.eq.881) then
         open(LUFILE,file='FinalCBbetE.dat',status='unknown')
      end if
C
      DO J = 1,M
         IC = 0
         MAX = 0
  100    CONTINUE
            MIN = MAX+1
            MAX = MAX+5
            IC = IC+1
            IF (MAX .GT. N) MAX = N
            MODJ  = MOD(J ,100 )
            MODIC = MOD(IC,1000)
            WRITE (LUFILE,9008) MODJ,MODIC,(V(I,J),I = MIN,MAX)
         IF (MAX.LT.N) GO TO 100
      ENDDO
      close(LUFILE)
      RETURN
C
 9008 FORMAT(I2,I3,1P,5E15.8)
      END

C======================================================================
      subroutine prt_lower_triangle(nbf,nbfLT,mat)
C======================================================================
      implicit none
      integer nbf
      integer nbfLT
      double precision mat(nbf,nbf)
      double precision mat_LT(nbfLT)
      integer i
      integer j
      integer ia

C Pack Fock into lower triangle
      ia=0
      do i=1,nbf
         do j=1,i
            ia=ia+1
            mat_LT(ia)=mat(i,j)
         end do
      end do

c     CALL PRTRIL(fockLT,nbf)
      CALL PRTRI(mat_LT,nbf)


      return
      end

C*MODULE MTHLIB  *DECK PRTRI
      SUBROUTINE PRTRI(D,N)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      DIMENSION D(*)
C
c     COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
c     COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
c     COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     ----- PRINT SYMMETRIC MATRIX -D- OF DIMENSION -N- -----
C
c     IF (MASWRK) THEN
      MAX = 5
c     IF (NPRINT .EQ. 6) MAX = 10
      MM1 = MAX-1
      DO 120 I0=1,N,MAX
         IL = MIN(N,I0+MM1)
         WRITE(*,9008)
         WRITE(*,9028) (I,I=I0,IL)
         WRITE(*,9008)
         IL = -1
         DO 100 I=I0,N
            IL=IL+1
            J0=I0+(I*I-I)/2
            JL=J0+MIN(IL,MM1)
            WRITE(*,9048) I,(D(J),J=J0,JL)
  100    CONTINUE
  120 CONTINUE
c     END IF
      RETURN
 9008 FORMAT(1X)
 9028 FORMAT(6X,10(4X,I4,4X))
 9048 FORMAT(I5,1X,10F12.7)
      END

c     SUBROUTINE PREVNU(V,E,M,LDV,ISTART,IEND)
      SUBROUTINE PREVNU(V,E,M,LDV,NBF)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
c     LOGICAL GOPARR,DSKWRK,MASWRK
c
c     parameter (MXAO=8192)
C
c     CHARACTER*8 QNUN,QNN
c     CHARACTER*10 PBFLAB
C
      DIMENSION V(LDV,M),E(M)
C
c     COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
c     COMMON /NUCMON/ QNUN(20),QNN(20),PBFLAB(MXAO)
c     COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
c     COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     ----- PRINT OUT EIGENDATA (VECTORS AND VALUES) -----
C     THE ROWS WILL BE LABELED WITH THE BASIS FUNCTION TAGS.
C     -V- IS N X M, WITH TRUE LEADING DIMENSION -LDV-
C
c     IF (MASWRK) THEN
      MAX = 5
c     IF (NPRINT .EQ. 6) MAX = 10
      IMAX = 0
C
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. M) IMAX = M
      WRITE (*,9008)
      WRITE (*,9028) (I,I = IMIN,IMAX)
      WRITE (*,9008)
      WRITE (*,9068) (E(I),I = IMIN,IMAX)
      WRITE (*,9008)
      J = 0
c     DO 120 IDX = ISTART,IEND
      DO 120 IDX = 1,NBF
         J = J + 1
c        WRITE (IW,9048) J,PBFLAB(IDX),(V(J,I),I = IMIN,IMAX)
         WRITE (*,9048) J,(V(J,I),I = IMIN,IMAX)
  120 CONTINUE
      IF (IMAX .LT. M) GO TO 100
c     ENDIF
      RETURN
C
 9008 FORMAT(1X)
 9028 FORMAT(17X,10(4X,I4,3X))
c9048 FORMAT(I5,2X,A10,10F11.6)
 9048 FORMAT(I5,2X,'LABEL     ',10F11.6)
 9068 FORMAT(17X,10F11.4)
      END

