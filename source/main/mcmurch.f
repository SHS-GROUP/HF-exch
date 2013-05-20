C ---------------- References ----------------
C
C    McMurchie & Davidson, J. Comput. Phys. 26, 218 (1978)
C
C    Persson & Taylor, Theor. Chem. Acc. 97, 240 (1997)
C
C    Mamedov, J. Math. Chem. 36, 301 (2004)
C
C --------------------------------------------

C*********************************************************************************
      subroutine TwoElecRepul(I1,J1,K1,alp1,A1,
     *                        I2,J2,K2,alp2,A2,
     *                        L1,M1,N1,beta1,B1,
     *                        L2,M2,N2,beta2,B2,
     *                        ans)
C
C Evaluates < f(1)f(2) | (1/r_12) | f(1)f(2) >
C
C*********************************************************************************

      implicit none

      integer, intent(in) :: I1,J1,K1,I2,J2,K2
      integer, intent(in) :: L1,M1,N1,L2,M2,N2
      real*8, intent(in)  :: alp1,alp2,beta1,beta2
      real*8, intent(in), dimension(3) :: A1,A2,B1,B2
      real*8, intent(out) :: ans

      integer              :: t1,u1,v1,t2,u2,v2
      integer              :: Nmax,Lmax,Mmax
      integer              :: Nboys,Lboys,Mboys
      integer              :: intbico
      real*8               :: p1,q1,p2,q2
      real*8               :: pref30,F0alp,Eherm
      real*8               :: aboys,bboys
      real*8               :: deriv,totsum
      real*8               :: pi
      real*8, dimension(3) :: P1mat,Q1mat,P2mat,Q2mat
      real*8, dimension(3) :: xKab1,xKab2
      real*8, dimension(3) :: RP1P2mat
C      real*8, allocatable  :: Ei1l1(:),Ej1m1(:),Ek1n1(:)
C      real*8, allocatable  :: Ei2l2(:),Ej2m2(:),Ek2n2(:)
C      real*8, allocatable  :: RNLM(:,:,:)
      real*8 :: Ei1l1(0:I1+L1),Ej1m1(0:J1+M1),Ek1n1(0:K1+N1)
      real*8 :: Ei2l2(0:I2+L2),Ej2m2(0:J2+M2),Ek2n2(0:K2+N2)
      real*8 :: RNLM(0:I1+L1+I2+L2,0:J1+M1+J2+M2,0:K1+N1+K2+N2)

      pi=4.0d0*datan(1.0d0)

C Generate overlap distributions for electron 1 (centred on P1==P)
C p1 == p, q1 == q
C      allocate(Ei1l1(0:I1+L1),Ej1m1(0:J1+M1),Ek1n1(0:K1+N1))
      call ovlapdistr(I1,L1,alp1,beta1,A1(1),B1(1),p1,q1,
     *             P1mat(1),Q1mat(1),xKab1(1),Ei1l1)
      call ovlapdistr(J1,M1,alp1,beta1,A1(2),B1(2),p1,q1,
     *             P1mat(2),Q1mat(2),xKab1(2),Ej1m1)
      call ovlapdistr(K1,N1,alp1,beta1,A1(3),B1(3),p1,q1,
     *             P1mat(3),Q1mat(3),xKab1(3),Ek1n1)

C Generate overlap distributions for electron 2 (centred on P2==P')
C p2 == p', q2 == q'
C      allocate(Ei2l2(0:I2+L2),Ej2m2(0:J2+M2),Ek2n2(0:K2+N2))
      call ovlapdistr(I2,L2,alp2,beta2,A2(1),B2(1),p2,q2,
     *             P2mat(1),Q2mat(1),xKab2(1),Ei2l2)
      call ovlapdistr(J2,M2,alp2,beta2,A2(2),B2(2),p2,q2,
     *             P2mat(2),Q2mat(2),xKab2(2),Ej2m2)
      call ovlapdistr(K2,N2,alp2,beta2,A2(3),B2(3),p2,q2,
     *             P2mat(3),Q2mat(3),xKab2(3),Ek2n2)

C Calculate numerical prefactor for integral over spherical Gaussians
C c.f. [McMurchie] Eqs. 3.30-31
      pref30=2.0d0*dsqrt(pi*pi*pi*pi*pi)/(p1*p2*dsqrt(p1+p2))

C Calculate parameters for Boys function (Fo)
C c.f. [McMurchie] Eqs. 3.30,32
      F0alp=p1*p2/(p1+p2)
      RP1P2mat(1)=P1mat(1)-P2mat(1)                       ! R_PP'
      RP1P2mat(2)=P1mat(2)-P2mat(2)
      RP1P2mat(3)=P1mat(3)-P2mat(3)

C Store table of Boys function (Fo) derivatives (R_{NLM})
C   N=1,...,i1+l1+i2+l2
C   L=1,...,j1+m1+j2+m2
C   M=1,...,k1+n1+k2+n2
      Nmax=i1+l1+i2+l2
      Lmax=j1+m1+j2+m2
      Mmax=k1+n1+k2+n2
C      allocate(RNLM(0:Nmax,0:Lmax,0:Mmax))
      call F0derivtab(Nmax,Lmax,Mmax,F0alp,
     *             RP1P2mat(1),RP1P2mat(2),RP1P2mat(3),RNLM)

      totsum=0.0d0

C Loop over hermite integrals
C [Persson] Eq. 19
      do t1=0,i1+l1
      do u1=0,j1+m1
      do v1=0,k1+n1

      do t2=0,i2+l2
      do u2=0,j2+m2
      do v2=0,k2+n2

C Calculate Hermite expansion coefficient appearing in [Persson] Eq. 19
        Eherm=Ei1l1(t1)*Ej1m1(u1)*Ek1n1(v1)
        Eherm=Eherm*Ei2l2(t2)*Ej2m2(u2)*Ek2n2(v2)

C Take derivative of F0 (Boys function) from R_{NLM} table with appropriate coefficients arising from chain rule
        Nboys=t1+t2
        Lboys=u1+u2
        Mboys=v1+v2
        aboys=1.0d0
        bboys=(-1.0d0)**(t2+u2+v2)
       
        deriv=aboys*bboys*RNLM(Nboys,Lboys,Mboys)

        totsum=totsum+Eherm*deriv

      end do
      end do
      end do

      end do
      end do
      end do

      ans=totsum*pref30

C      deallocate(RNLM)
C      deallocate(Ei1l1,Ej1m1,Ek1n1,Ei2l2,Ej2m2,Ek2n2)

      return
      end subroutine TwoElecRepul


      subroutine ovlapdistr(i,j,a,b,Acoord,Bcoord,
     *                      p,q,Pcoord,Qcoord,KABcoord,Eijt)

C Calculate 1-D overlap density of two Gaussian-type orbitals
C Gi(x)=x^i exp(-a*(x-Ax)^2) and Gj(x)=x^j exp(-b*(x-Bx)^2)
C
C See [Persson] Eqs. 1-17 and discussions in between

      implicit none

      integer, intent(in) :: i,j
      real*8, intent(in)  :: a,b,Acoord,Bcoord

      real*8, intent(out) :: p,q,Pcoord,Qcoord,KABcoord
      real*8, dimension(0:i+j), intent(out) :: Eijt

      p=a+b
      q=a*b/(a+b)
      Pcoord=(a*Acoord+b*Bcoord)/(a+b)
      Qcoord=Acoord-Bcoord
      KABcoord=dexp(-q*Qcoord*Qcoord)

      call calcEijt(i,j,a,b,p,q,Qcoord,KABcoord,Eijt)

      return
      end subroutine ovlapdistr


      subroutine calcEijt(imax,jmax,a,b,p,q,Qcoord,E0,Eijt)

C Calculates E_t^{ij} for t=0,...,i+j using recursion relations in [Persson]
C c.f. [Persson] Eqs. 16-17

      implicit none

      integer, intent(in) :: imax, jmax
      real*8, intent(in)  :: a,b,p,q,Qcoord,E0

      real*8, dimension(0:imax+jmax), intent(out) :: Eijt

      integer :: i, j, t
      real*8  :: coeff1,coeff2i,coeff2j,coeff3
      real*8, dimension(0:imax,0:jmax,0:imax+jmax) :: E

C Initialize since E_t^{ij} = 0 \forall t > i+j
      do i=0,imax
      do j=0,jmax
      do t=0,imax+jmax
       E(i,j,t)=0.0d0
      end do
      end do
      end do

C Calculate coefficients for recursion relations
      coeff1=1.0d0/(2.0d0*p)
      coeff2i=-q*Qcoord/a
      coeff2j=q*Qcoord/b

C Set initial value to K_AB
      E(0,0,0)=E0

C Calculate E(i,0,t) by incrementing i, then incrementing t to the maximum nonzero entry (t=i)
      j=0
      do i=0,imax-1
       E(i+1,j,0)=coeff2i*E(i,j,0)+E(i,j,1)
       do t=1,i
        coeff3=dble(t+1)
        E(i+1,j,t)=coeff1*E(i,j,t-1)+coeff2i*E(i,j,t)+coeff3*E(i,j,t+1)
       end do
       E(i+1,j,i+1)=coeff1*E(i,j,i)+coeff2i*E(i,j,i+1)
      end do

C Calculate E(i,j,t) for each i by incrementing j, then t to the maximum nonzero entry (t=i+j)
      do i=0,imax
       do j=0,jmax-1
        E(i,j+1,0)=coeff2j*E(i,j,0)+E(i,j,1)
        do t=1,i+j
         coeff3=dble(t+1)
         E(i,j+1,t)=coeff1*E(i,j,t-1)+coeff2j*E(i,j,t)+coeff3*E(i,j,t+1)
        end do
        E(i,j+1,i+j+1)=coeff1*E(i,j,i+j)+coeff2j*E(i,j,i+j+1)
       end do
      end do

C Generate output
      do t=0,imax+jmax
        Eijt(t)=E(imax,jmax,t)
      end do

      return
      end subroutine calcEijt


      subroutine F0derivtab(Nmax,Lmax,Mmax,alp,a,b,c,RNLM)

C Calculate table of Boys function (F0) derivatives (R_{NLM}) where
C   N=1,...,Nmax
C   L=1,...,Lmax
C   M=1,...,Mmax
C using recursion relations in [McMurchie] Eqs. 4.4-4.8

      implicit none

      integer, intent(in) :: Nmax, Lmax, Mmax
      real*8, intent(in)  :: alp, a, b, c
      real*8, intent(out) :: RNLM(0:Nmax,0:Lmax,0:Mmax)

      integer :: j, n, l, m, Jmax
      real*8  :: T, Fj
      real*8  :: RNLMj(0:Nmax,0:Lmax,0:Mmax,0:Nmax+Lmax+Mmax)

      T    = alp*((a*a)+(b*b)+(c*c))
      Jmax = Nmax + Lmax + Mmax 

C Initialize all {RNLM,Jmax} = 0 since unspecified by recursion relations in [McMurchie]
      do n=0,Nmax
      do l=0,Lmax
      do m=0,Mmax
       RNLMj(n,l,m,Jmax)=0.0d0
      end do
      end do
      end do

C Calculate all R000,j from Fj(T)
      do j=0,Jmax
       call iboysA(j,T,Fj)
       RNLMj(0,0,0,j) = ((-2.0d0*alp)**j)*Fj
      end do

C Calculate all other RNLM,j

      m=0
      if (Mmax.ge.1) then

C {R001,j}
       do j=0,Jmax-1
         RNLMj(0,0,m+1,j) = c*RNLMj(0,0,m,j+1)
       end do

      end if

      l=0
      if (Lmax.ge.1) then

C {R010,j} and {R011,j}
       do j=0,Jmax-1
        RNLMj(0,l+1,m,j) = b*RNLMj(0,l,m,j+1)
       end do

       if (Mmax.ge.1) then

        do j=0,Jmax-1
         RNLMj(0,l+1,m+1,j) = b*RNLMj(0,l,m+1,j+1)
        end do

       end if

      end if

      n=0
      if (Nmax.ge.1) then

C {R100,j}, {R110,j}, {R101,j} and {R111,j}
       do j=0,Jmax-1
        RNLMj(n+1,l,m,j) = a*RNLMj(n,l,m,j+1)
       end do

       if ((Lmax.ge.1).and.(Mmax.ge.1)) then

        do j=0,Jmax-1
         RNLMj(n+1,l+1,m+1,j) = a*RNLMj(n,l+1,m+1,j+1)
         RNLMj(n+1,l+1,m,j) = a*RNLMj(n,l+1,m,j+1)
         RNLMj(n+1,l,m+1,j) = a*RNLMj(n,l,m+1,j+1)
        end do

       else if ((Lmax.ge.1).and.(Mmax.eq.0)) then

        do j=0,Jmax-1
         RNLMj(n+1,l+1,m,j) = a*RNLMj(n,l+1,m,j+1)
        end do

       else if ((Lmax.eq.0).and.(Mmax.ge.1)) then

        do j=0,Jmax-1
         RNLMj(n+1,l,m+1,j) = a*RNLMj(n,l,m+1,j+1)
        end do

       end if

      end if


C {Rn00,j}, {Rn10,j}, {Rn01,j} and {Rn11,j} for [2 <= n <= Nmax]
      if ((Lmax.ge.1).and.(Mmax.ge.1)) then

       do n=1,Nmax-1

        do j=0,Jmax-1
         RNLMj(n+1,l,m,j)       = a*RNLMj(n,l,m,j+1)+
     *                            dble(n)*RNLMj(n-1,l,m,j+1)
         RNLMj(n+1,l+1,m,j)     = a*RNLMj(n,l+1,m,j+1)+
     *                            dble(n)*RNLMj(n-1,l+1,m,j+1)
         RNLMj(n+1,l,m+1,j)     = a*RNLMj(n,l,m+1,j+1)+
     *                            dble(n)*RNLMj(n-1,l,m+1,j+1)
         RNLMj(n+1,l+1,m+1,j)   = a*RNLMj(n,l+1,m+1,j+1)+
     *                            dble(n)*RNLMj(n-1,l+1,m+1,j+1)
        end do

       end do

      else if ((Lmax.ge.1).and.(Mmax.eq.0)) then

       do n=1,Nmax-1

        do j=0,Jmax-1
         RNLMj(n+1,l,m,j)       = a*RNLMj(n,l,m,j+1)+
     *                            dble(n)*RNLMj(n-1,l,m,j+1)
         RNLMj(n+1,l+1,m,j)     = a*RNLMj(n,l+1,m,j+1)+
     *                            dble(n)*RNLMj(n-1,l+1,m,j+1)
        end do

       end do


      else if ((Lmax.eq.0).and.(Mmax.ge.1)) then

       do n=1,Nmax-1

        do j=0,Jmax-1
         RNLMj(n+1,l,m,j)       = a*RNLMj(n,l,m,j+1)+
     *                            dble(n)*RNLMj(n-1,l,m,j+1)
         RNLMj(n+1,l,m+1,j)     = a*RNLMj(n,l,m+1,j+1)+
     *                            dble(n)*RNLMj(n-1,l,m+1,j+1)
        end do

       end do

      else

       do n=1,Nmax-1

        do j=0,Jmax-1
         RNLMj(n+1,l,m,j)       = a*RNLMj(n,l,m,j+1)+
     *                            dble(n)*RNLMj(n-1,l,m,j+1)
        end do

       end do

      end if


      if ((Nmax.ge.1).and.(Mmax.ge.1)) then

       do l=1,Lmax-1

C {R0l0,j} and {R0l1,j} for [2 <= l <= Lmax]
        do j=0,Jmax-1
         RNLMj(0,l+1,m,j)     = b*RNLMj(0,l,m,j+1)+
     *                          dble(l)*RNLMj(0,l-1,m,j+1)
         RNLMj(0,l+1,m+1,j)   = b*RNLMj(0,l,m+1,j+1)+
     *                          dble(l)*RNLMj(0,l-1,m+1,j+1)
        end do

        n=0
C {R1l0,j} and {R1l1,j} for [2 <= l <= Lmax]
        do j=0,Jmax-1
         RNLMj(n+1,l+1,m,j)     = a*RNLMj(n,l+1,m,j+1)
         RNLMj(n+1,l+1,m+1,j)   = a*RNLMj(n,l+1,m+1,j+1)
        end do

C {Rnl0,j} and {Rnl1,j} for [2 <= n <= Nmax] & [2 <= l <= Lmax]
        do n=1,Nmax-1

         do j=0,Jmax-1
          RNLMj(n+1,l+1,m,j)     = a*RNLMj(n,l+1,m,j+1)+
     *                             dble(n)*RNLMj(n-1,l+1,m,j+1)
          RNLMj(n+1,l+1,m+1,j)   = a*RNLMj(n,l+1,m+1,j+1)+
     *                             dble(n)*RNLMj(n-1,l+1,m+1,j+1)
         end do

        end do

       end do

      else if ((Nmax.ge.1).and.(Mmax.eq.0)) then

       do l=1,Lmax-1

C {R0l0,j} for [2 <= l <= Lmax]
        do j=0,Jmax-1
         RNLMj(0,l+1,m,j)     = b*RNLMj(0,l,m,j+1)+
     *                          dble(l)*RNLMj(0,l-1,m,j+1)
        end do

        n=0
C {R1l0,j} for [2 <= l <= Lmax]
        do j=0,Jmax-1
         RNLMj(n+1,l+1,m,j)     = a*RNLMj(n,l+1,m,j+1)
        end do

C {Rnl0,j} for [2 <= n <= Nmax] & [2 <= l <= Lmax]
        do n=1,Nmax-1

         do j=0,Jmax-1
          RNLMj(n+1,l+1,m,j)     = a*RNLMj(n,l+1,m,j+1)+
     *                             dble(n)*RNLMj(n-1,l+1,m,j+1)
         end do

        end do

       end do

      else if ((Nmax.eq.0).and.(Mmax.ge.1)) then

       do l=1,Lmax-1

C {R0l0,j} for [2 <= l <= Lmax]
        do j=0,Jmax-1
         RNLMj(0,l+1,m,j)     = b*RNLMj(0,l,m,j+1)+
     *                          dble(l)*RNLMj(0,l-1,m,j+1)
         RNLMj(0,l+1,m+1,j)   = b*RNLMj(0,l,m+1,j+1)+
     *                          dble(l)*RNLMj(0,l-1,m+1,j+1)
        end do

       end do


      else

       do l=1,Lmax-1

C {R0l0,j} for [2 <= l <= Lmax]
        do j=0,Jmax-1
         RNLMj(0,l+1,m,j)     = b*RNLMj(0,l,m,j+1)+
     *                          dble(l)*RNLMj(0,l-1,m,j+1)
        end do

       end do

      end if


      if ((Nmax.ge.1).and.(Lmax.ge.1)) then

       do m=1,Mmax-1

C {R00m,j} for [2 <= m <= Mmax]
        do j=0,Jmax-1
         RNLMj(0,0,m+1,j)   = c*RNLMj(0,0,m,j+1)+
     *                        dble(m)*RNLMj(0,0,m-1,j+1)
        end do

        l=0
C {R01m,j} for [2 <= m <= Mmax]
        do j=0,Jmax-1
         RNLMj(0,l+1,m+1,j)   = b*RNLMj(0,l,m+1,j+1)
        end do

        n=0
C {R10m,j} and {R11m,j} for [2 <= m <= Mmax]
        do j=0,Jmax-1
         RNLMj(n+1,l,m+1,j)     = a*RNLMj(n,l,m+1,j+1)
         RNLMj(n+1,l+1,m+1,j)   = a*RNLMj(n,l+1,m+1,j+1)
        end do

C {Rn0m,j} and {Rn1m,j} for [2 <= n <= Nmax] & [2 <= m <= Mmax]
        do n=1,Nmax-1

         do j=0,Jmax-1
          RNLMj(n+1,l,m+1,j)     = a*RNLMj(n,l,m+1,j+1)+
     *                             dble(n)*RNLMj(n-1,l,m+1,j+1)
          RNLMj(n+1,l+1,m+1,j)   = a*RNLMj(n,l+1,m+1,j+1)+
     *                             dble(n)*RNLMj(n-1,l+1,m+1,j+1)
         end do

        end do

        do l=1,Lmax-1

C {R0lm,j} for [2 <= l <= Lmax] & [2 <= m <= Mmax]
         do j=0,Jmax-1
          RNLMj(0,l+1,m+1,j)   = b*RNLMj(0,l,m+1,j+1)+
     *                           dble(l)*RNLMj(0,l-1,m+1,j+1)
         end do

         n=0
C {R1lm,j} for [2 <= l <= Lmax] & [2 <= m <= Mmax]
         do j=0,Jmax-1
          RNLMj(n+1,l+1,m+1,j)   = a*RNLMj(n,l+1,m+1,j+1)
         end do

C {Rnlm,j} for [2 <= n <= Nmax] & [2 <= l <= Lmax] & [2 <= m <= Mmax]
         do n=1,Nmax-1

          do j=0,Jmax-1
           RNLMj(n+1,l+1,m+1,j)   = a*RNLMj(n,l+1,m+1,j+1)+
     *                              dble(n)*RNLMj(n-1,l+1,m+1,j+1)
          end do

         end do

        end do

       end do

      else if ((Nmax.ge.1).and.(Lmax.eq.0)) then

       do m=1,Mmax-1

C {R00m,j} for [2 <= m <= Mmax]
        do j=0,Jmax-1
         RNLMj(0,0,m+1,j)   = c*RNLMj(0,0,m,j+1)+
     *                        dble(m)*RNLMj(0,0,m-1,j+1)
        end do

        l=0
        n=0
C {R10m,j} for [2 <= m <= Mmax]
        do j=0,Jmax-1
         RNLMj(n+1,l,m+1,j)     = a*RNLMj(n,l,m+1,j+1)
        end do

C {Rn0m,j} for [2 <= n <= Nmax] & [2 <= m <= Mmax]
        do n=1,Nmax-1

         do j=0,Jmax-1
          RNLMj(n+1,l,m+1,j)     = a*RNLMj(n,l,m+1,j+1)+
     *                             dble(n)*RNLMj(n-1,l,m+1,j+1)
         end do

        end do

       end do

      else if ((Nmax.eq.0).and.(Lmax.ge.1)) then

       do m=1,Mmax-1

C {R00m,j} for [2 <= m <= Mmax]
        do j=0,Jmax-1
         RNLMj(0,0,m+1,j)   = c*RNLMj(0,0,m,j+1)+
     *                        dble(m)*RNLMj(0,0,m-1,j+1)
        end do

        l=0
C {R01m,j} for [2 <= m <= Mmax]
        do j=0,Jmax-1
         RNLMj(0,l+1,m+1,j)   = b*RNLMj(0,l,m+1,j+1)
        end do

        do l=1,Lmax-1

C {R0lm,j} for [2 <= l <= Lmax] & [2 <= m <= Mmax]
         do j=0,Jmax-1
          RNLMj(0,l+1,m+1,j)   = b*RNLMj(0,l,m+1,j+1)+
     *                           dble(l)*RNLMj(0,l-1,m+1,j+1)
         end do

        end do

       end do

      else

       do m=1,Mmax-1

C {R00m,j} for [2 <= m <= Mmax]
        do j=0,Jmax-1
         RNLMj(0,0,m+1,j)   = c*RNLMj(0,0,m,j+1)+
     *                        dble(m)*RNLMj(0,0,m-1,j+1)
        end do

       end do

      end if


C Store RNLM,0 in RNLM
      do n=0,Nmax
      do l=0,Lmax
      do m=0,Mmax
        RNLM(n,l,m)=RNLMj(n,l,m,0)
      end do
      end do
      end do

      return
      end subroutine F0derivtab


      subroutine incompgamm(m,T,F)

C Computes incomplete gamma function Fm(T) using methods
C discussed in [Mamedov] along with downward recursion where
C
C               /1  2m (-Tx^2)   
C     Fm(T)  =  |   x  e       dx
C               /0               

      implicit none

      integer, intent(in) :: m
      real*8, intent(in)  :: T
      real*8, intent(out) :: F(0:m)

C Use suggested value for Tmax in [Mamedov]
      integer, parameter :: Tmax=35, itermax=1000
      real*8, parameter  :: accy=1.0d-12
      integer :: i
      real*8 :: pi, denomfac, term, summ

      pi=4.0d0*datan(1.0d0)

C Use [Mamedov] Eqs. 7-8 when T > Tmax
      if (T.ge.Tmax) then

       F(0)=0.5d0*dsqrt(pi/T)
       do i=1,m
C Use i-1/2 instead of 2i-1 to mimic what was done previously (checked with Mathematica)
C         F(i)=dble(2*i-1)*F(i-1)/T
         F(i)=dble(i-0.5d0)*F(i-1)/T
       end do

C Use [Mamedov] Eqs. 6&3 when T < Tmax
      else

       denomfac=dble(2*m+1)
       term=1.0d0/denomfac
       summ=term
       do i=1,itermax
         term=term*2.0d0*T/(denomfac+2.0d0*dble(i))
         summ=summ+term
         if (term.le.accy) goto 10
       end do

   10  continue
       F(m)=dexp(-T)*summ

       do i=m-1,0,-1
         F(i)=(2.0d0*T*F(i+1)+dexp(-T))/dble(2*i+1)
       end do

      end if

      return
      end subroutine incompgamm


      subroutine iboysA(m,x,Fmx)

C Dummy subroutine to calculate F_m(x) where m is integer

      implicit none

      integer, intent(in) :: m
      real*8, intent(in)  :: x
      real*8, intent(out) :: Fmx

      real*8 :: F(0:M)

      call incompgamm(m,x,F)
      Fmx = F(m)

      return
      end subroutine iboysA


      subroutine Gderivtab(Nmax,Lmax,Mmax,alp,a,b,c,GNLM)

C Calculate table of a Gaussian function (G) derivatives (G_{NLM}) where
C   N=1,...,Nmax
C   L=1,...,Lmax
C   M=1,...,Mmax
C for
C     G = Exp[-alp*[a^2+b^2+c^2]]
C using direct approach

      implicit none

      integer, intent(in) :: Nmax, Lmax, Mmax
      real*8, intent(in)  :: alp, a, b, c
      real*8, intent(out) :: GNLM(0:Nmax,0:Lmax,0:Mmax)

      real*8, parameter :: zero=0.0d0
      integer :: n, l, m
      real*8 :: nderiv, lderiv, mderiv

      do n=0,Nmax
        call gderiv_Ax(n,zero,alp,a,nderiv)

        do l=0,Lmax
          call gderiv_Ax(l,zero,alp,b,lderiv)

          do m=0,Mmax
            call gderiv_Ax(m,zero,alp,c,mderiv)

            GNLM(n,l,m)=nderiv*lderiv*mderiv

          end do
        end do
      end do

      return
      end subroutine Gderivtab


      subroutine gderiv_Ax(N,x,alp,Ax,gderAx)

C This routine computes the derivative of Gaussian func
C               d     
C     gd0   = (---)^N Exp[-alp*(x-Ax)^2]
C              dAx    
C
C using substitution
C     y = sqrt(alp)*(x-Ax)
C
C properties of HermiteA polynomials
C     (d^N/dy^N) Exp[-y^2] = (-1)^N * H_n(y) * Exp[-y^2]
C
C and chain rule since dy/dAx is constant

      implicit none

      integer, intent(in) :: N
      real*8, intent(in)  :: x, alp, Ax
      real*8, intent(out) :: gderAx

      real*8 :: y,dy,de,HermiteAH

C Leave out (-1)^N factors from dy and de since they cancel with each other

C y and its Nth derivative wrt Ax
      y=dsqrt(alp)*(x-Ax)
      dy=(dsqrt(alp))**N
     
C Nth derivative of Exp[-y^2] wrt y
      de=HermiteAH(N,y)*dexp(-y*y)

      gderAx=dy*de

      return
      end subroutine gderiv_Ax


      function HermiteAH(n,x)

C Returns H_n(x) where H_n is the nth HermiteA polynomial
C H_n(x) = (-1)^n * Exp[x^2] * (d^n/dx^n) Exp[-x^2]

      implicit none

      integer, intent(in) :: n
      real*8, intent(in)  :: x

      real*8 :: HermiteAH

      integer :: i
      real*8  :: H(0:n)

      H(0) = 1.0d0

      if (n.gt.0) then
       H(1) = 2.0d0*x

       if (n.gt.1) H(2) = (4.0d0*x*x)-2.0d0

       do i=2,n-1
          H(i+1)=(2.0d0*x*H(i))-(2.0d0*dble(i)*H(i-1))
       end do

      end if

      HermiteAH = H(n)

      return
      end function HermiteAH


      function intbico(n,k)
C Returns the binomial coefficient (n  k) as an integer using standard recursive relation
C (n+1  k+1) = (n+1)/(k+1) * (n  k)
C Assumes n >= k
      implicit none

      integer, intent(in) :: n,k
      integer :: i,ans,intbico

      ans=1.0d0
      do i=1,k
        ans=ans*(n-k+i)
        ans=ans/i
      end do
      intbico=ans

      return
      end function intbico

