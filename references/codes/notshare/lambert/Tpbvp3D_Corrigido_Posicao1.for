! This program uses a gradient method to solve the TPBVP in the 3-body problem in 3-D. Initial guess are in the rotating system. Final conditions are in the rotating System (SC=1) or in the fixed system (SC=0). Velocity in the initial guess is in the rotating frame!!

      PROGRAM Tpbvp3D 
      IMPLICIT REAL*8(A-H,O-Z)
      Real*8 ms, mj, mu, Pi, PerSJ, d_SJm, d_STm, d_SMm 
      COMMON /JAC/ U
	COMMON /FILTER/ FCC
	COMMON /START/ X0,Y0,Z0
	COMMON /FIM/ XFF,YF,ZF
	COMMON /SOL/ XRF
        EXTERNAL DERIVS,RKQC
      DIMENSION XF(6),XR(6),F(3),V(3),DELV(3),DV(3),XRF(6),XRI(6)
     *,XRFF(6),XFI(6),XFFF(6), Xfixo(10000,3), Xrota(10000,3)
      DIMENSION DelTT(10000), DelvTF(10000), DelvTR(10000)
      OPEN (UNIT=1, FILE='NEW_inp3d.dat') 
      OPEN (UNIT=2, FILE='FIX.DAT')
        OPEN (UNIT=3, FILE='ROT.DAT')
	OPEN (UNIT=5, FILE='JAC.DAT')
	OPEN (UNIT=6, FILE='ORBIT.DAT')
	OPEN (UNIT=7, FILE='DELTAV.DAT')

        OPEN (UNIT=4, FILE='MEUS_DADOS_DE_ENTRADA.DAT')
	OPEN (UNIT=8, FILE='DT_DV_Orbitas.DAT')
        OPEN (UNIT=10, FILE='xyz_Orbitas_SistFixo.DAT')
        OPEN (UNIT=11, FILE='xyz_Orbitas_SistGirante.DAT')
        OPEN (UNIT=15, FILE='DT_DELVTOTR_adimensionais.DAT')
!*******************************************************************
      ms = 1988500D+24    ! Sun mass (kg)
      mj = 1898.19D+24    ! Jupiter mass (kg)
      mu = mj/(ms+mj)
     
      Pi = (4.0D0)*Datan(1.0D0)
!!PerSJ = (4332.589D0*86400.0D0)  ! Período de Júpiter em segundos (Per. em Dias * Segs./Dia)
      PerSJ = 11.862615D0*(365.25636D0/1.0000174D0)*86400.0D0  ! Em segundos
! Periodo sideral orbital de Jupiter = 11.862615 anos
! Periodo sideral orbital da Terra = 365.25636 dias = 1.0000174 anos
! 24*60*60 = 86400 segundos por dia
      PerSJanos = 11.862615D0

      d_SJm = 7.7857D+8         ! Distância média Sol-Júpiter (km)
      d_STm = 1.495978707D+8    ! Distância média Sol-Terra (km)
      d_SMm = 2.27925D+8        ! Distância média Sol-Marte (km)
!*******************************************************************
44       FORMAT(I,3E9.2)
1010     FORMAT(6x,10F14.11)
33       FORMAT(5x,3F14.10)    
!*******************************************************************
        U = mu
C
C	Read: number of steps in the plot, Prec. for int., filter,
C	tolerance in final numbers.
C
      READ (1,44) NP,EPS,FCC,TOLF
C
C	Read: X,Y,Z initial (fixed S.); X,Y,Z final (rot. S.), Initial time, 
C	minimum transfer time, maximum transfer time, delta transfer time.
C
	READ(1,1010) X0,Y0,Z0,XFF,YF,ZF,T00,TF1,TF2,DTFF
C
C	Read: X,Y,Z initial guess for velocity (Fixed S.), regular. final time
C
	READ(1,33) XV0,YV0,ZV0
C
C	Read: X,Y,Z at the initial point (Fixed S.)
C

	READ(1,33) VXI,VYI,VZI
C
C	Read: X,Y,Z at the final point (Fixed S.)
C

	READ(1,33) VXF,VYF,VZF
c
c	 If the final point is in orbit around the Earth (circular)  
c
	VCIRC=dsqrt((1.0d0-u)/(DABS(xff+u)))-(-xff-u)
c
        write(*,*) "VCIRC ao redor do Sol", VCIRC !!!!!! Escrevi !ao redor do Sol, mas no ponto final que e a orbita de Marte

!*******************************************************************
! Tornando os valores de velocidades de entrada, dados em km/s, adimensionais:

      XV0 = XV0*PerSJ/((2.0D0*Pi)*d_SJm)
      YV0 = YV0*PerSJ/((2.0D0*Pi)*d_SJm)
      ZV0 = ZV0*PerSJ/((2.0D0*Pi)*d_SJm)
  
      VXI = VXI*PerSJ/((2.0D0*Pi)*d_SJm) 
      VYI = VYI*PerSJ/((2.0D0*Pi)*d_SJm) 
      VZI = VZI*PerSJ/((2.0D0*Pi)*d_SJm)
  
      VXF = VXF*PerSJ/((2.0D0*Pi)*d_SJm) 
      VYF = VYF*PerSJ/((2.0D0*Pi)*d_SJm) 
      VZF = VZF*PerSJ/((2.0D0*Pi)*d_SJm)

1001  FORMAT(F12.6)
1002  FORMAT(2F12.6)
1003  FORMAT(3F12.6)

      write(4,*) "Valores adimensionais só para conferir"
      write(4,*) "mu ="
      write(4,1001) mu
      write(4,*) "(X0,Y0,Z0) ="
      write(4,1003) - d_STm/d_SJm, 0, 0
      write(4,*) "(XFF,YF,ZF) ="
      write(4,1003)  d_SMm/d_SJm, 0, 0
      write(4,*) "(XV0,YV0,ZV0) ="
      write(4,1003) XV0, YV0, ZV0
      write(4,*) "(VXI,VYI,VZI) ="
      write(4,1003) VXI, VYI, VZI
      write(4,*) "(VXF,VYF,VZF) ="
      write(4,1003) VXF, VYF, VZF
      write(4,*) "(TF1,TF2) ="
      write(4,1002) TF1, TF2
 
!*******************************************************************
        i = 1
        j = 1
	PI=4.0D0*DATAN(1.0D0)
	DO 69 TFAU=TF1,TF2,DTFF
	IF (TFAU.NE.TF1)THEN
	  XV0=V(1)
	  YV0=V(2)
	  ZV0=V(3)
	ENDIF
	TF=TFAU
	T0=T00
	XRI(1)=X0
	XRI(2)=Y0
	XRI(3)=Z0
	XRI(4)=VXI
	XRI(5)=VYI
	XRI(6)=VZI
	IF (Z0.NE.0.0D0) THEN
	ANGH1=DATAN(Y0/Z0)
		ANGH1=ANGH1*180.0D0/PI
	  IF (Y0.LT.0.0D0) ANGH1=ANGH1+180.0D0
	  IF (ANGH1.LT.0.0D0) ANGH1=ANGH1+360.0D0
	ENDIF
	IF (Z0.EQ.0.0D0) THEN
		ANGH1=90.0d0
	  IF (Y0.LT.0.0D0) ANGH1=ANGH1+180.0D0
	ENDIF

	CALL ROTOFI(XRI,T0,XFI)
	XRFF(1)=XFF
 	XRFF(2)=YF
	XRFF(3)=ZF
	XRFF(4)=VXF
	XRFF(5)=VYF
	XRFF(6)=VZF
		IF (ZF.NE.0.0D0) THEN
 	ANGH2=DATAN(YF/ZF)
		ANGH2=ANGH2*180.0D0/PI
	  IF (YF.LT.0.0D0) ANGH2=ANGH2+180.0D0
	  IF (ANGH2.LT.0.0D0) ANGH2=ANGH2+360.0D0
	ENDIF
	IF (ZF.EQ.0.0D0) THEN
		ANGH2=90.0d0
	  IF (Y0.LT.0.0D0) ANGH2=ANGH2+180.0D0
	ENDIF

	U1=U	
	XR(1)=X0
	XR(2)=Y0
	XR(3)=Z0
	XR(4)=XV0
	XR(5)=YV0
	XR(6)=ZV0
	CALL ROTOFI(XR,T0,XF)
	V(1)=XV0
	V(2)=YV0
	V(3)=ZV0
1       FORMAT(F9.6,A1,F9.6,A1,F9.6)
9       FORMAT(I6,A1,F12.6,A1,F12.6,A1,F12.6)
	N2=3
	X1=T0
	X2=TF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!! Diminui H1 de 0.01 para 0.0001 e os resultados melhoraram!!
	!H1=0.01D0
        H1=0.000001D0 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	HMIN=0.00000000d0
!!! Diminui DELV(1,2,3) como teste, como ele altera o chute de velocidade inicial antes de integrar o sistema em SHOOT, mas os resultados permaneceram praticamente inalterados.
	DELV(1)=0.00001D0  
	DELV(2)=0.00001D0
	DELV(3)=0.00001D0
	NVAR=6
66	CONTINUE
        CALL SHOOT(NVAR,V,DELV,N2,X1,X2,EPS,H1,HMIN,F,DV)
	TDF=DSQRT(F(1)*F(1)+F(2)*F(2)+F(3)*F(3))
	IF (TDF-TOLF) 11,11,12
12	CONTINUE 
10	GOTO 66
11	CONTINUE
	DTF=(TF-T0)/NP
        !*******************************************************************
        deltaTT = (TF-T0)
        !*******************************************************************
	write(*,*) 'PASSOU',TFAU
        !*******************************************************************
        DelTT(i) = deltaTT        
        !*******************************************************************

	TF=DTF
	XR(1)=X0
 	XR(2)=Y0
	XR(3)=Z0
	XR(4)=V(1)
	XR(5)=V(2)
	XR(6)=V(3)
        !XR(4), XR(5), XR(6) = V(1),V(2) e V(3) ajustados por SHOOT, aqui.
	DELV1R=DSQRT((V(1)-VXI)**2+(V(2)-VYI)**2+(V(3)-VZI)**2)
        !O próximo valor de V(1), V(2) e V(3), da próxima iteração, será igual ao valor apresentado aqui acima, que é o novo chute de velocidade inicial.
	 CALL ROTOFI(XR,T0,XF)
	CALL ENEINERTIAL(XR,XF,U,EINE1)
	DELV1F=DSQRT((XF(4)-XFI(4))**2+(XF(5)-XFI(5))**2+(XF(6)-
     *XFI(6))**2)
        WRITE(2,1) XF(1),CHAR(9),XF(2),CHAR(9),XF(3)
        Xfixo(i,1) = XF(1)
        Xfixo(i,2) = XF(2)
        Xfixo(i,3) = XF(3)
        WRITE(10,9) i,CHAR(9),Xfixo(i,1),CHAR(9),Xfixo(i,2),CHAR(9),Xfixo(i,3)
        WRITE(3,1) XR(1),CHAR(9),XR(2),CHAR(9),XR(3)
        Xrota(i,1) = XR(1)
        Xrota(i,2) = XR(2)
        Xrota(i,3) = XR(3)
        WRITE(11,9) i,CHAR(9),Xrota(i,1),CHAR(9),Xrota(i,2),CHAR(9),Xrota(i,3)
       DO 18 II=1,NP
        CALL ODEINT(XR,6,T0,TF,EPS,h1,hmin,NOK,NBAD,DERIVS,RKQC)
	CALL JACOBI(XR,U,C1)
	WRITE(5,*)TF,C1
        T0=TF
        TF=TF+DTF
        CALL ROTOFI(XR,T0,XF)
        WRITE(2,1) XF(1),CHAR(9),XF(2),CHAR(9),XF(3)
        Xfixo(j,1) = XF(1)
        Xfixo(j,2) = XF(2)
        Xfixo(j,3) = XF(3)
        WRITE(10,9) j,CHAR(9),Xfixo(j,1),CHAR(9),Xfixo(j,2),CHAR(9),Xfixo(j,3)
        WRITE(3,1) XR(1),CHAR(9),XR(2),CHAR(9),XR(3)
        Xrota(j,1) = XR(1)
        Xrota(j,2) = XR(2)
        Xrota(j,3) = XR(3)
        WRITE(11,9) j,CHAR(9),Xrota(j,1),CHAR(9),Xrota(j,2),CHAR(9),Xrota(j,3)
18    CONTINUE
		CALL ROTOFI(XRFF,T0,XFFF)
		CALL ENEINERTIAL(XR,XF,U,EINE2)
	EINETOT=(EINE1+EINE2)/2.0d0
      DELV2R=DSQRT((XR(4)-VXF)**2+(XR(5)-VYF)**2+(XR(6)-VZF)**2)
	DELV2F=DSQRT((XF(4)-XFFF(4))**2+(XF(5)-XFFF(5))**2+(XF(6)-
     *XFFF(6))**2)
C
C	 Test for the sense of the orbits
C
	XRFF(5)=-VCIRC  !! Vcirc é a velocidade circular no ponto final ao redor do Sol, mas que é em uma órbita próxima de Marte, da forma como considerei  (eu escrevi) !!!!  
	CALL ROTOFI(XRFF,T0,XFFF)
	DELVSE=DSQRT((XF(4)-XFFF(4))**2+(XF(5)-XFFF(5))**2+(XF(6)-
     *XFFF(6))**2)


	  XHI=XR(1)
	 YHI=XR(2)
	VXH=XR(4)
	VYH=XR(5)
	VEHIR=XR(4)**2+XR(5)**2+XR(6)**2+2.0d0*(xhi*vyh-vxh*yhi)+xhi*xhi
     *+yhi*yhi+2.0d0*U*(xhi+vyh)+U*U
	 write(*,*) "vehir", vehir
	 DVEHIR=DSQRT(VEHIR)

	DELVTOTR=DELV1R+DELV2R
	DELVTOTF=DELV1F+DELV2F
        DelvTF(i) = DELVTOTF
        DelvTR(i) = DELVTOTR
100      FORMAT(4F14.9)
101      FORMAT(10F14.9)
102	 FORMAT(4F14.9,I6)
103	FORMAT(2F14.9,F4.1)
	FP=DATAN(V(2)/V(1))
	FP=FP*180.0D0/PI
	  IF (V(1).LT.0.0D0) FP=FP+180.0D0
	  IF (FP.LT.0.0D0) FP=FP+360.0D0
      WRITE(7,101)T0,FP,C1,EINETOT,ANGH1,ANGH2,DELVTOTF,
     *DELVTOTR,DELVSE
       WRITE(6,*)'Fixed System'
        WRITE(6,100) XF
        WRITE(6,*)'Rotating System'
        WRITE(6,100) XR
        WRITE(6,*)'Jacobian Constant, Mu, DT, EPS, NP'
        WRITE(6,102) C1,U,DTF,EPS,NP
	WRITE(6,*)'Energy, Initial flight path angle'
	WRITE(6,100)CB,FP
	WRITE(6,*)'Filter, Tolerance for f, System for constraints'
	WRITE(6,103) FCC,TOLF
        WRITE(6,*)
        WRITE(6,*)'Fixed System'
        WRITE(6,100) XF
        WRITE(6,*)'Rotating System'
        WRITE(6,100) XR
        WRITE(6,*)'Vector F (error)'
        WRITE(6,100) F
	WRITE(6,*)"DELTA-vs"
        WRITE(6,100)DELV1F,DELV2F,DELVTOTF
!*******************************************************************
! Escrevendo dados finais a serem plotados em um arquivo a parte:
       ! Passando para unidades dimensionais de anos e km/s 
       DelTT(i) = (DelTT(i)*PerSJanos)/(2.0D0*Pi)
       DelvTF(i) = DelvTF(i)*(2.0D0*Pi*d_SJm)/PerSJ
       DelvTR(i) = DelvTR(i)*(2.0D0*Pi*d_SJm)/PerSJ
       write(8,9) i,CHAR(9),DelTT(i),CHAR(9),DelvTF(i),CHAR(9),DelvTR(i)
       write(15,1002) deltaTT, DELVTOTR 
       norbits = i
       i = i + 1
       j = j + 1
!*******************************************************************       
69	CONTINUE        
         
       STOP
       END

C
C
	SUBROUTINE LOAD(X1,V,Y)
	IMPLICIT REAL*8(A-H,O-Z)
	COMMON /START/ X0,Y0,Z0
	DIMENSION V(3),Y(6)
	Y(1)=X0
	Y(2)=Y0
	Y(3)=Z0
	Y(4)=V(1)
	Y(5)=V(2)
	Y(6)=V(3)
	RETURN
	END
C
C
	SUBROUTINE SCORE(X2,Y,F)
	IMPLICIT REAL*8(A-H,O-Z)
	COMMON /FIM/ XFF,YF,ZF
	DIMENSION F(3),Y(6)
	F(1)=-XFF+Y(1)
	F(2)=-YF+Y(2)
	F(3)=-ZF+Y(3)
	RETURN
	END
C
C
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
        implicit real*8(a-h,o-z)
      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
C
C
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
        implicit real*8(a-h,o-z)
      PARAMETER (NMAX=100,TINY=1.0e-20)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      D=1.d0
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1.d0/AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1.d0/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END
C
C
      SUBROUTINE SHOOT(NVAR,V,DELV,N2,X1,X2,EPS,H1,HMIN,F,DV)
        implicit real*8(a-h,o-z)
	common /filter/ fcc
	COMMON /SOL/ YI
      PARAMETER (NP=6)
      DIMENSION V(N2),DELV(N2),F(N2),DV(N2),Y(NP),DFDV(NP,NP),INDX(NP)
     *,YI(6)
      EXTERNAL DERIVS,RKQC
      CALL LOAD(X1,V,Y)
	DO 15 II=1,6
	YI(II)=Y(II)
15	CONTINUE
      CALL ODEINT(Y,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS,RKQC)
      CALL SCORE(X2,Y,F)
      DO 12 IV=1,N2
        SAV=V(IV)
        V(IV)=V(IV)+DELV(IV)
        CALL LOAD(X1,V,Y)
        CALL ODEINT(Y,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS,RKQC)
        CALL SCORE(X2,Y,DV)
        DO 11 I=1,N2
          DFDV(I,IV)=(DV(I)-F(I))/DELV(IV)
11      CONTINUE
        V(IV)=SAV
12    CONTINUE
      DO 13 IV=1,N2
        DV(IV)=-F(IV)
13    CONTINUE
      CALL LUDCMP(DFDV,N2,NP,INDX,DET)
      CALL LUBKSB(DFDV,N2,NP,INDX,DV)
      DO 14 IV=1,N2
        V(IV)=V(IV)+DV(IV)*fcc
14    CONTINUE
      RETURN
      END
C
C
      SUBROUTINE ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS,RK
     *QC)
        IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (MAXSTP=10000,NMAX=6,TWO=2.D0,ZERO=0.0,TINY=1.D-30)
      COMMON /PATH/ KMAX,KOUNT,DXSAV,XP(200),YP(10,200)
        EXTERNAL DERIVS,RKQC
      DIMENSION YSTART(NVAR),YSCAL(NMAX),Y(NMAX),DYDX(NMAX)
      X=X1
      H=SIGN(H1,X2-X1)
      NOK=0
      NBAD=0
      KOUNT=0
      DO 11 I=1,NVAR
        Y(I)=YSTART(I)
11    CONTINUE     
      XSAV=X-DXSAV*TWO
      DO 16 NSTP=1,MAXSTP
        CALL DERIVS(X,Y,DYDX)
        DO 12 I=1,NVAR
          YSCAL(I)=ABS(Y(I))+ABS(H*DYDX(I))+TINY
12      CONTINUE
        IF(KMAX.GT.0)THEN
          IF(ABS(X-XSAV).GT.ABS(DXSAV)) THEN
            IF(KOUNT.LT.KMAX-1)THEN
              KOUNT=KOUNT+1
              XP(KOUNT)=X
              DO 13 I=1,NVAR
                YP(I,KOUNT)=Y(I)
13            CONTINUE
              XSAV=X
            ENDIF
          ENDIF
        ENDIF
        IF((X+H-X2)*(X+H-X1).GT.ZERO) H=X2-X
        CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        IF(HDID.EQ.H)THEN
          NOK=NOK+1
        ELSE
          NBAD=NBAD+1
        ENDIF
        IF((X-X2)*(X2-X1).GE.ZERO)THEN
          DO 14 I=1,NVAR
            YSTART(I)=Y(I)
14        CONTINUE
          IF(KMAX.NE.0)THEN
            KOUNT=KOUNT+1
            XP(KOUNT)=X
            DO 15 I=1,NVAR
              YP(I,KOUNT)=Y(I)
15          CONTINUE
          ENDIF
          RETURN
        ENDIF
        IF(ABS(HNEXT).LT.HMIN) THEN 
        WRITE(*,*) 'Stepsize smaller than minimum.'
	RETURN
	ENDIF
        H=HNEXT
16    CONTINUE
C      PAUSE 'Too many steps.'
	write(*,*) 'too many steps'
      RETURN
      END
C
C
      SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
        IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=6,FCOR=.0666666667D0,
     *    ONE=1.D0,SAFETY=0.9D0,ERRCON=6.D-4)
      EXTERNAL DERIVS
      DIMENSION Y(N),DYDX(N),YSCAL(N),YTEMP(NMAX),YSAV(NMAX),DYSAV(NMAX)
      PGROW=-0.2D0
      PSHRNK=-0.25D0
      XSAV=X
      DO 11 I=1,N
        YSAV(I)=Y(I)
        DYSAV(I)=DYDX(I)
11    CONTINUE
      H=HTRY
1     HH=0.5D0*H
      CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X=XSAV+HH
      CALL DERIVS(X,YTEMP,DYDX)
      CALL RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X=XSAV+H
      IF(X.EQ.XSAV)PAUSE 'Stepsize not significant in RKQC.'
      CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX=0.
      DO 12 I=1,N
        YTEMP(I)=Y(I)-YTEMP(I)
        ERRMAX=MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
12    CONTINUE
      ERRMAX=ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
        H=SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID=H
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT=SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT=4.D0*H
        ENDIF
      ENDIF
      DO 13 I=1,N
        Y(I)=Y(I)+YTEMP(I)*FCOR
13    CONTINUE
      RETURN
      END
C
C
        SUBROUTINE ROTOFI(XR,T,XF)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION XF(6),XR(6)
	CT=DCOS(T)
	ST=DSIN(T)
        XF(1)=XR(1)*CT-XR(2)*ST
        XF(2)=XR(2)*CT+XR(1)*ST
	  XF(3)=XR(3)
        XF(4)=-XR(1)*ST-XR(5)*ST+XR(4)*CT-XR(2)*CT
        XF(5)=XR(4)*ST-XR(2)*ST+XR(1)*CT+XR(5)*CT
        XF(6)=XR(6)
	  RETURN
        END
C
C
      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
        IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=6)
      DIMENSION Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
      HH=H*0.5D0
      H6=H/6.D0
      XH=X+HH
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
14    CONTINUE
      RETURN
      END
C
	SUBROUTINE JACOBI(XR,U,C1)

C	 This routine uses Broucke's system, like in
C	 Traveling Between the Lagrangian Points and the Moon
C
         IMPLICIT REAL*8(A-H,O-Z)
	 DIMENSION XR(6)
	 R2=DSQRT((-1.0D0+U+XR(1))**2+XR(2)*XR(2)+XR(3)**2)
	 R1=DSQRT((U+XR(1))**2+XR(2)*XR(2)+XR(3)**2)
	 C=XR(5)*XR(5)+XR(4)*XR(4)+XR(6)*XR(6)
	 OME=(XR(1)**2+XR(2)**2)/2.D0+(1.D0-U)/R1+U/R2
	 C1=C/2-OME
         RETURN
         END
C
C		This routine calculates the energy in the inertial frame
C
		SUBROUTINE ENEINERTIAL(XR,XF,U,EINE)

C	 This routine uses Broucke's system, like in
C	 Traveling Between the Lagrangian Points and the Moon
C
         IMPLICIT REAL*8(A-H,O-Z)
    	 DIMENSION XR(6),XF(6)
	  R2=DSQRT((-1.0D0+U+XR(1))**2+XR(2)*XR(2)+XR(3)**2)
	  R1=DSQRT((U+XR(1))**2+XR(2)*XR(2)+XR(3)**2)
	  VEL2=XF(5)*XF(5)+XF(4)*XF(4)+XF(6)*XF(6)
	  OME=(1.D0-U)/R1+U/R2
	  EINE=VEL2/2.0D0-OME
         RETURN
         END

c
c
	SUBROUTINE FITORO(XF,T,XR)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION XF(6),XR(6)
	XR(1)=XF(1)*DCOS(T)+XF(2)*DSIN(T)
	XR(2)=XF(2)*DCOS(T)-XF(1)*DSIN(T)
	XR(3)=XF(3)
	XR(5)=-XF(1)*DCOS(T)+XF(5)*DCOS(T)-XF(4)*DSIN(T)-XF(2)*DSIN(T)
	XR(4)=XF(4)*DCOS(T)+XF(2)*DCOS(T)-XF(1)*DSIN(T)+XF(5)*DSIN(T)
	XR(6)=XF(6)
	RETURN
	END
C           
       SUBROUTINE DERIVS (T,X,DX)
	 REAL*8 R1,R2,X(6),DX(6),T,MI,MISTAR,PR1,PR2
C       
       COMMON /JAC/ MI
	 MISTAR=1.0d0-MI
      
       R1=(X(1)+MI)**2+X(2)**2+X(3)**2
       PR1=R1**1.5
       R2=(X(1)-MISTAR)**2+X(2)**2+X(3)**2
       PR2=R2**1.5
       DX(1)=X(4)
       DX(2)=X(5)
       DX(3)=X(6)
       DX(4)=2*X(5)+X(1)-(MISTAR*(X(1)+MI))/PR1-(MI*(X(1)-MISTAR))/PR2
       DX(5)=-2*X(4)+X(2)-(MISTAR*X(2))/PR1-(MI*X(2))/PR2
       DX(6)=-(MISTAR*X(3))/PR1-(MI*X(3))/PR2
       RETURN 
       END 
C
  
