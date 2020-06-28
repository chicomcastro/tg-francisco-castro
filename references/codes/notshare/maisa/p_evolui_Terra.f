      program p_Earth_Venus_no_SunJupiter
c-------------------------
c     Este programa de analise preliminar investiga trajetorias geradas 
c     por CIs de uma orbita circular media da Terra em torno do SOl.
c     (i) Parte-se de CIs no Referencial Heliocêntrico Inercial,
c     (ii) Transforma-se estados para o Ref. Baricêntrico Girante do Sistema Sol-Jupiter.
c     Após acrescimo de velocidade, evolui no tempo pelo PR3C Sol-Jupiter até cruzar 
c     a órbita circular de Venus. Ajustando cruzamento com Venus. 
c 
c     Problema de 3 corpos, usando escolha de posicoes dos primarios da Caltech
c-------------------------      
c      use m_conserva
c---------------------------------------------------
      implicit real*8 (A-H,O-Z)
      real*8 y0(5),ysj(5),yorig(5),y1(5),bless(7)   ! net=5  ! ndim=net
      parameter (PI=acos(-1.0d0),TPI=2.0d0*PI,HPI=0.5d0*PI)
c-----------------------------------------------------------------------
c     Constantes importantes
c-----------------------------------------------------------------------
      amu_sol=1.32712d11 ! km^3 s^{-2} mu_do_Sol =(G*Msol)  (2BP)
      rc=149597870.7d0 ! km    Raio medio da orbita da Terra em torno do Sol (circular)   1.A.U

      dsj=778.54720d06 ! km fator de normalização de distâncias no Sol-Júpiter
      wsj=2.0d0*PI/(4332.589d0*24.0d0*60.0d0*60.d0) ! em segundos, 4332.589d0 days (periodo de jupiter em torno do Sol)
      amu_sj=9.53160991d-04   !(mu do Sol-Jupiter PR3C)

      r_venus=0.72333199d0*rc/dsj  ! Ja no Sistema Sol-Jupiter
      r_mars=1.52366231d0*rc/dsj   ! Ja no Sistema S0l-Jupiter
c      write(*,*)' raio da orb de venus=',r_venus; ! read(*,*)BOBO
c      write(*,*)' raio da orb de terra=',rc; ! read(*,*)BOBO
c      write(*,*)' raio da orb de marte=',r_mars; ! read(*,*)BOBO
c-----------------------------------------------------------------------
c     iniciando o programa 
c--------------

      phi_sj=0.0d0   ! PI/6.0d0    ! FASE INICIAL DO SOL-JUPITER  30 graus

      ti=0.0d0 
      Fat_day_to_sec=24.0d0*60.0d0*60.d0

      tfis=600.0d0; tf=tfis*Fat_day_to_sec*wsj    ! tfis em dias (600 dias=0.8701289 uat)   
      write(*,*)'tfinal (days),(uat) =',tfis,tf   

      vc=dsqrt(amu_sol/rc)  ! (km/s)  ! velocidade da orb circular da Terra em torno do Sol
c      write(*,*)'Razao de raios orb_Terra/orb_Jupiter=',rc/dsj
c=======================================================================
      theta=0.0d0; ntheta=0
      tol=TPI; dfasesj=PI/8.0d0; nphi=0
      do while (phi_sj.lt.tol)
         write(*,*)'phi_sj,theta=',phi_sj,theta
         nphi=nphi+1 

         y0(1)=rc*dcos(theta);  y0(2)=rc*dsin(theta)
         y0(3)=-vc*dsin(theta); y0(4)=vc*dcos(theta)
         y0(5)=ti;                              ! time=0.0d0
c-----------------------------------------------------------------------
c        Incremento no modulo do incremento de velocidade na partida da Terra
c        dv_maximo/n_deltav_max=10km/s / 100 valores=0.1km/s de passo em dv
c-----------------------------------------------------------------------
         n_deltav_max=500; d_dv=0.02d0                   ! Define dv max e passo em dv conforme acima
         n_csi_max=500; dcsi=TPI/(real(n_csi_max))       ! variação angular do incremento de velocidade na partida da Terra
         nble=0                                          ! nble=Numero total de traj analisadas
         ndv=0;  dv=0.0d0   ! dv em km/s 
         do while (ndv.le.n_deltav_max)
            ndv=ndv+1
c              write(*,*)'nble, ndv=',nble, ndv

            ncsi=0; csi=0.0d0        ! INICIALIZACAO DAS VARIAVEIS ANGULAR DA VARIACAO DA VELOCIDADE
            do while (ncsi.lt.n_csi_max)
               ncsi=ncsi+1; nble=nble+1
c=======================================================================
               bless(4)=dv   ! deltav na Terra em km/s
               dvx=dv*cos(csi)
               dvy=dv*sin(csi)
c              write(*,*)' dv fisico (km/s)=',dv
c              write(*,*)'dv calculado=',sqrt(dvx*dvx+dvy*dvy)
c              write(*,*)'dv sol-jup=',dv/(dsj*wsj) 
 
               y1(1)=y0(1);y1(2)=y0(2);y1(3)=y0(3)+dvx;y1(4)=y0(4)+dvy
               y1(5)=y0(5)       ! tempo inicial nulo

               call Transf_RHI_to_RBG (y1,ysj,phi_sj,theta)  ! (Uni.Fisicas para Adimensionais do PR3C SJ)
               Czero=CteJac(ysj)  ! Cte de energia da ci
               C_ini=Czero

               write(31,453)(ysj(k),k=1,4),Czero,ndv,ncsi,nble,nphi      ! outra opção (theta/PI)
c--------------
               call evolucao4D (ysj,amu_sj,Czero,tf,kopcao,tout,
     &                 kfase,bless,phi_sj,theta,nble,ndv,ncsi,nrp,nphi)   ! Tempo inicial no ysj(5)
100            continue
               csi=csi+dcsi
            enddo
            dv=dv+d_dv
         enddo
         write(*,*)'Numero total de traj analisadas:',nble
         phi_sj=phi_sj+dfasesj
      enddo
c-----------------------------------------------------------------------
405   format(5E16.6)      
453   format(5E16.6,2I5,I7,I5)      
c-----------------------------------------------------------------------
      STOP
      END
c=======================================================================
      subroutine Transf_RHI_to_RBG (xin,xout,phi_sj,theta)    ! Ref Heliocentrico Inercial para Baricentrico do SJ Girante (PR3C)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      dimension xin(5),xout(5),x1(5),x2(5),x3(5)   ! tempo na quinta posicao dos vetores
      parameter (PI=acos(-1.0d0),TPI=2.0d0*PI,HPI=0.5d0*PI)
c-----INPUTS------------------------------------------------------------
c     xin: pos and vel in Heliocentric Inertial (unidades fisicas), time in seconds
c-----OUTPUTS-----------------------------------------------------------
c     xout: pos and vel in Baricentric Rotating (unidades Sol-Jupiter), tsj: time in SJ dimensionless time units
c-----------------------------------------------------------------------
      do i=1,5; x1(i)=xin(i); enddo
c-----Constantes
      dsj=778.54720d06 ! km fator de normalização de distâncias no Sol-Júpiter 
      Per_Jup=4332.589d0*24.0d0*60.0d0*60.d0 ! em segundos, 4332.589d0 days (periodo de jupiter em torno do Sol)
      wsj=2.0d0*dacos(-1.0d0)/Per_Jup
      amu_sj=9.53160991d-04   !(mu do Sol-Jupiter PR3C)
      dP1=-amu_sj   ! Convencao Caltech (Nao BCN)
c-----Rescaling
      x1(1)=x1(1)/dsj
      x1(2)=x1(2)/dsj
      x1(3)=x1(3)/(dsj*wsj)
      x1(4)=x1(4)/(dsj*wsj)
      x1(5)=x1(5)*wsj            ! c tse=tfis*we; tfis=tem/wm (do programa velho)
      alfa=x1(5)+phi_sj
c    Translation of center from the Sun to the baricenter of the SJ-system
      x2(1)=x1(1)+dP1*cos(alfa)
      x2(2)=x1(2)+dP1*sin(alfa)
      x2(3)=x1(3)-dP1*sin(alfa)
      x2(4)=x1(4)+dP1*cos(alfa)
      x2(5)=x1(5)
c     Rotation around z of SE by alfa
      x3(1)= x2(1)*cos(alfa)+x2(2)*sin(alfa)
      x3(2)=-x2(1)*sin(alfa)+x2(2)*cos(alfa)
      x3(3)= x2(3)*cos(alfa)-x2(1)*sin(alfa)+
     &       x2(4)*sin(alfa)+x2(2)*cos(alfa)
      x3(4)=-x2(3)*sin(alfa)-x2(1)*cos(alfa)+
     &       x2(4)*cos(alfa)-x2(2)*sin(alfa)
      x3(5)=x2(5)  
c     Escrevendo vetor de estado de saida
      do i=1,5; xout(i)=x3(i); enddo
c-----
      return
405   format(5E16.6)
      end
c=======================================================================
      subroutine Transf_RBG_to_RHI (xin,xout,phi_sj,theta)    !  Do Ref Baricentrico do SJ Girante (PR3C) PARA Ref Heliocentrico Inercial
c=======================================================================
      implicit real*8 (A-H,O-Z)
      dimension xin(5),xout(5),x1(5),x2(5),x3(5),x4(5)   ! tempo na quinta posicao dos vetores
      parameter (PI=acos(-1.0d0),TPI=2.0d0*PI,HPI=0.5d0*PI)
c-----INPUTS------------------------------------------------------------
c     xin: pos and vel in Baricentric Rotating (unidades Sol-Jupiter), tsj: time in SJ dimensionless time units
c-----OUTPUTS-----------------------------------------------------------
c     xout: pos and vel in Heliocentric Inertial (unidades fisicas), time in seconds
c-----------------------------------------------------------------------
      do i=1,5; x1(i)=xin(i); enddo
c-----Constantes
      dsj=778.54720d06 ! km fator de normalização de distâncias no Sol-Júpiter 
      Per_Jup=4332.589d0*24.0d0*60.0d0*60.d0 ! em segundos, 4332.589d0 days (periodo de jupiter em torno do Sol)
      wsj=2.0d0*dacos(-1.0d0)/Per_Jup
      amu_sj=9.53160991d-04   !(mu do Sol-Jupiter PR3C)
      dP1=-amu_sj   ! Convencao Caltech (Nao BCN)
c     INVERSE Rotation around z of SE by alfa
      alfa=x1(5)+phi_sj
      x2(1)= x1(1)*cos(alfa)-x1(2)*sin(alfa)
      x2(2)= x1(1)*sin(alfa)+x1(2)*cos(alfa)
      x2(3)= x1(3)*cos(alfa)-x1(1)*sin(alfa)-
     &       x1(4)*sin(alfa)-x1(2)*cos(alfa)
      x2(4)= x1(3)*sin(alfa)+x1(1)*cos(alfa)+
     &       x1(4)*cos(alfa)-x1(2)*sin(alfa)
      x2(5)=x1(5)
c    Translation of center of RF from the baricenter of the SJ-system TO the Sun
      x3(1)=x2(1)-dP1*cos(alfa)
      x3(2)=x2(2)-dP1*sin(alfa)
      x3(3)=x2(3)+dP1*sin(alfa)
      x3(4)=x2(4)-dP1*cos(alfa)
      x3(5)=x2(5)
c-----Rescaling
      x4(1)=x3(1)*dsj
      x4(2)=x3(2)*dsj
      x4(3)=x3(3)*(dsj*wsj)
      x4(4)=x3(4)*(dsj*wsj)
      x4(5)=x3(5)/wsj            ! c tse=tfis*we; tfis=tem/wm (do programa velho)
c     Escrevendo vetor de estado de saida
      do i=1,5; xout(i)=x4(i); enddo
c-----------------------------------------------------------------------
c      write(42,405)(x4(k),k=1,4),(theta/PI)   ! dis-normalizados
c      write(43,405)(x3(k),k=1,4),(theta/PI)   ! dis-transladados
c      write(44,405)(x2(k),k=1,4),(theta/PI)   ! dis-girados
c-----
      return
405   format(5E16.6)
      end
c=======================================================================
      function CteJac (yin)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      real*8,parameter :: u=9.53160991d-04   ! amu_sj  (mu do Sol-Jupiter PR3C)
      real*8,parameter :: xP1=-u, xP2=1.0d0-u, yP1=0.0d0, yp2=0.0d0    ! CALTECH
      dimension yin(5)
      x1=yin(1); x2=yin(2); x3=yin(3); x4=yin(4) 
      
      R1=dsqrt((x1-xP1)**2.0d0+x2**2.0d0)
      R2=dsqrt((x1-xP2)**2.0d0+x2**2.0d0)
      
      Omega=0.5d0*(x1*x1+x2*x2)+((1.0d0-u)/R1)+
     &     (u/R2)+0.5d0*u*(1.0d0-u)
      CteJac=2.0d0*Omega-(x3*x3+x4*x4)
      return
      end
c=======================================================================
