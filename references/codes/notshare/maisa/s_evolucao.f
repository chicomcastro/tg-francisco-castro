c=======================================================================
      subroutine evolucao4D (yin,rmu,C0,tf,kopcao,t,kfase,bless,phi_sj,
     &                            theta_0_earth,nble,ndv,ncsi,nrp,nphi)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      dimension yin(5),y(4),yold(4),bless(7),yvenus(4),ynew(4),xfb(4)
      dimension vel_IN(3),z_venus(4),ymarte(4)
      dimension v_fake(3),z_marte(4),y_EN_Ine(4)
      real*8 N_v,N_e,N_m
      parameter (PI=acos(-1.0d0),TPI=2.0d0*PI,HPI=0.5d0*PI)

      dsj=778.54720d06 ! km fator de normalização de distâncias no Sol-Júpiter 
      raio_venus=6051.8d0/dsj                ! Raio medio do planeta Venus 
      wsj=2.0d0*PI/(4332.589d0*24.0d0*60.0d0*60.d0) ! em segundos, 4332.589d0 days (periodo de jupiter em torno do Sol)
      rc=149597870.7d0 ! 1 A.U.(km)= Raio medio da orbita da Terra em torno do Sol (circular)
      rterra=rc/dsj

      Fat_day_to_sec=24.0d0*60.0d0*60.d0

c     N_e, N_m em Unidade de Segundos (e nao em u.a.t do Sol-Jupiter)       
      Pv=224.701d0*Fat_day_to_sec; N_v=TPI/Pv     ! 0.32586474777749747        19.281574180800266  
      Pe=365.250d0*Fat_day_to_sec; N_e=TPI/Pe     ! 0.52969100772017463        11.861982203969884 
      Pm=686.980d0*Fat_day_to_sec; N_m=TPI/Pm     ! 0.99626866114607948        6.3067178083787025   
      PEM_Hohmman=2.2362d7                        ! segundos, Curtis pag 412 Example 8.2
      PEM_days=PEM_Hohmman/Fat_day_to_sec         ! write(*,*)'PEM_days (Hohmann time in days)=',PEM_days
       
      do i=1,4; y(i)=yin(i); yold(i)=y(i); enddo
      t=yin(5); told=t-1.0d0 ! dummy diferente de t
      hnext=0.002d0; dt=0.001d0

      kopcao=0; kfase=1 ! 1: entre Terra e Marte
      nrp=0 
      call testar(t,y,kopcao,C0,C,kfase,nrp)

      do 10 while (kopcao.eq.0); kopcao=0
         call evoluinorm(tf,t,told,y,yold,hnext,dt,kopcao,C0,C,rmu,
     &                                                 kfase,nrp)
 10   continue    
c-----------------------------------------------------------------------
c           Chegada em Marte ou atingiu tf?
c-----------------------------------------------------------------------
      if (kopcao.eq.7) then
         call Bissecao(t,told,y,yold,tmarte,ymarte,2) ! ymarte é o vetor da EN em Marte, não é o vetor de Marte!!
         bless(2)=(tmarte/wsj)/Fat_day_to_sec         ! Tempo Earth-Mars em dias

         call Transf_Fly_by_Venus (ymarte,phi_sj,tmarte,v_fake,
     &                                  z_marte,y_EN_Ine,theta_marte,2)      !pos da EN em y_EN_Ine
         bless(7)=theta_marte*180.0d0/PI                       

         del_Vx=z_marte(3)-y_EN_Ine(3)
         del_Vy=z_marte(4)-y_EN_Ine(4) 
         delV_Marte=sqrt(del_Vx*del_Vx+del_Vy*del_Vy)
         delV_Marte_fis=delV_Marte*(dsj*wsj)
         bless(5)=delV_Marte_fis
c-----------------------------------------------------------------------
         phi_good=-PI+N_e*PEM_Hohmman ! Fase da Terra - Fase de Marte na partida do retorno a Terra
c         write(*,*)(N_e*PEM_Hohmman),(-PI),phi_good
cc         read(*,*)bobo
         phi_good_deg=phi_good*180.0d0/PI
c--------------
         tmarte=tmarte/wsj
         Fase_0_Mars=theta_marte-N_m*tmarte      ! Fase_ONDE_(de quem), i.e., fase inicial de Marte (na partida da Terra)
         Fase_M_Mars=theta_marte                 ! Fase em Marte de Marte (ie, qdo EN em Marte)
         Fase_0_Eart=theta_0_earth               ! fase inicial da Terra (angulo Theta do Prog Principal)
         Fase_M_Eart=theta_0_earth+N_e*tmarte    ! fase em Marte da Terra
         r=180.0d0/PI 
         at=Fase_M_Eart*r
         call Mod_TWOPI (Fase_M_Mars,fmd,fmm_mod)
         call Mod_TWOPI (Fase_M_Eart,fed,fem_mod)
         call Mod_TWOPI (Fase_0_Mars,fzmd,fzm_mod)
         call Mod_TWOPI (Fase_0_Eart,fzed,fzm_mod)

         alfa=Fase_M_Mars-Fase_M_Eart          ! diferença de fase entre (M-E) em Marte
         call Mod_TWOPI (alfa,ald,alfa_mod)
         N=0 
15       continue
         b12=phi_good-alfa-(TPI*float(N))
         twait=b12/(N_m-N_e)
         if (twait.lt.0.0d0) then; N=N+1; goto 15; endif
         twait_days=twait/Fat_day_to_sec
c-----------------------------------------------------------------------
         bless(3)=twait_days ! Tempo de espera em Marte (days)
         stp=bless(2)+bless(3)
         svp=bless(4)+bless(5)
         st=bless(2)+bless(3)+258.8636d0
         sv=bless(4)+bless(5)+5.593447d0

c         if (bless(2).lt.500.0d0.AND.bless(3).lt.600.0d0.AND.
c     &                 bless(5).lt.16.0d0) then
          write(500,*)(bless(jj),jj=2,5),bless(7),fed,nphi,ndv,ncsi,nble

          write(200,*)nphi,fzed,fzmd,fed,fmd,twait_days,ald

c          write(21,*)Fase_0_Mars*r,theta_marte*r,N_m,tmarte
c          write(22,*)b12*r,phi_good*r,alfa*r,b12_mod*r
c         endif
c-----------------------------------------------------------------------
c        A partir de agora EN segue orbita circular de Marte até momento de retorno a Terra com Hohmann                   
      endif 
      return; end subroutine evolucao4D
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c     kopcao = 0 integracao sem regularizacao          R E V I S A R
c              1 integracao com regularizacao
c              2 colisao com o primario menor
c              3 colisao com o primario maior
c              4 cte de Jacobi nao conservada
c              5 tempo final atingido
c              6 alguma secao é criterio de parada
c=======================================================================
      subroutine testar(time,y,kopcao,C0,C,kfase,nrp)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      real*8,parameter :: rmu=9.53160991d-04   ! amu_sj  (mu do Sol-Jupiter PR3C)
      dimension y(4),xsecao(4)
      external sisdim
c-----
      kopcao=0 
      epsJac=1.0d-11
      d_max=5.0d0

      rc=149597870.7d0 ! km  Raio medio da orbita da Terra em torno do Sol (circular)   1.A.U
      dsj=778.54720d06 ! km fator de normalização de distâncias no Sol-Júpiter 

      r_venus=0.72333199d0*rc/dsj  ! Ja no Sistema Sol-Jupiter
      r_mars=1.52366231d0*rc/dsj   ! Raio medio da orbita de Marte Em A.U. (Semimajor axis) Mars Mean Orbital Elements (J2000)
c-----
      d1=sqrt((y(1)+rmu)**2.0d0+y(2)**2.0d0)
      d2=sqrt((y(1)+rmu-1.0d0)**2.0d0+y(2)**2.0d0)
c-----
      C=CteJac(y)
      dif=abs(C-C0); ! write(205,*)time,dif
      if(kopcao.ne.1)then
         if(abs(C-C0).gt.epsJac)then
            write(4,*)'Cte de Jacobi não conservada! t =',time,nci
            write(4,*)'C-C0 =',abs(C-C0); write(4,*)'d2= ',d2
c            write(*,*)'see log! read bobo'; read(*,*)bobo
            kopcao=4; return
         endif
      endif

      if (d1.gt.d_max) then    ! P3 se afastou demais do Sistema Sol_Jupiter
c         write(*,*)'time,d1,d_max=',time,d1,d_max
         kopcao=8; return
      endif
c      write(299,*)time,d1,r_mars,(r_mars-d1)

      if (kfase.eq.0.AND.d1.lt.r_venus) then
c         write(*,*)'Atingiu VENUS  time,d1,r_venus=',time,d1,r_venus
c         write(*,*)'distancia P3 a Venus=', abs(d1-r_venus) 
c         write(60,*)abs(d1-r_venus) ! distancia P3 a Venus
         kopcao=9
         RETURN
      endif

c      write(*,*)'rc,r_venus,r_mars=',(rc/dsj),r_venus,r_mars; STOP 
      if (kfase.eq.1.AND.d1.ge.r_mars) then
c         write(*,*)'Atingiu MARTE  time,d1,r_mars=',time,d1,r_mars
         kopcao=7
      endif
     
 20   format(9e25.10)
 21   format(7f20.15)
 22   format(4f20.15)
      RETURN; end subroutine testar
c-----------------------------------------------------------------------
c     kopcao=0 (evolui normal); 4 (não conservou C)
c            5 (chegou ao tempo final); 8 (afastamento demasiado do Sol)
c            9 (chegou a Venus); 7 (chegou a Marte)
c=======================================================================
      subroutine Mod_TWOPI (teta_in,tetadeg,tetamod)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      parameter (PI=acos(-1.0d0),TPI=2.0d0*PI,HPI=0.5d0*PI)
      kn=0; kp=0; ! write(*,*)'tetaIN=',teta_in
      teta=teta_in
      if (teta.lt.0.0d0) then
         do while (teta.lt.0.0d0)
            kn=kn+1
            teta=teta+TPI
         enddo
      else
         do while (teta.gt.TPI)
            kp=kp+1
            teta=teta-TPI
         enddo
      endif
      tetamod=teta   
      tetadeg=tetamod*180.0d0/PI
c      write(*,*)'kn,kp,tetamod,tetadeg=',kn,kp,tetamod,tetadeg
      RETURN; END subroutine Mod_TWOPI
c=======================================================================
      subroutine sisdim4D (t,y,yprime)
c=======================================================================
      implicit real*8 (a-h,o-z)
      dimension y(4),yprime(4)
c-----
      rmu=9.53160991d-04   ! (mu do Sol-Jupiter PR3C)
c-----
      x1=y(1); x2=y(2); x3=y(3); x4=y(4)
      R1=sqrt((x1+rmu)**2.0d0+x2**2.0d0)
      R2=sqrt((x1-1.0d0+rmu)**2.0d0+x2**2.0d0)
c-----
      Omega_x1=x1-(1.0d0-rmu)*(x1+rmu)/R1**3.0d0-
     &         rmu*(x1-1.0d0+rmu)/R2**3.0d0
      Omega_x2=x2-(1.0d0-rmu)*x2/R1**3.0d0-rmu*x2/R2**3.0d0
c-----
      yprime(1)=x3
      yprime(2)=x4
      yprime(3)=Omega_x1+2.0d0*x4
      yprime(4)=Omega_x2-2.0d0*x3
      RETURN; END subroutine sisdim4D
c=======================================================================
      subroutine evoluinorm (tfim,time,told,y,yold,hnext,dt,
     &                       kopcao,C0,C,rmu,kfase,nrp)
c=======================================================================
      implicit real*8 (a-h,o-z)
      dimension y(4),yold(4),relerr(4),abserr(4),work(56) ! 56=4*14=N*14
c-----
      external sisdim4D
c-----
      do i=1,4; relerr(i)=1.0D-14; abserr(i)=1.0D-15; enddo
      if(abs(dt).le.abs(hnext))then; tout=time+dt; else
      tout=time+hnext; endif; iflag=1
c-----------------------------------------------------------------------
      do 100 while (abs(time).le.abs(tfim).and.kopcao.eq.0)
         do i=1,4; yold(i)=y(i); enddo; told=time
         call rkf78(sisdim4D,4,y,time,tout,
     &              relerr,abserr,iflag,work,dt)  
         if(abs(time-tout).gt.abs(hnext))then
            write(4,*)'Passo de tempo excede hnext.',nci
            write(4,*)'Passo=',abs(time-tout),' time=',time,kopcao
            write(*,*)'see log. read bobo'; read(*,*)bobo
         endif
         if(iflag.ne.2)then
            write(4,*)'read bobo. iflag=',iflag,nci; read(*,*)bobo
         else ! time=tout (isso já é feito)
            if(abs(dt).lt.abs(hnext))then
               tout=tout+dt; else
               tout=tout+hnext
            endif
         endif
         call testar(time,y,kopcao,C0,C,kfase,nrp)
 100  continue
c-----------------------------------------------------------------------
      if(abs(time).gt.abs(tfim))kopcao=5
c-----------------------------------------------------------------------
c     kopcao=0 (evolui normal); 4 (não conservou C)
c            5 (chegou ao tempo final); 8 (afastamento demasiado do Sol)
c-----------------------------------------------------------------------    
      return; end subroutine evoluinorm
c=======================================================================

c     Sub abaixo nao usada

c=======================================================================
      subroutine sisdim_secao (t,y,yprime)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      real*8 x1,x2,x3,x4,t,Omega_x1,Omega_x2,R1,R2,aka
      real*8 y(4),yprime(4)

      u=9.53160991d-04        ! (mu do Sol-Jupiter PR3C, amu_sj)
      xP1=-u; xP2=1.0d0-u     ! Convencao Caltech (Nao BCN)

      x1=y(1); x2=y(2); x3=y(3); x4=y(4)
      R1=dsqrt((x1-xP1)*(x1-xP1)+x2*x2)
      R2=dsqrt((x1-xP2)*(x1-xP2)+x2*x2)

      Omega_x1=x1-(1.0d0-u)*(x1-xP1)/R1**3.0d0-u*(x1-xP2)/R2**3.0d0
      Omega_x2=x2-(1.0d0-u)*x2/R1**3.0d0-u*x2/R2**3.0d0
c-----sistema dinamico--------------------------------------------------
      yprime(1)=x3
      yprime(2)=x4
      yprime(3)=Omega_x1+2.0d0*x4
      yprime(4)=Omega_x2-2.0d0*x3
c-----------------------------------------------------------------------
      rdot=(x1-xP1)*x3+x2*x4/sqrt((x1-xP1)*(x1-xP1)+x2*x2) 
      aka=1.0d0/rdot
      do i=1,4; yprime(i)=yprime(i)*aka; enddo
c-----------------------------------------------------------------------
      return
      end subroutine sisdim_secao
c=======================================================================
