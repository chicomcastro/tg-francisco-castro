c=======================================================================
      subroutine Transf_Fly_by_Venus (xin,phi_sj,tin,vel_IN,z_venus,x3,
     &                                            theta_venus,kplaneta)    !  Do Ref Baricentrico do SJ Girante (PR3C) PARA Ref Heliocentrico Inercial
c=======================================================================
      implicit real*8 (A-H,O-Z)
      dimension xin(4),x1(4),x2(4),x3(4),x4(4),x5(4)
      dimension vel_IN(3),z_venus(4)
      parameter (PI=acos(-1.0d0),TPI=2.0d0*PI,HPI=0.5d0*PI)
c-----INPUTS------------------------------------------------------------
c     xin: pos and vel in Baricentric Rotating (unidades Sol-Jupiter), tsj: time in SJ dimensionless time units
c     Transforma para Inercial Heliocentrico (mantem unidades Sol-Jupiter)
c-----Constantes--------------------------------------------------------
      amu_sol=1.32712d11 ! km^3 s^{-2} mu_do_Sol =(G*Msol)  (2BP)
      dsj=778.54720d06 ! km fator de normalização de distâncias no Sol-Júpiter 
      wsj=2.0d0*PI/(4332.589d0*24.0d0*60.0d0*60.d0) ! em segundos, 4332.589d0 days (periodo de jupiter em torno do Sol)
      amu_sj=9.53160991d-04   !(mu do Sol-Jupiter PR3C)
      dP1=-amu_sj   ! Convencao Caltech (Nao BCN)

      rc=149597870.7d0 ! 1 A.U.(km)= Raio medio da orbita da Terra em torno do Sol (circular)
c      r_venus=0.72333199d0*rc  ! Raio medio da orbita de Venus Em km
c      r_mars=1.52366231d0*rc   ! Raio medio da orbita de Marte Em km  Semimajor axis (AU) Mars Mean Orbital Elements (J2000)

      if (kplaneta.eq.1) then;  r_venus=0.72333199d0*rc  ! (Planeta VENUS, em km)
      else if (kplaneta.eq.2) then;  r_venus=1.52366231d0*rc; endif  ! (Planeta PARA MARTE, em km)
               
      raio_venus=6051.8d0/dsj                ! Raio medio do planeta Venus
      Gm_venus=0.32486d06/((dsj**3.0d0)*(wsj**2.0d0))   ! km^3/s^2 normalizando pelo Sol-Jupiter (Venus Fact Sheet Nasa)
c-----------------------------------------------------------------------
      do i=1,4; x1(i)=xin(i); enddo; t=tin
c     INVERSE Rotation around z of SE by alfa
      alfa=t+phi_sj
      x2(1)= x1(1)*cos(alfa)-x1(2)*sin(alfa)
      x2(2)= x1(1)*sin(alfa)+x1(2)*cos(alfa)
      x2(3)= x1(3)*cos(alfa)-x1(1)*sin(alfa)-
     &       x1(4)*sin(alfa)-x1(2)*cos(alfa)
      x2(4)= x1(3)*sin(alfa)+x1(1)*cos(alfa)+
     &       x1(4)*cos(alfa)-x1(2)*sin(alfa)
c    Translation of center of RF from the baricenter of the SJ-system TO the Sun
      x3(1)=x2(1)-dP1*cos(alfa)
      x3(2)=x2(2)-dP1*sin(alfa)
      x3(3)=x2(3)+dP1*sin(alfa)    
      x3(4)=x2(4)-dP1*cos(alfa)
      vin_mod=sqrt(x3(3)**2.0d0+x3(4)**2.0d0)   ! modulo da velocidade da espaconave wrt the Sun (Ref Inercial baricentrico)
c-----
      xvenus=x3(1); yvenus=x3(2)    ! Posicao de Venus = Posicao da Espaconave
      call angle(xvenus,yvenus,theta_venus)
      theta_venus_deg=theta_venus*180.0d0/PI
      v_venus=dsqrt(amu_sol/r_venus)/(dsj*wsj)  ! (km/s) depois (escalonado no Sistema Sol-Jupiter)
      z_venus(1)=x3(1)
      z_venus(2)=x3(2)
      z_venus(3)=-v_venus*sin(theta_venus)  ! vx_venus
      z_venus(4)= v_venus*cos(theta_venus)  ! vy_venus
c-----------------------------------------------------------------------
c     Transladando origem do baricentro para VENUS: afeta só velocidade para abordagem pontual de Swing-by        
      vel_IN(1)=x3(3)-z_venus(3)  ! vx_infinito_IN
      vel_IN(2)=x3(4)-z_venus(4)  ! vy_infinito_IN
      vel_IN(3)=sqrt(vel_IN(1)**2.0d0+vel_IN(2)**2.0d0)
c-----
      RETURN; end subroutine Transf_Fly_by_Venus
c-----
c     OBS: Para continuacao do calculo do SB, posicoes em x3 e velocidades em vel_IN
c     posicoes (em x3) no ref inercial centrado no Sol, velocidades no inercial centrado em Venus
c     z_venus vetor do planeta venus (pos e vel)
c=======================================================================
      subroutine Rot_Swingby (rp,ksinal,alfa,vel_IN,z_venus,x3,xin,xout,
     &                                         bless,ang_deg,delta_deg)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      dimension vel_IN(3),z_venus(4),x3(4),x5(4),xout(4),xin(4),bless(7)
      parameter (PI=acos(-1.0d0),TPI=2.0d0*PI,HPI=0.5d0*PI)
c-----------------------------------------------------------------------
      dsj=778.54720d06 ! km fator de normalização de distâncias no Sol-Júpiter 
      wsj=2.0d0*PI/(4332.589d0*24.0d0*60.0d0*60.d0) ! em segundos, 4332.589d0 days (periodo de jupiter em torno do Sol)
      amu_sj=9.53160991d-04; dP1=-amu_sj    !(mu do Sol-Jupiter PR3C) Convencao Caltech
      Gm_venus=0.32486d06/((dsj**3.0d0)*(wsj**2.0d0))   ! km^3/s^2 normalizando pelo Sol-Jupiter (Venus Fact Sheet Nasa)
c-----
      rc=149597870.7d0 ! km  Raio medio da orbita da Terra em torno do Sol (circular)   1.A.U
      r_venus=0.72333199d0*rc  ! Raio medio da orbita de Venus Em km
      raio_venus=6051.8d0/dsj  ! Raio medio do planeta Venus
c----------------------------
      v_mod=vel_IN(3)
      fator=rp*(v_mod**2.0d0)/Gm_venus
      sen_delta=1.0d0/(1.0d0+fator)
      delta=asin(sen_delta)
      delta_deg=delta*180.0d0/PI
c      write(*,*)'******************************************' 
c      write(*,*)'Altitude_rp(km/s),fator,sen_delta,delta_deg,v_mod(Sj)='
c      write(*,*)((rp-raio_venus)*dsj),fator,sen_delta,delta_deg,v_mod; 
c      write(*,*)
c      write(68,*)rp,fator,sen_delta,delta_deg
c      read(*,*) BOBO
c--------
      beta=delta*2.0d0    ! Rotacao de FLY-BY , para beta positivo, sentido anti-horario
      if (ksinal.eq.2) beta=-beta
      vx_OUT= vel_IN(1)*cos(beta)+vel_IN(2)*sin(beta)
      vy_OUT=-vel_IN(1)*sin(beta)+vel_IN(2)*cos(beta)
      v_rela_out=sqrt(vx_OUT**2.0d0+vy_OUT**2.0d0)
c      write(*,*)'INFI OUT vx , vy, v_mod=',vx_OUT, vy_OUT,v_rela_out
c-------- 
c     Transladando origem de VENUS para o SOL  (SB afeta só velocidade para abordagem pontual de Swing-by) 
      x5(3)=vx_OUT+z_venus(3)  ! VETOR final no Ref Inercial Heliocentrico
      x5(4)=vy_OUT+z_venus(4) 
      v_out_mod=sqrt(x5(3)**2.0d0+x5(4)**2.0d0)   ! modulo da velocidade da espaconave wrt the Sun NA SAIDA
c      write(*,*)'v_mod,v_rela_out,abs(v_mod-v_rela_out)=' 
c      write(*,*) v_mod,v_rela_out,abs(v_mod-v_rela_out)  
c--------
c     Calculando variacao de velocidade no Swing-by no Ref Inercial Centrado no Sol  
      del_vx=x5(3)-x3(3); del_vy=x5(4)-x3(4)
      delv_sol=sqrt(del_vx**2.0d0+del_vy**2.0d0)
      delv_teorico=2.0d0*v_mod*sen_delta
      dif=delv_teorico-delv_sol
      fv=(dsj*wsj)
c      write(*,*)'Analise das velocidades: vinf, Delta v max,delta v efe'
c      write(*,*)(v_mod*fv),(2.0d0*v_mod*fv),delv_sol*fv
c      write(*,*)'ec = deltav efe/(2 * vinf)'
      ec=(2.0d0*v_mod)/delv_sol
c      write(*,*)'ecentricidade=',ec,(1.0d0/sen_delta)
c      read(*,*) BOBO

      bless(6)=delv_sol*(dsj*wsj)
c      write(*,*)'bless(6)=',bless(6); read(*,*) bobooooo   

      prod_esc=x5(3)*x3(3)+x5(4)*x3(4)
      v3m=sqrt(x3(3)*x3(3)+x3(4)*x3(4))
      v5m=sqrt(x5(3)*x5(3)+x5(4)*x5(4))
      angulo=acos(prod_esc/(v3m*v5m))
      if (prod_esc.lt.0.0d0) angulo=-angulo
      ang_deg=angulo*180.0d0/PI          ! variacao de direcao do vetor velocidade

c      write(90,*)rp,delta_deg,ksinal,del_mod,dif,ang_deg

      bcos=(x5(3)-x3(3))/(-2.0d0*dsin(delta)*v_mod)
      bsin=(x5(4)-x3(4))/(-2.0d0*dsin(delta)*v_mod)
      call angle_v2 (bcos,bsin,psi); psi_deg=psi*180.0d0/PI
      ddd=(PI-delta*2.0d0)*180.0d0/PI
c      write(91,*)delta_deg,psi_deg,ddd,ang_deg
c      write(*,*)'delta_deg,(psi*180.0d0/PI)=',delta_deg,(psi*180.0d0/PI)
c--------
c     Transladando origem de VENUS para o SOL e logo para o baricentro (afeta só velocidade para abordagem pontual de Swing-by) 
      x5(3)=vx_OUT+z_venus(3)-dP1*sin(alfa)  ! VETOR DE SAIDA
      x5(4)=vy_OUT+z_venus(4)+dP1*cos(alfa) 
      v_out_mod=sqrt(x5(3)**2.0d0+x5(4)**2.0d0)   ! modulo da velocidade da espaconave wrt the Sun NA SAIDA
c     Voltando ao Referencial Girante do PR3C   
      cs=cos(alfa); sn=sin(alfa)
      xout(1)=xin(1); xout(2)=xin(2)
      xout(3)= x5(3)*cs-xin(1)*sn+x5(4)*sn+xin(2)*cs
      xout(4)=-x5(3)*sn-xin(1)*cs+x5(4)*cs-xin(2)*sn
c-----
c      write(*,*)'Voltando ao Refgirante: vx,vy,vmod=',
c     &      xout(3),xout(4),sqrt(xout(3)**2.0d0+xout(4)**2.0d0)
c      read(*,*)bobo 
c-----
      return; end subroutine Rot_Swingby 
c=======================================================================
      subroutine angle (x,y,theta)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      pi=acos(-1.0d0); hpi=0.5d0*pi; tpi=2.0d0*pi
      zero=0.0d0;
      r=sqrt(x*x+y*y)
      xs=dasin(y/r); xc=dacos(x/r)
      if(xs.ge.zero)then
         if(xc.le.hpi)then; theta=xc    ! 1º quad
         else; theta=xc; endif          ! 2º quad
      else if(xs.lt.zero)then
         if(xc.ge.hpi)then; theta=pi-xs ! 3º quad
         else; theta=tpi+xs; endif      ! 4º quad
      endif
      return; end subroutine angle  
c=======================================================================
      subroutine angle_v2 (xc,xs,theta)
c=======================================================================
      implicit real*8 (A-H,O-Z)
      pi=acos(-1.0d0); hpi=0.5d0*pi; tpi=2.0d0*pi
      zero=0.0d0;
      if(xs.ge.zero)then
         if(xc.le.hpi)then; theta=xc    ! 1º quad
         else; theta=xc; endif          ! 2º quad
      else if(xs.lt.zero)then
         if(xc.ge.hpi)then; theta=pi-xs ! 3º quad
         else; theta=tpi+xs; endif      ! 4º quad
      endif
      return; end subroutine angle_v2  
c=======================================================================
c     G=6.674184d-20   ! kg⁻1 km^3 s^-2
c     amassa_venus=4.8675d24    !  Massa de Venus em kg
c     Gmm=G*amassa_venus/((dsj**3.0d0)*(wsj**2.0d0))
c     write(*,*)Gmm,Gm_venus,abs(Gm_venus-Gmm) ! dif=5.9062 (ou 4.44d-11 no SJ)
c     r_v=108.21d06; ddd=abs(r_venus-r_v); write(*,*)r_venus,r_v,ddd
