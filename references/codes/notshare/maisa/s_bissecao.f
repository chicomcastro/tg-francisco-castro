c=======================================================================
      subroutine Bissecao (t,told,y,yold,tvenus,yvenus,kplaneta)
c=======================================================================
      implicit real*8 (a-h,o-z)
      real*8 yold(4),y(4),x1(4),x2(4),xmeio(4),yvenus(4)
      real*8,parameter :: rmu=9.53160991d-04   ! amu_sj  (mu do Sol-Jupiter PR3C)
      external sisdim4D
c-----------------------------------------------------------------------
      rc=149597870.7d0 ! km  Raio medio da orbita da Terra em torno do Sol (circular)   1.A.U
      dsj=778.54720d06 ! km fator de normalização de distâncias no Sol-Júpiter
      if (kplaneta.eq.1) then;  r_venus=0.72333199d0*rc/dsj  ! Ja no Sistema Sol-Jupiter  (Planeta VENUS)
      else if (kplaneta.eq.2) then;  r_venus=1.52366231d0*rc/dsj; endif  ! PARA MARTE
      tol=1.0d-10
c-----
      do k=1,4; x1(k)=yold(k); x2(k)=y(k); enddo; t1=told; t2=t

      d1=sqrt((x1(1)+rmu)**2.0d0+x1(2)**2.0d0)
      d2=sqrt((x2(1)+rmu)**2.0d0+x2(2)**2.0d0)
      dif1=d1-r_venus; dif2=d2-r_venus
      if (dif1*dif2.gt.0.0d0) then
         write(*,*)'Problema na entrada da bissecao'; read(*,*)BOBO 
      endif
      nbis=0
10    continue
      passo=(t2-t1)/2.0d0; nbis=nbis+1; do k=1,4; xmeio(k)=x1(k); enddo
      call myRK4 (sisdim4D,4,xmeio,t1,passo); tmeio=t1+passo
      d_meio=sqrt((xmeio(1)+rmu)**2.0d0+xmeio(2)**2.0d0)
      dif_meio=d_meio-r_venus; prod=dif1*dif_meio

      if (abs(dif_meio).lt.tol) goto 20 

      if (prod.le.0.0d0) then; 
         do k=1,4; x2(k)=xmeio(k); enddo; t2=tmeio; dif2=dif_meio
      else
         do k=1,4; x1(k)=xmeio(k); enddo; t1=tmeio; dif1=dif_meio
      endif

      if (nbis.gt.200) then; 
      write(*,*)'Excedeu numero de passos na bis de Venus'
      STOP; endif
        
      if (abs(dif_meio).gt.tol.AND.nbis.lt.200) goto 10
 
20    continue
      tvenus=tmeio; do k=1,4; yvenus(k)=xmeio(k); enddo
cc      write(335,*)(yvenus(k),k=1,4),tvenus,dif_meio
cc      write(*,*)'Exito: n_iteracoes,dif=',nbis,dif_meio; write(*,*)
cc      write(*,*)'BISSECAO    t,told,tvenus=',t,told,tvenus;
c      write(*,*)'told,tvenus,t=',told,tvenus,t; write(*,*)
c      write(*,*)'yold(k) = ',(yold(k),k=1,4)
c      write(*,*)'yvenus(k)=',(yvenus(k),k=1,4)
c      write(*,*)'y(k) =    ',(y(k),k=1,4) 

      RETURN; END subroutine Bissecao
c=======================================================================
      subroutine myRK4 (f,nn,y,t,h)
c=======================================================================
      implicit none
      integer i,nn
      real*8 :: t,h,hh,h6,th
      real*8, dimension (nn) :: dydx,y,dym,dyt,yt

      hh=h*0.5d0
      h6=h/6.0d0
      th=t+hh

      CALL f(t,y,dydx)
      do 11 i=1,nn
         yt(i)=y(i)+hh*dydx(i)
 11   continue
      CALL f(th,yt,dyt)
      do 12 i=1,nn
         yt(i)=y(i)+hh*dyt(i)
 12   continue
      CALL f(th,yt,dym)
      do 13 i=1,nn
         yt(i)=y(i)+h*dym(i)
 13   continue
      do 14 i=1,nn
         dym(i)=dyt(i)+dym(i)
 14   continue
      call f(t+h,yt,dyt)
      do 15 i=1,nn
         y(i)=y(i)+h6*(dydx(i)+dyt(i)+(2.0d0*dym(i)))
 15   continue

      return; end subroutine myRK4
c=======================================================================
