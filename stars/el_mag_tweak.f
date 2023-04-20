	IMPLICIT NONE
	real pi,q,t,ksi0,f,fs,ksi,ksi1,omega0,R0cr,R0,beta,Te(100,200)
	real omega1,theta,lam(100),phi(100,200),nju(100,200)
	real r(100,200),x(100,200),y(100,200),z(100,200),gsr,Tef
	real r1k(100,200),r2k(100,200),dx(100,200),dy(100,200)
	real dz(100,200),g(100,200),cbeta(100,200),w(100,200),s,sg,angle
	real ph(100),gam_and,climb_and
	real a(100,200),dI_and,sI(100),res(100)
	real omega,modul,casat,limb,pi2,integ
	real tem(20),a1(20),a2(20),a3(20),a4(20)
	real wl(1000),wl5(1000),filtr(1000)
	real angle1, angle2, angle_step, q1, q2, q_step, ang,qq
	integer nn,mm,i,m(200),k,l,band,n,te_output
	common tem,a1,a2,a3,a4
	common /a_common/wl,wl5,filtr
	pi=3.1415926
	pi2=2*pi
	te_output=0
	print*, 'input filling factor R0/R0cr'
	read*, t
	print*, 'input q=Mg/Mh (min) (max) (step)'
	read*, q1, q2, q_step
	print*, 'input i (min) (max) (step)'
	read*, angle1, angle2, angle_step
	print*, 'input Tef'
	read*, Tef
	print*, 'input beta'
	read*, beta
	open(11,file='xyz.dat')
	open(13,file='Te.dat')
	open(14,file='J_mag.res')
	print*, 'input 1 - for J band,2 - for H, 3 - for K'
	read*, band
	open(21,file='about.txt')
	write(21,*)'Filling factor R0/R0cr=',t
	write(21,*)'q1,q2,q_step=',q1,q2,q_step
	write(21,*)'i1,i2,i_step=',angle1,angle2,angle_step
	write(21,*)'Teff=',Tef
	write(21,*)'betta=',beta
	write(21,*)'band=',band
	close(21)
	
	
	if(band.EQ.1) then 
	   open(16,file='T_35_J.dat')
	else if(band.EQ.2) then
	   open(16,file='T_35_H.dat')
	else 
	   open(16,file='T_35_K.dat')
	end if
	do 111 i=1,11
	   read(16,*) tem(i),a1(i),a2(i),a3(i),a4(i)
 111	continue
	close(16)
	
	n=1000
	if(band.EQ.1) then
	   do 910 i=1,n
	      wl(i)=1040.+0.4*i
	      wl5(i)=wl(i)**5.
	      filtr(i)=exp(((wl(i)-1238.)**2.)/1.69e4)
 910	   continue
	else if(band.EQ.2) then
	   do 911 i=1,n
	      wl(i)=1440.+0.4*i
	      wl5(i)=wl(i)**5.
	      filtr(i)=1.0/(-3.2e-5*wl(i)*wl(i)+0.104*wl(i)-84.1)
 911	   continue
	else
	   do 912 i=1,n
	      wl(i)=1900.+0.61*i
	      wl5(i)=wl(i)**5.
	      filtr(i)=1.0/(-1.14e-5*wl(i)*wl(i)+0.0503*wl(i)-54.1)
 912	   continue
	end if
	
	open(22,file='filtr.dat')
	do 588 i=1,1000
	   write(22,*)wl(i),' ',filtr(i)
 588	continue
	close(22)
	
cccc    Вычисление критического потенциала и радиуса R0cr
cccc    **************************************************
	do 920 qq=q1,q2,q_step
	   q=1.0/qq
	   print *,'q= ',q
	   ksi0=0.9
 11	   f=-1./(ksi0*ksi0)-q+q/((ksi0-1.)**2)+(1.+q)*ksi0
	   fs=2./(ksi0**3)-2.*q/((ksi0-1.)**3)+1.+q
	   ksi=(fs*ksi0-f)/fs
	   if(abs(ksi-ksi0).le.0.005*ksi) then
	      goto 10
	   else
	      ksi0=ksi
	      goto 11
	   end if
 10	   ksi1=ksi
	   print*, 'ksi1=', ksi1
	   omega0=omega(q,1.,0.,ksi1)
	   print*, 'omega0=', omega0
	   R0cr=casat(q,omega0,0.,1.)
	   print*, 'R0cr=', R0cr
	   R0=t*R0cr
	   omega1=omega(q,0.,1.,R0)             
cccc    ************************************************************
cccc    Вычисление геометрической формы звезды и силы тяжести на ее 
cccc    поверхности. ***********************************************
	   nn=100
	   mm=199
	   s=0.
	   sg=0.
	   do 1 i=2,nn
c       print*, 'i=',i
	      theta=(i-1)*pi/nn
	      m(i)=int(mm*sin(theta))
	      lam(i)=cos(theta)
	      do 2 k=1,m(i)+1
c       print*,'k=',k
		 phi(i,k)=pi2*(k-1)/m(i)
		 nju(i,k)=sin(theta)*cos(phi(i,k))
		 r(i,k)=casat(q,omega1,lam(i),nju(i,k))
c       print*, 'casat'
		 x(i,k)=r(i,k)*lam(i)
		 y(i,k)=r(i,k)*nju(i,k)*tan(phi(i,k))
		 z(i,k)=r(i,k)*nju(i,k)
		 r1k(i,k)=x(i,k)*x(i,k) + y(i,k)*y(i,k) + z(i,k)*z(i,k)
		 r2k(i,k)=(x(i,k)-1.)*(x(i,k)-1.)+y(i,k)*y(i,k)+z(i,k)*z(i,k)
		 dx(i,k)=q*(1.-x(i,k))/(r2k(i,k)**1.5) - x(i,k)/(r1k(i,k)**1.5)
     *            + (1.+q)*x(i,k) - q
		 dy(i,k)=-q*y(i,k)/(r2k(i,k)**1.5)-y(i,k)/(r1k(i,k)**1.5)
     *            +(1.+q)*y(i,k)
		 dz(i,k)=-q*z(i,k)/(r2k(i,k)**1.5)-z(i,k)/(r1k(i,k)**1.5)
		 g(i,k)=modul(dx(i,k),dy(i,k),dz(i,k))
		 cbeta(i,k)=-(x(i,k)*dx(i,k) + y(i,k)*dy(i,k) + z(i,k)*dz(i,k))/
     *           (g(i,k)*modul(x(i,k),y(i,k),z(i,k)))
		 w(i,k)=r(i,k)*r(i,k)*(sin(theta)/cbeta(i,k))*(pi/nn)*
     *           (pi2/m(i))
		 write(11,*) x(i,k),y(i,k),z(i,k)
c       write(12,*) g(i,k),cbeta(i,k),w(i,k)
		 s=s+w(i,k)
		 sg=sg+w(i,k)*g(i,k)
 2	      continue
 1	   continue	
	   gsr=s/sg
	   print*,'GSR', gsr, 1/gsr
	   do 3 i=2,nn
	      do 4 k=1,m(i)+1
		 Te(i,k)=Tef*((g(i,k)*gsr)**beta)
		 a(i,k)=integ(Te(i,k))
 4	      continue
 3	   continue
	   if(te_output.eq.0)then
	      te_output=1
	      do 913 i=2,nn
		 do 914 k=1,m(i)+1
 914		    if(abs(y(i,k)).lt.1.e-6)write(13,*) x(i,k),y(i,k),Te(i,k)
 913		 continue
	      endif
	      close(13)
	      
	      do 20 ang=angle1,angle2,angle_step
		 write(14,*)qq,' ',ang
		 print *, 'i=',ang
		 angle=ang*0.0174533
		 do 5 l=1,100
		    ph(l)=pi2*l*0.01
		    sI(l)=0.
		    do 6 i=2,nn
		       theta=(i-1)*pi/nn
		       m(i)=int(mm*sin(theta))
		       do 7 k=1,m(i)+1  
			  gam_and=(sin(angle)*cos(ph(l))*dx(i,k)+sin(angle)*
     *                    sin(ph(l))*dy(i,k)-cos(angle)*dz(i,k))/g(i,k)
			  if(gam_and.GE.0.) then
			     climb_and=limb(gam_and,Te(i,k))
     	                     dI_and = climb_and*gam_and*w(i,k)*a(i,k)
			     sI(l)=sI(l)+dI_and
			  end if	
 7		       continue
 6		    continue
		    res(l)=-2.5*log10(sI(l))-10.
		    write(14,*)ph(l)/pi2,' ',res(l)
 5		 continue
 20	      continue
 920	   continue
	   end
	
	
	real function integ(T)
	real h,kon,c,xx(1000),a(1000),abb(1000),f(1000),s, T
	real kont1,wl(1000),wl5(1000),filtr(1000)
	integer n,i
	common /a_common/wl,wl5,filtr
	h=6.626e-27
	kon=1.38e-16
	kont1=h/(kon*T)
	c=3.e10	
	n=1000
	do 10 i=1,n
	   xx(i)=c/(wl(i)*1.e-7)
	   a(i)=exp(xx(i)*konT1)-1.
	   abb(i)=1./(a(i)*wl5(i))
	   f(i)=abb(i)/filtr(i)
c       if(t.lt.3100)f(i)=f(i)/1.1
 10	continue
	s=0.
	do 1 i=1,n-1
	   s=s+0.5*(f(i)+f(i+1))*(wl(i+1)-wl(i))
 1	continue
	integ=s
	return
	end
	  
	real function omega(q,lam,nju,rr)
	real q,lam,nju,rr,rr2
	rr2=rr*rr
        omega=1./rr+q*((1.-2.*lam*rr+rr2)**(-0.5))-q*lam*rr+
     *  0.5*(1.+q)*(1.-nju*nju)*rr2
	return
	end




	real function casat(q,tt,lam,nju)
	real lam,nju,tt,q,r,f,r0,fs,omega,ff,r02
	r0=0.01
 1	ff=omega(q,lam,nju,r0)
	f=ff-tt
	r02=r0*r0
	fs=-1./(r02)-q*(r0-lam)*((1.-2.*lam*r0+r02)**(-1.5))-
     *     q*lam+r0*(1.+q)*(1.-nju*nju)
	r=(fs*r0-f)/fs
	if(abs(r-r0).le.0.001*r) then
	   goto 10
	else
	   r0=r
	   goto 1
	end if
 10     casat=r
	return
	end

	real function modul(x,y,z)
	real x,y,z
        modul=sqrt(x*x+y*y+z*z)
	return
	end
	
	real function limb(gam,Te)
	real gam,Te,tem(20),a1(20),a2(20),a3(20),a4(20),f,f1,fff
	integer k
	common tem,a1,a2,a3,a4
	k=int(1.+(Te-2000.)/200.)
	f=fff(gam,a1(k),a2(k),a3(k),a4(k))
	f1=fff(gam,a1(k+1),a2(k+1),a3(k+1),a4(k+1))
	limb=((f1-f)*Te+f*tem(k+1)-f1*tem(k))/(tem(k+1)-tem(k))
	return
	end

	real function fff(gam,a1,a2,a3,a4)
	real gam,a1,a2,a3,a4
	fff=1.-a1*(1.-sqrt(gam))-a2*(1.-gam)-a3*
     *  (1.-(gam**1.5))-a4*(1.-gam*gam)
	return
	end
