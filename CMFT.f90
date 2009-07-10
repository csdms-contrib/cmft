program scarp
!author: Giulio Mariotti, Boston University, 2008/2009, giuliom@bu.edu

implicit none

!constants parameters
double precision,parameter:: pi=3.14, g=9.81, rhow=1000, rhos=2300*0.7, rhoa=1.3 ![kg/m3]
!wave period
double precision,parameter:: T=2  ![s]
!Diffusion coefficient for the suspended sediment
double precision,parameter:: D=1 ![m2/s]
!erosion parameter
double precision, parameter::alpha=4.12/10000,beta=1.5/100000 !erosion parameter [],[]
double precision, parameter::teta_cr=0.7, P_cr=20 ! critical value for shear stress and for impact erosion [Pa],[W]
!settling parameter
double precision, parameter::ws=.1/1000,dp=0.001,teta_dep=.1  !settling velocity [m/s], particel diameter [m],shear threshold deposition [Pa]
!vegetation parameters
double precision, parameter::hmax=.9,hmin=.1,dh=hmax-hmin,B_max=2000,K_veg=5  !vegetation increasing effect on critical shear stress
double precision,parameter::k_b=0.009/12/30/24 !organogenic sedimentation rate [m/h]
double precision,parameter::att=0.03 !wave attenuation [%]
!tide setting
double precision, parameter:: h_tide=2,ho=0 !tide amplitude [m], mean sea level [m]
!wind setting
double precision, parameter::durata=12,cicli_calma=10 !lenght of the wind event in [h], number of calm cycle between two wind events
double precision,parameter::Umax=20.
!concentartion of sediment at the seaward boundary
double precision,parameter:: Ccost=.01![g/l]                               !USE THIS AS CONTROL PARAMETER
!Sea level rise
double precision, parameter:: RSLR=0*0.005 !relative sea level rise [m/year]  !USE THIS AS CONTROL PARAMETER
!space and time structure
double precision,parameter::dx=0.1, dxC=2, dt=0.5 ![m],[m],[h]
integer, parameter:: cell=5000,cellS=floor(cell*dx/dxC),tmax=10000! bottom cell size, water cell cell size, max # iter


!model variables
!mean water depth [m], amount of sediment in the water column, mean tidal current, concentration of sediemnt in the water coulumn
double precision:: ymean(cellS)=1,S(cellS)=0,qmean(cellS)=0,C(cellS)=0
!water depth [m],bottom elevation [m],tide current[m/s],wave energy,Wave Power[W], Wave length [m],shear stress [Pa], aboveground biomass []
double precision :: y(cell)=0,zg(cell),q(cell),E(cell),P(cell),L(cell),tau(cell),B(cell)=0
double precision:: U,h=h_tide/2,htemp,tempo=0,tempodt=0,tempold,rr
!wind velocity [m/s], water elevation [m],times...
integer ::tt=0,cicli=1 

!MAIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call random_seed()
call init()  !initialization 
do while  (tt<tmax)  !time loop
tt=tt+1
tempodt=tempodt+dt ! [h], lenght of the wind event
tempold=tempo
tempo=tempo+dt/24  ! [day], cumulative time

!main simulation: one period of calm plus one period of wind 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!calm event simulation
	if  (tempodt>durata)    then 
	tempodt=0		
	call calma(cicli)
	end if

	!wind event simulation
	call random_number(rr)
	U=rr*Umax
	htemp=h
	h=ho+h_tide/2*sin(tt*2*pi*dt/12)  !update water level with the tide
	
	!processes
	zg=zg-(tempo-tempold)/365*RSLR !relative sea level rise
	call biomass(B,zg,cicli) 
	call waterdepth(h,zg,y,ymean)
	call tidalcurrent(h,htemp,y,q,qmean)	
	call windwave(E,y,P,L,B,U)
	call shearstress(E,L,y,tau,q)
	call erosion(zg,S,P,y,tau,B,cicli)
	call advdiff(S,qmean,ymean)	
	call sedimentation(S,y,ymean,tau,zg,cicli)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	!drawing
	if (mod(tt,24*30+1)==0) call output(tempo,dx,h,zg,cell)
			
end do
  call outputfinal(tempo,dx,h,zg,cell)
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init()
integer:: i
double precision:: dumb

!read the initial bottom configurtation
open (1, FILE ='zgin.dat')
do i=1,cell
read(1,*) dumb,zg(i)
end do
close(1)

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calma(cicli)
integer:: cicli	
		U=0
		!The first cycle is run simply setting zero wind for 12 hours (a tidal cycle)
		call evo_nowind()
		!the second cycle is run for 12 hours, but moltypling the erosional and depositional results
		cicli=cicli_calma	
		call evo_nowind()
		tempo=tempo+dt*(1+cicli) 
		cicli=1
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine evo_nowind()  !evolution with zero wind (calm period)
integer:: i
	do i=1,24  
		htemp=h
	        h=ho+h_tide/2*sin(tt*2*pi*dt/12)
		call waterdepth(h,zg,y,ymean)
		call biomass(B,zg,cicli)
		call tidalcurrent(h,htemp,y,q,qmean)
		call shearstress0(y,tau,q)
		P=0
		E=0
		call erosion(zg,S,P,y,tau,B,cicli)
		call advdiff(S,qmean,ymean)
		call sedimentation(S,y,ymean,tau,zg,cicli)
		tt=tt+1
	end do
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine waterdepth(h,zg,y,ymean)  
double precision:: h
double precision:: y(cell), zg(cell)
double precision:: ymean(floor(cell*dx/dxC))
integer:: i
ymean=0
do i=1,cell
 if (h>zg(i)) then
 y(i)=(h-zg(i))   
 ymean(1+floor((i-0.5)*dx/dxC))=ymean(1+floor((i-0.5)*dx/dxC))+y(i)
 else
 y(i)=0
 end if
end do
ymean=ymean/dxC*dx
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tidalcurrent(h,htemp,y,q,qmean)
double precision :: qmean(floor(cell*dx/dxC)),h,htemp
double precision :: q (cell)
double precision:: y(cell)
integer:: i
q(1)=0.
do i=2,cell
 q(i)=q(i-1);   
 if (y(i)>0)    then
 q(i)=q(i)+(htemp-h)*dx/(3600*dt)  !current positive if going toward right
 else
 q(i)=0    
 end if
 qmean(1+floor((i-0.5)*dx/dxC))=qmean(1+floor((i-0.5)*dx/dxC))+q(i)
end do
qmean=qmean/dxC*dx;
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine windwave(E,y,P,L,B,U)
double precision::  E(cell),y(cell),P(cell),L(cell),B(cell),U
double precision:: cg, lengthdx!, U
integer:: i,j
!first loop, calcultaing the wind on the first seaward cell, propagating from a far point, 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
P=0
E=0
E(cell)=.01**2/8*rhow*g  !initial value for the wave height: 1 cm
j=cell
lengthdx=100  !length of the cell, in this case 100 times the cell. 300*100*dx=3000m= 3 km 
do i=1,300
  call wavepropagation(lengthdx,cg,j)    
 end do

!main loop for the wave energy propagation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
lengthdx=1
E(cell)=E(cell)
i=0
do while (i<cell ) 
j=cell-i
	if (y(j)<=0.01) then 
		P(j)=E(j)*cg/dx
		P(j+1)=P(j)  
		i=cell+9999 !to go out of the loop   
		E(j)=0
	else
		call wavepropagation(lengthdx,cg,j)
	end if 
	E(j-1)=E(j) 
	i=i+1 
end do
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wavepropagation(lengthdx,cg,j)
double precision:: LL,Swg,Sbrk,Swc,Sbf,H,cg,k,sigma
double precision:: Qb,Qb1,lengthdx!
integer:: j

		L(j)=1.56*T**2  !guess an initial wave length
		LL=0
		do while (abs(L(j)-LL)>L(j)*0.1)
		LL=L(j)
		L(j)=((9.81*T**2)/(2*pi))*tanh(2*pi*y(j)/L(j))  !dispersion equation
		end do
		k=2*pi/L(j) 	
		sigma=2*pi/T
		cg=L(j)/T*0.5*(1+2*k*y(j)/sinh(2*k*y(j)))  !group velocity
		H=sqrt(E(j)*8/rhow/g)
		Swg=80*(0.001*rhoa/rhow/g/k)**2*sigma*U**4 + E(j)*5*rhoa/rhow/T*(U/L(j)*T-0.9)
		Sbf=-4*0.015*pi*H/T*k/sinh(k*y(j))/sinh(2*k*y(j))*E(j)
		Swc=-3./100000/T*(E(j)/T**4/g**2/0.00457)**2*E(j)
        	Sbrk=0.
              		  if (E(j)>0) then
              	 	  Qb=0.2;
              	          Qb1=Qb;           
             	          Qb=exp((Qb1-1)/((H/(0.78*(y(j))))**2))
             	          do while (abs(Qb1-Qb)>0.001)
           	  	  Qb1=Qb   
           	          Qb=exp((Qb-1)/((H/(0.78*(y(j))))**2))
       	   		  end do
       	   		  Sbrk=-2/T*Qb*(0.78*(y(j))/H)**2*E(j)
			  end if
		E(j)=E(j)+(Swg+Sbf+Swc+Sbrk)/cg*dx*lengthdx  !wave propagate to the next cell
		!controls
		if (E(j)<0) E(j)=0.    !no wave less than zero
		if (B(j)>0)	E(j)=((1-att*B(j)*dx)*sqrt(E(j)))**2  !attenuation from the vegetation 
		if (sqrt(E(j)*8/rhow/g)>0.78*(y(j))) E(j)=(0.78*y(j))**2/8*rhow*g !breaking if wave too high
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine shearstress(E,L,y,tau,q)
double precision:: L(cell),E(cell),tau(cell),y(cell),q(cell)
double precision:: tau_w, tau_c,k
integer:: i
tau=0
do i=1,cell
if (y(i)>0.01) then 
k=2*pi/L(i)
tau_w=0.5*0.004*rhow*(pi*sqrt(E(i)*8/rhow/g)/T/sinh(k*y(i)))
tau_c=0.01*rhow*((q(i))/y(i))**2
if (tau_w>0 .or. tau_c>0) tau(i)=tau_w+tau_c*(1+1.2*(tau_w/(tau_c+tau_w))**1.5)
end if
end do
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine shearstress0(y,tau,q)
double precision:: tau(cell),y(cell),q(cell)
double precision:: tau_w, tau_c,k
integer:: i
tau=0
do i=1,cell
if (y(i)>0.01) then
tau(i)=0.01*rhow*((q(i))/y(i))**2
end if
end do
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine erosion(zg,S,P,y,tau,B,cicli)
double precision:: zg(cell),S(cellS),B(cell)
double precision:: y(cell),P(cell),tau(cell)
integer::i,j,cicli
double precision:: R
do j=1,cell
i=cell+1-j

!bottom shear stress erosion
if (y(i)>0.) then
R=cicli*dt*3600*alpha*(tau(i)-teta_cr*(1+B(i)*K_veg))
if (R<0) R=0.
zg(i)=zg(i)-R/rhos
S(1+floor((i-0.5)*dx/dxC))=S(1+floor((i-0.5)*dx/dxC))+R*dx
end if

!impact erosion
if ( P(i)>0 .and. zg(i)>zg(i+1)) then 
R=cicli*dt*3600*beta*(P(i)-P_cr*(1+B(i)*K_veg))
if (R<0) R=0
if (R>rhos*(zg(i)-zg(i+1))) then
R=rhos*(zg(i)-zg(i+1))
zg(i)=zg(i+1)
else
zg(i)=zg(i)-R/rhos
end if
if (i<cell-floor(dxC/dx)) S(2+floor((i-0.5)*dx/dxC))=S(2+floor((i-0.5)*dx/dxC))+R*dx
if (i>=cell-floor(dxC/dx)) S(1+floor((i-0.5)*dx/dxC))=S(1+floor((i-0.5)*dx/dxC))+R*dx
end if

end do
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine biomass (B,zg,cicli)
integer::cicli
double precision:: B(cell),zg(cell)
double precision::season
integer::j
season=(0.1+0.9*0.5*(1+sin(2*pi*tempo/365-pi*0.5)))
B=0
do j=1,cell
if (zg(j)>hmin+ho .and. zg(j)<hmax+ho) B(j)= (-4*((zg(j)-(hmin+ho))/dh)**2+4*(zg(j)-(hmin+ho))/dh)*season
end do
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sedimentation(S,y,ymean,tau,zg,cicli)
double precision:: S(cellS),ymean(cellS),C(cellS),y(cell),tau(cell),zg(cell)
double precision:: Sed
integer:: i,cicli
double precision:: uu,ds
C=0
do i=1,cellS
if (ymean(i)>0) C(i)=S(i)/(ymean(i)*dxC)
end do
do i=1,cell
	Sed=0
	if (y(i)>0 ) then  !you can setlle	
		if (tau(i)<teta_dep) Sed=C(1+floor((i-0.5)*dx/dxC))*2*ws*(1-tau(i)/teta_dep)*dt*3600   !sedimentation 
		if (B(i)>0) then
			uu=abs(q(i))/y(i)*10
			ds=0.0006*(B(i)*B_max)**0.3
			Sed=Sed+ds*uu*0.224*(uu*ds*1000000)**0.718*(dp/ds)**2.08*250*(B(i)*B_max)**0.3032*0.0609*(B(i)*B_max)**0.1876 
		end if
		if (Sed>y(i)*C(1+floor((i-0.5)*dx/dxC))) Sed=y(i)*C(1+floor((i-0.5)*dx/dxC))  !to much sedimentation
		zg(i)=zg(i) +Sed/rhos;
		S(1+floor((i-0.5)*dx/dxC))=S(1+floor((i-0.5)*dx/dxC))-Sed*dx; 
	end if
	if (B(i)>0) zg(i)=zg(i)+B(i)*k_b*dt*cicli 
end do   

!to sediment everything when the water depth become zero
do i=1,cellS-1
	if (ymean(i)==0 .and. S(i)>0) then
	S(i+1)=S(i+1)+S(i)
	S(i)=0
	end if
end do
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine advdiff(S,qmean,ymean)
double precision::S(cellS),qmean(cellS),ymean(cellS)
  integer, parameter :: itmax = 100, maxnz = 3,n = cellS,inw = 2*cellS,nw =20*cellS 
  external omin
  integer i,ier,ip(n),iparm(30), iwksp(inw),jcoef(n,maxnz), maxnzz,mdim,ndim,p(n)
  external mic1
  real ( kind = 8 ) coef(n,maxnz),rhs(n),rparm(30),u(n),ubar(n), wksp(nw)
  
!This is a subroutine to solve the advection diffusion equation of the suspendend sediment
!the input are the sediemnt distribution, the current velcoity and the diffusion coefficient
!the output is the concentration at the next time step
!The adv diff equation is solved witht centered difference for diffusion and upwind for advection
!Since the water flow direction change direction with time, the strucute of the matrix is different if the current goes to right or to left, because
!the discrete difference have always to be upwind!!
!the linear system is solve witth the biconjugate gradient method
!the coefficient matrix is handled as a sparse matrix
  
 if (qmean(cellS)<0)  then
 S(cellS)=S(cellS)+Ccost*  (-qmean(cellS)/dxC +D/(dxC**2))*dt*3600
 end if
 
 if (qmean(cellS)>0)  then
 S(cellS)=S(cellS) +  Ccost*D/(dxC**2)*dt*3600
 end if
    
  do i = 1, n
  ubar(i) = real ( i, kind = 8 )
  end do
  rhs=S
  call dfault ( iparm, rparm )
  iparm(2) = itmax
  iparm(12) = 1
  mdim = 3
  ndim = n
  maxnzz = maxnz
  
  coef=0
  jcoef(1,1:3) =  (/ 0, 1, 2 /)
  jcoef(n,1:3) = (/ n-1, n, 0 /) 
  
  
if (qmean(cellS)>0) then !current from left to right


    if (ymean(1)>0) coef(1,2)=-qmean(1)/ymean(1)/dxC 
    if (ymean(cellS-1)>0) coef(cellS,1)=+qmean(cellS-1)/ymean(cellS-1)/(dxC)
   
    !termine no cost
    if (ymean(cellS)>0)   coef(cellS,2)=-qmean(cellS)/ymean(cellS)/(dxC) !questo � quello che esce dall'ultima cella
    
    if (ymean(cellS)>0 .and. ymean(cellS-1)>0) then
    coef(cellS,2)=coef(cellS,2)-D/(dxC**2)/ymean(cellS) 
    coef(cellS,1)=coef(cellS,1)+D/(dxC**2)/ymean(cellS-1) 
    
    !termine aggiuntivo di diffusione
    coef(cellS,2)=coef(cellS,2)-D/(dxC**2)/ymean(cellS)
       
    end if
    if (ymean(1)>0 .and. ymean(2)>0) then
    coef(1,2)=coef(1,2)-D/(dxC**2)/ymean(1) 
    coef(1,3)=coef(1,3)+D/(dxC**2)/ymean(2) 
    end  if 
   
    do i = 2, n-1
    if (ymean(i-1)>0) coef(i,1)=+qmean(i-1)/ymean(i-1)/dxC ;
    if (ymean(i)>0)    coef(i,2)=-qmean(i)/ymean(i)/dxC;  
    if (ymean(i+1)>0 .and. ymean(i)>0) then
    coef(i,3)=coef(i,3)+D/(dxC**2)/ymean(i+1) ;
    coef(i,2)=coef(i,2)-D/(dxC**2)/ymean(i) ;
    end if
    if (ymean(i-1)>0 .and. ymean(i)>0) then
    coef(i,1)=coef(i,1)+D/(dxC**2)/ymean(i-1) ;
    coef(i,2)=coef(i,2)-D/(dxC**2)/ymean(i) ;
    end if
    jcoef(i,1:3) = (/ i-1, i, i+1 /)
    end do
    
        
else !current from right to left
 
  
    if (ymean(2)>0) coef(1,3)=-qmean(2)/ymean(2)/(dxC) 
    if (ymean(cellS)>0)  coef(cellS,2)=+qmean(cellS)/ymean(cellS)/(dxC)      
    if (ymean(cellS)>0 .and. ymean(cellS-1)>0) then
    coef(cellS,2)=coef(cellS,2)-D/(dxC**2)/ymean(cellS) 
    coef(cellS,1)=coef(cellS,1)+D/(dxC**2)/ymean(cellS-1) 
    !termine aggiuntivo di diffusione
    coef(cellS,2)=coef(cellS,2)-D/(dxC**2)/ymean(cellS)     
    
    end if
    if (ymean(1)>0 .and. ymean(2)>0) then
    coef(1,2)=coef(1,2)-D/(dxC**2)/ymean(1) 
    coef(1,3)=coef(1,3)+D/(dxC**2)/ymean(2) 
    end  if
 
	do i=2,cellS-1
    if (ymean(i)>0)  coef(i,2)=+qmean(i)/ymean(i)/(dxC)         
    if (ymean(i+1)>0)  coef(i,3)=-qmean(i+1)/ymean(i+1)/(dxC)      
    if (ymean(i+1)>0 .and. ymean(i)>0) then
    coef(i,3)=coef(i,3)+D/(dxC**2)/ymean(i+1) 
    coef(i,2)=coef(i,2)-D/(dxC**2)/ymean(i) 
    end if
    if (ymean(i-1)>0 .and. ymean(i)>0) then
    coef(i,1)=coef(i,1)+D/(dxC**2)/ymean(i-1) 
    coef(i,2)=coef(i,2)-D/(dxC**2)/ymean(i) 
    end if
    jcoef(i,1:3) = (/ i-1, i, i+1 /)
	end do
end if
   do i=1,cellS
   coef(i,2)=1-dt*3600*coef(i,2)
   coef(i,1)=-dt*3600*coef(i,1)
   coef(i,3)=-dt*3600*coef(i,3)
   end do

   !final call of the solver
    u=S 
    call nspcg ( mic1, omin, ndim, mdim, n, maxnzz, coef, jcoef, p, ip, u, ubar, &
    rhs, wksp, iwksp, nw, inw, iparm, rparm, ier ) 
    S=u

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output(tempo,dx,h,zg,cell)
double precision:: dx,tempo
integer::cell
double precision:: zg(cell),h
integer:: i

open (unit=2,file="zg.dat",status="replace")
do i=1,cell
write(2,*) i*dx,zg(i)
end do
close(2)


end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine outputfinal(tempo,dx,h,zg,cell)
double precision:: dx,tempo
integer::cell
double precision:: zg(cell),h
integer:: i

open (unit=2,file="zgfinal.dat",status="replace")
do i=1,cell
write(2,*) i*dx,zg(i)
end do
close(2)


end subroutine
end program scarp 
