module constants
  integer,parameter::Nh=10,nxn=100
  real*8,parameter::R_con=1.D0,a_con=1.2D0,b_con=-1.D0,c0_con=1.1D0,eta=1D-6,epsil=1D-12
end module constants

program solovev
  use constants
  implicit none
  real*8::z_X,zeta_X,psi_h,delt_r,delt_z,delt_tht,psi,r_min,r_max,z_min,z_max,tht,r,z,zeta
  real*8::Coe(1:4),x(3),psi_X,psi_p,multipol,AA(Nh,Nh),bb(Nh),c(1:Nh),temp1,temp2
  real*8::temp3,temp4,z_temp
  integer::i,j,k,ii,jj,start1,start2,finish1,finish2,finish3,start,finish
  real*8,dimension(1000,1000)::psi_t,RR,ZZ
  !double-null equilibria
  call system_clock(start)
  psi_X=(a_con-c0_con)*(b_con+c0_con)**2*R_con**4/(8*c0_con**2)
  zeta_X=-(b_con+c0_con)/(2.D0*c0_con)*R_con
  z_X=((b_con+c0_con)*(a_con-c0_con)/(2*c0_con**2))**(0.5)*R_con
  do i=1,Nh
     tht=(i-0.5D0)*4.D0*atan(1.D0)/(Nh*1.D0)!tht 选取需要满足tan(tht).nq.0,不趋近无穷
     Coe(1)=c0_con*R_con*(tan(tht))**2
     Coe(2)=0.50*(b_con+c0_con)*(tan(tht))**2*R_con**2+0.5*(a_con-c0_con)*R_con**2
     Coe(3)=0.D0
     Coe(4)=-(1.D0-eta)*psi_X
     call getCubicRoot(Coe,x)
     if(cos(tht)>=0.D0)then
        zeta=max(x(1),x(2),x(3));r=sqrt(2.D0*R_con*zeta+R_con**2);z=tan(tht)*zeta
     else if(cos(tht)<0.D0)then
        zeta=x(1)+x(2)+x(3)-max(x(1),x(2),x(3))-min(x(1),x(2),x(3));
        r=sqrt(2.D0*R_con*zeta+R_con**2);z=tan(tht)*zeta
     endif
     do j=1,Nh
        AA(i,j)=multipol(2*(j-1),r,z)
     enddo
     call calc_psi_p(r,z,psi_p)
     bb(i)=(1.D0-eta)*psi_X-psi_p
  enddo

  open(1,file='dat1')
    do i=1,Nh
       write(1,*)AA(i,1:Nh),bb(i)
    enddo
  close(1)

  call doolittle(AA,bb,c)
  !  call Gauss_Seidel(AA,bb,c)
  print*,c
  print*,(matmul(AA,c)-bb)/bb(1)
  !计算psi_total
  r_min=0.9D0;r_max=1.07D0;
  z_min=-0.09D0;z_max=0.09D0
!  !$OMP PARALLEL
  !  !$OMP DO
  call system_clock(start1)
        delt_r=(r_max-r_min)/100.D0
        delt_z=(z_max-z_min)/100.D0
!open(2,file='psito.txt')
 !       !$OMP PARALLEL
!        !$OMP DO
  do i=1,100
     do j=1,100
        r=r_min+i*delt_r;RR(i,j)=r
        z=z_min+j*delt_z;ZZ(i,j)=z
        zeta=(r**2-R_con**2)/(2.D0*R_con);
        temp1=0.5*(b_con+c0_con)*R_con**2*z**2+c0_con*R_con*zeta*z**2
        temp2=0.5*(a_con-c0_con)*R_con**2.D0*zeta**2.D0
        psi_t(i,j)=temp1+temp2
        if(psi_t(i,j)<(1-eta)*psi_X.and.z<z_X.and.z>-z_X)then
           temp3=(1.D0-eta)*psi_X-0.5D0*(a_con-c0_con)*R_con**2*zeta**2
           temp4=0.5D0*(b_con+c0_con)*R_con**2+c0_con*R_con*zeta
           z_temp=sqrt(temp3/temp4)
           if(z>-z_temp.and.z<z_temp)then
           else
              call calc_psi_p(r,z,psi_p)
              psi_h=0.D0
              do k=1,Nh
                 psi_h=psi_h+c(k)*multipol(2*(k-1),r,z)
              enddo
              psi_t(i,j)=psi_p+psi_h
           endif
        else
           call calc_psi_p(r,z,psi_p)
           psi_h=0.D0
           do k=1,Nh
              psi_h=psi_h+c(k)*multipol(2*(k-1),r,z)
           enddo
           psi_t(i,j)=psi_p+psi_h
        endif
        write(2,*)r,z,psi_t(i,j)
     enddo
  enddo
  !!$OMP end do
  !!$OMP end parallel
!close(2)
open(2,file='psito_my.txt')
do i=1,100
   do j=1,100
      write(2,*)RR(i,j),zz(i,j),psi_t(i,j)
   enddo
enddo
close(2)
call system_clock(finish1)
open(4,file='psitoso.txt')
  do i=1,100
     do j=1,100
        r=r_min+i*delt_r
        z=z_min+j*delt_z
        zeta=(r**2-R_con**2)/(2.D0*R_con);
        temp1=0.5*(b_con+c0_con)*R_con**2*z**2+c0_con*R_con*zeta*z**2
        temp2=0.5*(a_con-c0_con)*R_con**2.D0*zeta**2.D0
        psi_t(i,j)=temp1+temp2
        write(4,*)r,z,psi_t(i,j)
     enddo
  enddo
close(4)
call system_clock(finish2)
  !   !$OMP END DO
 ! !$OMP END PARALLEL
open(3,file='LCFS.txt')
  do i=1,200
     delt_tht=8.D0*atan(1.D0)/200.D0
     tht=(i-0.5D0)*delt_tht
     Coe(1)=c0_con*R_con*(tan(tht))**2
     Coe(2)=0.50*(b_con+c0_con)*(tan(tht))**2*R_con**2+0.5*(a_con-c0_con)*R_con**2
     Coe(3)=0.D0
     Coe(4)=-psi_X
     call getCubicRoot(Coe,x)
     if(cos(tht)>=0.D0)then
        zeta=max(x(1),x(2),x(3));r=sqrt(2.D0*R_con*zeta+R_con**2);z=tan(tht)*zeta
     else if(cos(tht)<0.D0)then
        zeta=x(1)+x(2)+x(3)-max(x(1),x(2),x(3))-min(x(1),x(2),x(3));
        r=sqrt(2.D0*R_con*zeta+R_con**2);z=tan(tht)*zeta
     endif
     call calc_psi_p(r,z,psi_p)
     psi_h=0.D0
     do k=1,Nh
        psi_h=psi_h+c(k)*multipol(2*(k-1),r,z)
     enddo
!     write(3,*)tht/(8*atan(1.D0)),r,z,psi_p
     write(3,*)tht/(8*atan(1.D0)),(psi_p+psi_h-psi_X)/psi_X
  enddo
close(3)
call system_clock(finish3)
z=0.D0
open(7,file='psitoz0.txt')
  do i=1,100
!     do j=1,100
        r=r_min+i*delt_r
!        z=z_min+j*delt_z
        zeta=(r**2-R_con**2)/(2.D0*R_con);
        temp1=0.5*(b_con+c0_con)*R_con**2*z**2+c0_con*R_con*zeta*z**2
        temp2=0.5*(a_con-c0_con)*R_con**2.D0*zeta**2.D0
        psi_t(i,j)=temp1+temp2
        if(psi_t(i,j)<psi_X.and.z<z_X.and.z>-z_X)then
           temp3=(1.D0-eta)*psi_X-0.5D0*(a_con-c0_con)*R_con**2*zeta**2
           temp4=0.5D0*(b_con+c0_con)*R_con**2+c0_con*R_con*zeta
           z_temp=sqrt(temp3/temp4)
           if(z>z_temp.or.z<-z_temp)then
              call calc_psi_p(r,z,psi_p)
              psi_h=0.D0
              do k=1,Nh
                 psi_h=psi_h+c(k)*multipol(2*(k-1),r,z)
              enddo
              psi_t(i,j)=psi_p+psi_h
           endif
        else
              call calc_psi_p(r,z,psi_p)
              psi_h=0.D0
              do k=1,Nh
                 psi_h=psi_h+c(k)*multipol(2*(k-1),r,z)
              enddo
              psi_t(i,j)=psi_p+psi_h
        endif
        write(7,*)r,z,psi_t(i,j)/psi_X
 !    enddo
  enddo
close(7)
  
open(8,file='psitosoz0.txt')
  do i=1,100
!     do j=1,100
        r=r_min+i*delt_r
!        z=z_min+j*delt_z
        zeta=(r**2-R_con**2)/(2.D0*R_con);
        temp1=0.5*(b_con+c0_con)*R_con**2*z**2+c0_con*R_con*zeta*z**2
        temp2=0.5*(a_con-c0_con)*R_con**2.D0*zeta**2.D0
        psi_t(i,j)=temp1+temp2
        write(8,*)r,z,psi_t(i,j)/psi_X
!     enddo
  enddo
close(8)

call system_clock(finish)
print*,'c time:',start1-start,'ms'
print*,'psito_my time:',finish1-start1,'ms'
print*,'psito_so time:',finish2-finish1,'ms'
print*,'LCFS time',finish3-finish2,'ms'
print*,'z0 time',finish-finish3,'ms'
print*,'all time:',finish-start,'ms'

end program solovev
subroutine doolittle(A,b,x)
  use constants
  implicit none
  real*8::A(Nh,Nh),b(Nh),x(Nh),y(Nh),L(Nh,Nh),U(Nh,Nh),sum1
  integer::i,j,k,r
  u(:,:)=0.D0
  l(:,:)=0.D0
  u(1,1:Nh)=A(1,1:Nh);
  do i=1,Nh
     l(i,i)=1.D0
  enddo
  do i=2,Nh
     l(i,1)=a(i,1)/u(1,1)
  enddo

  do k=2,Nh
     do j=k,Nh
        sum1=0.D0
        do r=1,k-1
           sum1=sum1+l(k,r)*u(r,j)
        enddo
        u(k,j)=a(k,j)-sum1
     enddo
     if(k<Nh)then
        do i=k+1,Nh
           sum1=0.D0
           do r=1,k-1
              sum1=sum1+l(i,r)*u(r,k)
           enddo
           l(i,k)=(a(i,k)-sum1)/u(k,k)
        enddo
     endif
  end do

  y(1)=b(1)
  do i=2,Nh
     sum1=0.D0
     do j=1,i-1
        sum1=sum1+l(i,j)*y(j)
     enddo
     y(i)=b(i)-sum1
  end do
  x(Nh)=y(Nh)/u(Nh,Nh)
  do k=2,Nh
     i=Nh-k+1
     sum1=0.D0
     do j=i+1,Nh
        sum1=sum1+u(i,j)*x(j)
     enddo
     x(i)=(y(i)-sum1)/u(i,i)
  enddo
end subroutine doolittle


subroutine Gauss_Seidel(A,b,x)
  implicit none
  integer,parameter::Nh=10
  real*8::A(Nh,Nh),b(Nh),x(Nh),temp_err(Nh),err,tol,temp1,sum1,sum2
  integer::i,j,k
!  x(:)=1.D0;
  tol=1.D-14
  err=2*tol
  k=0
  do while(err>tol)
     do i=1,Nh
        sum1=0.D0
        do j=1,Nh
           if(j.ne.i)then
              sum1=sum1+A(i,j)*x(j)
           endif
        enddo
        temp1=-(sum1-b(i))/A(i,i)
        temp_err(i)=abs(x(i)-temp1)
        x(i)=temp1
     enddo
     do i=1,Nh-1
        err=max(temp_err(i),temp_err(i+1))
     enddo
     print*,err,temp_err
     k=k+1
     if(k==5)then
        stop 50
     endif
  enddo
end subroutine Gauss_Seidel

function multipol(n,r,z)
  implicit none
  integer::n,i
  real*8::r,z,A,multipol
     A=1.D0;multipol=0.D0
if(n>0.D0)then
  do i=0,n/2-1
     if(i>0.D0)then
        A=-(n/2+1-i)*(n/2-i)/(i*(i-0.5D0))*A
     endif
     multipol=multipol+A*r**(n-2.D0*i)*z**(2.D0*i)
  enddo
else if(n==0)then
   multipol=1.D0
endif
  RETURN
end function multipol

subroutine calc_psi_p(r0,z0,psi_p)
  use constants
  implicit none
  integer::i,j,k
  real*8::psi_X,delt_psi,delt_tht,r0,z0,r,z,zeta,temp1,temp2,Kk,Ek,k2,psi_p,psi,tht,Coe(1:4),x(1:3),k22,psi_p2
  psi_X=(a_con-c0_con)*(b_con+c0_con)**2*R_con**4/(8.D0*c0_con**2);
  delt_psi=(1.D0-eta)*psi_X/(nxn*1.D0);delt_tht=8.D0*atan(1.D0)/(nxn*1.D0)
  psi_p=0.D0;psi_p2=0.D0;k=1;
!open(10,file='test.txt')
  do i=1,nxn
     do j=1,nxn
        !psi=int((i-0.5D0)*delt_psi*1.D12)/1.D12;
        psi=(i-0.5D0)*delt_psi;
        tht=(j*1.D0-0.5D0)*delt_tht
!        tht=8*atan(1.D0)-(j*1.D0-0.5D0)*delt_tht
     Coe(1)=c0_con*R_con*(tan(tht))**2
     !     Coe(2)=0.5*(b_con+c0_con)*(tan(tht))**2*R_con**2+0.5*(a_con-c0_con)*R_con**2
     Coe(2)=0.5*(b_con+c0_con)*(tan(tht))**2*R_con**2+0.5*(a_con-c0_con)*R_con**2
     Coe(3)=0
     Coe(4)=-psi
     call getCubicRoot(Coe,x)
     if(cos(tht)>=0.D0)then
        zeta=max(x(1),x(2),x(3));r=sqrt(2.D0*R_con*zeta+R_con**2);z=tan(tht)*zeta
     else if(cos(tht)<0.D0)then
        zeta=x(1)+x(2)+x(3)-max(x(1),x(2),x(3))-min(x(1),x(2),x(3));
        r=sqrt(2.D0*R_con*zeta+R_con**2);z=tan(tht)*zeta
     endif
     k2=(4.D0*r0*r)/((r+r0)**2.D0+(z0-z)**2.D0+epsil)
   !!!!!!!!!!!!call 椭圆积分
     call ellipke(k2,1,Kk,Ek)
     temp1=sqrt(r0*r)/(8.D0*atan(1.D0))*((2.D0-k2)*Kk-2.D0*Ek)/sqrt(k2)*(a_con+b_con*R_con**2/r**2.D0)
     temp2=(zeta**2.D0+z**2.D0)/(3.D0*c0_con*zeta*z**2.D0+(a_con-c0_con)*R_con*zeta**2+(b_con+c0_con)*R_con*z**2)
     psi_p=psi_p-temp1*temp2*delt_psi*delt_tht
 !    write(10,*)k,k2,Kk,Ek,r,z,psi_p
     enddo
  enddo
!  close(10)
end subroutine calc_psi_p

!http://blog.sciencenet.cn/blog-548663-627770.html
subroutine getCubicRoot(P,x)
  implicit none
  integer::root_num
  real*8,parameter::twopi=8.D0*atan(1.D0)
  real*8::P(1:4),a,b,c,d,Alph,Beta,Delt,R1,R2,tht,X(3)
  X=0.D0
  a=P(1)
  b=P(2)/(3.D0*a)
  c=P(3)/(6.D0*a)
  d=p(4)/(2.D0*a)

  Alph=-b**3.D0+3.D0*b*c-d
  Beta=b**2.D0-2.D0*c
  Delt=Alph**2-beta**3

  if(Delt>0.D0)then
     tht=Alph+sqrt(Delt);R1=sign(abs(tht)**(1.D0/3.D0),tht)
     tht=Alph-sqrt(Delt);R2=sign(abs(tht)**(1.D0/3.D0),tht)
     x(1)=-b+R1+R2
     stop 'one root'
  else if(Delt==0.D0)then
     R1=sign(abs(Alph)**(1.D0/3.D0),Alph)
     if(R1==0.D0)then
        x(1)=-b
     else
        x(1)=-b+2.D0*R1
        x(2)=-b-R1
     endif
     stop '2 roots'
  elseif(Delt<0.D0)then
     tht=acos(Alph/(sqrt(Beta)*Beta))
     x(1)=-b+2.D0*sqrt(Beta)*cos(tht/3.D0)
     x(2)=-b+2.D0*sqrt(Beta)*cos((tht+twopi)/3.D0)
     x(3)=-b+2.D0*sqrt(Beta)*cos((tht-twopi)/3.D0)
     root_num=3
  endif
!print*,root_num
end subroutine getCubicRoot

    subroutine ellipke(k2,nk,K,E)
!C******************************************************************                                             
!       calculate complete elliptic intergral of                                                                 
!      the first and second kinds ( K and E )                                                                    
!        note: the argument k2=kappa**2                                                                          
!            k2=k2(1:nk)                                                                                         
!     use ellf elle                                                                                              
!       by Y. Sun                                                                                                
!C******************************************************************                                             
        double precision k2(*),K(*),E(*)
        !real(8) :: k2(:), K(:), E(:)!！！！！！！是语言版本的问题嘛？？？？？？？？？？？？？                   
    integer i,nk 
        !double precision k0,ellf,elle                                                                           
    double precision k0
    real(8) :: ellf,elle
        do i=1,nk
           k0=sqrt(k2(i))       !  change the argument to k^2                                                    
           K(i)=ellf(3.1415926535897932d0/2.0,k0)
           E(i)=elle(3.1415926535897932d0/2.0,k0)
    end do
    end subroutine ellipke

    FUNCTION ellf(phi,ak)
        double precision ellf,ak,phi
!******************************************************************                                              
!       C USES rf                                                                                                
!       Legendre elliptic integral of the 1st kind F(.; k), evaluated using Carlson's function RF .              
!       The argument ranges are 0 <=phi<=pi/2, 0 <=k sin(phi)<=1.                                                
!******************************************************************                                              
        !double precision s,rf  ! for internal function rf, it is not necessary to define the data type          
    double precision s
    real(8) ::rf 
    !double precision s                                                                                          
        s=sin(phi)
        ellf=s*rf(cos(phi)**2.0d0,(1.0d0-s*ak)*(1.0d0+s*ak),1.0d0)
        return
        END FUNCTION ellf

        FUNCTION elle(phi,ak)
        double precision elle,ak,phi
!C******************************************************************                                             
!       C USES rd,rf                                                                                             
!       Legendre elliptic integral of the 2nd kind E(.; k), evaluated using Carlson's functions RD               
!       and RF . The argument ranges are  0 <=phi<=pi/2, 0 <=k sin(phi)<=1.                                      
!C******************************************************************                                             
        !double precision cc,q,s,rd,rf                                                                           
    double precision cc,q,s
    real(8) ::rd,rf 
    !double precision cc,q,s                                                                                     

        s=sin(phi)
        cc=cos(phi)**2
        q=(1.-s*ak)*(1.+s*ak)
        elle=s*(rf(cc,q,1.0d0)-((s*ak)**2)*rd(cc,q,1.0d0)/3.0d0)
        return
      end function elle

        FUNCTION rd(x,y,z)
!C******************************************************************                                             
!       Computes Carlson's elliptic integral of the second kind, RD(x; y; z). x and y must be                    
!       nonnegative, and at most one can be zero. z must be positive. TINY must be at least twice                
!       the negative 2/3 power of the machine overflow limit. BIG must be at most 0:1.ERRTOL                     
!       times the negative 2/3 power of the machine underflow limit.                                             
!C******************************************************************                                             
        double precision rd,x,y,z,ERRTOL,TINY,BIG
        double precision C1,C2,C3,C4,C5,C6
        PARAMETER (ERRTOL=.05,TINY=1.e-25,BIG=4.5E21,C1=3./14.,C2=1./6.,C3=9./22.,C4=3./26.,C5=.25*C3,C6=1.5*C4)

        double precision alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee
        double precision fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
        if(min(x,y).lt.0..or.min(x+y,z).lt.TINY.or.max(x,y,z).gt.BIG)pause 'invalid arguments in rd'
        xt=x
        yt=y
        zt=z
        sum=0.
        fac=1.
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=.25*fac
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=.2*(xt+yt+3.*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
        ea=delx*dely
        eb=delz*delz
        ec=ea-eb
        ed=ea-6.*eb
        ee=ed+ec+ec
        rd=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
    return
    END FUNCTION rd

    !rf(cos(phi)**2,(1.0d0-s*ak)*(1.0d0+s*ak),1.0d0)                                                             
        FUNCTION rf(x,y,z)
!C******************************************************************                                             
!       Computes Carlson's elliptic integral of the ^Lrst kind, RF (x; y; z). x, y, and z must be                
!       nonnegative, and at most one can be zero. TINY must be at least 5 times the machine                      
!       underflow limit, BIG at most one ^Lfth the machine overflow limit.                                       
!C******************************************************************                                             
        double precision rf,x,y,z,ERRTOL,TINY,BIG,THIRD
        double precision C1,C2,C3,C4
        PARAMETER (ERRTOL=.08,TINY=1.5e-38,BIG=3.E37,THIRD=1./3.,C1=1./24.,C2=.1,C3=3./44.,C4=1./14.)
        double precision alamb,ave,delx,dely,delz,e2,e3
        double precision sqrtx,sqrty,sqrtz,xt,yt,zt
        if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,z).gt.BIG)pause 'invalid arguments in rf'
        xt=x
        yt=y
        zt=z
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=THIRD*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
        e2=delx*dely-delz**2
        e3=delx*dely*delz
        rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
        return
    END FUNCTION rf
