program solovev
      implicit none
      integer,parameter::Nh=10
      integer::i,j
      real*8::A(Nh,Nh),b(Nh),x(Nh)
      open(1,file="dat1")
      do i=1,Nh
      read(1,*)A(i,1:Nh),b(i)
      end do
      close(1)
      call doolittle(A,b,x)
      print*,A
      print*,b
      print*,x
      print*,matmul(A,x)-b
!      call Gauss_Seidel(A,b,x)
!      print*,x
!      print*,matmul(A,x)-b
end program solovev

subroutine Gauss_Seidel(A,b,x)
  implicit none
  integer,parameter::Nh=10
  real*8::A(Nh,Nh),b(Nh),x(Nh),temp2,sum1,temp_err(Nh),err,tol,temp1(Nh),temp3(Nh)
  integer::i,j,k,order(Nh)
  !  x(:)=1.D0;
  tol=1.D-15
  err=2*tol
  do i=1,Nh
     do j=i,Nh-1
        temp2=max(A(i,j),A(i,j+1))
     enddo
     do j=1,Nh
        if(A(i,j)==temp2.and.j.ne.i)then
           temp1(1:Nh)=A(1:Nh,i)
           A(1:Nh,i)=A(1:Nh,j)
           A(1:Nh,j)=temp1(1:Nh)
           order(j)=i;order(i)=j
        endif
     enddo
  enddo
  
  do while(err>tol)
     do i=1,Nh
        sum1=0
        do j=1,Nh
           if(j.ne.i)then
              sum1=sum1+A(i,j)*x(j)
           endif
        enddo
        temp1(i)=-(sum1-b(i))/A(i,i)
        temp_err(i)=abs(x(i)-temp1(i))
        x(i)=temp1(i)
     enddo
     do i=1,Nh-1
        err=max(temp_err(i),temp_err(i+1))
     enddo
     print*,err,tol
  enddo
  do i=1,Nh
     temp3(i)=x(order(i))
  enddo
  x(:)=temp3(:)
end subroutine Gauss_Seidel

subroutine doolittle(A,b,x)
      implicit none
      integer,parameter::Nh=10
      real*8::A(Nh,Nh),b(Nh),x(Nh),y(Nh),L(Nh,Nh),U(Nh,Nh),sum1
      integer::i,j,k,r
      u(1,:)=A(1,:);
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
