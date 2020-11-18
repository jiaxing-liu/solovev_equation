program trans
  implicit none
  integer,parameter::n=100
  integer::i,j,k
  real*8::r(n,n),z(n,n),psi(n,n)
  open(1,file='psitoso.txt')
  do i=1,n
     do j=1,n
        read(1,*)r(i,j),z(i,j),psi(i,j)
     enddo
  enddo
  close(1)
  open(4,file='psitosog.txt')
  do i=1,n
     do j=1,n
        if(psi(i,j)>0.D0.and.psi(i,j)<1.D-3)then
           write(4,*)r(i,j),z(i,j),(psi(i,j))
        else if(psi(i,j)>=1.D-3)then
           write(4,*)r(i,j),z(i,j),1.D-3
        else
           write(4,*)r(i,j),z(i,j),1.D-6
        endif
     enddo
  enddo
  close(4)
  open(2,file='psirz.txt')
  do i=1,n
        write(2,*)r(i,1:n)
  enddo
  do i=1,n
        write(2,*)z(i,1:n)
  enddo
  do i=1,n
        write(2,*)psi(i,1:n)
  enddo
  close(2)

  open(11,file='psito.txt')
  do i=1,n
     do j=1,n
        read(11,*)r(i,j),z(i,j),psi(i,j)
     enddo
  enddo
  close(11)
 ! open(44,file='psitog_my.txt')
  !do i=1,n
   !  do j=1,n
!        if(psi(i,j)>0.D0.and.psi(i,j)<1.D-3)then
 !          write(44,*)r(i,j),z(i,j),(psi(i,j))
  !      else if(psi(i,j)>=1.D-3)then
   !        write(4,*)r(i,j),z(i,j),1.D-3
    !    else
     !      write(44,*)r(i,j),z(i,j),1.D-6
      !  endif
!     enddo
 ! enddo
 ! close(44)
  open(22,file='psirz_my.txt')
  do i=1,n
        write(22,*)r(i,1:n)
  enddo
  do i=1,n
        write(22,*)z(i,1:n)
  enddo
  do i=1,n
        write(22,*)psi(i,1:n)
  enddo
  close(22)
endprogram trans
