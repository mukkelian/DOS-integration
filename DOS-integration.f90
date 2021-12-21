! Density of State (DOS) integration program
! Copyright (C) 2021  Mukesh Kumar Sharma, Department of Physics,
! Indian Institute of Technology Roorkee, Uttrakhand, India, PIN code 247667

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

!       UPDATE :: Mar 12 2021 10:14 PM
!       integeration based on simpson's method
	implicit none
	real, allocatable :: x(:),y(:,:)
	real::u,l,ans,n1,h, Emax, Emin
	integer:: i,j,n, var, cl, n_cl, e_n, e_m, e_max, e_min
	character :: opt, tag(150), choice
	open(1,file='input.dat',status='unknown')
	open(2,file='output.dat',status='unknown')
	print*, "Only numerical datas are required, if there's string then &
        kindly remove it"
        print*,''
	print*,'Enter the no of coloumns in the data file'
	read*,n_cl
        call system('rm -rf file')
        open(11,file ='file', status = 'unknown' )
        call system ("echo $(wc -l < input.dat) >> file")
        read(11,*) var
	close(11)
        call system('rm -rf file')

	allocate(x(0:var-1), y(0:var-1,2:n_cl))
	x = 0; y = 0
        print"(' found ',i6,' lines')",var
	do i= 0,var-1
	read(1,*)x(i), (y(i,j),j =2, n_cl) 
	end do
        close (1)

101	print*,'Coloumn no. to integrate: LOOK into the input.dat file'
	read*,cl
        print*,''
        
	h = abs(x(0) - x(1))
	
	tag = '#'
	write(*,11) tag

100	print*,''
	print*,''
102	print*, 'Do you want to do for whole data set (y)'
	print*,' or for data range (n)'
	read*, opt

	if (opt.eq.'n')then

	print*, 'Here are the points on x-axis'
	write(*,10) x
10	format (12f12.8,/)
	write(*,11) tag

	print*, 'Enter the value of Xmax, Xmin!'

	print*, 'Xmin (lower limit)'
	read*, Emin
	print*, 'Xmax (upper limit)'
	read*, Emax

	do i = 0, var-1

	e_min = x(i)/Emin
        if(abs(x(i)-Emin).le.0.00001)then
!	if (e_min == 1)then
        print*, i,Emin,'(line no., X)'
	e_n = i
	end if

	e_max = x(i)/Emax
        if (abs(x(i)-Emax).le.0.00001)then
!	if (e_max == 1)then
        print*, i,Emax,'(line no., X)'
	e_m = i
	end if

	end do

	n = e_m+1-e_n

        elseif(opt.eq.'y')then

	n = var-1
	e_n = 0; e_m = var-1
	
	Emin = x(e_m); Emax =x(e_n)
	 
        else
        print*,''
        print*,'please choose the option wisely!'
        print*,'````````````````````````````````'
        print*,''
        goto 102
	end if

	write(2,3) n
3	format('Number of intervals is',1x,i4,/)
	write(2,4)x(e_m),x(e_n)
4	format('upper limit = ',f10.5,5x,'lower limit = ',f10.5,/)

	if(mod(n,2)==0.0)then
	call simpson_1_3 (h,n,ans,var,cl,n_cl,e_n,e_m)
	write(2,*)'value of integral is =', ans
	print*,''
	print*,''
	write(*,*)'value of INTEGRAL is =', ans
	else
	call simpson_3_8 (h,n,ans,var,cl,n_cl,e_n,e_m)
	write(2,*)'value of integral is =', ans
	print*,''
	print*,''
	write(*,*)'value of INTEGRAL is =', ans
	end if

        print*,''
        print*,''
        print*, 'for coloumn no.',cl
        print*, 'for the X range',Emin,'to', Emax
	print*,''
	print*,''
	print*, 'May the Force be with you!'
	print*,'~ Mukesh Kumar Sharma'
	print*,'e-mail :: msharma1@ph.iitr.ac.in'
	print*,''
	print*,''
11	format (150A1)
        print*,''
        print*,''
        print*,'Do you want to continue for the same coloumn'
        print*,"press 'y' for same coloumn or,"
        print*,"press 'n' for other coloumn"
        print*,'press any other key to exit'
        read*, choice
        choice_ : if(choice .eq.'y')then
                goto 100
        elseif(choice .eq.'n')then
                goto 101
        else
              print*,"Kindly see the 'output.dat', output file"
              print*,''
              print*,'Thank you!'
              print*,''
              exit choice_
        endif choice_

	stop
	end


	subroutine simpson_1_3 (h,n,ans,var,cl,n_cl,e_n,e_m)
	integer i,j,var, cl, n_cl, e_n, e_m
	real, allocatable :: x_(:),y_(:,:),x_1(:),y_1(:,:)
	real:: h,ans,x1,y_odd,y_even
	allocate(x_(0:var-1), y_(0:var-1,2:n_cl))
        allocate(x_1(0:var-1), y_1(0:var-1,2:n_cl))
	x_ = 0; y_ = 0
	print*, "Entering into the Simpson's 1/3 rd rule"
        open(1,file='input.dat',status='unknown')

        do i=0,var-1
        read(1,*)x_1(i), (y_1(i,j),j =2, n_cl)
        end do
        close (1)


	do i = e_n, e_m
	x_(i) = x_1(i); y_(i,cl) = y_1(i,cl)
	end do

	y_odd=0.
       	do i = e_n+1, e_m-1,2
	y_odd = y_(i,cl) + y_odd
	end do

	y_even=0.
	do i=e_n+2,e_m-2,2
	y_even = y_even + y_(i,cl)
	end do

	ans=(h/3.)*(y_(e_n,cl)+y_(e_m,cl)+(4*y_odd)+(2*y_even))
        deallocate(x_, y_, x_1, y_1)
	end subroutine simpson_1_3


	subroutine simpson_3_8 (h,n,ans,var,cl,n_cl,e_n,e_m)
	integer i, j, n, var, cl, n_cl, e_n, e_m
        real:: h,ans,x1,y_odd,y_even,y1,y2,y3
	real, allocatable :: x_(:),y_(:,:), x_1(:),y_1(:,:)
	allocate(x_(e_n:e_m), y_(e_n:e_m,2:n_cl))
        allocate(x_1(0:var-1), y_1(0:var-1,2:n_cl))

	print*, "Entering into the Simpson's 3/8 th rule"
	x_ = 0; y_ = 0; x_1 = 0; y_1 = 0

        open(1,file='input.dat',status='unknown')

        do i=0,var-1
        read(1,*)x_1(i), (y_1(i,j),j =2, n_cl)
        end do
        close (1)

	do i = e_n, e_m
	x_(i) = x_1(i); y_(i,cl) = y_1(i,cl)
	end do

	y1=0.
	y2=0.
	y3=0.

	do i=e_n+3,e_m-3,3!1,n-2,3
	y1=y1+y_(i,cl)
	end do

	do i=e_n+1,e_m-2,3!2,n-1,3
        y2=y2+y_(i,cl)
        end do

	do i=e_n+2,e_m-1,3!3,n-3,3
        y3=y3+y_(i,cl)
        end do

	ans=3.*(h/8.)*(y_(e_n,cl)+y_(e_m,cl)+(2*y1)+(3*y2)+(3*y3))
        deallocate(x_, y_, x_1, y_1)
	end subroutine simpson_3_8

