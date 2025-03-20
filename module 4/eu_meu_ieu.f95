module variables
    implicit none
 integer:: i,niter
 real*8:: y,x,Ymid,yend,f1,fend,f0,f
real*8,parameter::dx=0.001
    
end module variables

program name
    use variables
    implicit none
    open(10,file="euler.dat", status="unknown")
    open(20,file="modified_euler.dat", status="unknown")
    open(30,file="improved_euler.dat", status="unknown")
    call euler
    call modified_euler
    call  imp_euler
end program name


subroutine euler
    use variables
    implicit none
    y=0.0d0
    x=0.0d0
   
    niter=1550
    do i= 1,niter

    y=y+dx*(y**2+1)
    x=x+dx

    write(10,*)  y,x
    IF (i==niter) then
        print*,48.078-y
    end if
    end do
end subroutine

subroutine modified_euler
    use variables
    implicit none
    
    y=0.0d0
    x=0.0d0
    Ymid=0.0d0
    niter=1550
    
    do i= 1,niter
        Ymid=y+(dx/2.0d0)*((y**2)+1)
         
        y= y+ dx*((Ymid**2)+1)
         x=x+dx
        write(20,*)  y , x
        IF (i==niter) then
            print*,48.078-y
        end if
        end do
    end subroutine

subroutine imp_euler
    use variables
    implicit none
    y=0.0d0
    x=0.0d0
    Yend=0.0d0
    f0=0.0d0
    niter=1550
    do i= 1,niter
        f0=(y**2)+1.0d0
        Yend=y+dx*f0
        fend=((yend**2) +1.0d0)
        f=(f0+fend)/(2.0d0)
        y= y+ dx*(f)

         x=x+dx
        write(30,*)  y , x
        IF (i==niter) then
            print*,48.078-y
        end if
        end do
    end subroutine