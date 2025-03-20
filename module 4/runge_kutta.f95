program name
    implicit none
    real*8:: f1,f2,f3,f4,dx,y,ytemp1,ytemp2,ytemp3,x
    integer:: i,niter
    dx=0.01
    niter=155
open(10,file="RK4.dat",status="unknown")
    do i=1,niter
        f1=y**2+1
        ytemp1=y+(dx/2.0d0)*f1
       
        f2=ytemp1**2+1
        ytemp2=y+(dx/2.0d0)*f2
        
        f3=ytemp2**2+1
        ytemp3=y+(dx)*f3
 
        f4=ytemp3**2+1

        y=y+(dx/6.0d0)*(f1+2*f2+2*f3+f4)

        x=x+dx
        
        write(10,*) x    ,    y

        IF (i==niter) then
            print*,48.078-y
        end if

    end do
    print*, "done"
end program name