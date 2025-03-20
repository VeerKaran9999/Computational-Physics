program name
    implicit none
    integer:: i,niter,T=50
    real*8::h,fx1,fx2,fx3,fx4,fv1,fv2,fv3,fv4,x,v,time,xtemp1,xtemp2,xtemp3,vtemp1,vtemp2,vtemp3

    h=0.01
    niter=5000
    x=0.1d0
    v=1.9d0
    time=0.0d0
open(10, file="q5outputdata.dat" ,status="unknown")
    do i= 1, niter
        !position calculation
        fx1=v
        xtemp1=x+(h/2.0d0)*fx1
        fv1=-sin(x)
        vtemp1=v+(h/2.0d0)*fv1
       
        fx2=vtemp1
        xtemp2=x+(h/2.0d0)*fx2
        fv2=-sin(xtemp1)
        vtemp2=v+(h/2.0d0)*fv2
        
        fx3=vtemp2
        xtemp3=x+(h)*fx3
        fv3=-sin(xtemp2)
        vtemp3=v+(h)*fv3
 
        fx4=vtemp3
        fv4=-sin(xtemp3)

        v=v+(h/6.0d0)*(fv1+2*fv2+2*fv3+fv4) ; x=x+(h/6.0d0)*(fx1+2*fx2+2*fx3+fx4)

        time=time+h
write(10,*) time ,x     ,v
  end do
  print*,"your file is ready"
end program name