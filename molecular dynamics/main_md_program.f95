program main_md_program
    use sim_para
    implicit none
   real*8:: ran1,zz1

   zz1=ran1(zzz)
   llx=dfloat(lx) ;lly=dfloat(ly) ;llz=dfloat(lz)
    
   call intialize_postion
   call velocity_ani_force_intialization
    
   c=0

open(90,file="dist.dat", status="unknown")
open(82,file="energy.dat", status="unknown",form="formatted")

 do k= 1,niter

    if(mod(k,100)==0) write(*,*) k

    call postion_velocity_update
    write(82,"(4g25.15)") dfloat(k),(new_pot_energy)/dfloat(n_part), ke/dfloat(n_part),(ke+new_pot_energy)/dfloat(n_part)
  
end do

end program main_md_program


! postion ,force and velocity update

subroutine postion_velocity_update
    use sim_para
    implicit none
    
    real*8:: zz1,ran1
    real*8::  dt2by2 ,dt_by_2

    ! (1)-updating postion

    dt2by2 =0.50d0*delta_t**2

    do i=1, n_part
        pos(3*i-2)= pos(3*i-2)+ vel(3*i-2)*delta_t +dt2by2*force(3*i-2)     !x componet of postion
        pos(3*i-1)=pos(3*i-1)+ vel(3*i-1)*delta_t +dt2by2*force(3*i-1)      !y componet of postion
        pos(3*i)=pos(3*i)+ vel(3*i)*delta_t +dt2by2*force(3*i)              !z componet of postion
         
        pos(3*i-2) = modulo(pos(3*i-2),llx)
        pos(3*i-1) = modulo(pos(3*i-1),lly)    ! periodic boundary condition 
        pos(3*i) = modulo(pos(3*i),lly)

    end do

    !(2)-calculation of force

    new_force= 0.0d0 ; new_pot_energy= 0.0d0

    acc= force/dfloat(mass)

    do i=1, n_part-1
        x1=pos(3*1-2) ; y1=pos(3*1-1) ;  z1=pos(3*1)

    do j=i+1, n_part                                       !differce find kela ahe postion madhala
        c=c+1
        x2=pos(3*j-2) ; y1=pos(3*j-1) ;z2=pos(3*j)
         
        x=x1-x2 ; y=y1-y2 ; z=z1-z2  
     
        ! implamation of minimum image convention
        if(abs(x).gt.(dfloat(lx/2)))  x=(dfloat(lx)-abs(x))*((-1.0d0*x)/abs(x))
        if(abs(y).gt.(dfloat(ly/2)))   y=(dfloat(ly)-abs(y))*((-1.0d0*y)/abs(y))
        if(abs(z).gt.(dfloat(lz/2)))   z=(dfloat(lz)-abs(z))*((-1.0d0*z)/abs(z))

        r=dsqrt(x*x+y*y+z*z)
         
        !atta apan only rc yevdya distance varil particle la force calculation sathi genar ahe 
         
        if (r<=rc) then

            lj=eps*((sigma/r)**12-(sigma/r)**6)-ufc+fc*r !modified lenard jones potential

            new_pot_energy=new_pot_energy+lj
            lj_force=eps*(12.0d0*(sigma12/(r)**13)-6*(sigma6/(r)**7))-fc

            write(90,*) new_force(3*i-2),new_force(3*i-1),new_force(3*i)

            new_force(3*i-2) =new_force(3*i-2) + (lj_force)*(x/r)
            new_force(3*i-1) =new_force(3*i-1) + (lj_force)*(y/r)   !eka particle var bakichya particle mule lagnara force 
            new_force(3*i)   =new_force(3*i)   + (lj_force)*(z/r)
            new_force(3*j-2) =new_force(3*j-2) - (lj_force)*(x/r)   
            new_force(3*j-1) =new_force(3*j-1) - (lj_force)*(y/r)   ! equal and opposite force lagele dusarya particle var
            new_force(3*j)   =new_force(3*j)   - (lj_force)*(z/r)

        end if   
    end do
end do
new_acc=new_force/(dfloat(mass))


!(3)- updating velocity

ke = 0.0d0
avr_vel_x = 0.0d0 
avr_vel_y = 0.0d0 
avr_vel_z= 0.0d0 


dt_by_2 = 0.50d0*delta_t
do i= 1, n_part
 
    vel(3*i-2)=vel(3*i-2)+(dt_by_2*(force(3*i-2)+new_force(3*i-2)))
    vel(3*i-1)=vel(3*i-1)+(dt_by_2*(force(3*i-1)+new_force(3*i-1)))
    vel(3*i)=vel(3*i)+(dt_by_2*(force(3*i)+new_force(3*i)))

    avr_vel_x= avr_vel_x + vel(3*i-2)
    avr_vel_y= avr_vel_y + vel(3*i-1)
    avr_vel_z= avr_vel_z + vel(3*i)

    ke= ke + vel(3*i-2)*vel(3*i-2)
    ke= ke + vel(3*i-1)*vel(3*i-1)
    ke= ke + vel(3*i)*vel(3*i)

end do
ke= 0.50d0*dfloat(mass)*ke
 
end subroutine postion_velocity_update



! intialization velocity , energy and force

subroutine velocity_ani_force_intialization
    use sim_para
    implicit none

    real*8:: vel_const, avg_vx, avg_vy, avg_vz
    real*8:: zz1, ran1
    zz1 =ran1(zzz)
    vel_const=dsqrt(12.0d0)*dfloat(kb_T)

    !intialize velocity

    do i= 1,n_part
        vel(3*i-2)=vel_const*(ran1(zzzz)-0.5d0)
        vel(3*i-1)=vel_const*(ran1(zzzz)-0.5d0)
        vel(3*i)=vel_const*(ran1(zzzz)-0.5d0)
    end do
avg_vx=0.0d0 ;avg_vy=0.0d0 ;avg_vz=0.0d0
    do i=1,n_part
        avg_vx=vel(3*i-2)+ avg_vx 
        avg_vy=vel(3*i-1)+ avg_vy
        avg_vz=vel(3*i)+ avg_vz
    end do
    avg_vx = avg_vx/dfloat(n_part) 
    avg_vy = avg_vy/dfloat(n_part) 
    avg_vz = avg_vz/dfloat(n_part) 

    !write (*,*) avg_vx,avg_vy,avg_vz

    do i=1,n_part
        vel(3*i-2)=avg_vx-vel(3*i-2)
        vel(3*i-1)=avg_vx-vel(3*i-1)
        vel(3*i)=avg_vx-vel(3*i)
    end do

   ! avg_vx=0.0d0 ;avg_vy=0.0d0 ;avg_vz=0.0d0
    !do i=1,n_part
     !   avg_vx=vel(3*i-2)+ avg_vx 
      !  avg_vy=vel(3*i-1)+ avg_vy
       ! avg_vz=vel(3*i)+ avg_vz
    !end do
    !avg_vx = avg_vx/dfloat(n_part) 
    !avg_vy = avg_vy/dfloat(n_part) 
    !avg_vz = avg_vz/dfloat(n_part) 

    !write (*,*) avg_vx,avg_vy,avg_vz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    !intialize force

    force= 0.0d0 ; new_pot_energy= 0.0d0

    acc= force/dfloat(mass)

    do i=1, n_part-1
        x1=pos(3*1-2) ; y1=pos(3*1-1) ;  z1=pos(3*1)

    do j=i+1, n_part                                       !differance find kela ahe postion madhala
        c=c+1
        x2=pos(3*j-2) ; y1=pos(3*j-1) ;z2=pos(3*j)
         
        x=x1-x2 ; y=y1-y2 ; z=z1-z2  
     
        ! implamation of minimum image convention
        if(abs(x).gt.(dfloat(lx)))  x=(dfloat(lx)-abs(x))*((-1.0d0*x)/abs(x))
        if(abs(y).gt.(dfloat(ly)))   y=(dfloat(ly)-abs(y))*((-1.0d0*y)/abs(y))
        if(abs(z).gt.(dfloat(lz)))   z=(dfloat(lz)-abs(z))*((-1.0d0*z)/abs(z))

        r=dsqrt(x*x+y*y+z*z)
         
        !atta apan only rc yevdya distance varil particle la force calculation sathi genar ahe 
         
        if (r<=rc) then

            lj=eps*((sigma/r)**12-(sigma/r)**6)-ufc+fc*r !modified lenard jones potential

            new_pot_energy=new_pot_energy+lj
            lj_force=eps*(12.0d0*(sigma12/(r)**13)-6*(sigma6/(r)**7))-fc

            write(90,*) new_force(3*i-2),new_force(3*i-1),new_force(3*i)

            new_force(3*i-2) =new_force(3*i-2) + (lj_force)*(x/r)
            new_force(3*i-1) =new_force(3*i-1) + (lj_force)*(y/r)   !eka particle var bakichya particle mule lagnara force 
            new_force(3*i)   =new_force(3*i)   + (lj_force)*(z/r)
            new_force(3*j-2) =new_force(3*j-2) - (lj_force)*(x/r)   
            new_force(3*j-1) =new_force(3*j-1) - (lj_force)*(y/r)   ! equal and opposite force lagele dusarya particle var
            new_force(3*j)   =new_force(3*j)   - (lj_force)*(z/r)

        end if   
    end do
end do