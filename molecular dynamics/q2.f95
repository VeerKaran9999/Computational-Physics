program main_md_program
    use sim_para
    implicit none

   llx=dfloat(lx) ;lly=dfloat(ly) ;llz=dfloat(lz)
    
   call intialize_postion
   call velocity_ani_force_intialization
    
   c=0

open(90,file="r.dat", status="unknown")
open(82,file="energy.dat", status="unknown",form="formatted")
open(10,file="position.dat",status="unknown",form="formatted")
 
do k= 1,niter

  ! if(mod(k,100)==0) then
    !write(*,*) k

    call postion_velocity_update
    write(82,*) dfloat(k),(new_pot_energy)/dfloat(n_part), ke/dfloat(n_part),(ke+new_pot_energy)/dfloat(n_part)
   ! end if
end do

end program main_md_program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine intialize_postion
    use sim_para
    implicit none
    real :: rx, ry, rz
    integer :: accept
    
    call random_seed()
    pos=0.0d0
    do i = 1, n_part
        accept = 0
        do while (accept == 0)
            call random_number(rx)
            call random_number(ry)
            call random_number(rz)
            
            rx = rx * lx
            ry = ry * ly
            rz = rz * lz
            
            accept = 1
            do j = 1, i-1
                r = dsqrt((rx - pos(3*j-2))**2 + (ry - pos(3*j-1))**2 + (rz - pos(3*j))**2)
                if (r <1.2* sigma) then
                    accept = 0
                    exit
                end if
            end do
        end do
        
        pos(3*i-2) = rx
        pos(3*i-1) = ry
        pos(3*i) = rz
        
    end do
    
    end subroutine intialize_postion


    ! intialization velocity , energy and force

subroutine velocity_ani_force_intialization
    use sim_para
    implicit none

    real*8:: vel_const, avg_vx, avg_vy, avg_vz
    real*8::  ran_1,ran_2,ran_3
    
    vel_const= dsqrt(12.0d0)*dfloat(kb_T)

    !intialize velocity

    do i= 1,n_part
        call random_number(ran_1)
        call random_number(ran_2)
        call random_number(ran_3)
        vel(3*i-2)=vel_const*(ran_1-0.5d0)
        vel(3*i-1)=vel_const*(ran_2-0.5d0)
        vel(3*i)=vel_const*(ran_3-0.5d0)

    end do
    !avg_vx=0.0d0 ;avg_vy=0.0d0 ;avg_vz=0.0d0
    !do i=1,n_part
     !   avg_vx=vel(3*i-2)+ avg_vx 
      !  avg_vy=vel(3*i-1)+ avg_vy
      ! avg_vz=vel(3*i)+ avg_vz
    !end do
    !avg_vx = avg_vx/(dfloat(n_part))
    !avg_vy = avg_vy/(dfloat(n_part)) 
    !avg_vz = avg_vz/(dfloat(n_part)) 

    !write (*,*) avg_vx,avg_vy,avg_vz

   ! do i=1,n_part                                  !this modification done by karan
      !  vel(3*i-2)=avg_vx-vel(3*i-2)
       ! vel(3*i-1)=avg_vy-vel(3*i-1)
        !vel(3*i)=avg_vz-vel(3*i)
    !end do

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !intialize force

    force= 0.0d0 ; pot_energy= 0.0d0 ; lj=0.0d0

    !acc= force/dfloat(mass)

    do i=1, n_part-1
        x1=pos(3*i-2) ; y1=pos(3*i-1) ;  z1=pos(3*i)

    do j=i+1, n_part                                       !differance find kela ahe postion madhala
      !  c=c+1
        x2=pos(3*j-2) ; y2=pos(3*j-1) ;z2=pos(3*j)
         
        x=x1-x2 ; y=y1-y2 ; z=z1-z2  
     
        ! implamation of minimum image convention
        if(abs(x).gt.(dfloat(lx/2)))  x=(dfloat(lx)-abs(x))*((-1.0d0*x)/abs(x))
        if(abs(y).gt.(dfloat(ly/2)))   y=(dfloat(ly)-abs(y))*((-1.0d0*y)/abs(y))
        if(abs(z).gt.(dfloat(lz/2)))   z=(dfloat(lz)-abs(z))*((-1.0d0*z)/abs(z))

        r=dsqrt(x*x+y*y+z*z)
         
        !atta apan only rc yevdya rance varil particle la force calculation sathi genar ahe 
         
        if (r<=rc) then

            lj=eps*(((sigma/r)**12)-(sigma/r)**6)-ufc+fc*r !modified lenard jones potential

            pot_energy=pot_energy+lj
            print*,pot_energy
            lj_force=eps*(12.0d0*(sigma12/(r)**13)-6*(sigma6/(r)**7))-fc

          !  write(90,*) new_force(3*i-2),new_force(3*i-1),new_force(3*i)

            force(3*i-2) =force(3*i-2) + (lj_force)*(x/r)
            force(3*i-1) =force(3*i-1) + (lj_force)*(y/r)   !eka particle var bakichya particle mule lagnara force ani (x/r ne multiply kel ahe becouse direction sathi)
            force(3*i)   =force(3*i)   + (lj_force)*(z/r)
            force(3*j-2) =force(3*j-2) - (lj_force)*(x/r)   
            force(3*j-1) =force(3*j-1) - (lj_force)*(y/r)   ! equal and opposite force lagele dusarya particle var
            force(3*j)   =force(3*j)   - (lj_force)*(z/r)

        end if   
    end do
    end do

 end subroutine velocity_ani_force_intialization




 ! postion ,force and velocity update

subroutine postion_velocity_update
    use sim_para
    implicit none
    
    real*8::  dt2by2 ,dt_by_2

    ! (1)-updating postion

    dt2by2 =0.50d0*delta_t**2
write(10,*) n_part
write(10,*) 
    do i=1, n_part
        pos(3*i-2)= pos(3*i-2)+ vel(3*i-2)*delta_t +dt2by2*force(3*i-2)     !x componet of postion
        pos(3*i-1)=pos(3*i-1)+ vel(3*i-1)*delta_t +dt2by2*force(3*i-1)      !y componet of postion
        pos(3*i)=pos(3*i)+ vel(3*i)*delta_t +dt2by2*force(3*i)              !z componet of postion
         
        pos(3*i-2) = modulo(pos(3*i-2),llx)
        pos(3*i-1) = modulo(pos(3*i-1),lly)    ! periodic boundary condition 
        pos(3*i) = modulo(pos(3*i),llz)
write(10,*) pos(3*i-2), pos(3*i-1),  pos(3*i)
    end do

    !(2)-calculation of force

    new_force= 0.0d0 ; new_pot_energy= 0.0d0

   ! acc= force/dfloat(mass)

    do i=1, n_part-1
        x1=pos(3*i-2) ; y1=pos(3*i-1) ;  z1=pos(3*i)

    do j=i+1, n_part                                       !differce find kela ahe postion madhala
        c=c+1
        x2=pos(3*j-2) ; y2=pos(3*j-1) ;z2=pos(3*j)
         
        x=x1-x2 ; y=y1-y2 ; z=z1-z2  
     
        ! implamation of minimum image convention
        if(abs(x).gt.(dfloat(lx/2)))  x=(dfloat(lx)-abs(x))*((-1.0d0*x)/abs(x))
        if(abs(y).gt.(dfloat(ly/2)))   y=(dfloat(ly)-abs(y))*((-1.0d0*y)/abs(y))
        if(abs(z).gt.(dfloat(lz/2)))   z=(dfloat(lz)-abs(z))*((-1.0d0*z)/abs(z))

        r=dsqrt(x*x+y*y+z*z)
         
        !atta apan only rc yevdya range varil particle la force calculation sathi genar ahe 
         
        if (r<=rc) then

            lj=eps*(((sigma/r)**12)-(sigma/r)**6)-ufc+fc*r !modified lenard jones potential

            new_pot_energy=new_pot_energy+lj
            lj_force=eps*(12.0d0*(sigma12/(r)**13)-6*(sigma6/(r)**7))-fc

            write(90,*) new_force(3*i-2),new_force(3*i-1),new_force(3*i)

            new_force(3*i-2) =new_force(3*i-2) + (lj_force)*(x/r)
            new_force(3*i-1) =new_force(3*i-1) + (lj_force)*(y/r)   !eka particle var bakichya particle mule lagnara force ani (x/r ne multiply kel ahe becouse direction sathi 
            new_force(3*i)   =new_force(3*i)   + (lj_force)*(z/r)
            new_force(3*j-2) =new_force(3*j-2) - (lj_force)*(x/r)   
            new_force(3*j-1) =new_force(3*j-1) - (lj_force)*(y/r)   ! equal and opposite force lagele dusarya particle var
            new_force(3*j)   =new_force(3*j)   - (lj_force)*(z/r)

        end if   
    end do
 end do
  !new_acc=new_force/(dfloat(mass))


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

   ! avr_vel_x= avr_vel_x + vel(3*i-2)
    !avr_vel_y= avr_vel_y + vel(3*i-1)
    !avr_vel_z= avr_vel_z + vel(3*i)

    ke= ke + vel(3*i-2)*vel(3*i-2)
    ke= ke + vel(3*i-1)*vel(3*i-1)
    ke= ke + vel(3*i)*vel(3*i)
 
    force(3*i-2) = new_force(3*i-2) 
    force(3*i-1) =new_force(3*i-1)  
    force(3*i)   =new_force(3*i)              ! next itteration sathi old force ha equat karto apan new force shi
      
 end do
 ke= 0.50d0*dfloat(mass)*ke
 
 end subroutine postion_velocity_update



