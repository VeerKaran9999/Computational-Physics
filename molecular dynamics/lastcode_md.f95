module parameter
    implicit none
    integer :: zzz=4280145, zzzz=77777, i, j, k, c
    integer, parameter :: npart = 2197, niter=20000
    integer, parameter :: lx=20, ly=20, lz=20
    real*8, parameter :: mass=1.0d0, kb_T=1.0d0, rc=2.50d0, rs=-2.0d0*rc, sigma=1.0d0, eps=4.0d0, delta_t=0.005d0
    real*8, parameter :: sigma6=sigma**6, sigma12=sigma**12
    real*8, parameter :: fc=eps*((12.0d0*sigma12/(rc**13))-(6.0d0*sigma6/(rc**7)))
    real*8, parameter :: ufc=fc*rc + eps*(((sigma/rc)**12)-((sigma/rc)**6))
    real*8 :: avg_vx, avg_vy, avg_vz, PE, KE, gx, gy, gz
    real*8 :: x, y, z, x1, y1, z1, x2, y2, z2, dx, dy, dz,fs, r, lj, lj_force
    real*8 :: invr, ir2, ir6, theoryKE, scalef
    real*8, parameter :: llx=dfloat(lx), lly=dfloat(ly), llz=dfloat(lz)
    real*8, parameter :: llxby2=llx/2.0d0, llyby2=lly/2.0d0, llzby2=llz/2.0d0
    real*8, dimension(3*npart) :: pos, vel, force, old_force, acc, old_acc
 end module parameter
 
 !==========MAIN PROGRAM==========
 program md
    use parameter
    implicit none
    integer :: time
    
    call position_init
    call velocity_init
    call calc_force
    
    open(14,file='force.dat',status='unknown',form='formatted')
    do i=1,npart
       write(14,*) force(3*i-2), force(3*i-1), force(3*i)
    end do
    close(14)  
      
    open(10,file='Positiondata.dat',status='unknown',form='formatted')
    open(11,file='PE+KE+TE.dat',status='unknown',form='formatted')
    
    
    do k=1,niter
       call update_pos
       call calc_force
       call update_vel
       !call thermostat
       if(mod(k,100)==0) then
       write(11,*) dfloat(k),PE/dfloat(npart),kE/dfloat(npart),(KE+PE)/dfloat(npart)
       end if
    end do

    close(11)
           
 end program md
 
 !==========INITIATION OF POSITION==========  
 
 subroutine position_init
    use parameter
    implicit none
    integer :: n
    integer,allocatable :: seed(:)
 
    call random_seed(size=n)
    allocate(seed(n))
    seed = zzz  
    call random_seed(put=seed)
 
    do i=1,npart
 62    call random_number(x)
       x= 0.5d0 + (x*(lx-1))
       call random_number(y)
       y= 0.5d0 + (y*(ly-1))
       call random_number(z)
       z= 0.5d0 + (z*(lz-1))
       do j= 1, i-1
          if (sqrt((x - pos(3*j-2))**2 + (y - pos(3*j-1))**2 + (z - pos(3*j))**2) < sigma) then
             goto 62
          end if
       end do
      
       pos(3*i-2) = x
       pos(3*i-1) = y
       pos(3*i) = z
    end do      
    open(1,file='initial_positions.dat',status='unknown')
    do i=1,npart
       write(1,*) pos(3*i-2), pos(3*i-1), pos(3*i)  
    end do
    close(1)
    deallocate(seed)
 end subroutine position_init
 
 !==========INITIATION OF VELOCITY==========
 
 subroutine velocity_init
    use parameter
    implicit none
    real*8 :: vel_const, vx, vy, vz
    integer :: n
    integer,allocatable :: seed(:)
 
    call random_seed(size=n)
    allocate(seed(n))
    seed = zzzz  
    call random_seed(put=seed)
    
    vel_const = dsqrt(12.0d0*kb_T/mass)
    
    do i= 1, npart
       call random_number(vx)
       call random_number(vy)
       call random_number(vz)
       vel(3*i-2) = vel_const*(vx-0.5d0)
       vel(3*i-1) = vel_const*(vy-0.5d0)
       vel(3*i) = vel_const*(vz-0.5d0)
    end do  
    
    avg_vx = 0.0d0; avg_vy = 0.0d0; avg_vz = 0.0d0
    do i = 1, npart
       avg_vx = vel(3*i-2)+avg_vx; avg_vy = vel(3*i-1)+avg_vy; avg_vz = vel(3*i)+avg_vz
    end do
    avg_vx = avg_vx/dfloat(npart); avg_vy = avg_vy/dfloat(npart); avg_vz = avg_vz/dfloat(npart)
    
    do i= 1, npart
       vel(3*i-2) = avg_vx - vel(3*i-2)
       vel(3*i-1) = avg_vy - vel(3*i-1)
       vel(3*i) = avg_vz - vel(3*i)
       KE = KE + (vel(3*i-2)*vel(3*i-2)+vel(3*i-1)*vel(3*i-1)+vel(3*i)*vel(3*i))
    end do
    KE = 0.50d0*mass*KE
    open(2,file='initial_velocities.dat',status='unknown')
    do i=1,npart
       write(2,*) vel(3*i-2), vel(3*i-1), vel(3*i)  
    end do
    close(2)
    deallocate(seed)
 end subroutine velocity_init
 
 !==========CALCULATION OF FORCE==========
 
 subroutine calc_force
    use parameter
    implicit none
    PE=0.0d0
    force=0.0d0
    lj_force=0.0d0
    
    do i=1,npart-1
       x1=pos(3*i-2); y1=pos(3*i-1); z1=pos(3*i)
       do j=i+1,npart
          x2=pos(3*j-2); y2=pos(3*j-1); z2=pos(3*j)
          dx=x1-x2; dy=y1-y2; dz=z1-z2
          
          if(abs(dx).ge.llxby2) dx=(llx-abs(dx))*((-1.0d0*dx)/abs(dx))
          if(abs(dy).ge.llyby2) dy=(lly-abs(dy))*((-1.0d0*dy)/abs(dy))
          if(abs(dz).ge.llzby2) dz=(llz-abs(dz))*((-1.0d0*dz)/abs(dz))
          r=dsqrt(dx*dx+dy*dy+dz*dz)
          
          if(r<=rc) then
             lj=eps*(((sigma/r)**12)-((sigma/r)**6))-ufc+fc*r
             PE=PE+lj
             lj_force=eps*((12.0d0*sigma12/(r**13))-(6.0d0*sigma6/(r**7)))-fc
             force(3*i-2)=force(3*i-2)+(lj_force*dx/r)
             force(3*i-1)=force(3*i-1)+(lj_force*dy/r)
             force(3*i)=force(3*i)+(lj_force*dz/r)
             force(3*j-2)=force(3*j-2)-(lj_force*dx/r)
             force(3*j-1)=force(3*j-1)-(lj_force*dy/r)
             force(3*j)=force(3*j)-(lj_force*dz/r)
          end if
       end do
    end do
    acc=force/mass  
 end subroutine calc_force
 
 !==========POSITION UPDATION==========
 
 subroutine update_pos
    use parameter
    implicit none
    real*8 :: dt2by2
    old_force=0.0d0
    
    dt2by2 = 0.50d0*delta_t**2    
    write(10,*) npart 
    write(10,*)
    do i=1,npart
       pos(3*i-2)=pos(3*i-2) + (vel(3*i-2)*delta_t) + (dt2by2*force(3*i-2))
       pos(3*i-1)=pos(3*i-1) + (vel(3*i-1)*delta_t) + (dt2by2*force(3*i-1))
       pos(3*i)=pos(3*i) + (vel(3*i)*delta_t) + (dt2by2*force(3*i))
      
       pos(3*i-2)=modulo(pos(3*i-2),llx)
       pos(3*i-1)=modulo(pos(3*i-1),lly)
       pos(3*i)=modulo(pos(3*i),llz)
      
       old_force(3*i-2)=force(3*i-2)
       old_force(3*i-1)=force(3*i-1)
       old_force(3*i)=force(3*i)  
       
       write(10,*) pos(3*i-2),pos(3*i-1),pos(3*i)
    end do
    old_acc=acc  
 end subroutine update_pos
 
 !==========VELOCITY UPDATION==========
 
 subroutine update_vel
    use parameter
    implicit none
    KE=0.0d0
    fs = 0.5d0*delta_t/mass
    
    do i=1,npart
       vel(3*i-2) = vel(3*i-2) +(fs*(old_force(3*i-2)+force(3*i-2)))
       vel(3*i-1) = vel(3*i-1) +(fs*(old_force(3*i-1)+force(3*i-1)))
       vel(3*i) = vel(3*i) +(fs*(old_force(3*i)+force(3*i)))
       KE = KE + (vel(3*i-2)*vel(3*i-2)+vel(3*i-1)*vel(3*i-1)+vel(3*i)*vel(3*i))
    end do
    KE = 0.50d0*mass*KE  
 end subroutine update_vel  
 
 !==========VELOCITY RESCALING(THERMOSTAT)=========
 
 subroutine thermostat
    use parameter
    implicit none
 
    if(mod(k,100)==0) then
       theoryKE=1.5d0*dfloat(npart)*kb_T
       scalef=dsqrt(theoryKE/KE)
       vel = vel*scalef
      
       KE = 0.0d0
       do i=1,npart
          KE = KE + (vel(3*i-2)*vel(3*i-2)+vel(3*i-1)*vel(3*i-1)+vel(3*i)*vel(3*i))
       end do  
       KE = 0.50d0*mass*KE
    end if
 end subroutine thermostat 