module sim_para
    implicit none
    integer:: zzz=-4280145, zzzz=77777
     integer,parameter:: lx=20,ly=20,lz=20,n_part=600
     integer,parameter::kb_T=1, niter=1000,mass=1,epsilon=1

     real*8,parameter:: rc=2.5d0, rs=2.0d0*rc, sigma=1.0d0
     real*8,parameter:: sigma6=sigma**6, sigma12=sigma**12 ,eps=4.0d0*epsilon
     real*8,parameter:: fc=eps*((12.0d0*sigma12/(rc**13))-(6.0d0*sigma6/(rc**7)))
     real*8,parameter:: ufc=eps*((sigma12/(rc**12))-(sigma6/(rc**6)))

     real*8::pos(3*n_part),vel(3*n_part),force(3*n_part),acc(3*n_part),new_acc(3*n_part)
     real*8::avr_vel_x,avr_vel_y,avr_vel_z
     integer:: i,j,k,c
     real*8::x1,y1,z1,x2,y2,z2,x,y,z,dx,dy,dz,r,lj,pot_energy,ke,gx,gy,gz,delta_t=0.005d0
     real*8:: lj_force
     real*8::new_pot_energy,new_force(3*n_part)
     real*8::invr, ir2, ir6
     real*8::llx,lly,llz,theoryke,scalef
end module sim_para