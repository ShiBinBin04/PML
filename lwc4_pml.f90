module lwc4v
  implicit none
  integer,parameter :: rp=4!selected_real_kind(16)
  integer,parameter :: ip=selected_int_kind(6)   
  contains
  subroutine expand(c,np,nx,nz,nx_pml,nz_pml)
    implicit none
    real   (kind=rp),dimension(:,:)::c
    integer(kind=ip),intent(in)    ::nx,nz,np,nx_pml,nz_pml
    integer(kind=ip)               ::ix,iz
       do ix=np+nx_pml+1,np+nx_pml+nx;do iz=1,np+nz_pml
        c(ix,iz)=c(ix,np+nz_pml+1)
       enddo;enddo
       do ix=np+nx_pml+1,np+nx_pml+nx;do iz=np+nz_pml+nz+1,np+nz_pml+nz+nz_pml+np
       c(ix,iz)=c(ix,np+nz_pml+nz)
       enddo;enddo
       do ix=1,np+nx_pml;do iz=np+nz_pml+1,np+nz_pml+nz
       c(ix,iz)=c(np+nx_pml+1,iz)
       enddo;enddo
       do ix=np+nx_pml+nx+1,np+nx_pml+nx+nx_pml+np;do iz=np+nz_pml+1,np+nz_pml+nz
       c(ix,iz)=c(np+nx_pml+nx,iz)
       enddo;enddo
       do ix=1,np+nx_pml;do iz=1,np+nz_pml
       c(ix,iz)=c(np+nx_pml+1,np+nz_pml+1)
       enddo;enddo
       do ix=np+nx_pml+nx+1,np+nx_pml+nx+nx_pml+np;do iz=1,np+nz_pml
       c(ix,iz)=c(np+nx_pml+nx,np+nz_pml+1)
       enddo;enddo
       do ix=np+nx_pml+nx+1,np+nx_pml+nx+nx_pml+np;do iz=np+nz_pml+nz+1,np+nz_pml+nz+np+nz_pml
       c(ix,iz)=c(np+nx_pml+nx,np+nz_pml+nz)
       enddo;enddo
    do ix=1,np+nx_pml;do iz=np+nz_pml+nz+1,np+nz_pml+nz+np+nz_pml
    c(ix,iz)=c(np+nx_pml+1,np+nz_pml+nz)
    enddo;enddo
    ! end of expland
    end subroutine expand

    subroutine spmlt(np,nx_pml,nz_pml,nx,nz,u1,u,s1,s,s0,c,cdt2,dt,dx,dz)
      !this subroutine aims to perform pml for top boundary
      !np for half operator length
      !nx_pml, nz_pml, nx and nz denotes for
      !u1
      implicit none
      integer(kind=ip)                :: ix,iz,iix,iiz
      integer(kind=ip),intent(in)     :: np,nx_pml,nz_pml,nx,nz
      real(kind=rp),dimension(:,:)    :: u,u1,c,cdt2
      real(kind=rp),dimension(:,:,:)  :: s1,s,s0
      real(kind=rp)       :: temp,d_x,dd_x,d_z,dd_z,dt,dx,dz,u2x,u2z,ux,uz
      
      do ix=np+Nx_pml+1,np+nx_pml+Nx;do iz=np+1,np+nz_pml    
      temp=-1.50*log(0.0001)*c(ix,iz)
      d_z=temp*(iz-np-nz_pml-1)*(iz-np-nz_pml-1)/(nz_pml*nz_pml*nz_pml)/dz
      dd_z=d_z*2/(dz*(iz-np-nz_pml-1))
      iix=ix-np-nx_pml
      iiz=iz-np
      u2x=(-u(ix+2,iz)+16.0*u(ix+1,iz)-30.0*u(ix,iz)+16.0*u(ix-1,iz)-u(ix-2,iz))/12.0
      u2x=u2x/dx/dx
      u2z=(-u(ix,iz+2)+16.0*u(ix,iz+1)-30.0*u(ix,iz)+16.0*u(ix,iz-1)-u(ix,iz-2))/12.0
      u2z=u2z/dz/dz
      ux=(u(ix+1,iz)-u(ix-1,iz))/(2.0*dx)
      uz=(u(ix,iz+1)-u(ix,iz-1))/(2.0*dz)
      s1(1,iix,iiz)=(      cdt2(ix,iz)*u2z+(2.-d_z*d_z*dt*dt)*s(1,iix,iiz)+(d_z*dt-1)*s0(1,iix,iiz))/(d_z*dt+1)
      s1(4,iix,iiz)=(-dd_z*cdt2(ix,iz)*uz +(2.-d_z*d_z*dt*dt)*s(4,iix,iiz)+(d_z*dt-1)*s0(4,iix,iiz))/(d_z*dt+1)
      s1(3,iix,iiz)=cdt2(ix,iz)*u2x+2.*s(3,iix,iiz)-s0(3,iix,iiz)
      s1(2,iix,iiz)=s(4,iix,iiz)*dt+s(2,iix,iiz)*(1-d_z*dt)
      u1(ix,iz)=s1(1,iix,iiz)+s1(2,iix,iiz)+s1(3,iix,iiz)
      enddo;enddo
      s0=s 
      s =s1
    end subroutine spmlt    
    subroutine spmlb(np,nx_pml,nz_pml,nx,nz,u1,u,s1,s,s0,c,cdt2,dt,dx,dz)
      implicit none
      integer(kind=4)                :: ix,iz,iix,iiz,np,nx_pml,nz_pml,nx,nz
      real(kind=rp),dimension(:,:)    :: u,u1,c,cdt2
      real(kind=rp),dimension(:,:,:)  :: s1,s,s0
      real(kind=rp)       :: temp,d_z,dd_z,dt,dx,dz,u2x,u2z,ux,uz
      
      do ix=np+nx_pml+1,np+nx_pml+Nx;do iz=np+nz_pml+Nz+1,np+nz_pml+Nz+nz_pml
        temp=-1.50*log(0.0001)*c(ix,iz)
        d_z=temp*((iz-nz-nz_pml-np)*dz)**2/(nz_pml*dz)**3
        dd_z=d_z*2/(dz*(iz-nz-nz_pml-np))
        iix=ix   -Nx_pml-Np
        iiz=iz-Nz-Nz_pml-Np
        u2x=(-u(ix+2,iz)+16.0*u(ix+1,iz)-30.0*u(ix,iz)+16.0*u(ix-1,iz)-u(ix-2,iz))/12.0
        u2x=u2x/dx/dx
        u2z=(-u(ix,iz+2)+16.0*u(ix,iz+1)-30.0*u(ix,iz)+16.0*u(ix,iz-1)-u(ix,iz-2))/12.0
        u2z=u2z/dz/dz
        ux=(u(ix+1,iz)-u(ix-1,iz))/(2.0*dx)
        uz=(u(ix,iz+1)-u(ix,iz-1))/(2.0*dz)
        s1(1,iix,iiz)=(      cdt2(ix,iz)*u2z+(2.-d_z*d_z*dt*dt)*s(1,iix,iiz)&
                       &+(d_z*dt-1)*s0(1,iix,iiz))/(d_z*dt+1)
        s1(4,iix,iiz)=(-dd_z*cdt2(ix,iz)*uz +(2.-d_z*d_z*dt*dt)*s(4,iix,iiz)&
                       &+(d_z*dt-1)*s0(4,iix,iiz))/(d_z*dt+1)
        s1(3,iix,iiz)=       cdt2(ix,iz)*u2x+2.*s(3,iix,iiz)-s0(3,iix,iiz)
        s1(2,iix,iiz)=s(4,iix,iiz)*dt+s(2,iix,iiz)*(1-d_z*dt)
        u1(ix,iz)=s1(1,iix,iiz)+s1(2,iix,iiz)+s1(3,iix,iiz)
      enddo;enddo
      s0=s; s=s1
    end subroutine spmlb
    
    subroutine spmll(np,nx_pml,nz_pml,nx,nz,u1,u,r1,r,r0,c,cdt2,dt,dx,dz)
      implicit none
      integer(kind=4)                :: ix,iz,iix,iiz,np,nx_pml,nz_pml,nx,nz
      real(kind=rp),dimension(:,:)    :: u,u1,c,cdt2
      real(kind=rp),dimension(:,:,:)  :: r1,r,r0
      real(kind=rp)       :: temp,d_x,dd_x,dt,dx,dz,u2x,u2z,ux,uz
      
      do ix=np+1,np+nx_pml ;do iz=np+nz_pml+1,np+nz_pml+Nz
      temp=-1.50*log(0.0001)*c(ix,iz)
      d_x=temp*(ix-np-nx_pml-1)*(ix-np-nx_pml-1)/(nx_pml*nx_pml*nx_pml)/dx
      dd_x=d_x*2/(dx*(ix-np-nx_pml-1))
      iix=ix-np
      iiz=iz-nx_pml-np
      u2x=(-u(ix+2,iz)+16.0*u(ix+1,iz)-30.0*u(ix,iz)+16.0*u(ix-1,iz)-u(ix-2,iz))/12.0
      u2x=u2x/dx/dx
      u2z=(-u(ix,iz+2)+16.0*u(ix,iz+1)-30.0*u(ix,iz)+16.0*u(ix,iz-1)-u(ix,iz-2))/12.0
      u2z=u2z/dz/dz
      ux=(u(ix+1,iz)-u(ix-1,iz))/(2.0*dx)
      uz=(u(ix,iz+1)-u(ix,iz-1))/(2.0*dz)
      r1(1,iix,iiz)=(      cdt2(ix,iz)*u2x+(2.-d_x*d_x*dt*dt)*r(1,iix,iiz)+(d_x*dt-1)*r0(1,iix,iiz))/(d_x*dt+1)
      r1(4,iix,iiz)=(-dd_x*cdt2(ix,iz)*ux +(2.-d_x*d_x*dt*dt)*r(4,iix,iiz)+(d_x*dt-1)*r0(4,iix,iiz))/(d_x*dt+1)
      r1(3,iix,iiz)=       cdt2(ix,iz)*u2z+2.*r(3,iix,iiz)-r0(3,iix,iiz)
      r1(2,iix,iiz)=r(4,iix,iiz)*dt+r(2,iix,iiz)*(1-d_x*dt)
      u1(ix,iz)=r1(1,iix,iiz)+r1(2,iix,iiz)+r1(3,iix,iiz)
      enddo;enddo
      r0=r; r=r1
    end subroutine spmll
    
    subroutine spmlr(np,nx_pml,nz_pml,nx,nz,u1,u,r1,r,r0,c,cdt2,dt,dx,dz)
      implicit none
      integer(kind=4)                :: ix,iz,iix,iiz,np,nx_pml,nz_pml,nx,nz
      real(kind=rp),dimension(:,:)    :: u,u1,c,cdt2
      real(kind=rp),dimension(:,:,:)  :: r1,r,r0
      real(kind=rp)       :: temp,d_x,dd_x,dt,dx,dz,u2x,u2z,ux,uz
      
      do ix=np+nx_pml+Nx+1,np+nx_pml+Nx+nx_pml ;do iz=np+nz_pml+1,np+nz_pml+Nz
      temp=-1.50*log(0.0001)*c(ix,iz)
      d_x=temp*(ix-np-nx_pml-Nx)*(ix-np-nx_pml-Nx)/(nx_pml*nx_pml*nx_pml)/dx
      dd_x=d_x*2/(dx*(ix-np-nx_pml-Nx))
      iix=ix-np-nx_pml-Nx
      iiz=iz-np-nz_pml
      u2x=(-u(ix+2,iz)+16.0*u(ix+1,iz)-30.0*u(ix,iz)+16.0*u(ix-1,iz)-u(ix-2,iz))/12.0
      u2x=u2x/dx/dx
      u2z=(-u(ix,iz+2)+16.0*u(ix,iz+1)-30.0*u(ix,iz)+16.0*u(ix,iz-1)-u(ix,iz-2))/12.0
      u2z=u2z/dz/dz
      ux=(u(ix+1,iz)-u(ix-1,iz))/(2.0*dx)
      uz=(u(ix,iz+1)-u(ix,iz-1))/(2.0*dz)
      r1(1,iix,iiz)=(cdt2(ix,iz)*u2x+(2.-d_x*d_x*dt*dt)*r(1,iix,iiz)+(d_x*dt-1)*r0(1,iix,iiz))/(d_x*dt+1)
      r1(4,iix,iiz)=(-dd_x*cdt2(ix,iz)*ux+(2.-d_x*d_x*dt*dt)*r(4,iix,iiz)+(d_x*dt-1)*r0(4,iix,iiz))/(d_x*dt+1)
      r1(3,iix,iiz)=cdt2(ix,iz)*u2z+2.*r(3,iix,iiz)-r0(3,iix,iiz)
      r1(2,iix,iiz)=r(4,iix,iiz)*dt+r(2,iix,iiz)*(1-d_x*dt)
      u1(ix,iz)=r1(1,iix,iiz)+r1(2,iix,iiz)+r1(3,iix,iiz)
      enddo;enddo
      r0=r; r=r1
    end subroutine spmlr

    subroutine spml_tl(np,nx_pml,nz_pml,nx,nz,u1,u,s1,s,s0,r1,r,r0,c,cdt2,dt,dx,dz)
      implicit none
      integer(kind=4)                :: ix,iz,iix,iiz,np,nx_pml,nz_pml,nx,nz
      real(kind=rp),dimension(:,:)    :: u,u1,c,cdt2
      real(kind=rp),dimension(:,:,:)  :: s1,s,s0,r1,r,r0
      real(kind=rp)       :: temp,d_x,dd_x,d_z,dd_z,dt,dx,dz,u2x,u2z,ux,uz
      
      do ix=np+1,np+nx_pml;do iz=np+1,np+nz_pml    
      temp=-1.50*log(0.0001)*c(ix,iz)
      d_x=temp*(ix-np-nx_pml-1)*(ix-np-nx_pml-1)/(nx_pml*nx_pml*nx_pml)/dx
      dd_x=d_x*2/(dx*(ix-np-nx_pml-1))
      d_z=temp*(iz-np-nz_pml-1)*(iz-np-nz_pml-1)/(nz_pml*nz_pml*nz_pml)/dz
      dd_z=d_z*2/(dz*(iz-np-nz_pml-1))
      iix=ix-np
      iiz=iz-np
      u2x=(-u(ix+2,iz)+16.0*u(ix+1,iz)-30.0*u(ix,iz)+16.0*u(ix-1,iz)-u(ix-2,iz))/12.0
      u2x=u2x/dx/dx
      u2z=(-u(ix,iz+2)+16.0*u(ix,iz+1)-30.0*u(ix,iz)+16.0*u(ix,iz-1)-u(ix,iz-2))/12.0
      u2z=u2z/dz/dz
      ux=(u(ix+1,iz)-u(ix-1,iz))/(2.0*dx)
      uz=(u(ix,iz+1)-u(ix,iz-1))/(2.0*dz)
      
      s1(1,iix,iiz)=(      cdt2(ix,iz)*u2z+(2.-d_z*d_z*dt*dt)*s(1,iix,iiz)+(d_z*dt-1)*s0(1,iix,iiz))/(d_z*dt+1)
      s1(3,iix,iiz)=(-dd_z*cdt2(ix,iz)*uz +(2.-d_z*d_z*dt*dt)*s(3,iix,iiz)+(d_z*dt-1)*s0(3,iix,iiz))/(d_z*dt+1)
      s1(2,iix,iiz)=s(3,iix,iiz)*dt+s(2,iix,iiz)*(1-d_z*dt)
      
      r1(1,iix,iiz)=(      cdt2(ix,iz)*u2x+(2.-d_x*d_x*dt*dt)*r(1,iix,iiz)+(d_x*dt-1)*r0(1,iix,iiz))/(d_x*dt+1)
      r1(3,iix,iiz)=(-dd_x*cdt2(ix,iz)*ux +(2.-d_x*d_x*dt*dt)*r(3,iix,iiz)+(d_x*dt-1)*r0(3,iix,iiz))/(d_x*dt+1)
      r1(2,iix,iiz)=r(3,iix,iiz)*dt+r(2,iix,iiz)*(1-d_x*dt)
      
      u1(ix,iz)=s1(1,iix,iiz)+s1(2,iix,iiz)+r1(1,iix,iiz)+r1(2,iix,iiz)
      enddo;enddo
      s0=s; s=s1; r0=r; r=r1
  end subroutine spml_tl

  subroutine spml_tr(np,nx_pml,nz_pml,nx,nz,u1,u,s1,s,s0,r1,r,r0,c,cdt2,dt,dx,dz)
      implicit none
      integer(kind=4)                :: ix,iz,iix,iiz,np,nx_pml,nz_pml,nx,nz
      real(kind=rp),dimension(:,:)    :: u,u1,c,cdt2
      real(kind=rp),dimension(:,:,:)  :: s1,s,s0,r1,r,r0
      real(kind=rp)       :: temp,d_x,dd_x,d_z,dd_z,dt,dx,dz,u2x,u2z,ux,uz
      
      do ix=np+nx_pml+nx+1,np+nx_pml+nx+nx_pml;do iz=np+1,np+nz_pml    
      temp=-1.50*log(0.0001)*c(ix,iz)
      d_x=temp*(ix-np-nx_pml-nx)*(ix-np-nx_pml-nx)/(nx_pml*nx_pml*nx_pml)/dx
      dd_x=d_x*2/(dx*(ix-np-nx_pml-nx))
      d_z=temp*(iz-np-nz_pml-1)*(iz-np-nz_pml-1)/(nz_pml*nz_pml*nz_pml)/dz
      dd_z=d_z*2/(dz*(iz-np-nz_pml-1))
      iix=ix-nx-nx_pml-np
      iiz=iz-np
      u2x=(-u(ix+2,iz)+16.0*u(ix+1,iz)-30.0*u(ix,iz)+16.0*u(ix-1,iz)-u(ix-2,iz))/12.0
      u2x=u2x/dx/dx
      u2z=(-u(ix,iz+2)+16.0*u(ix,iz+1)-30.0*u(ix,iz)+16.0*u(ix,iz-1)-u(ix,iz-2))/12.0
      u2z=u2z/dz/dz
      ux=(u(ix+1,iz)-u(ix-1,iz))/(2.0*dx)
      uz=(u(ix,iz+1)-u(ix,iz-1))/(2.0*dz)
      
      s1(1,iix,iiz)=(      cdt2(ix,iz)*u2z+(2.-d_z*d_z*dt*dt)*s(1,iix,iiz)+(d_z*dt-1)*s0(1,iix,iiz))/(d_z*dt+1)
      s1(3,iix,iiz)=(-dd_z*cdt2(ix,iz)*uz +(2.-d_z*d_z*dt*dt)*s(3,iix,iiz)+(d_z*dt-1)*s0(3,iix,iiz))/(d_z*dt+1)
      s1(2,iix,iiz)=s(3,iix,iiz)*dt+s(2,iix,iiz)*(1-d_z*dt)
      
      r1(1,iix,iiz)=(      cdt2(ix,iz)*u2x+(2.-d_x*d_x*dt*dt)*r(1,iix,iiz)+(d_x*dt-1)*r0(1,iix,iiz))/(d_x*dt+1)
      r1(3,iix,iiz)=(-dd_x*cdt2(ix,iz)*ux +(2.-d_x*d_x*dt*dt)*r(3,iix,iiz)+(d_x*dt-1)*r0(3,iix,iiz))/(d_x*dt+1)
      r1(2,iix,iiz)=r(3,iix,iiz)*dt+r(2,iix,iiz)*(1-d_x*dt)
      
      u1(ix,iz)=s1(1,iix,iiz)+s1(2,iix,iiz)+r1(1,iix,iiz)+r1(2,iix,iiz)
      enddo;enddo
      s0=s; s=s1; r0=r; r=r1
  end subroutine spml_tr

  subroutine spml_bl(np,nx_pml,nz_pml,nx,nz,u1,u,s1,s,s0,r1,r,r0,c,cdt2,dt,dx,dz)
      implicit none
      integer(kind=4)                :: ix,iz,iix,iiz,np,nx_pml,nz_pml,nx,nz
      real(kind=rp),dimension(:,:)    :: u,u1,c,cdt2
      real(kind=rp),dimension(:,:,:)  :: s1,s,s0,r1,r,r0
      real(kind=rp)       :: temp,d_x,dd_x,d_z,dd_z,dt,dx,dz,u2x,u2z,ux,uz
      
      do ix=np+1,np+nx_pml;do iz=np+nz_pml+nz+1,np+nz_pml+nz+nz_pml    
      temp=-1.50*log(0.0001)*c(ix,iz)
      d_x=temp*(ix-np-nx_pml-1)*(ix-np-nx_pml-1)/(nx_pml*nx_pml*nx_pml)/dx
      dd_x=d_x*2/(dx*(ix-np-nx_pml-1))
      d_z=temp*(iz-np-nz_pml-nz)*(iz-np-nz_pml-nz)/(nz_pml*nz_pml*nz_pml)/dz
      dd_z=d_z*2/(dz*(iz-np-nz_pml-nz))
      iix=ix-np
      iiz=iz-nz-nz_pml-np
      u2x=(-u(ix+2,iz)+16.0*u(ix+1,iz)-30.0*u(ix,iz)+16.0*u(ix-1,iz)-u(ix-2,iz))/12.0
      u2x=u2x/dx/dx
      u2z=(-u(ix,iz+2)+16.0*u(ix,iz+1)-30.0*u(ix,iz)+16.0*u(ix,iz-1)-u(ix,iz-2))/12.0
      u2z=u2z/dz/dz
      ux=(u(ix+1,iz)-u(ix-1,iz))/(2.0*dx)
      uz=(u(ix,iz+1)-u(ix,iz-1))/(2.0*dz)
      
      s1(1,iix,iiz)=(      cdt2(ix,iz)*u2z+(2.-d_z*d_z*dt*dt)*s(1,iix,iiz)+(d_z*dt-1)*s0(1,iix,iiz))/(d_z*dt+1)
      s1(3,iix,iiz)=(-dd_z*cdt2(ix,iz)*uz +(2.-d_z*d_z*dt*dt)*s(3,iix,iiz)+(d_z*dt-1)*s0(3,iix,iiz))/(d_z*dt+1)
      s1(2,iix,iiz)=s(3,iix,iiz)*dt+s(2,iix,iiz)*(1-d_z*dt)
      
      r1(1,iix,iiz)=(      cdt2(ix,iz)*u2x+(2.-d_x*d_x*dt*dt)*r(1,iix,iiz)+(d_x*dt-1)*r0(1,iix,iiz))/(d_x*dt+1)
      r1(3,iix,iiz)=(-dd_x*cdt2(ix,iz)*ux +(2.-d_x*d_x*dt*dt)*r(3,iix,iiz)+(d_x*dt-1)*r0(3,iix,iiz))/(d_x*dt+1)
      r1(2,iix,iiz)=r(3,iix,iiz)*dt+r(2,iix,iiz)*(1-d_x*dt)
      
      u1(ix,iz)=s1(1,iix,iiz)+s1(2,iix,iiz)+r1(1,iix,iiz)+r1(2,iix,iiz)
      enddo;enddo
      s0=s; s=s1; r0=r; r=r1
  end subroutine spml_bl

  subroutine spml_br(np,nx_pml,nz_pml,nx,nz,u1,u,s1,s,s0,r1,r,r0,c,cdt2,dt,dx,dz)
      implicit none
      integer(kind=4)                :: ix,iz,iix,iiz,np,nx_pml,nz_pml,nx,nz
      real(kind=rp),dimension(:,:)    :: u,u1,c,cdt2
      real(kind=rp),dimension(:,:,:)  :: s1,s,s0,r1,r,r0
      real(kind=rp)       :: temp,d_x,dd_x,d_z,dd_z,dt,dx,dz,u2x,u2z,ux,uz
      
      do ix=np+nx_pml+nx+1,np+nx_pml+nx+nx_pml;do iz=np+nz_pml+nz+1,np+nz_pml+nz+nz_pml    
      temp=-1.50*log(0.0001)*c(ix,iz)
      d_x=temp*(ix-np-nx_pml-nx)*(ix-np-nx_pml-nx)/(nx_pml*nx_pml*nx_pml)/dx
      dd_x=d_x*2/(dx*(ix-np-nx_pml-nx))
      d_z=temp*(iz-np-nz_pml-nz)*(iz-np-nz_pml-nz)/(nz_pml*nz_pml*nz_pml)/dz
      dd_z=d_z*2/(dz*(iz-np-nz_pml-nz))
      iix=ix-nx-nx_pml-np
      iiz=iz-nz-nz_pml-np
      u2x=(-u(ix+2,iz)+16.0*u(ix+1,iz)-30.0*u(ix,iz)+16.0*u(ix-1,iz)-u(ix-2,iz))/12.0
      u2x=u2x/dx/dx
      u2z=(-u(ix,iz+2)+16.0*u(ix,iz+1)-30.0*u(ix,iz)+16.0*u(ix,iz-1)-u(ix,iz-2))/12.0
      u2z=u2z/dz/dz
      ux=(u(ix+1,iz)-u(ix-1,iz))/(2.0*dx)
      uz=(u(ix,iz+1)-u(ix,iz-1))/(2.0*dz)
      
      s1(1,iix,iiz)=(      cdt2(ix,iz)*u2z+(2.-d_z*d_z*dt*dt)*s(1,iix,iiz)+(d_z*dt-1)*s0(1,iix,iiz))/(d_z*dt+1)
      s1(3,iix,iiz)=(-dd_z*cdt2(ix,iz)*uz +(2.-d_z*d_z*dt*dt)*s(3,iix,iiz)+(d_z*dt-1)*s0(3,iix,iiz))/(d_z*dt+1)
      s1(2,iix,iiz)=s(3,iix,iiz)*dt+s(2,iix,iiz)*(1-d_z*dt)
      
      r1(1,iix,iiz)=(      cdt2(ix,iz)*u2x+(2.-d_x*d_x*dt*dt)*r(1,iix,iiz)+(d_x*dt-1)*r0(1,iix,iiz))/(d_x*dt+1)
      r1(3,iix,iiz)=(-dd_x*cdt2(ix,iz)*ux +(2.-d_x*d_x*dt*dt)*r(3,iix,iiz)+(d_x*dt-1)*r0(3,iix,iiz))/(d_x*dt+1)
      r1(2,iix,iiz)=r(3,iix,iiz)*dt+r(2,iix,iiz)*(1-d_x*dt)
      
      u1(ix,iz)=s1(1,iix,iiz)+s1(2,iix,iiz)+r1(1,iix,iiz)+r1(2,iix,iiz)
      enddo;enddo
      s0=s; s=s1; r0=r; r=r1
  end subroutine spml_br

end module lwc4v
 
program lwc4
    ! the 4th order for accoustic wave equation
    ! with PML
    use lwc4v
    implicit none
    integer                   :: Nx, Nz,Nt
    integer                   :: ix, iz, it,sx,sz,recz,recx,ke
    integer                   :: np=2,nx_pml=20,nz_pml=20
    real(kind=rp),parameter    :: pi=3.1415926
    real(kind=rp),allocatable  :: u(:,:),u0(:,:),u1(:,:)! u0 at n-1, u at n, u1 at (n+1)
    real(kind=rp),allocatable  :: ts1(:,:,:),ts(:,:,:),ts0(:,:,:),lr1(:,:,:),lr0(:,:,:),lr(:,:,:)
    real(kind=rp),allocatable  :: bs1(:,:,:),bs(:,:,:),bs0(:,:,:),rr1(:,:,:),rr0(:,:,:),rr(:,:,:)
    real(kind=rp),allocatable  :: tls1(:,:,:),tls(:,:,:),tls0(:,:,:),tlr1(:,:,:),tlr0(:,:,:),tlr(:,:,:)
    real(kind=rp),allocatable  :: trs1(:,:,:),trs(:,:,:),trs0(:,:,:),trr1(:,:,:),trr0(:,:,:),trr(:,:,:)
    real(kind=rp),allocatable  :: bls1(:,:,:),bls(:,:,:),bls0(:,:,:),blr1(:,:,:),blr0(:,:,:),blr(:,:,:)
    real(kind=rp),allocatable  :: brs1(:,:,:),brs(:,:,:),brs0(:,:,:),brr1(:,:,:),brr0(:,:,:),brr(:,:,:)
    real(kind=rp),allocatable  :: c(:,:),cc(:,:),cdt2(:,:),srec1(:,:),srec2(:,:),vrec(:,:)
    real(kind=rp)              :: dx,dz,dt,u2t
    real(kind=rp)              :: t,td,f0,ffs,fft,ff2t,tt
    real(kind=rp)              :: x,z
    real(kind=rp)              :: u2x,u2z,u2t,u4t,u2x0,u2x1,u2z0,u2z1,u2x2z,zz1,zz2,zz3,zz4,u2x2t,u2z2t,u4x,u4z
    integer(kind=4) :: iix,iiz
    real   (kind=rp) :: d_z,dd_z,d_x,dd_x,temp

    !x=3.0; z=3.0;  
    tt=2.0
    dx=0.01;   dz=0.01;   dt=0.001
    !x=7.62; z=2.44;  tt=1.5
    !dx=0.02;   dz=0.02;   dt=0.0005
    Nx=floor(x/dx)+1; Nz=floor(z/dz)+1 ; Nt=floor(tt/dt)+1
    Nx=299; Nz=199
    write(*,*) "Nx=",nx,"Nz=",nz,"nt=",nt
    nx_pml=20; nz_pml=20; 
    recz=Np+nz_pml+1;
    recx=Np+nx_pml+20;
    f0=15.0;td=1.0/f0
    
    sx=np+nx_pml+Nx/2+1; 
    sz=recz;
    
    allocate(c   (Nx+2*np+2*nx_pml, Nz+2*np+2*nx_pml))
    write(*,*) size(c)   
    open(unit=10,file='model.txt',status='unknown')
    do ix=np+nx_pml+1,np+nx_pml+nx;do iz=np+nz_pml+1,np+nz_pml+nz
    read(10,*) c(ix,iz)
    !c(ix,iz)=c(ix,iz)/1000
    enddo;enddo
    close(10)
    call expand(c,np,nx,nz,nx_pml,nz_pml)
    allocate(u   (Nx+2*np+2*nx_pml, Nz+2*np+2*nz_pml))
    allocate(u0  (Nx+2*np+2*nx_pml, Nz+2*np+2*nz_pml))
    allocate(u1  (Nx+2*np+2*nx_pml, Nz+2*np+2*nz_pml))    
    allocate(cdt2(Nx+2*np+2*Nx_pml, Nz+2*np+2*Nz_pml))
    allocate(cc (Nx+2*np+2*Nx_pml, Nz+2*np+2*Nz_pml))
    allocate(srec1(Nx+2*np+2*Nx_pml, Nt              ))
    allocate(srec2(Nx+2*np+2*Nx_pml, Nt              ))
    allocate(vrec(Nz+2*np+2*Nz_pml, Nt              ))
    allocate(ts1(4,Nx,Nz_pml));allocate(ts (4,Nx,Nz_pml))
    allocate(ts0(4,Nx,Nz_pml));allocate(bs1(4,Nx,Nz_pml))
    allocate(bs (4,Nx,Nz_pml));allocate(bs0(4,Nx,Nz_pml))
    allocate(lr1(4,Nx_pml,Nz));allocate(lr (4,Nx_pml,Nz))
    allocate(lr0(4,Nx_pml,Nz));allocate(rr1(4,Nx_pml,Nz))
    allocate(rr (4,Nx_pml,Nz));allocate(rr0(4,Nx_pml,Nz))

    allocate(tls0(3,Nx_pml,Nz_pml));allocate(tls (3,Nx_pml,Nz_pml))
    allocate(tls1(3,Nx_pml,Nz_pml));allocate(tlr0(3,Nx_pml,Nz_pml))
    allocate(tlr (3,Nx_pml,Nz_pml));allocate(tlr1(3,Nx_pml,Nz_pml))

    allocate(trs0(3,Nx_pml,Nz_pml));allocate(trs (3,Nx_pml,Nz_pml))
    allocate(trs1(3,Nx_pml,Nz_pml));allocate(trr0(3,Nx_pml,Nz_pml))
    allocate(trr (3,Nx_pml,Nz_pml));allocate(trr1(3,Nx_pml,Nz_pml))

    allocate(bls0(3,Nx_pml,Nz_pml));allocate(bls (3,Nx_pml,Nz_pml))
    allocate(bls1(3,Nx_pml,Nz_pml));allocate(blr0(3,Nx_pml,Nz_pml))
    allocate(blr (3,Nx_pml,Nz_pml));allocate(blr1(3,Nx_pml,Nz_pml))
    
    allocate(brs0(3,Nx_pml,Nz_pml));allocate(brs (3,Nx_pml,Nz_pml))
    allocate(brs1(3,Nx_pml,Nz_pml));allocate(brr0(3,Nx_pml,Nz_pml))
    allocate(brr (3,Nx_pml,Nz_pml));allocate(brr1(3,Nx_pml,Nz_pml))
    
    u=0.0; u0=0.0; u1=0.0; 
    srec1=0.0;srec2=0;vrec=0.0
    ts=0.0; ts1=0.0; ts0=0;ts=0; bs1=0; bs0=0;
    lr=0.0; lr1=0.0; lr0=0;rr=0; rr1=0; rr0=0;
    tls=0; tls1=0.0; tls0=0; tlr=0;  tlr1=0; tlr0=0;
    trs=0; trs1=0; trs0=0; trr1=0; trr=0;  trr0=0;
    bls0=0;bls1=0; bls=0;  blr1=0; blr=0;  blr0=0;
    brs0=0;brs1=0; brs=0;  blr1=0; blr=0;  blr0=0;    
    cc=c*c; cdt2=(c*dt)**2
    open(unit=30,file="us.bin",form="unformatted",access="direct",status="replace",recl=Nx*Nz)
    do it=1, Nt
        ffs=exp(-(pi*f0*(t-td))**2)*(1.0-2.0*(pi*f0*(t-td))**2)
        t=(it-1)*dt
        !ffs=-5.76*f0**2*(1.0-16.0*(0.6*f0*t-1.0)**2)*exp(-8.0*(0.6*f0*t-1.0)**2)
        !fft=0.6*16.0*5.76*f0**3*(0.6*f0*t-1.0)*(3.0-16.0*(0.6*f0*t-1.0)**2)*exp(-8.0*(0.6*f0*t-1.0)**2)
        !ff2t=0.36*16.0*5.76*f0**4*(3.0-96.0*(0.6*f0*t-1.0)**2+16.0*16.0*(0.6*f0*t-1.0)**4)*&
        !     &exp(-8.0*(0.6*f0*t-1.0)**2)
        write(*,*) t+dt,ffs
       !omp parrallel default(private) shared(u,u1,u0)
       !omp do
        do ix=np+nx_pml+1,np+nx_pml+Nx; do iz=np+nz_pml+1,np+nz_pml+Nz
          ! computing partial 
          u2x =(-u(ix+2,iz)+16.0*u(ix+1,iz)-30.0*u(ix,iz)+16.0*u(ix-1,iz)-u(ix-2,iz))/(12.*dx**2)
          !u2x1=(16.*(u(ix+2,iz)+u(ix,iz))-(u(ix+3,iz)+u(ix-1,iz))-30.*u(ix+1,iz))/(12.*dx**2)
	  ! u2x0=(16.*(u(ix,iz)+u(ix-2,iz))-(u(ix+1,iz)+u(ix-3,iz))-30.*u(ix-1,iz))/(12.*dx**2)
         
          u2z=(-u(ix,iz+2)+16.0*u(ix,iz+1)-30.0*u(ix,iz)+16.0*u(ix,iz-1)-u(ix,iz-2))/(12.0*dz**2)
          !u2z1=(16.0*(u(ix,iz+2)+u(ix,iz))-(u(ix,iz+3)+u(ix,iz-1))-30.0*u(ix,iz+1))/(12.0*dz**2)
	  ! u2z0=(16.*(u(ix,iz)+u(ix,iz-2))-(u(ix,iz+1)+u(ix,iz-3))-30.0*u(ix,iz-1))/(12*dz**2)

          !u4x=(u2x1-2*u2x+u2x0)/dx**2
	   !u4z=(u2z1-2*u2z+u2z0)/dz**2
          
          !ZZ1=(16.*(U(ix+1,iz+1)+U(ix-1,iz+1))-(U(ix+2,iz+1)+U(ix-2,iz+1))-30.*U(ix,iz+1))/(12.*dx**2)
	   !ZZ2=(16.*(U(ix+1,iz-1)+U(ix-1,iz-1))-(U(ix+2,iz-1)+U(ix-2,iz-1))-30.*U(ix,iz-1))/(12.*dx**2)
	  ! ZZ3=(16.*(U(ix+1,iz+2)+U(ix-1,iz+2))-(U(ix+2,iz+2)+U(ix-2,iz+2))-30.*U(ix,iz+2))/(12.*dx**2)
	  ! ZZ4=(16.*(U(ix+1,iz-2)+U(ix-1,ix-2))-(U(ix+2,iz-2)+U(ix-2,iz-2))-30.*U(ix,iz-2))/(12.*dx**2)

	   !U2x2z=(16.*(ZZ1+ZZ2)-(ZZ3+ZZ4)-30.*U2x)/(12*dz**2)

          u2t= cc(ix,iz)*(u2x+u2z)
	   !U2x2t= U4x+U2x2z
	   !U2z2t= U2x2z+U4z
	   !U4t  = U2x2t+U2z2t
          ! applying source
          !if(ix.ge.(sx-ke).and.ix.le.(sx+ke).and.iz.ge.(sz-ke).and.iz.le.(sz+ke))then
          !if(ix.eq.sx.and.iz.eq.sz)then
          ! u2t=u2t+ffs
          !else
            !u2t=u2t+ffs/sqrt(((ix-sx)*dx)**2+((iz-sz)*dz)**2)
          !endif
          !endif
          !u1(ix,iz)=2.0*u(ix,iz)-u0(ix,iz)+cdt2(ix,iz)*u2t
          !if(((ix.eq.(sx-2)).or.(ix.eq.(sx+2))).and.(iz.ge.(sz-2)).and.(iz.le.(sz+2)))then
	   !  U2t=U2t+ffs/6.0	  
	     !U4t=U4t+ff2t/6.0
          !endif
          !if((ix.ge.(sx-1)).and.(ix.le.(sx+1)).and.((iz.eq.(sz-2)).or.(iz.eq.(sz+2))))then
	    !U2t=U2t+ffs/6.0	  
	     !U4t=U4t+ff2t/6.0
          !endif
	   !if((ix.ge.(sx-1)).and.(ix.le.(sx+1)).and.(iz.ge.(sz-1)).and.(iz.le.(sz+1)))then
	     if((ix.eq.sx).and.(iz.eq.sz))then
	        U2t=U2t+ffs/1.0
	        !U4t=U4t+ff2t/1.0
            !else
	        !U2t=U2t+ffs/2.0	  
	        !U4t=U4t+ff2t/2.0
	     endif
	   !endif
          u1(ix,iz)=2.0*u(ix,iz)-u0(ix,iz)+dt**2*u2t!+dt**2/12.0*cdt2(ix,iz)*U4t
        enddo; enddo
        !----------------- PML-----------------
        call spmlt  (np,nx_pml,nz_pml,nx,nz,u1,u,ts1,ts,ts0,c,cdt2,dt,dx,dz)
        call spmlb  (np,nx_pml,nz_pml,nx,nz,u1,u,bs1,bs,bs0,c,cdt2,dt,dx,dz)
        call spmll  (np,nx_pml,nz_pml,nx,nz,u1,u,lr1,lr,lr0,c,cdt2,dt,dx,dz)
        call spmlr  (np,nx_pml,nz_pml,nx,nz,u1,u,rr1,rr,rr0,c,cdt2,dt,dx,dz)
        call spml_tl(np,nx_pml,nz_pml,nx,nz,u1,u,tls1,tls,tls0,tlr1,tlr,tlr0,c,cdt2,dt,dx,dz)
        call spml_tr(np,nx_pml,nz_pml,nx,nz,u1,u,trs1,trs,trs0,trr1,trr,trr0,c,cdt2,dt,dx,dz)
        call spml_bl(np,nx_pml,nz_pml,nx,nz,u1,u,bls1,bls,bls0,blr1,blr,blr0,c,cdt2,dt,dx,dz)
        call spml_br(np,nx_pml,nz_pml,nx,nz,u1,u,brs1,brs,brs0,brr1,brr,brr0,c,cdt2,dt,dx,dz)
        !-----------------PML-------------
        srec1(:,it)= u(:, recz)
        vrec (:,it)= u(recx,:)
        write(30,rec=it) ((u(ix,iz),iz=np+nz_pml+1,np+nz_pml+Nz),ix=np+nx_pml+1,np+nx_pml+Nx)
        ! swap u0, u,u1
        u0=u; u=u1;
        !omp end do
    !omp end parrallel        
    enddo
    close(30)
    ! srec1(sx,:)=0.0
    open(unit=20,file='trec.txt')
    open(unit=30,file='lrec.txt')
    do ix=np+nx_pml+1,np+nx_pml+Nx;do it=1,Nt 
       write(20,*) srec1(ix,it)
    enddo;enddo
    do iz=np+nx_pml+1,np+nx_pml+Nx;do it=1,Nt 
       write(30,*) vrec(iz,it)
    enddo;enddo
    close(20);close(30)
    open(unit=10,file='usnap.txt')
    do ix=np+nx_pml+1,np+nx_pml+Nx; do iz=np+nz_pml+1,np+Nz_pml+Nz
    write(10,*) u(ix,iz)
    enddo; enddo
    open(unit=10,file='us.txt')
    do ix=1,2*np+2*nx_pml+Nx; do iz=1,2*np+2*Nz_pml+Nz
    write(10,*) u(ix,iz)
    enddo; enddo
    close(10)
end program
