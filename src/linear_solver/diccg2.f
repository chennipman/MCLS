C-------------------------------------------------------------------------------
C Deflated CG algorithm Variant (Iterative Inner Solver)
C x: Ax = y
C-------------------------------------------------------------------------------
C     subroutines called
C	calltime        <- in this file
C	makeElemAZ2     <- in this file
C	fill_mat        <- external, but added at end
C	computeAZandE2  <- in this file
C	fill_vec        <- external, but added at end
C	make_u2         <- in this file
C	make_cholinc    <- external, but added at end
C	copy_vec        <- external, but added at end
C	matvec_prod     <- external
C	lin_comb        <- external, but added at end
C	make_Px2        <- in this file
C	precond         <- external, but added at end
C	dot_prod        <- external, but added at end
C	make_Px2        <- in this file
C	make_Ptx2       <- in this file
C
C
C	

      subroutine diccg2(Nx, Ny, Nz, A, y, tol, x,iter, 
     1                  norm1,norm2, rx, epsIn)
      implicit none
      
      integer Nx, Ny, Nz, prec, it, NoVec, reson
      double precision A(Nx,Ny,Nz,4), y(Nx,Ny,Nz),tol,x(Nx,Ny,Nz)
      double precision timedcg
      
      integer Nt,i,j,MAXITER, iter
      double precision norm1, norm2
      integer Nr, Nb,elemAZ, rx
      
      double precision alpha, temp, tol2, rz, rz_old
      double precision time1, time2, time3, time4,timedef
      double precision timeit,realres, res(1000), epsIn
      double precision norm_res2, norm_res2_0, normr0

      double precision, dimension(:,:,:), 
     1        allocatable :: p, b, r, rh, rh_old, r_old, temp1,r0,
     1                       D,temp_vec,temp_vecd, z, z0, Axran
      double precision, dimension(:,:), 
     1        allocatable ::   E, AZ
      double precision, dimension(:), 
     1        allocatable :: Px, Ptx, u,xAZ1, xAZ2, xran
     
      INTEGER ISEED,k 
      REAL RANF
      
      character*128 FNAM
       
      parameter(MAXITER = 1000)

      NoVec=rx**3 

      allocate(p(Nx,Ny,Nz))
      allocate(b(Nx,Ny,Nz))
      allocate(r(Nx,Ny,Nz))
      allocate(r0(Nx,Ny,Nz))
      allocate(r_old(Nx,Ny,Nz))
      allocate(rh(Nx,Ny,Nz))
      
      allocate(rh_old(Nx,Ny,Nz))
      allocate(D(Nx,Ny,Nz))
      allocate(temp_vec(Nx,Ny,Nz))
      allocate(temp_vecd(Nx,Ny,Nz))
      allocate(z(Nx,Ny,Nz))
      
      allocate(z0(Nx,Ny,Nz))
      allocate(temp1(Nx,Ny,Nz))
      allocate(Px(Nx*Ny*Nz))
      allocate(Ptx(Nx*Ny*Nz)) 
      allocate(u(Nx*Ny*Nz))     
      allocate(E(NoVec,4))
      allocate(xAZ1(Nx*Ny*Nz))
      allocate(xAZ2(NoVec))
      allocate(Axran(Nx,Ny,Nz))
      allocate(xran(Nx*Ny*Nz))
      
      Nt = Nx*Ny*Nz
      tol2 = tol*tol
      
C--------------------------------------------------------------
C Preparation of the Deflated PCG method   
C--------------------------------------------------------------
      
      call calltime(time1) 
      
      call makeElemAZ2(elemAZ,Nb,Nr,rx,NoVec,Nt,Nx)
      
      allocate(AZ(elemAZ,3))
      call fill_mat(elemAZ,3,AZ)
      call fill_mat(NoVec,4,E)
      call computeAZandE2(NoVec,Nx,Ny,Nz,A,AZ,E,elemAZ,
     1                   Nb,rx) 

      call fill_vec(1000, 0d0, res)
      call fill_vec(Nt,   0d0, Px)
      call fill_vec(Nt,   0d0, Ptx)
      call fill_vec(Nt,   0d0, u)

      call make_u2(NoVec,Nt,Nb,rx, E,y,u,Nx,Ny,Nz,epsIn)
      
      call make_cholinc(Nx, Ny, Nz, A, D)
      call fill_vec(Nt, 0d0, x)
      call fill_vec(Nt, 0d0, rh)  
      
      call copy_vec(Nt, x, xran) 
      
      call matvec_prod(Nx, Ny, Nz, A, xran, Axran)      
      call lin_comb(Nt, y, Axran, -1d0, r)
      
      call copy_vec(Nt, r,r0)
    
      call make_Px2(r,NoVec, Nt, AZ, E,rh,elemAZ,Nb,
     1             rx,Nx,Ny,Nz,epsIn)
      
      call precond(Nx, Ny, Nz, A, D, rh, z)

      call copy_vec(Nt, z, p)
      call dot_prod(Nt, rh, z, rz)
      
      call precond(Nx, Ny, Nz, A, D, r, z0)
      call dot_prod(Nt, z0, z0, norm_res2_0)

      
      call calltime(time3)
      timedef=time3-time1
      
C-------------------------------------------------------------
C Deflated PCG method...    
C-------------------------------------------------------------

      do i=1,MAXITER
      
	call matvec_prod(Nx, Ny, Nz, A, p, temp_vec)

	call make_Px2(temp_vec,NoVec, Nt, AZ, E,
     1           temp_vecd,elemAZ,Nb,rx,Nx,Ny,Nz,epsIn)

	call dot_prod(Nt, p, temp_vecd, temp)

	alpha = rz/temp
	call lin_comb(Nt, x, p, alpha, x)
	call lin_comb(Nt, rh, temp_vecd, -alpha, rh)
	call dot_prod(Nt, z, z, norm_res2)
	norm_res2 = norm_res2 / norm_res2_0

	res(i) =  sqrt(norm_res2)

C	if (reson .eq. 1) then
C	  write(6,'(I10,2X,E9.3)') i, sqrt(norm_res2)
C	endif
	
	if (norm_res2 .le. tol2) goto 100
	
        call precond(Nx, Ny, Nz, A, D, rh, z)
	rz_old = rz
	call dot_prod(Nt, rh, z, rz)
	call lin_comb(Nt, z, p, rz/rz_old, p)
      end do
      
      write(6,"('DICCG2 didn''t converge')")
      stop
      
  100 continue
  
      iter = i
      norm1 = sqrt(norm_res2)
      
C-------------------------------------------------------------
C Correction steps of DPCG    
C-------------------------------------------------------------     
      call fill_vec(Nt,0d0,xAZ1)
      call fill_vec(NoVec,0d0,xAZ2)
       
      call make_Ptx2(NoVec, Nt, AZ,E,x, Ptx,u,x,
     1              elemAZ,Nb,rx,Nx,Ny,Nz,epsIn) 

      call matvec_prod(Nx, Ny, Nz, A, x, temp1)
      call lin_comb(Nt, temp1, y, -1d0, temp1)
      call dot_prod(Nt, temp1, temp1, realres)
      call dot_prod(Nt, r0, r0, normr0)
      realres =  sqrt(realres)/sqrt(normr0)
      norm2=realres

      call calltime(time2)
      timeit=time2-time3
      timedcg=time2-time1
      
      it=i

  200 continue
  
      deallocate(p)
      deallocate(b)
      deallocate(r)
      deallocate(r0)
      deallocate(r_old)
      deallocate(rh)
      
      deallocate(rh_old)
      deallocate(D)
      deallocate(temp_vec)
      deallocate(temp_vecd)
      deallocate(z)
      
      deallocate(z0)
      deallocate(Px)
      deallocate(Ptx)
      deallocate(u)
      deallocate(AZ)
      
      deallocate(E)
      deallocate(xAZ1)
      deallocate(xAZ2)
      deallocate(Axran)
      deallocate(xran)
      
      deallocate(temp1)
  
      return
      end
C---------------------------------------------------------------
C Compute number of elements of AZ and other related parameters
C---------------------------------------------------------------
      subroutine makeElemAZ2(elemAZ,Nb,Nr,rx,NoVec, Nt,Nx)
      implicit none

      integer Nr,Nb,nb1,nc,ni,nib,elemAZ,rx, Nt,NoVec,
     1        correc, Nx 

      Nr = Nt/NoVec
      Nb = Nx/rx
      
      nc=48*Nb**2-24*Nb+8
      ni=(rx-2)**3*(12*Nb**2-12*Nb+8)
      nb1=12*(rx-2)*(8*Nb**2-5*Nb+2)
      nib=6*(rx-2)**2*(10*Nb**2-8*Nb+4)
      correc=1
      elemAZ=ni+nb1+nc+nib

      return
      end
C-------------------------------------------------------------------------------
C   First Deflation Steps: making u - iterative E!
C   u=Z*inv(E)Zt*b
C-------------------------------------------------------------------------------
      subroutine make_u2(NoVec,Nt,Nb,rx,E,y,u,Nx,
     1                   Ny,Nz, epsIn)
      implicit none  
       
      integer Nt, NoVec,Nx,Ny,Nz, Nb,rx,i,j,k
 
C------Extra      
      integer itE
      real norm1, norm2
      double precision epsIn
C-----------
      
      double precision, dimension(:), 
     1        allocatable :: y1,y2
     
      double precision y(Nt),u(Nt),E(NoVec,4)
      allocate(y1(NoVec))
      allocate(y2(NoVec))     

  
      call fill_vec(NoVec,0d0,y1)
      call fill_vec(NoVec,0d0,y2)
      call fill_vec(NoVec,0d0,u)
      
      call Ztx(y,y1,Nx,Ny,Nz,NoVec,Nb,rx)
      
C----------- Changed!
      norm2=0d0
      call cg(rx, rx, rx,E,y1,epsIn,y2,itE,norm1,norm2)
     
C      call cg(imax, jmax, kmax, A, y, tolmax, p, iter, norm1, norm2)
C      write(*,*) 'Intern CG:', itE


C      call band_chol_solve(NoVec,E,y2,y1,rx**2) 

C-------------------


      call Zx(y2,u,Nx,Ny,Nz,NoVec,Nb,rx)
      
      
      deallocate(y1)
      deallocate(y2) 
      
      return
      end
C-------------------------------------------------------------------------------
C  End Deflation Steps - Iterative Solve E
C  Pt*x
C-------------------------------------------------------------------------------
      subroutine make_Ptx2(NoVec, Nt,AZ,E,xtd,Ptx,u,xdef,elemAZ,Nb,
     1                    rx,Nx,Ny,Nz,epsIn)
      implicit none

C------Extra      
      integer itE
      real norm1, norm2
      double precision epsIn
C-----------

      integer Nt, NoVec
      integer elemAZ,Nb,rx,Nx,Ny,Nz
      double precision AZ(elemAZ,3), E(NoVec,4),xdef(Nt)
      double precision xtd(Nt), Ptx(Nt), alfa, u(Nt), alfa1
      
      double precision, dimension(:), 
     1        allocatable :: xt1,xt2,xt3
     
      allocate(xt1(NoVec))
      allocate(xt2(NoVec)) 
      allocate(xt3(Nt))   
      
      call fill_vec(NoVec,0d0,xt1)
      call fill_vec(NoVec,0d0,xt2)
      call fill_vec(Nt,0d0,xt3)
      
      alfa=-1d0
      alfa1=1d0
      
      call AZtx(AZ,xtd,xt1,Nx,Ny,Nz,NoVec,elemAZ)

C----------- Changed!
      norm2=0d0
  
      call cg(rx, rx, rx,E,xt1,epsIn,xt2,itE,norm1,norm2)
     
C      call cg(imax, jmax, kmax, A, y, tolmax, p, iter, norm1, norm2)
C      write(*,*) 'Intern CG:', itE

C      call band_chol_solve(NoVec,E,xt2,xt1,rx**2)
C      call band_chol_solve(NoVec,E,y2,y1,rx**2) 

C-------------------
      

      
      call Zx(xt2,xt3,Nx,Ny,Nz,NoVec,Nb,rx)
      call lin_comb(Nt,xtd,xt3,alfa,Ptx)      
      call lin_comb(Nt,u,Ptx,alfa1,xdef) 

      deallocate(xt1)
      deallocate(xt2) 
      deallocate(xt3)   

      return
      end   
C-------------------------------------------------------------------------------
C  Main Deflation Steps - Iterative solve E!
C  Px
C-------------------------------------------------------------------------------
      subroutine make_Px2(xd, NoVec, Nt,AZ,E,Px1, elemAZ,Nb,
     1                   rx,Nx,Ny,Nz, epsIn)
      implicit none
      
C------Extra      
      integer itE
      real norm1, norm2
      double precision epsIn
C-----------
      
      integer Nt, NoVec, elemAZ,Nb,rx,Nx,Ny,Nz
      double precision AZ(elemAZ,3),E(NoVec,4)
      double precision xd(Nt),Px1(Nt)
      double precision, dimension(:), 
     1        allocatable :: x1,x2,x3
      
      allocate(x1(NoVec))
      allocate(x2(NoVec)) 
      allocate(x3(Nt))  
      
            
      call fill_vec(NoVec,0d0,x1)
      call fill_vec(NoVec,0d0,x2)
      call fill_vec(Nt,0d0,x3)

      call Ztx(xd,x1,Nx,Ny,Nz,NoVec,Nb,rx)
      
C----------- Changed!
      norm2=0d0
      
      call cg(rx, rx, rx,E,x1,epsIn,x2,itE,norm1,norm2)
     
C      call cg(imax, jmax, kmax, A, y, tolmax, p, iter, norm1, norm2)
C      write(*,*) 'Intern CG:', itE

C      call band_chol_solve(NoVec,E,x2,x1,rx**2)
C      call band_chol_solve(NoVec,E,y2,y1,rx**2) 

C-------------------
     
      
      

      
      
      
      
      
      
      call AZx(AZ,x2,x3,Nx,Ny,Nz,NoVec,elemAZ)  
      call lin_comb(Nt, xd, x3, -1d0, Px1)
  
      deallocate(x1)
      deallocate(x2) 
      deallocate(x3)  
       
      return
      end      
C-------------------------------------------------------------------------------
C  Case 1
C-------------------------------------------------------------------------------
      subroutine case2_1(leftright,A,AZ,E,Nb,noB,Nx,NxNy,k,j,
     1                 h,vert,pos,elemAZ,Ny,Nz,rx,NoVec)
      implicit none
      
      character leftright
      integer i,x,rx,pos, NoVec,h, NxNy
      integer Nb, j, Nx, Ny, Nz, elemAZ, k, noB, vert
      double precision A(Nx*Ny*Nz,4), AZ(elemAZ,3)
      double precision E(NoVec,4)

      if (leftright .EQ.'r') then
        do i=1,vert
	  x=pos+(j-1)*Nx+(h-1)*NxNy
	  
          AZ(k,1)=x
    	  AZ(k,2)=noB
          AZ(k,3)=-A(x,2)
          E(noB,1)=AZ(k,3)+E(noB,1)
          k=k+1
          x=x+1
	  
	  AZ(k,1)=x
          AZ(k,2)=noB
          AZ(k,3)=A(x-1,2)
          E(noB,2)=AZ(k,3)+E(noB,2)
          j=j+1
          k=k+1
        end do
      elseif (leftright .EQ. 'l') then
        do i=1,vert      
          x=pos+(j-1)*Nx+(h-1)*NxNy
	  
          AZ(k,1)=x
          AZ(k,2)=noB
          AZ(k,3)=A(x,2)
          k=k+1
          x=x+1
    
          AZ(k,1)=x 
          AZ(k,2)=noB 
          AZ(k,3)=-A(x-1,2)
          E(noB,1)=AZ(k,3)+E(noB,1)
          j=j+1
          k=k+1
        end do
      elseif (leftright .EQ. 'm') then
        do i=1,vert
          x=pos+(j-1)*Nx+(h-1)*NxNy
	  
          AZ(k,1)=x 
          AZ(k,2)=noB 
          AZ(k,3)=A(x,2)
          k=k+1
          x=x+1
    
          AZ(k,1)=x 
          AZ(k,2)=noB 
          AZ(k,3)=-A(x-1,2)
          E(noB,1)=AZ(k,3)+E(noB,1)
          k=k+1
          x=x+Nb-1
    
          AZ(k,1)=x 
          AZ(k,2)=noB
          AZ(k,3)=-A(x,2)
          E(noB,1)=AZ(k,3)+E(noB,1)
          k=k+1
          x=x+1
    
          AZ(k,1)=x 
          AZ(k,2)=noB
          AZ(k,3)=A(x-1,2)
          E(noB,2)=AZ(k,3)+E(noB,2) 
          j=j+1
          k=k+1 
        end do 
      end if
 
      return
      end
C-------------------------------------------------------------------------------
C  Case 2
C------------------------------------------------------------------------------- 
      subroutine case2_2(updown,A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,pos,
     1                 elemAZ,Ny,Nz,rx,NoVec)
      implicit none
      
      character updown
      integer Nb, j, Nx, Ny, Nz, elemAZ, k, noB, pos, NoVec,rx
      double precision A(Nx*Ny*Nz,4),AZ(elemAZ,3)
      double precision E(NoVec,4)
      integer i,x
      integer h, NxNy

      if (updown .EQ.'u') then
        do i=1,Nb
          x=pos-1+i+(j-1)*Nx+(h-1)*NxNy 
	  
          AZ(k,1)=x 
          AZ(k,2)=noB 
          AZ(k,3)=A(x-Nx,3) 
          E(noB,3)=AZ(k,3)+E(noB,3) 
          k=k+1 
        end do   
        j=j+1 
      elseif (updown .EQ. 'd') then
        do i=1,Nb
          x=pos-1+i+(j-1)*Nx+(h-1)*NxNy 
	  
          AZ(k,1)=x 
          AZ(k,2)=noB 
          AZ(k,3)=A(x,3) 
          k=k+1 
        end do  
        j=j+1   
      end if
      
      return
      end
C-------------------------------------------------------------------------------
C  Case 3
C------------------------------------------------------------------------------- 
      subroutine case2_3(leftright,updown,A,AZ,E,Nb,noB,Nx,NxNy,k,
     1                 j,h,pos,elemAZ,Ny,Nz,rx,NoVec)
      implicit none
      
      character leftright, updown
      integer Nb, j, Nx, Ny, Nz, elemAZ, k, noB, pos,i
      integer rx,NoVec,x
      double precision A(Nx*Ny*Nz,4), AZ(elemAZ,3)
      double precision E(NoVec,4)
      integer h, NxNy

      if (leftright .EQ. 'l') then
        if (updown .EQ. 'u') then
          x=pos+(j-1)*Nx+(h-1)*NxNy 
    
          do i=1,(Nb-1)
            AZ(k,1)=x 
            AZ(k,2)=noB 
            AZ(k,3)=-A(x,3) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do

          AZ(k,1)=x 
          AZ(k,2)=noB 
          AZ(k,3)=-A(x,3)-A(x,2) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          AZ(k,1)=x 
          AZ(k,2)=noB 
          AZ(k,3)=A(x-1,2) 
          E(noB,2)=AZ(k,3)+E(noB,2) 
          k=k+1 
          j=j+1 
    
        elseif (updown .EQ. 'd') then
          x=pos+(j-1)*Nx+(h-1)*NxNy   
      
          do i=1,(Nb-1)
            AZ(k,1)=x 
            AZ(k,2)=noB 
            AZ(k,3)=-A(x-Nx,3) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do
    
          AZ(k,1)=x 
          AZ(k,2)=noB 
          AZ(k,3)=-A(x-Nx,3)-A(x,2) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 
      
          AZ(k,1)=x 
          AZ(k,2)=noB 
          AZ(k,3)=A(x-1,2) 
          E(noB,2)=AZ(k,3)+E(noB,2) 
          k=k+1 
          j=j+1 
        end if
  
      elseif (leftright .EQ. 'r') then
        if (updown .EQ. 'u') then
          x=pos+(j-1)*Nx+(h-1)*NxNy 
	  
          AZ(k,1)=x 
          AZ(k,2)=noB 
          AZ(k,3)=A(x,2) 
          k=k+1 
          x=x+1 
   
          AZ(k,1)=x 
          AZ(k,2)=noB 
          AZ(k,3)=-A(x,3)-A(x-1,2) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          do i=1,(Nb-1)
            AZ(k,1)=x 
            AZ(k,2)=noB 
            AZ(k,3)=-A(x,3) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do
          j=j+1 
    
        elseif (updown .EQ. 'd') then
          x=pos+(j-1)*Nx+(h-1)*NxNy 
	  
          AZ(k,1)=x 
          AZ(k,2)=noB 
          AZ(k,3)=A(x,2) 
          k=k+1 
          x=x+1 

          AZ(k,1)=x 
          AZ(k,2)=noB 
          AZ(k,3)=-A(x-1,2)-A(x-Nx,3) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 
  
          do i=1,(Nb-1)
            AZ(k,1)=x 
            AZ(k,2)=noB 
            AZ(k,3)=-A(x-Nx,3) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do   
          j=j+1 

        end if
      end if
 
      return
      end 
C-------------------------------------------------------------------------------
C  Case 4
C-------------------------------------------------------------------------------
      subroutine case2_4(updown,A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,pos,
     1                 elemAZ,Ny,Nz,rx,NoVec)
      implicit none
      
      character  updown
      integer Nb, j, Nx, Ny, Nz, elemAZ, k, noB,pos,i
      integer rx,NoVec,x
      double precision A(Nx*Ny*Nz,4),AZ(elemAZ,3),
     1                 E(NoVec,4)  
      integer h, NxNy 

      if (updown .EQ. 'u') then
        x=pos+(j-1)*Nx+(h-1)*NxNy 
	
        AZ(k,1)=x  
        AZ(k,2)=noB  
        AZ(k,3)=A(x,2) 
        k=k+1 
        x=x+1 

        AZ(k,1)=x  
        AZ(k,2)=noB  
        AZ(k,3)=-A(x,3)-A(x-1,2) 
        E(noB,1)=AZ(k,3)+E(noB,1) 
        k=k+1 
        x=x+1 

        do i=1,Nb-2
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x,3) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 
        end do
  
        AZ(k,1)=x  
        AZ(k,2)=noB  
        AZ(k,3)=-A(x,3)-A(x,2) 
        E(noB,1)=AZ(k,3)+E(noB,1) 
        k=k+1 
        x=x+1 
  
        AZ(k,1)=x  
        AZ(k,2)=noB  
        AZ(k,3)=A(x-1,2) 
        E(noB,2)=AZ(k,3)+E(noB,2) 
        k=k+1 
        j=j+1 
      elseif (updown .EQ. 'd') then
        x=pos+(j-1)*Nx+(h-1)*NxNy 
	
        AZ(k,1)=x  
        AZ(k,2)=noB  
        AZ(k,3)=A(x,2) 
        k=k+1 
        x=x+1 
  
        AZ(k,1)=x  
        AZ(k,2)=noB  
        AZ(k,3)=-A(x-1,2)-A(x-Nx,3) 
        E(noB,1)=AZ(k,3)+E(noB,1) 
        k=k+1 
        x=x+1 

        do i=1,Nb-2
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x-Nx,3)   
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 
        end do
  
        AZ(k,1)=x  
        AZ(k,2)=noB  
        AZ(k,3)=-A(x-Nx,3)-A(x,2)   
        E(noB,1)=AZ(k,3)+E(noB,1) 
        k=k+1 
        x=x+1 

        AZ(k,1)=x  
        AZ(k,2)=noB  
        AZ(k,3)=A(x-1,2)   
        E(noB,2)=AZ(k,3)+E(noB,2) 
        k=k+1 
        j=j+1 

      end if
      
      return
      end
C-------------------------------------------------------------------------------
C  Case 5
C-------------------------------------------------------------------------------  
      subroutine case2_5(updown,A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,pos,
     1                 elemAZ,Ny,Nz,rx,NoVec)
      implicit none
      
      character  updown
      integer Nb, j, Nx, Ny, Nz, elemAZ, k, noB,pos,i
      integer rx,NoVec,x,m,pos0
      double precision A(Nx*Ny*Nz,4),AZ(elemAZ,3),
     1                 E(NoVec,4)  
      integer h, NxNy 
      
      pos0=pos
      if (updown .EQ. 'd') then
        
        do m=1,Nb
	  x=pos+(j-1)*Nx+(h-1)*NxNy
          do i=1,Nb
       	    AZ(k,1)=x 
	    AZ(k,2)=noB 
	    AZ(k,3)=A(x-NxNy,4) 
	    E(noB,4)=AZ(k,3)+E(noB,4) 
            k=k+1  
            x=x+1 
          end do
          pos=pos+Nx 
        end do
         
      elseif (updown .EQ. 'u') then
        do m=1,Nb
          x=pos+(j-1)*Nx+(h-1)*NxNy 
	  
          do i=1,Nb
      	    AZ(k,1)=x 
	    AZ(k,2)=noB 
	    AZ(k,3)=A(x,4) 
            k=k+1 
            x=x+1 
          end do
          pos=pos+Nx 
        end do
      end if
      
      pos=pos0
      
      return
      end
C-------------------------------------------------------------------------------
C  Case 6n
C-------------------------------------------------------------------------------  
      subroutine case2_6n(leftright,updown,A,AZ,E,Nb,noB,Nx,NxNy,
     1                  k,j,h,vert,pos,elemAZ,Ny,Nz,rx,NoVec)
      implicit none
      
      character leftright, updown
      integer Nb, j, Nx, Ny, Nz, elemAZ, k, noB,pos,i
      integer rx,NoVec,x,m,vert,m0
      double precision A(Nx*Ny*Nz,4),AZ(elemAZ,3),
     1                 E(NoVec,4)  
      integer h, NxNy 

      if (leftright .EQ. 'l') then
        if (updown .EQ. 'u') then

          do m=1,(vert-1)
   
            x=pos+(j-1)*Nx+(h-1)*NxNy 

            do i=1,(Nb-1)
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x-NxNy,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x,2)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x-1,2) 
            E(noB,2)=AZ(k,3)+E(noB,2) 
            j=j+1 
            k=k+1 
          end do
          x=pos+(j-1)*Nx+(h-1)*NxNy 
	  
          do i=1,(Nb-1)
             AZ(k,1)=x  
             AZ(k,2)=noB  
             AZ(k,3)=-A(x,3)-A(x-NxNy,4) 
             E(noB,1)=AZ(k,3)+E(noB,1) 
             k=k+1 
             x=x+1 
          end do

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x,3)-A(x,2)-A(x-NxNy,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)= A(x-1,2) 
          E(noB,2)=AZ(k,3)+E(noB,2) 
          k=k+1 
          j=j+1 
      
        elseif (updown .EQ. 'd') then

          x=pos+(j-1)*Nx+(h-1)*NxNy 
          do i=1,(Nb-1)
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-Nx,3)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x-Nx,3)-A(x,2)-A(x-NxNy,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=A(x-1,2) 
          E(noB,2)=AZ(k,3)+E(noB,2) 
          k=k+1 
          j=j+1 
          x=x+1 

          do m=1,(vert-1)
            x=pos+(j-1)*Nx+(h-1)*NxNy 
	    
            do i=1,(Nb-1)
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x-NxNy,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x,2)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x-1,2) 
            E(noB,2)=AZ(k,3)+E(noB,2) 
            j=j+1 
            k=k+1 
          end do
    
        elseif (updown .EQ. 'm') then 
          x=pos+(j-1)*Nx+(h-1)*NxNy   
    
          do i=1,(Nb-1)
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-Nx,3)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x-Nx,3)-A(x,2)-A(x-NxNy,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=A(x-1,2) 
          E(noB,2)=AZ(k,3)+E(noB,2) 
          k=k+1 
          j=j+1 
          x=x+1 

          do m=2,(vert-1)
            x=pos+(j-1)*Nx+(h-1)*NxNy  
       
     	    do i=1,(Nb-1)
	      AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x-NxNy,4) 
	      E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x,2)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x-1,2) 
            E(noB,2)=AZ(k,3)+E(noB,2) 
            j=j+1 
            k=k+1  
          end do
          x=pos+(j-1)*Nx+(h-1)*NxNy 
    
          do i=1,(Nb-1)
             AZ(k,1)=x  
             AZ(k,2)=noB  
             AZ(k,3)=-A(x,3)-A(x-NxNy,4) 
             E(noB,1)=AZ(k,3)+E(noB,1) 
             k=k+1 
             x=x+1 
          end do

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x,3)-A(x,2)-A(x-NxNy,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)= A(x-1,2) 
          E(noB,2)=AZ(k,3)+E(noB,2) 
          k=k+1 
          j=j+1 
          x=x+1     
        end if
  
      elseif (leftright .EQ.'r') then
        if (updown .EQ. 'u') then

          do m=j,(j+vert-2)
            x=pos+(m-1)*Nx+(h-1)*NxNy
       
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x,2) 
            k=k+1 
            x=x+1 
       
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-1,2)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1)         
            k=k+1 
            x=x+1 
    
            do i=1,(Nb-1)
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x-NxNy,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do
          end do

          j=m 
          x=pos+(j-1)*Nx+(h-1)*NxNy 
    
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=A(x,2) 
          k=k+1 
          x=x+1 
   
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x,3)-A(x-1,2)-A(x-NxNy,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          do i=1,(Nb-1)
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x,3)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do
          j=j+1 
	  
        elseif (updown .EQ. 'd') then
          x=pos+(j-1)*Nx+(h-1)*NxNy 
    
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)= A(x,2) 
          k=k+1 
          x=x+1 
 
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x-1,2)-A(x-Nx,3)-A(x-NxNy,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 
  
          do i=1,(Nb-1)
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-Nx,3)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do   
          j=j+1 
    
          do m=j,j+vert-2
            x=pos+(m-1)*Nx+(h-1)*NxNy   
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x,2) 
            k=k+1 
            x=x+1 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
C            AZ(k,3)=-A(x-1, 2)-A(x-1-NxNy,4)
	    AZ(k,3)=-A(x-1, 2)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
    
            do i=1,(Nb-1)
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x-NxNy,4) 
	      E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do
          end do
    
        elseif (updown .EQ. 'm') then
    
          x=pos+(j-1)*Nx+(h-1)*NxNy 
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)= A(x,2) 
          k=k+1 
          x=x+1 
 
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x-1,2)-A(x-Nx,3)-A(x-NxNy,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1   
    
          do i=1,(Nb-1)
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-Nx,3)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do  
          j=j+1 
	  m0=j-1
    
          do m=j,(j+vert-3)
            x=pos+(m-1)*Nx+(h-1)*NxNy
       
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x,2) 
            k=k+1 
            x=x+1 
       
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-1,2)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1)         
            k=k+1 
            x=x+1 
    
            do i=1,(Nb-1)
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x-NxNy,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do

	    m0=m
          end do
          
	  m=m0
          j=m+1
	  
          x=pos+(j-1)*Nx+(h-1)*NxNy 
    
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=A(x,2) 
          k=k+1 
          x=x+1 
   
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x,3)-A(x-1,2)-A(x-NxNy,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          do i=1,(Nb-1)
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x,3)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do
          j=j+1  
   
        end if
  
      elseif (leftright .EQ.'m') then
          if (updown .EQ. 'u') then
         
            do m=1,(vert-1)
               x=pos+(m-1)*Nx+(h-1)*NxNy 
	 
               AZ(k,1)=x  
               AZ(k,2)=noB  
               AZ(k,3)=A(x,2) 
               k=k+1 
               x=x+1 
          
               AZ(k,1)=x  
               AZ(k,2)=noB  
               AZ(k,3)=-A(x-1,2)-A(x-NxNy,4) 
               E(noB,1)=AZ(k,3)+E(noB,1)         
               k=k+1 
               x=x+1  
          
               do i=1,(Nb-2)
	          AZ(k,1)=x  
                  AZ(k,2)=noB  
                  AZ(k,3)=-A(x-NxNy,4) 
	          E(noB,1)=AZ(k,3)+E(noB,1) 
                  k=k+1 
                  x=x+1 
               end do

               AZ(k,1)=x  
               AZ(k,2)=noB  
               AZ(k,3)=-A(x,2)-A(x-NxNy,4) 
               E(noB,1)=AZ(k,3)+E(noB,1) 
               k=k+1 
               x=x+1 
      
               AZ(k,1)=x  
               AZ(k,2)=noB  
               AZ(k,3)=A(x-1,2) 
               E(noB,2)=AZ(k,3)+E(noB,2) 
               j=j+1 
               k=k+1 
             end do
             x=pos+(j-1)*Nx+(h-1)*NxNy 
       
             AZ(k,1)=x  
             AZ(k,2)=noB  
             AZ(k,3)=A(x,2) 
             k=k+1 
             x=x+1 

             AZ(k,1)=x  
             AZ(k,2)=noB  
             AZ(k,3)=-A(x,3)-A(x-NxNy,4)-A(x-1,2) 
             E(noB,1)=AZ(k,3)+E(noB,1) 
             k=k+1   
             x=x+1 
       
             do i=3,Nb
                AZ(k,1)=x  
                AZ(k,2)=noB  
                AZ(k,3)=-A(x,3)-A(x-NxNy,4) 
                E(noB,1)=AZ(k,3)+E(noB,1) 
                k=k+1 
                x=x+1 
             end do
       
             AZ(k,1)=x  
             AZ(k,2)=noB  
             AZ(k,3)= -A(x,3)-A(x,2)-A(x-NxNy,4) 
             E(noB,1)=AZ(k,3)+E(noB,1) 
             k=k+1 
             x=x+1 

             AZ(k,1)=x  
             AZ(k,2)=noB  
             AZ(k,3)= A(x-1,2) 
             E(noB,2)=AZ(k,3)+E(noB,2) 
             k=k+1 
             j=j+1 
       
          elseif (updown .EQ. 'd') then
            x=pos+(j-1)*Nx+(h-1)*NxNy
             
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)= A(x,2) 
            k=k+1 
            x=x+1 
 
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-1,2)-A(x-Nx,3)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 

            do i=1,(Nb-2)
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x-Nx,3)-A(x-NxNy,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-Nx,3)-A(x,2)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x-1,2) 
            E(noB,2)=AZ(k,3)+E(noB,2) 
            k=k+1 
            j=j+1 
            x=x+1 

            do m=j,(j+vert-2)
              x=pos+(m-1)*Nx+(h-1)*NxNy  
	 
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=A(x,2) 
              k=k+1 
              x=x+1 
      
              AZ(k,1)=x  
              AZ(k,2)=noB  
C              AZ(k,3)=-A(x-1,2)-A(x-1-NxNy,4) 
              AZ(k,3)=-A(x-1,2)-A(x-NxNy,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
    
              do i=1,(Nb-2)
	        AZ(k,1)=x  
                AZ(k,2)=noB  
                AZ(k,3)=-A(x-NxNy,4) 
            	E(noB,1)=AZ(k,3)+E(noB,1) 
                k=k+1 
                x=x+1 
             end do

              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x,2)-A(x-NxNy,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
      
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=A(x-1,2) 
              E(noB,2)=AZ(k,3)+E(noB,2) 
              j=j+1 
              k=k+1 
            end do
    
        
          elseif (updown .EQ. 'm') then
            x=pos+(j-1)*Nx+(h-1)*NxNy 
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x,2) 
            k=k+1 
            x=x+1 

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-1,2)-A(x-Nx,3)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
      
            do i=3,Nb
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x-Nx,3)-A(x-NxNy,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-Nx,3)-A(x,2)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x-1,2) 
            E(noB,2)=AZ(k,3)+E(noB,2) 
            k=k+1 
            j=j+1 
            x=x+1 

            do m=j,(j+vert-3)
              x=pos+(m-1)*Nx+(h-1)*NxNy 
	
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=A(x,2) 
              k=k+1 
              x=x+1 
        
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x-1,2)-A(x-NxNy,4) 
              E(noB,1)=AZ(k,3)+E(noB,1)         
              k=k+1 
              x=x+1 
          
              do i=3,Nb
	        AZ(k,1)=x  
                AZ(k,2)=noB  
                AZ(k,3)=-A(x-NxNy,4) 
	        E(noB,1)=AZ(k,3)+E(noB,1) 
                k=k+1 
                x=x+1 
              end do

              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x,2)-A(x-NxNy,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
      
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=A(x-1,2) 
              E(noB,2)=AZ(k,3)+E(noB,2) 
              j=j+1 
              k=k+1 
            end do
            x=pos+(j-1)*Nx+(h-1)*NxNy 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x,2) 
            k=k+1 
            x=x+1 

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-1,2)-A(x,3)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
      
            do i=3,Nb
               AZ(k,1)=x  
               AZ(k,2)=noB  
               AZ(k,3)=-A(x,3)-A(x-NxNy,4) 
               E(noB,1)=AZ(k,3)+E(noB,1) 
               k=k+1 
               x=x+1 
            end do

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x,3)-A(x,2)-A(x-NxNy,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)= A(x-1,2) 
            E(noB,2)=AZ(k,3)+E(noB,2) 
            k=k+1 
            j=j+1 
          end if
      end if 
      
      return
      end
C-------------------------------------------------------------------------------
C  Case 6p
C-------------------------------------------------------------------------------  
      subroutine case2_6p(leftright,updown,A,AZ,E,Nb,noB,Nx,NxNy,
     1                  k,j,h,vert,pos,elemAZ,Ny,Nz,rx,NoVec)
      implicit none
      
      character  leftright, updown
      integer Nb, j, Nx, Ny, Nz, elemAZ, k, noB,pos,i
      integer rx,NoVec,x,m, vert,m0
      double precision A(Nx*Ny*Nz,4),AZ(elemAZ,3),
     1                 E(NoVec,4)  
      integer h, NxNy 

      if (leftright .EQ. 'l') then
        if (updown .EQ. 'u') then
    
          do m=1,(vert-1)
            x=pos+(j-1)*Nx+(h-1)*NxNy 
      
            do i=1,(Nb-1)
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x,4) 
	      E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x,2)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x-1,2) 
            E(noB,2)=AZ(k,3)+E(noB,2) 
            j=j+1 
            k=k+1  
          end do
          x=pos+(j-1)*Nx+(h-1)*NxNy 
    
          do i=1,(Nb-1)
             AZ(k,1)=x  
             AZ(k,2)=noB  
             AZ(k,3)=-A(x,3)-A(x,4) 
             E(noB,1)=AZ(k,3)+E(noB,1) 
             k=k+1 
             x=x+1 
          end do

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x,3)-A(x,2)-A(x,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=A(x-1,2) 
          E(noB,2)=AZ(k,3)+E(noB,2) 
          k=k+1 
          j=j+1 
       
        elseif (updown .EQ. 'd') then
          x=pos+(j-1)*Nx+(h-1)*NxNy 
    
          do i=1,(Nb-1)
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-Nx,3)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x-Nx,3)-A(x,2)-A(x,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=A(x-1,2) 
          E(noB,2)=AZ(k,3)+E(noB,2) 
          k=k+1 
          j=j+1 
          x=x+1 

          do m=1,(vert-1)
            x=pos+(j-1)*Nx+(h-1)*NxNy 
      
            do i=1,(Nb-1)
	      AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x,4) 
	      E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x,2)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x-1,2) 
            E(noB,2)=AZ(k,3)+E(noB,2) 
            j=j+1 
            k=k+1 
          end do
    
        elseif (updown .EQ. 'm') then
          x=pos+(j-1)*Nx+(h-1)*NxNy  
     
          do i=1,(Nb-1)
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-Nx,3)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x-Nx,3)-A(x,2)-A(x,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=A(x-1,2) 
          E(noB,2)=AZ(k,3)+E(noB,2) 
          k=k+1 
          j=j+1 
          x=x+1 
      
          do m=2,(vert-1)
            x=pos+(j-1)*Nx+(h-1)*NxNy 
        
            do i=1,(Nb-1)
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x,2)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x-1,2) 
            E(noB,2)=AZ(k,3)+E(noB,2) 
            j=j+1 
            k=k+1 
          end do
          x=pos+(j-1)*Nx+(h-1)*NxNy 
    
          do i=1,(Nb-1)
             AZ(k,1)=x  
             AZ(k,2)=noB  
             AZ(k,3)=-A(x,3)-A(x,4) 
             E(noB,1)=AZ(k,3)+E(noB,1) 
             k=k+1 
             x=x+1 
          end do

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x,3)-A(x,2)-A(x,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)= A(x-1,2) 
          E(noB,2)=AZ(k,3)+E(noB,2) 
          k=k+1 
          j=j+1 
          x=x+1  
        end if
  
      elseif (leftright .EQ. 'r') then
        if (updown .EQ. 'u') then
    
          do m=j,(j+vert-2)
            x=pos+(m-1)*Nx+(h-1)*NxNy 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x,2) 
            k=k+1 
            x=x+1 
       
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-1,2)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1)         
            k=k+1 
            x=x+1 
    
            do i=1,(Nb-1)
	      AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x,4) 
	      E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do
          end do
          j=m
          x=pos+(j-1)*Nx+(h-1)*NxNy 
	  
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=A(x,2) 
          k=k+1 
          x=x+1 
   
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x,3)-A(x-1,2)-A(x,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          do i=1,(Nb-1)
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x,3)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do
          j=j+1 
    
        elseif (updown .EQ. 'd') then
          x=pos+(j-1)*Nx+(h-1)*NxNy
     
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)= A(x,2) 
          k=k+1 
          x=x+1 
 
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x-1,2)-A(x-Nx,3)-A(x,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 
  
          do i=1,(Nb-1)
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-Nx,3)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do  
          j=j+1 
    
          do m=j,(j+vert-2)
            x=pos+(m-1)*Nx+(h-1)*NxNy 
        
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x,2) 
            k=k+1 
            x=x+1 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
C            AZ(k,3)=-A(x-1,2)-A(x-1,4) 
            AZ(k,3)=-A(x-1,2)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
    
            do i=1,(Nb-1)
	      AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x,4) 
	      E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do

          end do
    
        elseif (updown .EQ. 'm') then
          x=pos+(j-1)*Nx+(h-1)*NxNy 
    
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)= A(x,2) 
          k=k+1 
          x=x+1 
 
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x-1,2)-A(x-Nx,3)-A(x,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1   
    
          do i=1,(Nb-1)
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-Nx,3)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do   
          j=j+1 
	  
C Extra line	%same line as in case2_6n 
	  m0=j-1	  
	  
C          write(*,*) j
C          write(*,*) vert
          do m=j,(j+vert-3)
            x=pos+(m-1)*Nx+(h-1)*NxNy 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x,2) 
            k=k+1 
            x=x+1 
       
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-1,2)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1)         
            k=k+1 
            x=x+1 
    
            do i=1,(Nb-1)
	      AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do
	    
	    m0=m
          end do
          
	  m=m0
          j=m+1

          x=pos+(j-1)*Nx+(h-1)*NxNy 

C          write(*,*) h
C          write(*,*) pos
C          write(*,*) m0
C          write(*,*) m
C          write(*,*) j
C          write(*,*) x, A(x,2)
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=A(x,2) 
          k=k+1 
          x=x+1 
   
          AZ(k,1)=x  
          AZ(k,2)=noB  
          AZ(k,3)=-A(x,3)-A(x-1,2)-A(x,4) 
          E(noB,1)=AZ(k,3)+E(noB,1) 
          k=k+1 
          x=x+1 

          do i=1,(Nb-1)
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x,3)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
          end do
          j=j+1  
   
        end if
  
      elseif (leftright .EQ.'m') then
          if (updown .EQ. 'u') then
            do m=1,(vert-1)
               x=pos+(m-1)*Nx+(h-1)*NxNy 
	 
               AZ(k,1)=x  
               AZ(k,2)=noB  
               AZ(k,3)=A(x,2) 
               k=k+1 
               x=x+1 
          
               AZ(k,1)=x  
               AZ(k,2)=noB  
               AZ(k,3)=-A(x-1,2)-A(x,4) 
               E(noB,1)=AZ(k,3)+E(noB,1)         
               k=k+1 
               x=x+1  
                
               do i=1,(Nb-2)
	          AZ(k,1)=x  
                  AZ(k,2)=noB  
                  AZ(k,3)=-A(x,4) 
	          E(noB,1)=AZ(k,3)+E(noB,1) 
                  k=k+1 
                  x=x+1 
               end do

               AZ(k,1)=x  
               AZ(k,2)=noB  
               AZ(k,3)=-A(x,2)-A(x,4) 
               E(noB,1)=AZ(k,3)+E(noB,1) 
               k=k+1 
               x=x+1 
      
               AZ(k,1)=x  
               AZ(k,2)=noB  
               AZ(k,3)=A(x-1,2) 
               E(noB,2)=AZ(k,3)+E(noB,2) 
               j=j+1 
               k=k+1 
             end do
             x=pos+(j-1)*Nx+(h-1)*NxNy 
       
             AZ(k,1)=x  
             AZ(k,2)=noB  
             AZ(k,3)=A(x,2) 
             k=k+1 
             x=x+1 
      
             AZ(k,1)=x  
             AZ(k,2)=noB  
             AZ(k,3)=-A(x,3)-A(x,4)-A(x-1,2) 
             E(noB,1)=AZ(k,3)+E(noB,1) 
             k=k+1   
             x=x+1 
       
             do i=3,Nb
                AZ(k,1)=x  
                AZ(k,2)=noB  
                AZ(k,3)=-A(x,3)-A(x,4) 
                E(noB,1)=AZ(k,3)+E(noB,1) 
                k=k+1 
                x=x+1 
             end do
       
             AZ(k,1)=x  
             AZ(k,2)=noB  
             AZ(k,3)= -A(x,3)-A(x,2)-A(x,4) 
             E(noB,1)=AZ(k,3)+E(noB,1) 
             k=k+1 
             x=x+1 

             AZ(k,1)=x  
             AZ(k,2)=noB  
             AZ(k,3)= A(x-1,2) 
             E(noB,2)=AZ(k,3)+E(noB,2) 
             k=k+1 
             j=j+1 
          elseif (updown .EQ. 'd') then
            x=pos+(j-1)*Nx+(h-1)*NxNy 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)= A(x,2) 
            k=k+1 
            x=x+1 
 
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-1,2)-A(x-Nx,3)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 

            do i=1,(Nb-2)
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x-Nx,3)-A(x,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-Nx,3)-A(x,2)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x-1,2) 
            E(noB,2)=AZ(k,3)+E(noB,2) 
            k=k+1 
            j=j+1 
            x=x+1 

            do m=j,(j+vert-2)
              x=pos+(m-1)*Nx+(h-1)*NxNy 
	  
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=A(x,2) 
              k=k+1 
              x=x+1 
      
              AZ(k,1)=x  
              AZ(k,2)=noB  
C              AZ(k,3)=-A(x-1,2)-A(x-1,4)
              AZ(k,3)=-A(x-1,2)-A(x,4)
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
    
              do i=1,(Nb-2)
	        AZ(k,1)=x  
                AZ(k,2)=noB  
                AZ(k,3)=-A(x,4) 
	        E(noB,1)=AZ(k,3)+E(noB,1) 
                k=k+1 
                x=x+1 
              end do

              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x,2)-A(x,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
      
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=A(x-1,2) 
              E(noB,2)=AZ(k,3)+E(noB,2) 
              j=j+1 
              k=k+1 
            end do
    
          elseif (updown .EQ. 'm') then
            x=pos+(j-1)*Nx+(h-1)*NxNy 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)= A(x,2) 
            k=k+1 
            x=x+1 

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-1,2)-A(x-Nx,3)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
      
            do i=3,Nb
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x-Nx,3)-A(x,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            end do

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-Nx,3)-A(x,2)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x-1,2) 
            E(noB,2)=AZ(k,3)+E(noB,2) 
            k=k+1 
            j=j+1 
            x=x+1 

            do m=j,(j+vert-3)
              x=pos+(m-1)*Nx+(h-1)*NxNy 
	
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=A(x,2) 
              k=k+1 
              x=x+1 
        
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x-1,2)-A(x,4) 
              E(noB,1)=AZ(k,3)+E(noB,1)         
              k=k+1 
              x=x+1 
          
              do i=3,Nb
	        AZ(k,1)=x  
                AZ(k,2)=noB  
                AZ(k,3)=-A(x,4) 
	        E(noB,1)=AZ(k,3)+E(noB,1) 
                k=k+1 
                x=x+1 
              end do

              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=-A(x,2)-A(x,4) 
              E(noB,1)=AZ(k,3)+E(noB,1) 
              k=k+1 
              x=x+1 
            
              AZ(k,1)=x  
              AZ(k,2)=noB  
              AZ(k,3)=A(x-1,2) 
              E(noB,2)=AZ(k,3)+E(noB,2) 
              j=j+1 
              k=k+1 
            end do
            x=pos+(j-1)*Nx+(h-1)*NxNy 
      
            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=A(x,2) 
            k=k+1 
            x=x+1 

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x-1,2)-A(x,3)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 
      
            do i=3,Nb
               AZ(k,1)=x  
               AZ(k,2)=noB  
               AZ(k,3)=-A(x,3)-A(x,4) 
               E(noB,1)=AZ(k,3)+E(noB,1) 
               k=k+1 
               x=x+1 
            end do

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)=-A(x,3)-A(x,2)-A(x,4) 
            E(noB,1)=AZ(k,3)+E(noB,1) 
            k=k+1 
            x=x+1 

            AZ(k,1)=x  
            AZ(k,2)=noB  
            AZ(k,3)= A(x-1,2) 
            E(noB,2)=AZ(k,3)+E(noB,2) 
            k=k+1 
            j=j+1 
    
          end if
      end if
      
      return
      end
C-------------------------------------------------------------------------------
C  Compute AZ and E efficiently (using Cases 1-6)
C-------------------------------------------------------------------------------
      subroutine computeAZandE2(NoVec,Nx,Ny,Nz,A,AZ,E,elemAZ,Nb, 
     1                     rx) 
      implicit none
      
      integer NoVec, Nx, Ny, Nz, elemAZ, noB, k,j, Nb, p,rx,q,i
      integer h, xpos,r, NxNy
      double precision A(Nx*Ny*Nz,4),AZ(elemAZ,3)
      double precision E(NoVec,4) 
      
      NxNy=Nx*Ny
      k=1 
C =========================================================================
      noB=1 
      xpos=1 

      do h=1,(Nb-1)
        j=1 
        call case2_1('r',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-1,
     1             xpos+Nb-1,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('l','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1             xpos,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1             elemAZ,Ny,Nz,rx,NoVec) 
      end do

      h=Nb 
      j=1 
      call case2_6p('l','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1            xpos,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1            elemAZ,Ny,Nz,rx,NoVec) 

      h=Nb+1 
      j=1 
      call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1           elemAZ,Ny,Nz,rx,NoVec) 

      noB=noB+1 
C=========================================================================
C noB=2 
      do p=2,(rx-1)
          xpos=(p-1)*Nb 
          do h=1,(Nb-1)
             j=1 
             call case2_1('m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-1,
     1                  xpos,elemAZ,Ny,Nz,rx,NoVec) 
             call case2_4('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1                  elemAZ,Ny,Nz,rx,NoVec) 
             call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1                  elemAZ,Ny,Nz,rx,NoVec) 
          end do
    
          xpos=(p-1)*Nb 
    

          h=Nb 
          j=1 
          call case2_6p('m','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1                xpos,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1                elemAZ,Ny,Nz,rx,NoVec) 

          h=Nb+1 
          j=1 
          call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1               elemAZ,Ny,Nz,rx,NoVec) 
    
          noB=noB+1 
      end do
C =========================================================================
C noB=3 
      xpos=(rx-1)*Nb 

      do h=1,(Nb-1)
        j=1 
        call case2_1('l',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-1,xpos,
     1             elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('r','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1             elemAZ,Ny,Nz,rx,NoVec) 
        call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1             elemAZ,Ny,Nz,rx,NoVec) 
      end do

      h=Nb 
      j=1 
      call case2_6p('r','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,xpos,
     1            elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,elemAZ,
     1            Ny,Nz,rx,NoVec) 

      h=Nb+1 
      j=1 
      call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,elemAZ,
     1           Ny,Nz,rx,NoVec) 

      noB=noB+1 
C =========================================================================
      do q=2,(rx-1)
C =========================================================================
C noB=4 
      xpos=1 

      do h=1,(Nb-1) 
        j=(q-1)*Nb 
        call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1             elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('l','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1             elemAZ,Ny,Nz,rx,NoVec) 
        call case2_1('r',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-2,
     1             xpos+Nb-1,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('l','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1             elemAZ,Ny,Nz,rx,NoVec) 
        call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1             elemAZ,Ny,Nz,rx,NoVec) 
      end do

      h=Nb 
      j=(q-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1           elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6p('l','m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,xpos,
     1            elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,elemAZ,
     1            Ny,Nz,rx,NoVec) 

      h=Nb+1 
      j=(q-1)*Nb+1 
      call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,elemAZ,
     1            Ny,Nz,rx,NoVec) 

      noB=noB+1 
C =========================================================================
C noB=5 
      do p=2,(rx-1)
          do h=1,(Nb-1)
              xpos=(p-1)*Nb 
              j=(q-1)*Nb   
              call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1                    elemAZ,Ny,Nz,rx,NoVec) 
              call case2_4('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1                    elemAZ,Ny,Nz,rx,NoVec) 
              call case2_1('m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-2,
     1                    xpos,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_4('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1                    elemAZ,Ny,Nz,rx,NoVec) 
              call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1                    elemAZ,Ny,Nz,rx,NoVec) 
          end do
    
          xpos=(p-1)*Nb 
    
          h=Nb 
          j=(q-1)*Nb 
          call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1               elemAZ,Ny,Nz,rx,NoVec) 
          call case2_6p('m','m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,xpos,
     1               elemAZ,Ny,Nz,rx,NoVec) 
          call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1               elemAZ,Ny,Nz,rx,NoVec) 

          h=Nb+1 
          j=(q-1)*Nb+1 
          call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1               elemAZ,Ny,Nz,rx,NoVec) 

          noB=noB+1 
      end do
C =========================================================================
C noB=6 
      xpos=(rx-1)*Nb 

      do h=1,(Nb-1)
          j=(q-1)*Nb 
          call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1               elemAZ,Ny,Nz,rx,NoVec) 
          call case2_3('r','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1               elemAZ,Ny,Nz,rx,NoVec) 
          call case2_1('l',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-2,xpos,
     1               elemAZ,Ny,Nz,rx,NoVec) 
          call case2_3('r','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1               elemAZ,Ny,Nz,rx,NoVec) 
          call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1               elemAZ,Ny,Nz,rx,NoVec) 
      end do

      h=Nb 
      j=(q-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,elemAZ,
     1          Ny,Nz,rx,NoVec) 
C      write(*,*) A(458,2)
C      write(*,*) Nx,Ny,Nz
      call case2_6p('r','m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,elemAZ,
     1          Ny,Nz,rx,NoVec) 

      h=Nb+1 
      j=(q-1)*Nb+1 
      call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,elemAZ,
     1          Ny,Nz,rx,NoVec) 

      noB=noB+1 
C =========================================================================
      end do
C =========================================================================
C noB=7 
      xpos=1 

      do h=1,(Nb-1)
        j=(rx-1)*Nb 
        call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,elemAZ,
     1          Ny,Nz,rx,NoVec) 
        call case2_3('l','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
        call case2_1('r',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-1,
     1          xpos+Nb-1,elemAZ,Ny,Nz,rx,NoVec) 
      end do

      h=Nb 
      j=(rx-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,elemAZ,
     1          Ny,Nz,rx,NoVec) 
      call case2_6p('l','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      
      h=Nb+1 
      j=(rx-1)*Nb+1 
      call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,elemAZ,
     1          Ny,Nz,rx,NoVec) 

      noB=noB+1 
C =========================================================================
C noB=8 
      do p=2,(rx-1)
          xpos=(p-1)*Nb 
    
          do h=1,(Nb-1)
              j=(rx-1)*Nb 
              call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1                    elemAZ,Ny,Nz,rx,NoVec) 
              call case2_4('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1                    elemAZ,Ny,Nz,rx,NoVec) 
              call case2_1('m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-1,
     1                    xpos,elemAZ,Ny,Nz,rx,NoVec) 
          end do

          xpos=(p-1)*Nb 

          h=Nb 
          j=(rx-1)*Nb 
          call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1                  elemAZ,Ny,Nz,rx,NoVec) 
          call case2_6p('m','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,xpos,
     1                  elemAZ,Ny,Nz,rx,NoVec) 

          h=Nb+1 
          j=(rx-1)*Nb+1 
          call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1               elemAZ,Ny,Nz,rx,NoVec) 

          noB=noB+1 
      end do
C =========================================================================
C noB=9 
      xpos=(rx-1)*Nb 

      do h=1,(Nb-1)
        j=(rx-1)*Nb 
        call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('r','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
        call case2_1('l',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-1,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      end do

      h=Nb 
      j=(rx-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6p('r','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 

      h=Nb+1 
      j=(rx-1)*Nb+1 
      call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      noB=noB+1 
C =========================================================================
C =========================================================================
C =========================================================================
C =========================================================================
C =========================================================================
C =========================================================================
      do r=2,(rx-1)
C =========================================================================
C noB=10 
      xpos=1 

      h=(r-1)*Nb 
      j=1 
      call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      h=(r-1)*Nb+1 
      j=1 
      call case2_6n('l','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      do h=((r-1)*Nb+2),(r*Nb-1)
        j=1 
        call case2_1('r',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-1,
     1          xpos+Nb-1,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('l','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      end do

      h=r*Nb 
      j=1 
      call case2_6p('l','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      h=r*Nb+1 
      j=1 
      call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      noB=noB+1 
C =========================================================================
C noB=11 
      do p=2,(rx-1)
          xpos=(p-1)*Nb 

          h=(r-1)*Nb 
          j=1 
          call  case5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 

          h=(r-1)*Nb+1 
          j=1 
          call case2_6n('m','u',A,AZ,E,Nb,noB,Nx,NxNy,
     1          k,j,h,Nb,xpos,elemAZ,Ny,Nz,rx,NoVec) 

     
          call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 

          do h=((r-1)*Nb+2),(r*Nb-1)
             j=1 
             call case2_1('m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb-1,xpos,elemAZ,Ny,Nz,rx,NoVec) 
             call case2_4('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
             call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
          end do

          h=r*Nb 
          j=1 
          call case2_6p('m','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

          h=r*Nb+1 
          j=1 
          call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

          noB=noB+1 
      end do
C =========================================================================
C noB=12 
      xpos=(rx-1)*Nb 

      h=(r-1)*Nb 
      j=1 
      call  case5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      h=(r-1)*Nb+1 
      j=1 
      call case2_6n('r','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      do h=((r-1)*Nb+2),(r*Nb-1)
        j=1 
        call case2_1('l',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-1,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('r','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
        call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      end do

      h=r*Nb 
      j=1 
      call case2_6p('r','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      h=r*Nb+1 
      j=1 
      call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      noB=noB+1 
C =========================================================================
      do q=2,(rx-1)
C =========================================================================
C noB=13 
      xpos=1 

      h=(r-1)*Nb 
      j=(q-1)*Nb+1 
      call  case5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      h=(r-1)*Nb+1 
      j=(q-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6n('l','m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      do h=((r-1)*Nb+2),(r*Nb-1)
        j=(q-1)*Nb 
        call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('l','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
        call case2_1('r',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-2,
     1          xpos+Nb-1,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('l','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
        call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      end do

      h=r*Nb 
      j=(q-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6p('l','m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,elemAZ,
     1          Ny,Nz,rx,NoVec) 

      h=r*Nb+1 
      j=(q-1)*Nb+1 
      call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      noB=noB+1 
C =========================================================================
C noB=14 
      do p=2,(rx-1)
          xpos=(p-1)*Nb 

          h=(r-1)*Nb 
          j=(q-1)*Nb+1 
          call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

          h=(r-1)*Nb+1 
          j=(q-1)*Nb 
          call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 
          call case2_6n('m','m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

          do h=((r-1)*Nb+2),(r*Nb-1)
              j=(q-1)*Nb   
              call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_4('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_1('m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-2,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_4('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
              call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 
          end do

          h=r*Nb 
          j=(q-1)*Nb 
          call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 
          call case2_6p('m','m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

          h=r*Nb+1 
          j=(q-1)*Nb+1 
          call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

          noB=noB+1 
      end do
C =========================================================================
C noB=15 
      xpos=(rx-1)*Nb 

      h=(r-1)*Nb 
      j=(q-1)*Nb+1 
      call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      h=(r-1)*Nb+1 
      j=(q-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6n('r','m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      do h=((r-1)*Nb+2),(r*Nb-1)
          j=(q-1)*Nb 
          call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_3('r','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_1('l',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-2,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_3('r','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
      end do

      h=r*Nb 
      j=(q-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6p('r','m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      h=r*Nb+1 
      j=(q-1)*Nb+1 
      call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      noB=noB+1 
C =========================================================================
      end do
C =========================================================================
C noB=16 
      xpos=1 

      h=(r-1)*Nb 
      j=(rx-1)*Nb+1 
      call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      h=(r-1)*Nb+1 
      j=(rx-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6n('l','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 

      do h=((r-1)*Nb+2),(r*Nb-1)
        j=(rx-1)*Nb 
        call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('l','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_1('r',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-1,
     1          xpos+Nb-1,elemAZ,Ny,Nz,rx,NoVec) 
      end do
      
      h=r*Nb 
      j=(rx-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6p('l','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 

      h=r*Nb+1 
      j=(rx-1)*Nb+1 
      call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      noB=noB+1 
C =========================================================================
C noB=17 
      do p=2,rx-1
          xpos=(p-1)*Nb 

          h=(r-1)*Nb 
          j=(rx-1)*Nb+1 
          call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 

          h=(r-1)*Nb+1 
          j=(rx-1)*Nb 
          call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_6n('m','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,
     1          h,Nb,xpos,elemAZ,Ny,Nz,rx,NoVec) 

          do h=((r-1)*Nb+2),(r*Nb-1)
              j=(rx-1)*Nb 
              call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_4('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_1('m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb-1,xpos,elemAZ,Ny,Nz,rx,NoVec) 
          end do

          h=r*Nb 
          j=(rx-1)*Nb 
          call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 
          call case2_6p('m','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 

          h=r*Nb+1 
          j=(rx-1)*Nb+1 
          call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

          noB=noB+1 
      end do
C =========================================================================
C noB=18 
      xpos=(rx-1)*Nb 

      h=(r-1)*Nb 
      j=(rx-1)*Nb+1 
      call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      h=(r-1)*Nb+1 
      j=(rx-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6n('r','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 

      do h=((r-1)*Nb+2),(r*Nb-1)
        j=(rx-1)*Nb 
        call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('r','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_1('l',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-1,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
      end do

      h=r*Nb 
      j=(rx-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6p('r','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 

      h=r*Nb+1 
      j=(rx-1)*Nb+1 
      call case2_5('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos+1,
     1          elemAZ,Ny,Nz,rx,NoVec) 

      noB=noB+1 
C =========================================================================
      end do
C =========================================================================
C =========================================================================
C =========================================================================
C =========================================================================
C =========================================================================
C =========================================================================
C =========================================================================
C noB=19 
      xpos=1 

      h=(rx-1)*Nb 
      j=1 
      call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 

      h=(rx-1)*Nb+1 
      j=1 
      call case2_6n('l','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb,xpos,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 

      do h=((rx-1)*Nb+2),(rx*Nb)
        j=1 
        call case2_1('r',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-1,
     1          xpos+Nb-1,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('l','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,xpos,
     1          elemAZ,Ny,Nz,rx,NoVec) 
      end do

      noB=noB+1 
C =========================================================================
C noB=20 
      do p=2,(rx-1)
          xpos=(p-1)*Nb 
    
          h=(rx-1)*Nb 
          j=1 
          call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 

          h=(rx-1)*Nb+1 
          j=1 
          call case2_6n('m','u',A,AZ,E,Nb,noB,Nx,NxNy,k,
     1          j,h,Nb,xpos,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 

          do h=((rx-1)*Nb+2),(rx*Nb)
              j=1 
              call case2_1('m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb-1,xpos,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_4('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
          end do
      
          noB=noB+1 
      end do
C =========================================================================
C noB=21 
      xpos=(rx-1)*Nb 

      h=(rx-1)*Nb 
      j=1 
      call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 

      h=(rx-1)*Nb+1 
      j=1 
      call case2_6n('r','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb,xpos,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 

      do h=((rx-1)*Nb+2),(rx*Nb)
        j=1 
        call case2_1('l',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb-1,xpos,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('r','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
      end do

      noB=noB+1 
C =========================================================================
      do q=2,(rx-1)
C =========================================================================
C noB=22 
      xpos=1 

      h=(rx-1)*Nb 
      j=(q-1)*Nb+1 
      call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 

      h=(rx-1)*Nb+1 
      j=(q-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6n('l','m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb,xpos,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 


      do h=((rx-1)*Nb+2),(rx*Nb)
        j=(q-1)*Nb 
        call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('l','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_1('r',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb-2,xpos+Nb-1,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('l','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
      end do

      noB=noB+1 
C =========================================================================
C noB=23 
      do p=2,(rx-1)
          xpos=(p-1)*Nb 

          h=(rx-1)*Nb 
          j=(q-1)*Nb+1 
          call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 

          h=(rx-1)*Nb+1 
          j=(q-1)*Nb 
          call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_6n('m','m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb,xpos,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 


          do h=((rx-1)*Nb+2),(rx*Nb)
              j=(q-1)*Nb   
              call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_4('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_1('m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb-2,xpos,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_4('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
          end do

          noB=noB+1 
      end do
C =========================================================================
C noB=24 
      xpos=(rx-1)*Nb 

      h=(rx-1)*Nb 
      j=(q-1)*Nb+1 
      call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 

      h=(rx-1)*Nb+1 
      j=(q-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6n('r','m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb,xpos,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 

      do h=((rx-1)*Nb+2),(rx*Nb)
          j=(q-1)*Nb 
          call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_3('r','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_1('l',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb-2,xpos,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_3('r','u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_2('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
      end do

      noB=noB+1 
C =========================================================================
      end do
C =========================================================================
C noB=25 
      xpos=1 

      h=(rx-1)*Nb 
      j=(rx-1)*Nb+1 
      call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 

      h=(rx-1)*Nb+1 
      j=(rx-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6n('l','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb,xpos,elemAZ,Ny,Nz,rx,NoVec) 

      do h=((rx-1)*Nb+2),(rx*Nb)
        j=(rx-1)*Nb 
        call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('l','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_1('r',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb-1,xpos+Nb-1,elemAZ,Ny,Nz,rx,NoVec) 
      end do

      noB=noB+1 
C =========================================================================
C noB=26 
      do p=2,(rx-1)
          xpos=(p-1)*Nb 

          h=(rx-1)*Nb 
          j=(rx-1)*Nb+1 
          call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 

          h=(rx-1)*Nb+1 
          j=(rx-1)*Nb 
          call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
          call case2_6n('m','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb,xpos,elemAZ,Ny,Nz,rx,NoVec) 

          do h=((rx-1)*Nb+2),(rx*Nb)
              j=(rx-1)*Nb 
              call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_4('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
              call case2_1('m',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb-1,xpos,elemAZ,Ny,Nz,rx,NoVec) 
          end do

          noB=noB+1 
      end do
      

      
C =========================================================================
C noB=27 
      xpos=(rx-1)*Nb 

      h=(rx-1)*Nb 
      j=(rx-1)*Nb+1 
      call case2_5('u',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 

      h=(rx-1)*Nb+1 
      j=(rx-1)*Nb 
      call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
      call case2_6n('r','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          Nb,xpos,elemAZ,Ny,Nz,rx,NoVec) 

      do h=((rx-1)*Nb+2),(rx*Nb)
        j=(rx-1)*Nb 
        call case2_2('d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos+1,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_3('r','d',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
        call case2_1('l',A,AZ,E,Nb,noB,Nx,NxNy,k,j,h,Nb-1,
     1          xpos,elemAZ,Ny,Nz,rx,NoVec) 
      end do
         
      return
      end



C-------------------------------------------------------------------------------
C  x(i) = alpha
C-------------------------------------------------------------------------------
      subroutine fill_vec(N, alpha, x)
      implicit none
      
      integer N
      
      double precision alpha, x(N)
      
      integer i
      
      do i=1,N                
        x(i) = alpha
      end do                  
      
      return
      end
C-------------------------------------------------------------------------------
C  z = x + alpha y
C-------------------------------------------------------------------------------
      subroutine lin_comb(N, x, y, alpha, z)
      implicit none
      
      integer N
      
      double precision x(N), y(N), z(N), alpha
      
      integer i
      
      do i=1,N                
        z(i) = x(i) + alpha * y(i)
      end do                  
      
      return
      end

C------------------------------------------------------------------------------- C  X(i) = 0
C-------------------------------------------------------------------------------
      subroutine fill_mat(M, N, X)
      implicit none
      
      integer M,N, i, j
      double precision alpha, X(M,N)
      
      do i=1,M 
        do j=1,N              
          X(i,j) = 0
        end do
      end do                  
      
      return
      end

C-------------------------------------------------------------------------------
C  x: Mx=y
C-------------------------------------------------------------------------------
      subroutine precond(Nx, Ny, Nz, A, D, y, x)
      implicit none
      
      integer Nx, Ny, Nz
      
      double precision A(Nx*Ny*Nz,4), D(Nx*Ny*Nz), y(Nx*Ny*Nz)
      double precision x(Nx*Ny*Nz)
      
      double precision, dimension(:), allocatable  :: q
      
      integer NxNy, Nt, i
      
      allocate(q(Nx*Ny*Nz))
      
      NxNy = Nx*Ny
      Nt   = Nx*Ny*Nz
      
C--------------------
C  Diagonal scaling
C--------------------
C      do i=1,Nt
C        x(i) = y(i) / A(i,1)
C      end do
C-----------------------
C  Incomplete Choleski
C-----------------------
      do i=1,1
        q(i) = y(i) / D(i)
      end do
      do i=2,Nx
        q(i) = (y(i) - A(i-1,2)*q(i-1)) / D(i)
      end do
      do i=Nx+1,NxNy
        q(i) = (y(i) - A(i-1,2)*q(i-1) - A(i-Nx,3)*q(i-Nx)) / D(i)
      end do
      do i=NxNy+1,Nt
        q(i) = (y(i) - A(i-1,2)*q(i-1) - A(i-Nx,3)*q(i-Nx)
     1                                 - A(i-NxNy,4)*q(i-NxNy)) / D(i)
      end do
      
      do i=Nt,Nt,-1
        x(i) = q(i)
      end do
      do i=Nt-1,Nt-Nx+1,-1
        x(i) = q(i) - A(i,2)*x(i+1) / D(i)
      end do
      do i=Nt-Nx,Nt-NxNy+1,-1
        x(i) = q(i) - (A(i,2)*x(i+1) + A(i,3)*x(i+Nx)) / D(i)
      end do
      do i=Nt-NxNy,1,-1
        x(i) = q(i) - (A(i,2)*x(i+1) + A(i,3)*x(i+Nx)
     1                               + A(i,4)*x(i+NxNy)) / D(i)
      end do
      
C      write(6,"(E25.15)") (x(i), i=1,Nt)
  
      deallocate(q)
      
      return
      end
C-------------------------------------------------------------------------------
C  
C-------------------------------------------------------------------------------
      subroutine make_cholinc(Nx, Ny, Nz, A, D)
      implicit none
      
      integer Nx, Ny, Nz
      
      double precision A(Nx*Ny*Nz,4)
      double precision D(Nx*Ny*Nz)
      
      integer NxNy, Nt, i
      
      NxNy = Nx*Ny
      Nt   = Nx*Ny*Nz
      
      do i=1,1
        D(i) = A(i,1)
      end do
      
      do i=2,Nx
        D(i) = A(i,1) - A(i-1,2)**2/D(i-1)
      end do
      
      do i=Nx+1,NxNy
        D(i) = A(i,1) - A(i-1,2)**2/D(i-1) - A(i-Nx,3)**2/D(i-Nx)
      end do
      
      do i=NxNy+1,Nt
        D(i) = A(i,1) - A(i-1,2)**2/D(i-1) - A(i-Nx,3)**2/D(i-Nx)
     1                                   - A(i-NxNy,4)**2/D(i-NxNy)
      end do
      
C      write(6,"(E25.15)") (D(i), i=1,Nx*Ny*Nz)
  
      return
      end


C-------------------------------------------------------------------------------
C  z = (x,y)
C-------------------------------------------------------------------------------
      subroutine dot_prod(N, x, y, z)
      implicit none
      
      integer N
      
      double precision x(N), y(N), z
      
      integer i
      
      z = 0.
      do i=1,N                
        z = z + x(i) * y(i)    
      end do                  
      
      return
      end

C-------------------------------------------------------------------------------
C  y = x
C-------------------------------------------------------------------------------
      subroutine copy_vec(N, x, y)
      implicit none
      
      integer N
      
      double precision x(N), y(N)
      
      integer i
      
      do i=1,N                
        y(i) = x(i)
      end do                  
      
      return
      end

C-------------------------------------------------------------------------------
C  y = Ax
C-------------------------------------------------------------------------------
      subroutine matvec_prod(Nx, Ny, Nz, A, x, y)
      implicit none
      
      integer Nx, Ny, Nz

      double precision A(Nx*Ny*Nz,4), x(Nx*Ny*Nz), y(Nx*Ny*Nz)

      integer NxNy, Nt, i

      NxNy = Nx*Ny
      Nt   = Nx*Ny*Nz

      do i=1,Nt
        y(i) = A(i,1)*x(i)
      end do
      
      do i=1,Nt-1
        y(i)   = y(i)   + A(i,2)*x(i+1)
        y(i+1) = y(i+1) + A(i,2)*x(i)
      end do
      
      do i=1,Nt-Nx
        y(i)    = y(i) +    A(i,3)*x(i+Nx)
        y(i+Nx) = y(i+Nx) + A(i,3)*x(i)
      end do
      
      do i=1,Nt-NxNy
        y(i)      = y(i)      + A(i,4)*x(i+NxNy)
        y(i+NxNy) = y(i+NxNy) + A(i,4)*x(i)
      end do
     
C      call output_cg(Nx,Ny,Nz, 'out.m', A, y, x)
C      stop
      
      return
      end
C-------------------------------------------------------------------------------
C  x: Ax = y
C    x should be filled
C-------------------------------------------------------------------------------
      subroutine cg(Nx, Ny, Nz, A, y, tol, x, iter, norm1, norm2)
      implicit none
      
      integer Nx, Ny, Nz
      
      double precision A(Nx*Ny*Nz,4), y(Nx*Ny*Nz), tol
      double precision x(Nx*Ny*Nz)
      
      integer iter
      
      double precision norm1, norm2
      
      double precision, dimension(:,:,:), allocatable :: b, r_old,
     1                                                   D, z,temp1 
      double precision, dimension(:), allocatable ::r,p,r0,temp_vec                                       
      double precision alpha, temp
      double precision tol2
      double precision rz, rz_old
      double precision norm_res2, norm_res2_0
      
      double precision maxx, minx, maxres
      
      double precision realres
      double precision normr0
      
C      double precision D(Nx,Ny,Nz)
      
      integer Nt, i, MAXITER
      
      integer j
      
      parameter(MAXITER = 1000)
      
      allocate(p(Nx*Ny*Nz))
      allocate(b(Nx,Ny,Nz))
      allocate(r(Nx*Ny*Nz))
      allocate(r0(Nx*Ny*Nz))
      allocate(r_old(Nx,Ny,Nz))
      allocate(D(Nx,Ny,Nz))
      allocate(temp_vec(Nx*Ny*Nz))
      allocate(z(Nx,Ny,Nz))
      
      allocate(temp1(Nx,Ny,Nz))
      
      
      Nt = Nx*Ny*Nz
      
      tol2 = tol*tol
      
C      write(*,*) x(1:10)

      call make_cholinc(Nx, Ny, Nz, A, D)
      
C      call fill_vec(Nt, 0d0, x)
C      call copy_vec(Nt, y, r)
      call matvec_prod(Nx, Ny, Nz, A, x, temp_vec)
      call lin_comb(Nt, y, temp_vec, -1d0, r)
      call copy_vec(Nt, r,r0)

      call precond(Nx, Ny, Nz, A, D, r, z)
      
      call copy_vec(Nt, z, p)
      
C      call dot_prod(Nt, r, r, norm_res2_0)
C      norm_res2_0 = -1d0
      call dot_prod(Nt, y, y, norm_res2_0)


C      write(*,*) 'Initial residual norm is',sqrt(norm_res2_0)
      
      if (norm_res2_0 .eq. 0d0) then
        iter = 0
        norm1 = 0.
        norm2 = 0.
        goto 110
      end if
      
       
      
      call dot_prod(Nt, r, z, rz)
      
C      write(*,*) r(1:10)
C      write(*,*) z(1:10,1,1)
C      write(*,*) rz
      
      

      do i=1,MAXITER
        
        call matvec_prod(Nx, Ny, Nz, A, p, temp_vec)
        call dot_prod(Nt, p, temp_vec, temp)

        alpha = rz/temp
        
        call lin_comb(Nt, x, p, alpha, x)
        call lin_comb(Nt, r, temp_vec, -alpha, r)
        
        call dot_prod(Nt, r, r, norm_res2)
        
        norm_res2 = norm_res2 / norm_res2_0
      
        maxx = 0d0
        minx = 0d0
        maxres = 0d0
        do j=1,Nt
          maxx = max(maxx, x(j))
          minx = min(minx, x(j))
          maxres = max(maxres, abs(r(j)))
        end do
        
C       write(6,*) i, tol, tol2
C       write(6,*) i, tol2, norm_res2
C       write(6,*) i, sqrt(norm_res2)
        
C, maxres/(maxx-minx)
C       write(6,'(I10,2X,E9.3)') i, sqrt(norm_res2)
        
        if (norm_res2 .le. tol2) goto 100
        
C       if ((maxres .le. tol*(maxx-minx)) .and.
C     1      (norm_res2.le.tol2)) goto 100
        
        call precond(Nx, Ny, Nz, A, D, r, z)
        
        rz_old = rz
        call dot_prod(Nt, r, z, rz)
        
        call lin_comb(Nt, z, p, rz/rz_old, p)
      end do
      
      write(6,"('CG didn''t converge?!?!')")
      stop
      
  100 continue
      
      
C      write(6,"(I5, 2E15.5, $)") i, sqrt(norm_res2), maxres/(maxx-minx)
      iter = i
      norm1 = sqrt(norm_res2)
C      norm2 = maxres/(maxx-minx)
      
      call matvec_prod(Nx, Ny, Nz, A, x, temp1)
      call lin_comb(Nt, temp1, y, -1d0, temp1)
      call dot_prod(Nt, temp1, temp1, realres)
C      call dot_prod(Nt, y, y, normr0)
      call dot_prod(Nt, r0, r0, normr0)
      realres =  sqrt(realres)/sqrt(normr0)
      norm2=realres
      
 110  continue
 
      deallocate(p)
      deallocate(b)
      deallocate(r)
      deallocate(r0)
      deallocate(r_old)
      deallocate(D)
      deallocate(temp_vec)
      deallocate(z)
      
      deallocate(temp1)
      
      return
      end