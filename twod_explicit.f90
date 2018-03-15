subroutine twod_explicit(u0,nx,ny,numt,dx,dy,dt,gamma,D)
	implicit none
	real(8), intent(inout), dimension(0:nx-1,0:ny-1) :: u0
!f2p intent(in,out) :: u0
	integer, intent(in) :: nx, ny, numt
!f2p intent(in) :: nx, ny, numt
	real(8), intent(in) :: dx,dy,dt,gamma,D
!f2p intent(in) :: dx,dy,gamma,dt,D
	real(8) :: unew(0:nx+3,0:nx+3) , uold(0:nx+3,0:nx+3)
	real(8) :: term1, term2, term3, term4, term5, term6, RHS
	integer :: n,j,k
!apply the boundary conditions
	unew = 0
	!print *, SHAPE(unew)
	!print *, SHAPE(u0)
	unew(2:nx+1,2:nx+1) = u0
	unew(0,2:nx+1) = u0(1,:)
	unew(2:nx+1,0) = u0(:,1)
	unew(1,2:nx+1) = u0(0,:)
	unew(2:nx+1,1) = u0(:,0)
	unew(nx+3,2:nx+1) = u0(nx-2,:)
	unew(2:nx+1,nx+3) = u0(:,nx-2)
	unew(nx+2,2:nx+1) = u0(nx-1,:)
	unew(2:nx+1,nx+2) = u0(:,nx-1)
	uold = unew
	do n = 1,numt-1
		do j = 2,nx+1
			do k = 2,nx+1
				term1 = uold(k,j+1)**3 - 2 * uold(k,j)**3 + uold(k,j-1)**3
				term2 = uold(k+1,j)**3 - 2 * uold(k,j)**3 + uold(k-1,j)**3
				term3 = uold(k,j+1) - 2*uold(k,j) + uold(k,j-1)
				term4 = uold(k+1,j) - 2*uold(k,j) + uold(k-1,j)
				term5 = uold(k,j-2) -4*uold(k,j-1) + 6*uold(k,j) -4*uold(k,j+1) + uold(k,j+2)
				term6 = uold(k-2,j) - 4*uold(k-1,j) + 6*uold(k,j) - 4*uold(k+1,j) + uold(k+2,j)
				RHS = term1/(dx**2) + term2/(dy**2) - term3/(dx**2) - term4/(dy**2) 
				RHS = RHS - gamma**2 * term5 / (dx**4) - gamma**2 * term6 / (dy**4)
				unew(k,j) = D*RHS*dt + uold(k,j)
			end do
		end do
		uold = unew;
		!apply boundary conditions again	unew(2:nx+1,2:nx+1) = u0
		uold(0,2:nx+1) = unew(3,2:nx+1)
		uold(2:nx+1,0) = unew(2:nx+1,3)
		uold(1,2:nx+1) = unew(2,2:nx+1)
		uold(2:nx+1,1) = unew(2:nx+1,2)
		uold(nx+3,2:nx+1) = unew(nx,2:nx+1)
		uold(2:nx+1,nx+3) = unew(2:nx+1,nx)
		uold(nx+2,2:nx+1) = unew(nx+1,2:nx+1)
		uold(2:nx+1,nx+2) = unew(2:nx+1,nx+1)
		unew = uold
	end do
	u0 = unew(2:nx+1,2:nx+1)
end subroutine