subroutine twod_ADI(u0,nx,ny,numt,dx,dy,dt,gamma,D)
	implicit none
	real*8, intent(inout), dimension(nx,ny) :: u0
!f2p intent(in,out) :: u0
	integer, intent(in) :: nx, ny, numt
!f2p intent(in) :: nx, ny, numt
	real*8, intent(in) :: dx,dy,dt,gamma,D
!f2p intent(in) :: dx,dy,gamma,dt,D
	real*8 :: unew(nx,nx) , uold(nx,nx), A(nx,nx)
	real*8:: b(nx), Diag(nx)
	real*8 :: rx, ry, rxx, ryy, uoldkp1, uoldkp2, uoldkm1, uoldkm2, uoldjp1, uoldjp2, uoldjm2, uoldjm1
	real*8 :: DL(nx-1), DU(nx-1)
	integer :: n,j,k, info, lower, upper, ii, jj, ku, kl
	integer :: IPIV(nx+4)
	CHARACTER(LEN=30) :: rowfmt
!apply the boundary conditions
	unew = 0
	uold = u0
	!first do x directions
	A = 0
	b = 0
	rx = DBLE(.5)*D*dt/(dx**2)
	ry = DBLE(.5)*D*dt/(dy**2)
	rxx = DBLE(.5)*D*(gamma**DBLE(2))*dt/(dx**4.00)
	ryy = DBLE(.5)*D*(gamma**DBLE(2))*dt/(dy**4.00)
	!print *, rx, ry, rxx, ryy
	!stop
	do n = 1,numt-1
		do k = 1,nx
			do j = 1,nx
				IF (j == 2) THEN
					uoldjm2 = uold(k,1)
					uoldjm1 = uold(k,j-1)
				ELSE IF (j == 1) THEN
					uoldjm2 = uold(k,2)
					uoldjm1 = uold(k,1)
				ELSE
					uoldjm1 = uold(k,j-1)
					uoldjm2 = uold(k,j-2)
				END IF
				IF (k == 2) THEN
					uoldkm2 = uold(1,j)
					uoldkm1 = uold(k-1,j)
				ELSE IF (k == 1) THEN
					uoldkm2 = uold(2,j)
					uoldkm1 = uold(1,j)
				ELSE
					uoldkm1 = uold(k-1,j)
					uoldkm2 = uold(k-2,j)
				END IF
				IF (j == nx) THEN
					uoldjp1 = uold(k,nx)
					uoldjp2 = uold(k,nx-1)
				ELSE IF (j == nx-1) THEN
					uoldjp1 = uold(k,j+1)
					uoldjp2 = uold(k,nx)
				ELSE
					uoldjp1 = uold(k,j+1)
					uoldjp2 = uold(k,j+2)
				END IF
				IF (k == nx) THEN
					uoldkp1 = uold(nx,j)
					uoldkp2 = uold(nx-1,j)
				ELSE IF (k == nx-1) THEN
					uoldkp1 = uold(k+1,j)
					uoldkp2 = uold(nx,j)
				ELSE
					uoldkp1 = uold(k+1,j)
					uoldkp2 = uold(k+2,j)
				END IF
				IF (j > 1 .AND. j < nx ) THEN
					A(j,j-1) = rx - rx*(uoldjm1**2)
					A(j,j) = DBLE(-2)*rx + DBLE(2)*rx*(uold(k,j)**2) + DBLE(1)
					A(j,j+1) = rx - rx*(uoldjp1**2)
				ELSE IF (j == 1) THEN
					A(j,j) = -rx + rx*(uold(k,j)**2) + DBLE(1)
					A(j,j+1) = rx - rx*(uoldjp1**2)
				ELSE IF (j == nx) THEN
					A(j,j) = -rx + rx*(uold(k,j)**2) + DBLE(1)
					A(j,j-1) = rx - rx*(uoldjm1**2)
				END IF
				b(j) = ry*(uoldkp1**3 - DBLE(2)*uold(k,j)**3 + uoldkm1**3)
				b(j) = b(j)  - ry*(uoldkp1 - DBLE(2)*uold(k,j) + uoldkm1)				
				b(j) = b(j) -rxx*(uoldjp2-DBLE(4)*uoldjp1 + DBLE(6)*uold(k,j) -DBLE(4)*uoldjm1+uoldjm2 )
				b(j) = b(j) - ryy*(uoldkp2-DBLE(4)*uoldkp1 + DBLE(6)*uold(k,j) -DBLE(4)*uoldkm1+uoldkm2 )
				b(j) = b(j) + uold(k,j)
			end do
			!sub ,super, digonal
			DO ii = 1,nx-1
				DL(ii) = A(ii+1,ii)
				Diag(ii) = A(ii,ii)
				DU(ii) = A(ii,ii+1)
			END DO
			Diag(nx) = A(nx,nx)
			! open (unit = 2, file = 'DL.txt')
			! do ii = 1,nx-1
				! write(2,'(ES)') DL(ii)
			! end do
			! close(2)
			!WRITE(rowfmt,'(A,I4,A)') '(',nx,'(ES))'
			! open (unit = 2, file = 'DU.txt')
			! do ii = 1,nx-1
				! write(2,'(ES)') DU(ii)
			! end do
			! close(2)
			! open (unit = 2, file = 'c0.txt')
			! do ii = 1,nx
				! write(2,FMT=rowfmt) (u0(ii,jj), jj = 1,nx)
			! end do
			! close(2)
			! open (unit = 2, file = 'D.txt')
			! do ii = 1,nx
				! write(2,'(ES)') Diag(ii)
			! end do
			! close(2)
			! open (unit = 2, file = 'b_b4.txt')
			! do ii = 1,nx
				! write(2,'(ES)') b(ii)
			! end do
			! close(2)
			CALL DGTSV(nx,1,DL,Diag,DU,b,nx,info)
			! open (unit = 2, file = 'b_after.txt')
			! do ii = 1,nx
				! write(2,'(ES)') b(ii)
			! end do
			! close(2)
			! open (unit = 2, file = 'A.txt')
			! do ii = 1,nx
				! write(2,FMT=rowfmt) (A(ii,jj), jj = 1,nx)
			! end do
			! close(2)
			! stop
			!print *, info, j, k
			!stop
			DO ii = 1,nx
				unew(k,ii) = b(ii)
			END DO
			!print *, unew
		end do
		! open (unit = 2, file = 'unew.txt')
		! do ii = 1,nx
			! write(2,FMT=rowfmt) (unew(ii,jj), jj = 1,nx)
		! end do
		! close(2)
		! stop
		!now do y directions
		uold = unew
		A = 0
		b = 0
		do j = 1,nx
			do k = 1,nx
				IF (j == 2) THEN
					uoldjm2 = uold(k,1)
					uoldjm1 = uold(k,j-1)
				ELSE IF (j == 1) THEN
					uoldjm2 = uold(k,2)
					uoldjm1 = uold(k,1)
				ELSE
					uoldjm1 = uold(k,j-1)
					uoldjm2 = uold(k,j-2)
				END IF
				IF (k == 2) THEN
					uoldkm2 = uold(1,j)
					uoldkm1 = uold(k-1,j)
				ELSE IF (k == 1) THEN
					uoldkm2 = uold(2,j)
					uoldkm1 = uold(1,j)
				ELSE
					uoldkm1 = uold(k-1,j)
					uoldkm2 = uold(k-2,j)
				END IF
				IF (j == nx) THEN
					uoldjp1 = uold(k,nx)
					uoldjp2 = uold(k,nx-1)
				ELSE IF (j == nx-1) THEN
					uoldjp1 = uold(k,j+1)
					uoldjp2 = uold(k,nx)
				ELSE
					uoldjp1 = uold(k,j+1)
					uoldjp2 = uold(k,j+2)
				END IF
				IF (k == nx) THEN
					uoldkp1 = uold(nx,j)
					uoldkp2 = uold(nx-1,j)
				ELSE IF (k == nx-1) THEN
					uoldkp1 = uold(k+1,j)
					uoldkp2 = uold(nx,j)
				ELSE
					uoldkp1 = uold(k+1,j)
					uoldkp2 = uold(k+2,j)
				END IF
				IF (k > 1 .AND. k < nx ) THEN
					A(k,k-1) = ry - ry*(uoldkm1**2)
					A(k,k) = DBLE(-2)*ry + DBLE(2)*ry*(uold(k,j)**2) + DBLE(1)
					A(k,k+1) = ry - ry*(uoldkp1**2)
				ELSE IF (k == 1) THEN
					A(k,k) = -ry + ry*(uold(k,j)**2) + DBLE(1)
					A(k,k+1) = ry - ry*(uoldkp1**2)
				ELSE IF (k == nx) THEN
					A(k,k) = -ry + ry*(uold(k,j)**2) + DBLE(1)
					A(k,k-1) = ry - ry*(uoldkm1**2)
				END IF
			
				b(k) = ry*(uoldjp1**3 - DBLE(2)*uold(k,j)**3 + uoldjm1**3)
				b(k) = b(k) - rx*(uoldjp1 - DBLE(2)*uold(k,j) + uoldjm1)
				b(k) = b(k) - rxx*(uoldjp2 -DBLE(4)*uoldjp1 + DBLE(6)*uold(k,j) -DBLE(4)*uoldjm1 + uoldjm2 )
				b(k) = b(k) - ryy*(uoldkp2 - DBLE(4)*uoldkp1+ DBLE(6)*uold(k,j) -DBLE(4)*uoldkm1 + uoldkm2 )
				b(k) = b(k) + uold(k,j)
			end do
			!sub ,super, digonal
			DO ii = 1,nx-1
				DL(ii) = A(ii+1,ii)
				Diag(ii) = A(ii,ii)
				DU(ii) = A(ii,ii+1)
			END DO			
			Diag(nx) = A(nx,nx)
			! open (unit = 2, file = 'b_b42.txt')
			! do ii = 1,nx
				! write(2,'(ES)') b(ii)
			! end do
			! close(2)
			! WRITE(rowfmt,'(A,I4,A)') '(',nx,'(ES))'
			! open (unit = 2, file = 'A2.txt')
			! do ii = 1,nx
				! write(2,FMT=rowfmt) (A(ii,jj), jj = 1,nx)
			! end do
			! close(2)
			CALL DGTSV(nx,1,DL,Diag,DU,b,nx,info)
			!print *, info, j, k
			DO ii = 1,nx
				unew(ii,j) = b(ii)
			END DO
			! open (unit = 2, file = 'b_after2.txt')
			! do ii = 1,nx
				! write(2,'(ES)') b(ii)
			! end do
			! close(2)
		end do
		uold = unew;

	end do
	!stop
	u0 = unew
end subroutine