Module symmetron_integrator

CONTAINS

SUBROUTINE solver(phi, res, P, y, a, b, EPS, omega, N)
	  !
		!		AUTHORS: Valeri Vardanyan <vardanyan@lorentz.leidenuniv.nl>
		!						Omar Contigiani <contigiani@lorentz.leidenuniv.nl>
		!
	  !   Solves Laplacian phi(y) = a * ((b*P - 1)*phi(y) + phi(y)**3) for a
		!		spherically symmetric rho using a Newton-Gauss-Seidel method.
		! 	Fourth order numerical derivatives are employed.
		!		The initial guess for phi(y) is refined over multiple iterations
		!		until the precision EPS is reached.
		!
		!
	  !		Arguments:
		! 		phi - The initial guess for the function phi at points y_i
		!			P - The
		!			y - ***Equispaced*** positions y_i
		!			a - Eq. parameter
		!			b - Eq. parameter
		!			EPS - threshold for the error between the LHS and RHS of the equation.
		!						Quantified as the square root of the quadratic sum of
		!						the differences for every y.
		!			omega - The final result of each iteration is phi_new. The next
		!							iteration is initialized with
		!							phi_(n+1) = phi_(n)*(omega-1) + phi_new*(omega-1)
		!			N  - size of the phi and P arrays.

		INTEGER MAXITS, imax
	  INTEGER N
	  DOUBLE PRECISION a, b
		DOUBLE PRECISION A_16, A_8, L_eqn(1:N)
		DOUBLE PRECISION, INTENT(INOUT) :: phi(1:N)
		DOUBLE PRECISION P(1:N), res(1:N), y(1:N), EPS
		PARAMETER (MAXITS=1000000)
		INTEGER i, n_iter
		DOUBLE PRECISION anorm, anorm_old, omega, phi_new, h
		DOUBLE PRECISION Delta_phi_sum

		! Initialize vars:
		anorm = 0.d0
		imax = N
		h = y(2)-y(1)
		DO i = 1, imax
				L_eqn(i) = 0.d0
		END DO

		!Boundary Conditions:
		phi(1) = 0.d0 !bnd cnd
		phi(2) = 0.d0 !bnd cnd
		phi(imax-1) = phi(imax - 2)
		phi(imax) = phi(imax - 1)
		!phi(imax) = phi(imax - 1)
		!phi(imax) = 1.d0 !bnd cnd
		!phi(imax - 1) = 1.d0 !bnd cnd

		!Zero-th iteration:
		DO i = 3, imax - 2
				A_16 = -phi(i + 2) + 16.d0*phi(i + 1) + 16.d0*phi(i - 1) - phi(i - 2)
				A_8 = -phi(i + 2) + 8.d0*phi(i + 1) - 8.d0*phi(i - 1) + phi(i - 2)

				L_eqn(i) = (A_16 - 3.d1*phi(i))/(12.d0*h*h) + (2.d0/y(i))*A_8/(12.d0*h) - a*((b*P(i) - 1.d0)*phi(i) + phi(i)**3.d0)
		END DO

		DO n_iter=1, MAXITS
			anorm_old = anorm
			anorm = 0.d0
	    Delta_phi_sum = 0.d0

			DO i = 3, imax - 2
					res(i) = L_eqn(i)
			    anorm = anorm + abs(L_eqn(i)*L_eqn(i))

			    phi_new = phi(i) - L_eqn(i)/(-30.d0/(12.d0*h*h) - a*((b*P(i) - 1.d0) + 3.d0*phi(i)**2.d0))
	        Delta_phi_sum = Delta_phi_sum + (phi(i) - phi_new)*(phi(i) - phi_new)
			    phi(i) = omega*phi_new + (1.d0 - omega)*phi(i)

			    A_16 = -phi(i + 2) + 16.d0*phi(i + 1) + 16.d0*phi(i - 1) - phi(i - 2)
			    A_8 = -phi(i + 2) + 8.d0*phi(i + 1) - 8.d0*phi(i - 1) + phi(i - 2)
			    L_eqn(i) = (A_16 - 30.d0*phi(i))/(12.d0*h*h) + (2.d0/y(i))*A_8/(12.d0*h) - a*((b*P(i) - 1.d0)*phi(i) + phi(i)**3.d0)
			END DO
			phi(imax-1) = phi(imax - 2)
			phi(imax) = phi(imax - 1)
			phi(1) = 0.d0
			phi(2) = 0.d0

			anorm = sqrt(anorm)
			!write(*,*) n_iter, h, sqrt(anorm), sqrt(Delta_phi_sum)


			IF (anorm.lt.EPS) THEN
				write(*, *) "Converged in"
				write(*, *) n_iter
				RETURN
			ENDIF
			IF (sqrt(Delta_phi_sum).lt.1.d-14) THEN
				write(*, *) "Not converged - machine precision reached"
				RETURN
			ENDIF
		END DO

		write(*,*) 'Not converged - MAXITER reached'

END SUBROUTINE solver

END MODULE symmetron_integrator
