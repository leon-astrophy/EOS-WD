PROGRAM WDEOS
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Integer !
INTEGER :: left, right, m

! Real !
REAL (DP) :: beta, re, dummy

! Open 
OPEN (UNIT = 9, FILE = 'HarrisonWheeler.eos', STATUS = 'REPLACE')

! Loop for Harrison-Wheeler EOS !
DO i = 1, n_space

	! mass number !
	abar = abarstart + (DBLE(i) - 1.0D0)*dabar
	
	! Assign !
	zbar = SQRT(0.5D0*b2/b5*abar)

	! Calculate neutron fractions !
	lhs_n = b1 + (2.0D0/3.0D0)*(b2/abar**(1.0D0/3.0D0)) + b4*(0.25D0 - (zbar/abar)**2) - (1.0D0/3.0D0)*b5*zbar**2/abar**(4.0D0/3.0D0)
	lhs_n = lhs_n*(m_u/m_n)
	lhs_n = lhs_n**2 - 1.0D0

	! For neutron !
	IF(lhs_n > 0.0D0) THEN
		x_n = SQRT(lhs_n)
		n_n = b_n*x_n**3
	ELSE
		x_n = 0.0D0
		n_n = 0.0D0
	END IF

	! Calculate electron fractions !
	lhs_e = b3 + b4*(1.0D0 - 2.0D0*zbar/abar) - 2.0D0*b5*zbar/abar**(1.0D0/3.0D0)
	lhs_e = lhs_e*(m_u/m_e)
	lhs_e = (lhs_e + 1.0D0)**2 - 1.0D0
	
	! Electron fractions !
	x_e = SQRT(lhs_e)
	n_e = b_e*x_e**3

	! Find nuclei fractions !
	n_atom = (n_e/zbar)

	! Find density !
	rho_e = m_e*n_e
	rho_n = m_n*n_n
	rho_atom = n_atom*liquid_drop(abar,zbar)
	rhomax2 = rho_e + rho_n + rho_atom

	! Print out the effective electron fractions !
	IF(i == 1) THEN
		ye = (m_u/((rhomax2/b_e/x_e**3) - m_e))
		WRITE (9,*) ye, zbar/abar
	END IF	

	! Pressure and epsilon !
	IF(x_e < 3.0D-3) THEN
		p_e = a_e*small_pressure(x_e)
		eps_e = a_e*small_energy(x_e)
	ELSE
		p_e = a_e*large_pressure(x_e)
		eps_e = a_e*large_energy(x_e)
	END IF
	IF(x_n < 3.0D-3) THEN
		p_n = a_n*small_pressure(x_n)
		eps_n = a_n*small_energy(x_n)
	ELSE
		p_n = a_n*large_pressure(x_n)
		eps_n = a_n*large_energy(x_n)
	END IF

	! Enthalpy !
	IF((rho_e + rho_atom) > 0.0D0) THEN
		h_e = (p_e/(rho_e + rho_atom)) + ((eps_e)/(rho_e + rho_atom))
	ELSE
		h_e = 0.0D0
	END IF
	IF(rho_n > 0.0D0) THEN
		h_n = (p_n/rho_n) + ((eps_n)/rho_n)
	ELSE
		h_n = 0.0D0
	END IF

	! Sum them up !
	p_out = p_e + p_n
	h_out = h_e + h_n
	eps_out = eps_e + eps_n

	! Print out !
	WRITE (9,*) log10(rhomax2), x_e, log10(h_out)

END DO

! Close !
CLOSE (9)

! Open !
OPEN (UNIT = 99, FILE = 'Neutronization.eos', STATUS = 'REPLACE')

! Loop for parameterized electron fractions EOS !
DO i = 1, n_space

	! Fermi momentum !
	x_e = xstart + (DBLE(i) - 1.0D0)*dx
	x_e = 10.0D0**(x_e)

	! Get Ye !
	ye = 1.0D0/(2.0D0 + 1.255D-2*x_e + 1.755D-5*x_e**2 + 1.376D-6*x_e**3)

	! Print out the effective electron fractions !
	IF(i == 1) THEN
		WRITE (99,*) ye, ye
	END IF	

	! Get density !
	rhomax2 = (m_u/ye + m_e)*b_e*x_e**3

	! Pressure and epsilon !
	IF(x_e < 3.0D-3) THEN
		p_e = a_e*small_pressure(x_e)
		eps_e = a_e*small_energy(x_e)
	ELSE
		p_e = a_e*large_pressure(x_e)
		eps_e = a_e*large_energy(x_e)
	END IF

	! Enthalpy !
	h_e = (p_e/(rhomax2)) + ((eps_e)/(rhomax2))

	! Print !
	WRITE (99,*) log10(rhomax2), x_e, log10(h_e)

END DO

! Close !
CLOSE (99)

! Open file !
OPEN (UNIT = 999, FILE = 'LiebendorferN13.eos', STATUS = 'REPLACE')

! For parameterized electron fractions !
rho1 = 2.0D7
rho2 = 2.0D13	
y1 = 0.5D0
y2 = 0.285D0
yc = 0.035D0

! Loop for parameterized electron fractions EOS !
DO i = 1, n_space

	! Density !
	rhomax2 = rhostart + (DBLE(i) - 1.0D0)*drho
	rhomax2 = 10.0D0**(rhomax2)

	! Evaluate !
	xrho = max(-1.0D0, min(1.0D0, ((2.0D0*log10(rhomax2) - log10(rho2) - log10(rho1))/(log10(rho2) - log10(rho1)))))
	ye = 0.5D0*(y2 + y1) + 0.5D0*xrho*(y2 - y1) + yc*(1.0D0 - abs(xrho) + 4.0D0*abs(xrho)*(abs(xrho) - 0.5D0)*(abs(xrho) - 1.0D0))

	! Print out the effective electron fractions !
	IF(i == 1) THEN
		WRITE (999,*) ye, ye
	END IF	

	! Find Fermi momentum !	
	x_e = (rhomax2/b_e/(m_u/ye + m_e))**(1.0D0/3.0D0)

	! Pressure and epsilon !
	IF(x_e < 3.0D-3) THEN
		p_e = a_e*small_pressure(x_e)
		eps_e = a_e*small_energy(x_e)
	ELSE
		p_e = a_e*large_pressure(x_e)
		eps_e = a_e*large_energy(x_e)
	END IF

	! Enthalpy !
	h_e = (p_e/(rhomax2)) + ((eps_e)/(rhomax2))

	! Print !
	WRITE (999,*) log10(rhomax2), x_e, log10(h_e)

END DO

CLOSE (999)

! Open file !
OPEN (UNIT = 999, FILE = 'LiebendorferG15.eos', STATUS = 'REPLACE')

! For parameterized electron fractions !
rho1 = 3.0D7
rho2 = 2.0D13	
y1 = 0.5D0
y2 = 0.278D0
yc = 0.035D0

! Loop for parameterized electron fractions EOS !
DO i = 1, n_space

	! Density !
	rhomax2 = rhostart + (DBLE(i) - 1.0D0)*drho
	rhomax2 = 10.0D0**(rhomax2)

	! Evaluate !
	xrho = max(-1.0D0, min(1.0D0, ((2.0D0*log10(rhomax2) - log10(rho2) - log10(rho1))/(log10(rho2) - log10(rho1)))))
	ye = 0.5D0*(y2 + y1) + 0.5D0*xrho*(y2 - y1) + yc*(1.0D0 - abs(xrho) + 4.0D0*abs(xrho)*(abs(xrho) - 0.5D0)*(abs(xrho) - 1.0D0))

	! Print out the effective electron fractions !
	IF(i == 1) THEN
		WRITE (999,*) ye, ye
	END IF	

	! Find Fermi momentum !	
	x_e = (rhomax2/b_e/(m_u/ye + m_e))**(1.0D0/3.0D0)

	! Pressure and epsilon !
	IF(x_e < 3.0D-3) THEN
		p_e = a_e*small_pressure(x_e)
		eps_e = a_e*small_energy(x_e)
	ELSE
		p_e = a_e*large_pressure(x_e)
		eps_e = a_e*large_energy(x_e)
	END IF

	! Enthalpy !
	h_e = (p_e/(rhomax2)) + ((eps_e)/(rhomax2))
	
	! Print !
	WRITE (999,*) log10(rhomax2), x_e, log10(h_e)

END DO

CLOSE (999)

! Open file !
OPEN (UNIT = 999, FILE = 'LiebendorferS15.eos', STATUS = 'REPLACE')

! For parameterized electron fractions !
rho1 = 2.2D8
rho2 = 9.5D12
y1 = 0.5D0
y2 = 0.279D0
yc = 0.022D0

! Loop for parameterized electron fractions EOS !
DO i = 1, n_space

	! Density !
	rhomax2 = rhostart + (DBLE(i) - 1.0D0)*drho
	rhomax2 = 10.0D0**(rhomax2)

	! Evaluate !
	xrho = max(-1.0D0, min(1.0D0, ((2.0D0*log10(rhomax2) - log10(rho2) - log10(rho1))/(log10(rho2) - log10(rho1)))))
	ye = 0.5D0*(y2 + y1) + 0.5D0*xrho*(y2 - y1) + yc*(1.0D0 - abs(xrho) + 4.0D0*abs(xrho)*(abs(xrho) - 0.5D0)*(abs(xrho) - 1.0D0))
	
	! Print out the effective electron fractions !
	IF(i == 1) THEN
		WRITE (999,*) ye, ye
	END IF	

	! Find Fermi momentum !	
	x_e = (rhomax2/b_e/(m_u/ye + m_e))**(1.0D0/3.0D0)

	! Pressure and epsilon !
	IF(x_e < 3.0D-3) THEN
		p_e = a_e*small_pressure(x_e)
		eps_e = a_e*small_energy(x_e)
	ELSE
		p_e = a_e*large_pressure(x_e)
		eps_e = a_e*large_energy(x_e)
	END IF

	! Enthalpy !
	h_e = (p_e/(rhomax2)) + ((eps_e)/(rhomax2))
	
	! Print !
	WRITE (999,*) log10(rhomax2), x_e, log10(h_e)

END DO

CLOSE (999)

! Read the number of lines in the file !
nlines = 0 
OPEN (999, file = 'y_e_vs_rho_at_7ms_Tc5.0d09_from_VULCAN2D.dat') 
DO 
    READ (999,*, END=10) 
    nlines = nlines + 1 
END DO 
10 CLOSE (999) 

! Allocate arrays !
ALLOCATE(rhotable(0:nlines-1))
ALLOCATE(yetable(0:nlines-1))

! Read !
OPEN(UNIT=999, FILE = 'y_e_vs_rho_at_7ms_Tc5.0d09_from_VULCAN2D.dat', ACTION='READ')
DO i = 1, nlines 
	READ(999,*) rhotable(nlines-i), yetable(nlines-i)
	rhotable(nlines-i) = log10(rhotable(nlines-i))
ENDDO
CLOSE(999)

! Open file !
OPEN (UNIT = 999, FILE = 'VULCAN2D.eos', STATUS = 'REPLACE')

! Loop for parameterized electron fractions EOS !
DO i = 1, n_space

	! Density !
	rhomax2 = rhostart2 + (DBLE(i) - 1.0D0)*drho2

	! Case by case !
	IF(rhomax2 < rhotable(0)) THEN

		! Table minimum !
		ye = yetable(0)

	ELSEIF(rhomax2 == rhotable(0)) THEN

		! Table minimum !
		ye = yetable(0)

	ELSE

		! Binary search !
		left = 0
		right = nlines - 1
		DO
			IF(left > right) THEN
				EXIT
			END IF
			m = floor(REAL((right + left)/2))
			IF(rhotable(m) < rhomax2) THEN
				left = m + 1	
			ELSEIF(rhotable(m) > rhomax2) THEN
				right = m - 1
			ELSEIF(rhotable(m) == rhomax2) THEN
				EXIT
			END IF
		END DO
		IF(rhotable(m) == rhomax2) THEN
			ye = yetable(m)
		ELSE
			m = m - 1
			IF(m > 1 .AND. m < nlines - 3) THEN
				CALL AKIMA(rhotable(m-2), rhotable(m-1), rhotable(m), rhotable(m+1), rhotable(m+2), rhotable(m+3), & 
					yetable(m-2), yetable(m-1), yetable(m), yetable(m+1), yetable(m+2), yetable(m+3), rhomax2, ye)
			ELSE
				! Linear !
				CALL LINEAR(rhotable(m), rhotable(m+1), yetable(m), yetable(m+1), rhomax2, ye)
			END IF
		END IF

	END IF

	! Power up the density !
	rhomax2 = 10.0D0**(rhomax2)

	! Print out the effective electron fractions !
	IF(i == 1) THEN
		WRITE (999,*) ye, ye
	END IF	

	! Find Fermi momentum !	
	x_e = (rhomax2/b_e/(m_u/ye + m_e))**(1.0D0/3.0D0)

	! Pressure and epsilon !
	IF(x_e < 3.0D-3) THEN
		p_e = a_e*small_pressure(x_e)
		eps_e = a_e*small_energy(x_e)
	ELSE
		p_e = a_e*large_pressure(x_e)
		eps_e = a_e*large_energy(x_e)
	END IF

	! Enthalpy !
	h_e = (p_e/(rhomax2)) + ((eps_e)/(rhomax2))

	! Print !
	WRITE (999,*) log10(rhomax2), x_e, log10(h_e)

END DO

CLOSE (999)

contains

	real(DP) function liquid_drop(a,z)
	implicit none
	real(DP) :: a, z
	liquid_drop = m_u*(b1*a + b2*a**(2.0D0/3.0D0) - b3*z + b4*a*(0.5D0 - z/a)**2 + b5*z**2/a**(1.0D0/3.0D0)) - z*m_e
	end function

	real(DP) function large_pressure(x)
	implicit none
	real(DP) :: x
	large_pressure = x*SQRT(x**2 + 1.0D0)*(2.0D0*x**2 - 3.0D0) + 3.0D0*log(x + SQRT(x**2 + 1.0D0))
	end function

	real(DP) function small_pressure(x)
	implicit none
	real(DP) :: x
	small_pressure = 1.6D0*x**5 - (4.0D0/7.0D0)*x**7 + (1.0D0/3.0D0)*x**9 - (5.0D0/2.2D1)*x**11 & 
			+ (3.5D1/2.08D2)*x**13 - (2.1D1/1.6D2)*x**15 + (2.31D2/2.176D3)*x**17
	end function

	real(DP) function large_energy(x)
	implicit none
	real(DP) :: x
	large_energy = 3.0D0*x*SQRT(x**2 + 1.0D0)*(1.0D0 + 2.0D0*x**2) - 3.0D0*log(x + SQRT(x**2 + 1.0D0)) - 8.0D0*x**3
	end function

	real(DP) function small_energy(x)
	implicit none
	real(DP) :: x
	small_energy = 8.0D0*x**3 + (1.2D1/5.0D0)*x**5 - (3.0D0/7.0D0)*x**7 + (1.0D0/6.0D0)*x**9 - (1.5D1/1.76D2)*x**11 & 
			+ (2.1D1/4.16D2)*x**13 - (2.1D1/6.40D2)*x**15 + (9.9D1/4.352D3)*x**17 - 8.0D0*x**3
	end function

	real(DP) function exchange_pressure(x, b)
	implicit none
	real(DP) :: x, b
	exchange_pressure = 0.03125D0*(b**4 + 1.0D0/b**4) + 0.25D0*(b**2 + 1.0D0/b**2) - 0.5625D0 - & 
			0.75D0*(b**2 - 1.0D0/b**2)*log(b) + 1.5D0*(log(b))**2 - & 
			x*(1.0D0 + x/sqrt(1.0D0 + x**2))/3.0D0*(0.125D0*(b**3 - 1.0D0/b**5) - 0.25D0*(b - 1.0D0/b**3) & 
			- 1.5D0*(b + 1.0D0/b**3)*log(b) + 3.0D0*log(b)/b)
	end function

	real(DP) function exchange_epsilon(x, b)
	implicit none
	real(DP) :: x, b
	exchange_epsilon = 0.25D0*(2.25D0 + 3.0D0*(b**2 - 1.0D0/b**2)*log(b) - 6.0D0*(log(b))**2 - &
			 (b**2 + 1.0D0/b**2) - 0.125D0*(b**4 + 1.0D0/b**2))/x**3

	end function

END PROGRAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine do 2-point interpolation for a given point of existing !
!datum and output the quantity that the user wished to interpolate 	 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LINEAR(x0, x1, y0, y1, x_in, y_out)
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! Input parameter !
REAL (DP), INTENT(IN) :: x0, x1, y0, y1, x_in

! Input parameter !
REAL (DP), INTENT(OUT) :: y_out

! Intermediate polynominal !
REAL (DP) :: nu0, nu1
REAL (DP) :: de0, de1
REAL (DP) :: l0, l1

! Assign numerator !
nu0 = (x_in - x1)
nu1 = (x_in - x0)

! Assign denominator !
de0 = (x0 - x1)
de1 = (x1 - x0)

! Assign polynominal coefficient !
l0 = nu0/de0
l1 = nu1/de1

! Compute the output !
y_out = l0*y0 + l1*y1

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Akima spline interpolation. See Hiroshi Akima 1970 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AKIMA(xm2, xm1, x0, xp1, xp2, xp3, ym2, ym1, y0, yp1, yp2, yp3, x_in, y_out)
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! Input parameter !
REAL (DP), INTENT(IN) :: xm2, xm1, x0, xp1, xp2, xp3, ym2, ym1, y0, yp1, yp2, yp3, x_in

! Input parameter !
REAL (DP), INTENT(OUT) :: y_out

! Intermediate polynominal !
REAL (DP) :: dm2, dm1, d0, dp1, dp2
REAL (DP) :: s0, s1

! Weights !
REAL (DP) :: w1, w2, w3, w4

! Coefficient of polynominal !
REAL (DP) :: p0, p1, p2, p3

! Temporal arrays !
REAL (DP) :: diff, temp

! Assign slopes !
dm2 = (ym1 - ym2)/(xm1 - xm2)
dm1 = (y0 - ym1)/(x0 - xm1)
d0 = (yp1 - y0)/(xp1 - x0)
dp1 = (yp2 - yp1)/(xp2 - xp1)
dp2 = (yp3 - yp2)/(xp3 - xp2)

! Assign weights (modified) !
w1 = abs(dp1 - d0) + 0.5D0*abs(dp1 + d0)
w2 = abs(dm1 - dm2) + 0.5D0*abs(dm1 + dm2)
w3 = abs(dp2 - dp1) + 0.5D0*abs(dp2 + dp1)
w4 = abs(d0 - dm1) + 0.5D0*abs(d0 + dm1)

! assign slopes !
IF(w1 == 0.0D0 .AND. w2 == 0.0D0) THEN
	s0 = 0.5D0*(dm1 + d0)
ELSE
	s0 = (w1*dm1 + w2*d0)/(w1 + w2)
END IF
IF(w3 == 0.0D0 .AND. w4 == 0.0D0) THEN
	s1 = 0.5D0*(d0 + dm1)
ELSE
	s1 = (w3*d0 + w4*dp1)/(w3 + w4)
END IF

! Assign temp !
diff = xp1 - x0
temp = x_in - x0

! assign coefficients !
p0 = y0
p1 = s0
p2 = (3.0D0*d0 - 2.0D0*s0 - s1)/diff
p3 = (s0 + s1 - 2.0D0*d0)/diff**2

! Output the interpolation !
y_out = p0 + p1*temp + p2*temp**2 + p3*temp**3

END SUBROUTINE