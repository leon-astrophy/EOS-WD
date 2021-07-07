! Define double precision !
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND (15, 307)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Mathematical constants and physical constants !
REAL (DP), PARAMETER :: pi = 3.1415926535897932384626433832795E0_DP
REAL (DP), PARAMETER :: hbar = 1.05457266D-27
REAL (DP), PARAMETER :: clight = 2.99792458D10
REAL (DP), PARAMETER :: gconst = 6.6743D-8
REAL (DP), PARAMETER :: solar = 1.98847D33
REAL (DP), PARAMETER :: fine = 7.2973525693D-3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

! Atomic mass unit !		
REAL (DP), PARAMETER :: m_u = 1.66053906660D-24

! Fermionic mass for electrons !		
REAL (DP), PARAMETER :: m_e = 9.109389699D-28

! Neutron mass !
REAL (DP), PARAMETER :: m_n = 1.67492749804D-24

! Proton mass !
REAL (DP), PARAMETER :: m_p = 1.67262192369D-24

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the scaling constant !
REAL (DP), PARAMETER :: a_e = (m_e**4*clight**5)/(2.4D1*pi**2*hbar**3)
REAL (DP), PARAMETER :: b_e = (m_e**3*clight**3)/(3.0D0*pi**2*hbar**3)

! Find the scaling constant !
REAL (DP), PARAMETER :: a_n = (m_n**4*clight**5)/(2.4D1*pi**2*hbar**3)
REAL (DP), PARAMETER :: b_n = (m_n**3*clight**3)/(3.0D0*pi**2*hbar**3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For HW equation of state !
REAL (DP), PARAMETER :: b1 = 0.991749D0
REAL (DP), PARAMETER :: b2 = 0.01911D0
REAL (DP), PARAMETER :: b3 = (m_n - m_p - m_e)/(m_u)
REAL (DP), PARAMETER :: b4 = 0.10175D0
REAL (DP), PARAMETER :: b5 = 0.000763D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Atomic number !
REAL (DP), PARAMETER :: zatom = 6.0D0
REAL (DP), PARAMETER :: aatom = 12.0D0

! For Salpeter EOS !
REAL (DP), PARAMETER :: pctf = -(m_e*clight**2)*(m_e*clight/hbar)**3
REAL (DP), PARAMETER :: a_c = pctf*(fine*zatom**(2.0D0/3.0D0)/10.0D0/pi**2)*(4.0D0/9.0D0/pi)**(1.0D0/3.0D0)
REAL (DP), PARAMETER :: a_tf = pctf*(1.62D2/1.75D2*(fine*zatom**(2.0D0/3.0D0))**2/9.0D0/pi**2)*(4.0D0/9.0D0/pi)**(2.0D0/3.0D0)
REAL (DP), PARAMETER :: a_ex = pctf*(fine/4.0D0/pi**3)
REAL (DP), PARAMETER :: a_cor = pctf*(fine**2*0.0311D0/9.0D0/pi**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For salpeter EOS !
REAL (DP), PARAMETER :: ry = 0.5D0*fine**2*m_e*clight**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For Salpeter EOS !
REAL (DP), PARAMETER :: e_c = -ry*1.8D0*(zatom)**(2.0D0/3.0D0)
REAL (DP), PARAMETER :: e_tf = -ry*(324.0D0/175.0D0)*(4.0D0/9.0D0/pi)**(2.0D0/3.0D0)*zatom**(4.0D0/3.0D0)
REAL (DP), PARAMETER :: e_ex = -(0.75D0/pi)*fine*m_e*clight**2
REAL (DP), PARAMETER :: e_cor = ry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Number of models per maximum density !
INTEGER, PARAMETER :: n_space = 50001

! Starting log maximum density !
REAL (DP), PARAMETER :: abarstart = 61.5D0

! Ending log maximum density !
REAL (DP), PARAMETER :: abarend = 96.6D0

! Step size !
REAL (DP), PARAMETER :: dabar = (abarend - abarstart)/(DBLE(n_space) - 1.0D0)

! Starting log maximum density !
REAL (DP), PARAMETER :: xstart = -2.0D0

! Ending log maximum density !
REAL (DP), PARAMETER :: xend = 1.537D0

! Step size !
REAL (DP), PARAMETER :: dx = (xend - xstart)/(DBLE(n_space) - 1.0D0)

! Starting log maximum density !
REAL (DP), PARAMETER :: rhostart = 7.3D0

! Ending log maximum density !
REAL (DP), PARAMETER :: rhoend = 11.0D0

! Step size !
REAL (DP), PARAMETER :: drho = (rhoend - rhostart)/(DBLE(n_space) - 1.0D0)

! Starting log maximum density !
REAL (DP), PARAMETER :: rhostart2 = 8.65D0

! Ending log maximum density !
REAL (DP), PARAMETER :: rhoend2 = 10.73D0

! Step size !
REAL (DP), PARAMETER :: drho2 = (rhoend2 - rhostart2)/(DBLE(n_space) - 1.0D0)

! Starting log maximum density !
REAL (DP), PARAMETER :: rhostart3 = 3.6D0

! Ending log maximum density !
REAL (DP), PARAMETER :: rhoend3 = 11.0D0

! Step size !
REAL (DP), PARAMETER :: drho3 = (rhoend3 - rhostart3)/(DBLE(n_space) - 1.0D0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Maxmimum iteration !
INTEGER, PARAMETER :: nmax = 1000

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!