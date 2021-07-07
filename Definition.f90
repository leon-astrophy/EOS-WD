MODULE DEFINITION
IMPLICIT NONE
INCLUDE "Parameter.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Density !
REAL (DP) :: rho_e
REAL (DP) :: rho_n
REAL (DP) :: rho_atom

! Pressure !
REAL (DP) :: p_c
REAL (DP) :: p_tf
REAL (DP) :: p_cor
REAL (DP) :: p_ex
REAL (DP) :: p_e
REAL (DP) :: p_n
REAL (DP) :: p_out
REAL (DP) :: p_ideal
REAL (DP) :: p_total

! Epsilon !
REAL (DP) :: eps_c
REAL (DP) :: eps_tf
REAL (DP) :: eps_cor
REAL (DP) :: eps_ex
REAL (DP) :: eps_e
REAL (DP) :: eps_n
REAL (DP) :: eps_out
REAL (DP) :: eps_ideal
REAL (DP) :: eps_total

! Enthalpy !
REAL (DP) :: h_e
REAL (DP) :: h_n
REAL (DP) :: h_out
REAL (DP) :: h_ideal
REAL (DP) :: h_total

! Electron fermi momentum !
REAL (DP) :: dlfmmo

! Electron fermi momentum !
REAL (DP) :: x_e
REAL (DP) :: x_n

! Mass number and atomic number
REAL (DP) :: abar
REAL (DP) :: zbar

! Number density !
REAL (DP) :: n_e
REAL (DP) :: n_n
REAL (DP) :: n_atom

! Dummy !
REAL (DP) :: lhs_e
REAL (DP) :: lhs_n

! Electron fractions !
REAL (DP) :: ye
REAL (DP) :: ye_ccsn

! Dimensionless variables !
REAL (DP) :: xrho

! For parameterized electron fractions !
REAL (DP) :: rho1
REAL (DP) :: rho2
REAL (DP) :: y1
REAL (DP) :: y2
REAL (DP) :: yc

! Number of lines in Ye table !
INTEGER :: nlines

! For EOS Table !
REAL (DP), ALLOCATABLE, DIMENSION(:) :: rhotable
REAL (DP), ALLOCATABLE, DIMENSION(:) :: yetable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Density !
REAL (DP) :: rhoscale2
REAL (DP) :: rhomax2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE