!==============================================================================!
! MODULE Projection                                                            !
!                                                                              !
! This module contains the variables and routines related to quantum-number    !
! projections.                                                                 !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_projection_angles                                           !
! - subroutine set_projection_matelem                                          !
! - subroutine determine_angle                                                 !
! - subroutine check_symmetries                                                !
! - subroutine check_compatibilities                                           !
! - subroutine reset_rotmatelem                                                !
! - subroutine integrate_rotmatelem                                            !
! - subroutine open_file_rotmatelem                                            !
! - subroutine readwrite_rotmatelem                                            !
! - subroutine print_rotmatelem                                                !
! - subroutine generate_rotation_spatial                                       !
! - subroutine generate_rotation_parity                                        !
! - subroutine generate_rotation_gauge                                         !
! - subroutine rotate_wavefunction                                             !
! - subroutine calculate_thouless                                              !
! - subroutine print_projmatelem_states                                        !
! - subroutine write_projmatelem_states                                        !
! - subroutine write_projmatelem_occnumb                                       !
! - subroutine write_projmatelem_electromagnetic                               !
! - subroutine print_files_complementary                                       !
! - subroutine reduce_projmatelem                                              !
!==============================================================================!
MODULE Projection 

!cmpi use MPI            
!cmpi use Parallelization
use Fields
use Operators

implicit none

!!! Miscellaneous      
integer :: misc_phys, & ! physics case computed
           misc_part, & ! full or partial calculations
           misc_frot    ! writes/reads the rotated matrix elements
real(r64) :: misc_cutrot ! cutoff on norm overlaps

!!! Rotated matrix elements
real(r64) :: frot_cutrot ! cutoff on norm overlaps (written in file)
complex(r64) :: rot_over,  & ! overlap
                rot_pari,  & ! parity
                rot_ener,  & ! energy
                rot_prot,  & ! proton  number
                rot_neut,  & ! neutron   "
                rot_prot2, & ! proton    "    squared
                rot_neut2, & ! neutron   "       "
                rot_amjx,  & ! angular momentum x axis
                rot_amjy,  & !    "       "     y  "
                rot_amjz,  & !    "       "     z  "
                rot_amjx2, & !    "       "     x  "   squared
                rot_amjy2, & !    "       "     y  "      "
                rot_amjz2, & !    "       "     z  "      "
                rot_amj2,  & !    "       "     total     "
                rot_amsop, & ! spin-orbit for protons
                rot_amson, & !    "   "    "  neutrons
                rot_istx,  & ! isospin x axis
                rot_isty,  & !    "    y  "
                rot_istz,  & !    "    z  "
                rot_istx2, & !    "    x  "   squared
                rot_isty2, & !    "    y  "      "
                rot_istz2, & !    "    z  "      "
                rot_ist2,  & !    "    total     "
                rot_ra2p,  & ! radius protons  squared
                rot_ra2n,  & !    "   neutrons    "
                rot_ra2m     !    "   matter      "
complex(r64), dimension(-1:1) :: rot_Q1mp, & ! Q_1m for protons
                                 rot_Q1mn, & !  "    "  neutrons
                                 rot_M1ma    ! M_1m  "  nucleons
complex(r64), dimension(-2:2) :: rot_Q2mp, & ! Q_2m for protons
                                 rot_Q2mn, & !  "    "  neutrons
                                 rot_M2ma    ! M_2m  "  nucleons
complex(r64), dimension(-3:3) :: rot_Q3mp, & ! Q_3m for protons
                                 rot_Q3mn    !  "    "  neutrons
complex(r64), dimension(:,:), allocatable :: rot_occn, & ! occupation numbers
                                             Ubar, & ! U matrix 
                                             Vbar, & !
                                             Util, & !
                                             Vtil, & !
                                             ROT, &  ! rotation matrix
                                             Rid     ! identity   "
character(len=:), allocatable :: frot_hamil ! hamiltonian name (written in file)

!!! Projected matrix elements 
integer :: tabid_J_dim, & ! dimension of table of indices for J,MJ,KJ
           tabid_P_dim    !     "     "    "   "     "     "  P
integer, dimension(:), allocatable :: tabid_P ! table of indices for J,MJ,KJ
integer, dimension(:,:,:), allocatable :: tabid_J !"  "     "     "  P           
complex(r64), dimension(:,:), allocatable :: proj_over,  & ! same as 
                                             proj_ener,  & ! rotated 
                                             proj_pari,  & ! quantities
                                             proj_prot,  &
                                             proj_neut,  &
                                             proj_prot2, &
                                             proj_neut2, &
                                             proj_amjz,  &
                                             proj_amjz2, &
                                             proj_amj2,  &
                                             proj_amsop, &
                                             proj_amson, &
                                             proj_istz,  &
                                             proj_istz2, &
                                             proj_ist2,  &
                                             proj_ra2p,  &
                                             proj_ra2n,  &
                                             proj_ra2m
complex(r64), dimension(:,:,:), allocatable :: proj_Q1mp, &
                                               proj_Q1mn, &
                                               proj_Q2mp, &
                                               proj_Q2mn, &
                                               proj_Q3mp, &
                                               proj_Q3mn, &
                                               proj_M1ma, &
                                               proj_M2ma
complex(r64), dimension(:,:,:,:), allocatable :: proj_occn

!!! Integrals
integer :: phiZ, & ! angle phiA
           phiN, & !   "   phiZ
           phiA, & !   "   phiN
           alpJ, & !   "   alphaJ
           betJ, & !   "   betaJ
           gamJ, & !   "   gammaJ
           parP, & !   "   parityP
           phiZ_dim, & ! number of angles for phiA  
           phiN_dim, & !   "    "    "     "  phiZ  
           phiA_dim, & !   "    "    "     "  phiN  
           alpJ_dim, & !   "    "    "     "  alpJ 
           betJ_dim, & !   "    "    "     "  betJ 
           gamJ_dim, & !   "    "    "     "  gamJ 
           parP_dim, & !   "    "    "     "  parP  
           pnp_facpi=2, & ! factor integral [0,X*pi] in PNP
           amp_2jmin, & ! smallest value of 2*j   
           amp_2jmax, & ! largest    "   "  2*j 
           pap_pmin,  & ! smallest value "  p
           pap_pmax,  & ! largest    "   "  p
           pnp_nosimp, & ! disable simplifications for PNP
           amp_nosimp, & !    "           "         "  AMP
           pap_nosimp    !    "           "         "  PAP
integer(i64) :: loop_ini, & ! start of the loop 
                loop_fin, & ! end   "   "   " 
                loop_dim    ! total number of angles in the loop
integer, dimension(:), allocatable :: angle_parP ! angles for parity 
real(r64), dimension(:), allocatable :: angle_phiZ,  & ! angles  for phiZ   
                                        angle_phiN,  & !   "      "  phiN   
                                        angle_phiA,  & !   "      "  phiA   
                                        angle_alpJ,  & !   "      "  alpJ 
                                        angle_betJ,  & !   "      "  betJ  
                                        angle_gamJ,  & !   "      "  gamJ 
                                        weight_phiZ, & ! weights  "  phiA   
                                        weight_phiN, & !    "     "  phiZ   
                                        weight_phiA, & !    "     "  phiN   
                                        weight_alpJ, & !    "     "  alpJ   
                                        weight_betJ, & !    "     "  betJ   
                                        weight_gamJ, & !    "     "  gamJ   
                                        weight_parP    !    "     "  parP   
logical :: is_pnpZ,  & ! switch projection Z
           is_pnpN,  & !    "       "      N
           is_pnpA,  & !    "       "      A
           is_ampJ,  & !    "       "      J
           is_ampMJ, & !    "       "      MJ
           is_ampKJ, & !    "       "      KJ 
           is_papP,  & !    "       "      P
           is_phiZ,  & !    "   integral  phiZ
           is_phiN,  & !    "       "     phiN
           is_phiA,  & !    "       "     phiA
           is_alpJ,  & !    "       "     alpJ
           is_betJ,  & !    "       "     betJ
           is_gamJ,  & !    "       "     gamJ
           is_parP     !    "       "     parP

!!! Good quantum numbers
integer, dimension(2) :: good_Z,   & ! good number of protons  Z
                         good_N,   & !  "     "    "  neutrons N
                         good_A,   & !  "     "    "  nucleons A
                         good_J,   & !  "     "    "  angular momentum 2*J 
                         good_MKJ, & !  "     "    "  z-component 2*MJ/2*KJ
                         good_P,   & !  "     "    "  parity  P
                         good_T,   & !  "     "    "  isospin 2*T
                         good_MKT    !  "     "    "  z-component 2*MT/2*KT
logical, dimension(2) :: is_good_1,   & ! switch good precision overlap
                         is_good_Z,   & !   "    "    Z    
                         is_good_N,   & !   "    "    N
                         is_good_A,   & !   "    "    A
                         is_good_ZN,  & !   "    "    factorization Z x N
                         is_good_J,   & !   "    "    J
                         is_good_MKJ, & !   "    "    MJ,KJ
                         is_good_P,   & !   "    "    A
                         is_good_T,   & !   "    "    T
                         is_good_MKT    !   "    "    MT,KT

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_projection_angles                                             !
!                                                                              ! 
! Set the discretization for rotations over Euler angles and other quantities  ! 
! related to the symmetry projections.                                         ! 
!------------------------------------------------------------------------------!
subroutine set_projection_angles

integer :: i, ialloc=0
real(r64), dimension(:), allocatable :: cosbetJ
!cmpi integer(i64) :: divide, rest

!!! Matrices related to rotations
allocate( Ubar(HOsp_dim,HOsp_dim), Vbar(HOsp_dim,HOsp_dim), &
          Util(HOsp_dim,HOsp_dim), Vtil(HOsp_dim,HOsp_dim), &
          ROT(HOsp_dim,HOsp_dim), Rid(HOsp_dim,HOsp_dim), &
          stat=ialloc ) 
if ( ialloc /= 0 ) stop 'Error during allocation of rotation matrices'

Rid = zzero
do i = 1, HOsp_dim
  Rid(i,i) = one
enddo

!!! Angles and weights
allocate( angle_phiZ(phiZ_dim), angle_phiN(phiN_dim), angle_phiA(phiA_dim), &
          angle_alpJ(alpJ_dim), angle_betJ(betJ_dim), angle_gamJ(gamJ_dim), & 
          angle_parP(parP_dim), cosbetJ(betJ_dim), stat=ialloc ) 
if ( ialloc /= 0 ) stop 'Error during allocation of angles'

angle_phiZ = zero
angle_phiN = zero
angle_phiA = zero
angle_alpJ = zero
angle_betJ = zero
angle_gamJ = zero
angle_parP = 1    
cosbetJ = one 

allocate( weight_phiZ(phiZ_dim), weight_phiN(phiN_dim), weight_phiA(phiA_dim), &
          weight_alpJ(alpJ_dim), weight_betJ(betJ_dim), weight_gamJ(gamJ_dim), &
          weight_parP(parP_dim), stat=ialloc ) 
if ( ialloc /= 0 ) stop 'Error during allocation of weights'

weight_phiZ = one 
weight_phiN = one 
weight_phiA = one 
weight_alpJ = one 
weight_betJ = one 
weight_gamJ = one 
weight_parP = one  

!!! Particle number
if ( is_phiZ ) then 
  if ( modulo(phiZ_dim,2) == 1 ) then
    do i = 1, phiZ_dim  
      angle_phiZ(i) = pnp_facpi * pi * (i-1) / phiZ_dim
    enddo 
  else
    do i = 1, phiZ_dim  
      angle_phiZ(i) = pnp_facpi * pi * (i-0.5) / phiZ_dim
    enddo 
  endif
  weight_phiZ = one / phiZ_dim
endif

if ( is_phiN ) then 
  if ( modulo(phiN_dim,2) == 1 ) then
    do i = 1, phiN_dim  
      angle_phiN(i) = pnp_facpi * pi * (i-1) / phiN_dim
    enddo 
  else
    do i = 1, phiN_dim  
      angle_phiN(i) = pnp_facpi * pi * (i-0.5) / phiN_dim
    enddo 
  endif
  weight_phiN = one / phiN_dim
endif

if ( is_phiA ) then 
  do i = 1, phiA_dim 
    angle_phiA(i) = pi * (i-1) / phiA_dim
  enddo 
  weight_phiA = one / phiA_dim
endif

!!! Angular momentum
if ( is_alpJ ) then 
  do i = 1, alpJ_dim
    angle_alpJ(i) = 2 * pi * (i-0.5) / alpJ_dim 
  enddo 
  weight_alpJ = one / alpJ_dim 
endif

if ( is_betJ ) then 
  call GaussLegendre(-1.0d+00,1.0d+00,cosbetJ,weight_betJ,betJ_dim)
  do i = 1, betJ_dim
    angle_betJ(i) = dacos(cosbetJ(i))
  enddo   
endif

if ( is_gamJ ) then 
  do i = 1, gamJ_dim
    angle_gamJ(i) = 2 * pi * (i-0.5) / gamJ_dim 
  enddo 
    weight_gamJ = one / gamJ_dim 
endif

!!! Parity
if ( is_parP .and. is_papP ) then 
  do i = 1, parP_dim
    angle_parP(i) = i
  enddo
  weight_parP = one / 2            
endif

!!! Loop dimensions
loop_dim = phiZ_dim * phiN_dim * phiA_dim * &
           alpJ_dim * betJ_dim * gamJ_dim * &
           parP_dim

loop_ini = 1        
loop_fin = loop_dim

!!! MPI distribution of angles
!cmpi divide = loop_dim / paral_teams
!cmpi rest = modulo(loop_dim,paral_teams)
!cmpi 
!cmpi paral_myangles = divide
!cmpi paral_myoffset = paral_myangles * paral_myteam
!cmpi 
!cmpi if ( paral_myteam < rest ) then
!cmpi   paral_myangles = paral_myangles + 1
!cmpi   paral_myoffset = paral_myoffset + paral_myteam
!cmpi else
!cmpi   paral_myoffset = paral_myoffset + rest
!cmpi endif
!cmpi 
!cmpi loop_ini = paral_myoffset + 1
!cmpi loop_fin = paral_myoffset + paral_myangles

end subroutine set_projection_angles

!------------------------------------------------------------------------------!
! subroutine set_projection_matelem                                            !
!                                                                              ! 
! Set the discretization for rotations over Euler angles and other quantities  ! 
! related to the symmetry projections.                                         ! 
!------------------------------------------------------------------------------!
subroutine set_projection_matelem

integer :: j, m, k, ialloc=0

!!! Sets the indices for parity projection
allocate( tabid_P(-parP_dim:parP_dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation tab_P in set_projection_&
                         &matelem'

tabid_P = 0
tabid_P_dim = 0

do j = pap_pmin, pap_pmax, 2
  tabid_P_dim = tabid_P_dim + 1
  tabid_P(j) = tabid_P_dim
enddo

!!! Sets the indices for angular momentum projection
allocate( tabid_J(-amp_2jmax:amp_2jmax,-amp_2jmax:amp_2jmax, &
                 & amp_2jmin:amp_2jmax), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation tab_J in set_projection_&
                         &matelem'
tabid_J = 0
tabid_J_dim = 0

do j = amp_2jmin, amp_2jmax, 2
  do m = -j, j, 2
    do k = -j, j, 2
      tabid_J_dim = tabid_J_dim + 1
      tabid_J(k,m,j) = tabid_J_dim
    enddo
  enddo
enddo

!!! Sets the arrays with the reduced indices
if ( misc_phys == 0 ) then 
  allocate( proj_over(tabid_P_dim,tabid_J_dim),  &
            proj_pari(tabid_P_dim,tabid_J_dim),  &
            proj_ener(tabid_P_dim,tabid_J_dim),  &
            proj_prot(tabid_P_dim,tabid_J_dim),  &
            proj_neut(tabid_P_dim,tabid_J_dim),  &
            proj_prot2(tabid_P_dim,tabid_J_dim), &
            proj_neut2(tabid_P_dim,tabid_J_dim), &
            proj_amjz(tabid_P_dim,tabid_J_dim),  &
            proj_amjz2(tabid_P_dim,tabid_J_dim), &
            proj_amj2(tabid_P_dim,tabid_J_dim),  &
            proj_amsop(tabid_P_dim,tabid_J_dim), &
            proj_amson(tabid_P_dim,tabid_J_dim), &
            proj_istz(tabid_P_dim,tabid_J_dim),  &
            proj_istz2(tabid_P_dim,tabid_J_dim), &
            proj_ist2(tabid_P_dim,tabid_J_dim),  &
            proj_ra2p(tabid_P_dim,tabid_J_dim),  &
            proj_ra2n(tabid_P_dim,tabid_J_dim),  &
            proj_ra2m(tabid_P_dim,tabid_J_dim),  &
            proj_Q1mp(-1:1,tabid_P_dim,tabid_J_dim), &
            proj_Q1mn(-1:1,tabid_P_dim,tabid_J_dim), &
            proj_Q2mp(-2:2,tabid_P_dim,tabid_J_dim), &
            proj_Q2mn(-2:2,tabid_P_dim,tabid_J_dim), &
            proj_Q3mp(-3:3,tabid_P_dim,tabid_J_dim), &
            proj_Q3mn(-3:3,tabid_P_dim,tabid_J_dim), &
            proj_M1ma(-1:1,tabid_P_dim,tabid_J_dim), &
            proj_M2ma(-2:2,tabid_P_dim,tabid_J_dim), &
            proj_occn(HOsh_dim,2,tabid_P_dim,tabid_J_dim), &
            rot_occn(HOsh_dim,2), &
            source=zzero, stat=ialloc ) 
endif
if ( ialloc /= 0 ) stop 'Error during allocation matelem in set_projection_&
                         &matelem!'

end subroutine set_projection_matelem

!------------------------------------------------------------------------------!
! subroutine determine_angle                                                   !
!                                                                              ! 
! Determines the position in the full integral over all angles. The order for  ! 
! the integral is the following:                                               ! 
! 1 = parP                                                                     ! 
! 2 = phiN                                                                     ! 
! 3 = phiZ                                                                     ! 
! 4 = phiA                                                                     ! 
! 5 = gamJ                                                                     ! 
! 6 = betJ                                                                     ! 
! 7 = alpJ                                                                     ! 
!                                                                              ! 
! Input: iloop = position in the loop over all angles                          ! 
!------------------------------------------------------------------------------!
subroutine determine_angle(iloop)

integer(i64), intent(in) :: iloop
integer(i64) :: i1, i2, i3, i4, i5, i6, i7, j1, j2, j3, j4, j5, j6, j7, &
                k1, k2, k3, k4, k5, k6, k7, m1, m2, m3, m4, m5, m6, m7, &
                iphiZ, iphiN, iphiA, ialpJ, ibetJ, igamJ, iparP

i1 = parP_dim
i2 = gamJ_dim
i3 = betJ_dim
i4 = alpJ_dim
i5 = phiA_dim
i6 = phiN_dim
i7 = phiZ_dim

j1 = i1 
j2 = j1 * i2
j3 = j2 * i3
j4 = j3 * i4
j5 = j4 * i5
j6 = j5 * i6
j7 = j6 * i7

k1 = iloop / j1
k2 = iloop / j2
k3 = iloop / j3
k4 = iloop / j4
k5 = iloop / j5
k6 = iloop / j6
k7 = iloop / j7

m1 = 0
m2 = 0
m3 = 0
m4 = 0
m5 = 0
m6 = 0
m7 = 0
if ( modulo(iloop,j1) == 0 ) m1 = 1
if ( modulo(iloop,j2) == 0 ) m2 = 1
if ( modulo(iloop,j3) == 0 ) m3 = 1
if ( modulo(iloop,j4) == 0 ) m4 = 1
if ( modulo(iloop,j5) == 0 ) m5 = 1
if ( modulo(iloop,j6) == 0 ) m6 = 1
if ( modulo(iloop,j7) == 0 ) m7 = 1

iparP = modulo(iloop,j1) + i1 * m1
igamJ = k1 + 1 - m1 - i2 * (k2 - m2)
ibetJ = k2 + 1 - m2 - i3 * (k3 - m3)
ialpJ = k3 + 1 - m3 - i4 * (k4 - m4)
iphiA = k4 + 1 - m4 - i5 * (k5 - m5)
iphiN = k5 + 1 - m5 - i6 * (k6 - m6)
iphiZ = k6 + 1 - m6 - i6 * (k7 - m7)

parP = int(iparP,i32)
gamJ = int(igamJ,i32)
betJ = int(ibetJ,i32)
alpJ = int(ialpJ,i32)
phiA = int(iphiA,i32)
phiN = int(iphiN,i32)
phiZ = int(iphiZ,i32)

end subroutine determine_angle  

!------------------------------------------------------------------------------!  
! subroutine check_symmetries                                                  !  
!                                                                              !  
! Check if the left or right states have a good number of protons, neutrons,   !  
! parity or third-component of th angular momentum. If it is the case, the     !  
! corresponding projection is skipped.                                         !  
!                                                                              !  
! Input: ndim = dimension of the sp basis                                      !
!------------------------------------------------------------------------------!  
subroutine check_symmetries(ndim)

integer, intent(in) :: ndim
integer :: i, j, t, occ, emp
real(r64) :: xtover, xtpari, xtprot, xsprot, xtneut, xsneut, xtamjz, xtamjv, &
             xsamjv, xtistz, xtistv, xsistv, xtnucl, xsnucl, ovva, eps=1.0d-8
real(r64), dimension(:), allocatable :: vou
real(r64), dimension(2) :: xj, xt, xpn
complex(r64), dimension(2) :: over, pari, & 
                              prot, neut, nucl, prot2, neut2, nucl2, &
                              amjx, amjy, amjz, amjx2, amjy2, amjz2, &
                              istx, isty, istz, istx2, isty2, istz2, &
                              amj2, ist2, amjv2, istv2, amsop, amson
complex(r64), dimension(ndim,ndim) :: ROT, U, V, D, rhoLR, kappaLR, kappaRL
character(3), dimension(2) :: ch_good_1, ch_good_P, ch_good_ZN, ch_good_Z,  &
                              ch_good_N, ch_good_A, ch_good_J, ch_good_MKJ, &
                              ch_good_T, ch_good_MKT
character(len=:), allocatable :: action_1, action_P, action_ZN, action_Z,  &
                                 action_N, action_A, action_J, action_MKJ, &
                                 action_T, action_MKT
character(len=*), parameter :: format1 = "(1a12,4x,1a3,1x,1f14.8,2x,1a3,1x, &
                                           &1f14.8,4x,1a)"

ch_good_1   = ' no'
ch_good_Z   = ' no'
ch_good_N   = ' no'
ch_good_A   = ' no'
ch_good_ZN  = ' no'
ch_good_J   = ' no'
ch_good_P   = ' no'
ch_good_MKJ = ' no'
ch_good_T   = ' no'
ch_good_MKT = ' no'

do t = 1, 2

  !!! Copies states 
  if ( t == 1 ) then ! Left states
    U = bogo_zUL 
    V = bogo_zVL
    D = bogo_zDL
    ovva = bogo_ovvaL
    occ = bogo_occL  
    emp = bogo_empL  
    allocate(vou(HOsp_dim-bogo_occL))
    vou = bogo_vouL
  else               ! Right states
    U = bogo_zUR 
    V = bogo_zVR
    D = bogo_zDR
    ovva = bogo_ovvaR
    occ = bogo_occR  
    emp = bogo_empR  
    allocate(vou(HOsp_dim-bogo_occR))
    vou = bogo_vouR
  endif

  !!! Evaluates operators  
  ROT = zzero
  do i = 1,ndim
    ROT(i,i) = zone
  enddo
                                                                                   
  call calculate_overlap(occ,emp,ovva,vou,D,occ,emp,ovva,vou,D,ROT,over(t), &
                         ndim)
  call calculate_densities(U,V,U,V,rhoLR,kappaLR,kappaRL,ndim)
  call calculate_particlenumber(0,rhoLR,kappaLR,kappaRL,prot(t),neut(t), &
                                prot2(t),neut2(t),ndim)
  call calculate_angularmomentum(0,rhoLR,kappaLR,kappaRL,amjx(t),amjy(t), & 
                                 amjz(t),amjx2(t),amjy2(t),amjz2(t),amsop(t), &
                                 amson(t),ndim)
  call calculate_isospin(0,rhoLR,kappaLR,kappaRL,istx(t),isty(t),istz(t), & 
                         istx2(t),isty2(t),istz2(t),ndim)
  call generate_rotation_parity(2,zone*ROT,ROT,ndim)
  call calculate_overlap(occ,emp,ovva,vou,D,occ,emp,ovva,vou,D,ROT,pari(t), &
                         ndim)

  nucl(t)  = prot(t) + neut(t)
  nucl2(t) = 2.0d0 * (prot2(t) + neut2(t) - 2.0d0*istz2(t))

  amj2(t) = amjx2(t) + amjy2(t) + amjz2(t)
  amjv2(t) = amjx(t)**2 + amjy(t)**2 + amjz(t)**2

  ist2(t) = istx2(t) + isty2(t) + istz2(t)
  istv2(t) = istx(t)**2 + isty(t)**2 + istz(t)**2

  pari(t) = pari(t) / over(t)

  !!! Performs the tests
  ! overlap 1
  xtover = abs(over(t) - one)         

  if ( xtover < eps ) then
    is_good_1(t) = .true.
    ch_good_1(t) = 'yes'
  endif

  ! particle number Z
  xtprot = abs(prot2(t) - prot(t)**2)
  xsprot = abs(prot(t) - valence_Z*zone)

  if ( (xtprot < eps) .and. (xsprot < eps) ) then  
    is_good_Z(t) = .true.
    ch_good_Z(t) = 'yes'
    good_Z(t) = int(abs(prot(t))+0.1d0) 
    if ( pnp_nosimp == 0 ) then
       phiZ_dim = 1
       is_phiZ = .false.
    endif
  endif

  ! particle number N
  xtneut = abs(neut2(t) - neut(t)**2)
  xsneut = abs(neut(t) - valence_N*zone)
  if ( (xtneut < eps) .and. (xsneut < eps) ) then  
    is_good_N(t) = .true.
    ch_good_N(t) = 'yes'
    good_N(t) = int(abs(neut(t))+0.1d0) 
    if ( pnp_nosimp == 0 ) then
       phiN_dim = 1
       is_phiN = .false.
    endif
  endif

  ! particle number A
  xtnucl = abs(nucl2(t) - nucl(t)**2)
  xsnucl = abs(nucl(t) - valence_A*zone)

  if ( (xtnucl < eps) .and. (xsnucl < eps) ) then  
    is_good_A(t) = .true.
    ch_good_A(t) = 'yes'
    good_A(t) = int(abs(nucl(t))+0.1d0) 
    if ( pnp_nosimp == 0 ) then
       phiA_dim = 1
       is_phiA = .false.
    endif
  endif

  ! separate Z and N
  xpn(t) = zero
  do j = 1, ndim/2
    do i = 1, ndim/2
      xpn(t) = xpn(t) + abs(kappaLR(i,j+ndim/2)) + abs(kappaLR(i+ndim/2,j))
    enddo
  enddo
  
  if ( xpn(t) < eps ) then
    is_good_ZN(t) = .true.
    ch_good_ZN(t) = 'yes'
    if ( pnp_nosimp == 0 ) pnp_facpi = 1 
  endif

  ! angular momentum J^2
  xj(t) = sqrt(abs(real(amjv2(t))))
  xtamjv = abs(amj2(t) - xj(t)*(xj(t)+1)) 
  xsamjv = abs(one*int(2*xj(t)+0.001) - 2*xj(t))
  if ( (aimag(amjv2(t)) > eps) .or. real(amjv2(t)) < zero ) xtamjv = one
  xj(t) = 0.5d0 * (-1 + sqrt(1+4*abs(amj2(t))))

  if ( (xtamjv < eps) .and. (xsamjv < eps) ) then  
    is_good_J(t) = .true.
    ch_good_J(t) = 'yes'
    good_J(t) = int(2*xj(t)+0.1d0) 
    if ( amp_nosimp == 0 ) then 
      betJ_dim = 1
      is_betJ = .false.
    endif
  endif

  ! angular momentum Jz              
  xtamjz = abs(amjz2(t) - amjz(t)**2)

  if ( xtamjz < eps ) then  
    is_good_MKJ(t) = .true.
    ch_good_MKJ(t) = 'yes'
    good_MKJ(t) = int(sign(one,real(amjz(t)))) * 2 * int(abs(amjz(t))+0.1d0) 
    if ( amp_nosimp == 0 ) then 
      if ( t == 1 ) then 
        alpJ_dim = 1
        is_alpJ = .false.
      else 
        gamJ_dim = 1
        is_gamJ = .false.
      endif 
    endif
  endif

  ! parity P
  xtpari = abs(abs(pari(t)) - one) 

  if ( xtpari < eps ) then
    is_good_P(t) = .true.
    ch_good_P(t) = 'yes'
    good_P(t) = int(sign(one,real(pari(t)))) * int(abs(pari(t))+0.1d0) 
    if ( pap_nosimp == 0 ) then 
      parP_dim = 1
      is_parP = .false.
    endif
  endif

  ! isospin T^2
  xt(t) = sqrt(abs(real(istv2(t))))
  xtistv = abs(ist2(t) - xt(t)*(xt(t)+1)) 
  xsistv = abs(one*int(2*xt(t)+0.001) - 2*xt(t))
  if ( (aimag(istv2(t)) > eps) .or. real(istv2(t)) < zero ) xtistv = one
  xt(t) = 0.5d0 * (-1 + sqrt(1+4*abs(ist2(t))))

  if ( (xtistv < eps) .and. (xsistv < eps) ) then  
    is_good_T(t) = .true.
    ch_good_T(t) = 'yes'
    good_T(t) = int(2*xt(t)+0.1d0) 
  endif

  ! isospin Tz              
  xtistz = abs(istz2(t) - istz(t)**2)

  if ( xtistz < eps ) then  
    is_good_MKT(t) = .true.
    ch_good_MKT(t) = 'yes'
    good_MKT(t) = int(sign(one,real(istz(t)))) * 2 * int(abs(istz(t))+0.1d0) 
  endif

  deallocate(vou)
enddo

!!! Printing 
! overlap 1
if ( .not. any(is_good_1)  ) then 
  action_1 = "none (but be careful)"   
else
  action_1 = "none"
endif 

! particle number Z
if ( any(is_good_Z) ) then 
  if ( pnp_nosimp == 0 ) then
    action_Z = "no projection on Z"   
  else
    action_Z = "none (simpl. disabled)"   
  endif
else
  action_Z = "none"
endif 

! particle number N
if ( any(is_good_N) ) then 
  if ( pnp_nosimp == 0 ) then
    action_N = "no projection on N"   
  else
    action_N = "none (simpl. disabled)"   
  endif
else
  action_N = "none"
endif 

! particule number A
if ( any(is_good_A) ) then 
  if ( pnp_nosimp == 0 ) then
    action_A = "no projection on A"   
  else
    action_A = "none (simpl. disabled)"   
  endif
else
  action_A = "none"
endif 

! separate Z and N
if ( any(is_good_ZN) ) then 
  if ( pnp_nosimp == 0 ) then
    action_ZN = "integral [1,pi] for Z,N"   
  else
    action_ZN = "none (simpl. disabled)"   
  endif
else
  action_ZN = "none"
endif 

! angular momtenum J^2
if ( any(is_good_J) ) then 
  if ( amp_nosimp == 0 ) then
    action_J = "no projection on J (beta)"   
  else
    action_J = "none (simpl. disabled)"   
  endif
else
  action_J = "none"
endif 

! angular momtenum Jz  
if ( any(is_good_MKJ) ) then 
  if ( amp_nosimp == 0 ) then
    if ( all(is_good_MKJ) ) then 
      action_MKJ = "no projection on MJ,KJ"   
    elseif ( is_good_MKJ(1) ) then 
      action_MKJ = "no projection on MJ"   
    elseif ( is_good_MKJ(2) ) then 
      action_MKJ = "no projection on KJ"   
    endif
  else
    action_MKJ = "none (simpl. disabled)"   
  endif
else
  action_MKJ = "none"
endif 

! parity P
if ( any(is_good_P) ) then 
  if ( pap_nosimp == 0 ) then
    action_P = "no projection on P"   
  else
    action_P = "none (simpl. disabled)"   
  endif
else
  action_P = "none"
endif 

! isospin T^2 and Tz
action_T = "none"
action_MKT = "none"

!cmpi if ( paral_myrank == 0 ) then 
print '(/,23x,"Left",15x,"Right",/,16x,38("-"),/2x,"Symmetry",7x,"?",7x, &
        &"Average",5x,"?",7x,"Average",5x,"Action",/,63("-"))'
print format1, 'Good overlap', ch_good_1(1), real(over(1)), &
                               ch_good_1(2), real(over(2)), action_1
print format1, 'Good Z      ', ch_good_Z(1), real(prot(1)), &
                               ch_good_Z(2), real(prot(2)), action_Z
print format1, 'Good N      ', ch_good_N(1), real(neut(1)), &
                               ch_good_N(2), real(neut(2)), action_N
print format1, 'Good A      ', ch_good_A(1), real(nucl(1)), &
                               ch_good_A(2), real(nucl(2)), action_A
print format1, 'Separate N/Z', ch_good_ZN(1), xpn(1)/ndim, & 
                               ch_good_ZN(2), xpn(2)/ndim, action_ZN
print format1, 'Good J      ', ch_good_J(1), xj(1), &
                               ch_good_J(2), xj(2), action_J
print format1, 'Good MJ/KJ  ', ch_good_MKJ(1), real(amjz(1)), &
                               ch_good_MKJ(2), real(amjz(2)), action_MKJ
print format1, 'Good P      ', ch_good_P(1), real(pari(1)), &
                               ch_good_P(2), real(pari(2)), action_P
print format1, 'Good T      ', ch_good_T(1), xt(1), & 
                               ch_good_T(2), xt(2), action_T
print format1, 'Good MT/KT  ', ch_good_MKT(1), real(istz(1)), &
                               ch_good_MKT(2), real(istz(2)), action_MKT
!cmpi endif

end subroutine check_symmetries

!------------------------------------------------------------------------------!  
! subroutine check_compatibilities                                             !  
!                                                                              !  
! Checks the compatibility of the input parameters with the possible good      !  
! quantum numbers found in check_symmetries. Stops the calculation if they are !  
! not compatible.                                                              !  
!------------------------------------------------------------------------------!  
subroutine check_compatibilities

integer :: ierror, l, jmin0, jmin1, jmin3, jmax0, jmax1, jmax3

!!! Counter for the number of errors (should be 0 at the end)
ierror = 0

!!!
!!! Particle number
!!!

if ( bogo_npaR /= bogo_npaL ) then
  ierror = ierror + 1
  print '(1a)', 'The number parity of the left state is not compatible with &
                &the number parity of the right state.'
endif

if ( bogo_npaR /= (-1)**(valence_A) ) then
  ierror = ierror + 1
  print '(1a)', 'The number parity of the right state is not compatible with &
                &the number of nucleons determined from the inputs.'
endif

if ( is_good_Z(1) ) then
  if ( ((misc_phys == 0) .and. ((good_Z(1) - valence_Z) /= 0)) .or. &
       ((misc_phys == 1) .and. (abs(good_Z(1) - valence_Z) > 1)) .or. &
       ((misc_phys == 2) .and. ((good_Z(1) - valence_Z) /= 2)) ) then
    ierror = ierror + 1
    print '(1a)', 'The number of protons of the left state is not compatible &
                  &with the number of protons determined from the inputs.'
  endif
endif

if ( is_good_Z(2) ) then
  if ( (good_Z(2) - valence_Z) /= 0 ) then
    ierror = ierror + 1
    print '(1a)', 'The number of protons of the right state is not compatible &
                  &with the number of protons determined from the inputs.'
  endif
endif

if ( is_good_N(1) ) then
  if ( ((misc_phys == 0) .and. ((good_N(1) - valence_N) /= 0)) .or. &
       ((misc_phys == 1) .and. (abs(good_N(1) - valence_N) > 1)) .or. &
       ((misc_phys == 2) .and. ((good_N(1) - valence_N) /= 2)) ) then
    ierror = ierror + 1
    print '(1a)', 'The number of neutrons of the left state is not compatible &
                  &with the number of neutrons determined from the inputs.'
  endif
endif

if ( is_good_N(2) ) then
  if ( (good_N(2) - valence_N) /= 0 ) then
    ierror = ierror + 1
    print '(1a)', 'The number of neutrons of the right state is not compatible &
                  &with the number of neutrons determined from the inputs.'
  endif
endif

if ( is_good_A(1) ) then
  if ( ((misc_phys == 0) .and. ((good_A(1) - valence_A) /= 0)) .or. &
       ((misc_phys == 1) .and. (abs(good_A(1) - valence_A) > 1)) .or. &
       ((misc_phys == 2) .and. ((good_A(1) - valence_A) /= 2)) ) then
    ierror = ierror + 1
    print '(1a)', 'The number of nucleons of the left state is not compatible &
                  &with the number of nucleons determined from the inputs.'
  endif
endif

if ( is_good_A(2) ) then
  if ( (good_A(2) - valence_A) /= 0 ) then
    ierror = ierror + 1
    print '(1a)', 'The number of nucleons of the right state is not compatible &
                  &with the number of nucleons determined from the inputs.'
  endif
endif

!!!
!!! Angular momentum
!!!

jmin0 = amp_2jmin
jmax0 = amp_2jmax

jmin1 = min(amp_2jmin,abs(amp_2jmin-2),abs(amp_2jmax-2))
jmax1 = amp_2jmax + 2

jmin3 = amp_2jmin
do l = 2, 6, 2
  jmin3 = min(jmin3,abs(amp_2jmin-l),abs(amp_2jmax-l))
enddo
jmax3 = amp_2jmax + 6

if ( is_good_J(1) ) then
  if ( ((misc_phys == 0) .and. (((good_J(1) - jmin3) < 0) .or. &
                                ((good_J(1) - jmax3) > 0))) .or. &
       ((misc_phys == 1) .and. (((good_J(1) - jmin1) < 0) .or. &
                                ((good_J(1) - jmax1) > 0))) .or. &
       ((misc_phys == 2) .and. (((good_J(1) - jmin0) < 0) .or. &
                                ((good_J(1) - jmax0) > 0))) ) then
    ierror = ierror + 1
    print '(1a)', 'The angular momentum of the left state is not compatible &
                  &with the angular momentum determined from the inputs.'
  endif
endif

if ( is_good_J(2) ) then
  if ( ((good_J(2) - amp_2jmin) < 0) .or. ((good_J(2) - amp_2jmax) > 0) ) then
    ierror = ierror + 1
    print '(1a)', 'The angular momentum of the right state is not compatible &
                  &with the angular momentum determined from the inputs.'
  endif
endif

!!!
!!! Parity
!!!

if ( is_good_P(1) ) then
  if ( (misc_phys == 2) .and. ((good_P(1) - 1) /= 0) ) then
    ierror = ierror + 1
    print '(1a)', 'The parity of the left state is not compatible with the &
                  &parity determined from the inputs.'
  endif
endif

if ( is_good_P(2) ) then
  if ( (misc_phys == 2) .and. ((good_P(2) - 1) /= 0) ) then
    ierror = ierror + 1
    print '(1a)', 'The parity of the right state is not compatible with the &
                  &parity determined from the inputs.'
  endif
endif

!!!
!!! Stops the code if an error has been found in the input file
!!!

if ( ierror /= 0 ) then
  !cmpi if ( paral_myrank == 0 ) then
  print "(a,1i2,a)", "The code has dectected ",ierror," incompatibily(ies) &
        &between the good quantum numbers and the input parameters."
  !cmpi endif 
  stop 
endif

end subroutine check_compatibilities

!------------------------------------------------------------------------------!
! subroutine reset_rotmatelem                                                  !
!                                                                              ! 
! Sets the rotated expectations values at the beginning of the iteration.      ! 
! Not really needed but is safer.                                              ! 
!------------------------------------------------------------------------------!
subroutine reset_rotmatelem   

if ( misc_phys == 0 ) then
  rot_over  = zzero  
  rot_pari  = zzero  
  rot_ener  = zzero  
  rot_prot  = zzero  
  rot_neut  = zzero  
  rot_prot2 = zzero  
  rot_neut2 = zzero  
  rot_amjx  = zzero  
  rot_amjy  = zzero  
  rot_amjz  = zzero  
  rot_amjx2 = zzero  
  rot_amjy2 = zzero  
  rot_amjz2 = zzero  
  rot_amsop = zzero  
  rot_amson = zzero  
  rot_istx  = zzero  
  rot_isty  = zzero  
  rot_istz  = zzero  
  rot_istx2 = zzero  
  rot_isty2 = zzero  
  rot_istz2 = zzero  
  rot_ra2p  = zzero  
  rot_ra2n  = zzero  
  rot_ra2m  = zzero  
  rot_Q1mp  = zzero  
  rot_Q1mn  = zzero  
  rot_Q2mp  = zzero  
  rot_Q2mn  = zzero  
  rot_Q3mp  = zzero  
  rot_Q3mn  = zzero  
  rot_M1ma  = zzero  
  rot_M2ma  = zzero  
  rot_occn  = zzero  
endif

end subroutine reset_rotmatelem   

!------------------------------------------------------------------------------!
! subroutine integrate_rotmatelem                                              !
!                                                                              ! 
! Sets the rotated expectations values to zero in the internal loop            ! 
!------------------------------------------------------------------------------!
subroutine integrate_rotmatelem

integer :: j, kj, mj, ikmj, p, ip, lambda, mu, xparP
real(r64) :: xphiZ, xphiN, xphiA, xalpJ, xbetJ, xgamJ, &
             wphiZ, wphiN, wphiA, walpJ, wbetJ, wgamJ, wparP,  &
             dwigbetJ
complex(r64) :: fac1, fac2, fac3, fac4, facmu
complex(r64), dimension(-3:3) :: fac5, fac6

!!! Initialization
lambda = 0
fac1 = zone
fac2 = zone
fac3 = zone
fac4 = zone
fac5 = zone
fac6 = zone

rot_amj2 = rot_amjx2 + rot_amjy2 + rot_amjz2
rot_ist2 = rot_istx2 + rot_isty2 + rot_istz2

if ( misc_phys == 0 ) then
  lambda = 6
elseif ( misc_phys == 1 ) then
  lambda = 2
endif

xphiZ = angle_phiZ(phiZ)
xphiN = angle_phiN(phiN)
xphiA = angle_phiA(phiA)
xalpJ = angle_alpJ(alpJ)
xbetJ = angle_betJ(betJ)
xgamJ = angle_gamJ(gamJ)
xparP = angle_parP(parP) 

wphiZ = weight_phiZ(phiZ)
wphiN = weight_phiN(phiN)
wphiA = weight_phiA(phiA)
walpJ = weight_alpJ(alpJ)
wbetJ = weight_betJ(betJ)
wgamJ = weight_gamJ(gamJ)
wparP = weight_parP(parP)

fac1 = wphiZ * exp(zimag * xphiZ * valence_Z) * &  
       wphiN * exp(zimag * xphiN * valence_N) * &
       wphiA * exp(zimag * xphiA * valence_A) 

!!! Loop on J,MJ,KJ,P (ifs for good quantum numbers)
do j = amp_2jmin, amp_2jmax, 2
  if ( (is_ampJ .neqv. is_betJ) .and. &
       ((is_good_J(1) .and. (abs(good_J(1) - j) > lambda)) .or. &
        (is_good_J(2) .and. (good_J(2) - j /= 0))) ) cycle
  fac2 = fac1 
  if ( is_betJ ) fac2 = fac1 * (j + one) / 2

  do mj = -j, j, 2
    if ( (is_ampMJ .neqv. is_alpJ) .and. is_good_MKJ(1) .and. & 
         (abs(good_MKJ(1) - mj) > lambda) ) cycle  
    fac3 = fac2 * wgamJ * exp(zimag * xalpJ * mj / 2.0d0)

    do kj = -j, j, 2
      if ( (is_ampKJ .neqv. is_gamJ) .and. is_good_MKJ(2) .and. & 
           (good_MKJ(2) /= kj) ) cycle 
      ikmj = tabid_J(kj,mj,j) 
      call dWigner(j,mj,kj,xbetJ,dwigbetJ) 
      fac4 = fac3 * wbetJ * dwigbetJ * &
                          wgamJ * exp(zimag * xgamJ * kj / 2.0d0)

      do p = pap_pmin, pap_pmax, 2
        if ( (is_papP .neqv. is_parP) .and. is_good_P(2) .and. & 
             (good_P(2) /= p) ) cycle
        ip = tabid_P(p) 

        do mu = -3, 3
          if ( (is_ampMJ .neqv. is_alpJ) .and. is_good_MKJ(1) ) then
            facmu = zone * kdelta(good_MKJ(1),2*mu+mj)
          else 
            facmu = zone
          endif

          fac5(mu) = facmu * fac4 * wparP * p**(xparP-1)
        enddo
        fac6 = fac5 * rot_over

        ! Remark: I could use additional symmetry restrictions but they make
        ! the code more difficult to read because are different for different
        ! operators. Only small gains expected, if any.
        ! if ( is_parP .and. is_good_P(1) .and. (good_P(1) /= +p) ) cycle
        ! if ( is_parP .and. is_good_P(1) .and. (good_P(1) /= -p) ) cycle
        ! if ( (is_ampMJ .neqv. is_alpJ) .and. is_good_MKJ(1) .and. &
        !      (good_MKJ(1) /= 2*mu+mj) ) cycle

        if ( misc_phys == 0 ) then 
          proj_over(ip,ikmj)  = proj_over(ip,ikmj)  + fac6(0) 
          proj_pari(ip,ikmj)  = proj_pari(ip,ikmj)  + fac5(0) * rot_pari
          proj_ener(ip,ikmj)  = proj_ener(ip,ikmj)  + fac6(0) * rot_ener
          proj_prot(ip,ikmj)  = proj_prot(ip,ikmj)  + fac6(0) * rot_prot
          proj_neut(ip,ikmj)  = proj_neut(ip,ikmj)  + fac6(0) * rot_neut
          proj_prot2(ip,ikmj) = proj_prot2(ip,ikmj) + fac6(0) * rot_prot2
          proj_neut2(ip,ikmj) = proj_neut2(ip,ikmj) + fac6(0) * rot_neut2
          proj_amjz(ip,ikmj)  = proj_amjz(ip,ikmj)  + fac6(0) * rot_amjz
          proj_amjz2(ip,ikmj) = proj_amjz2(ip,ikmj) + fac6(0) * rot_amjz2
          proj_amj2(ip,ikmj)  = proj_amj2(ip,ikmj)  + fac6(0) * rot_amj2
          proj_amsop(ip,ikmj) = proj_amsop(ip,ikmj) + fac6(0) * rot_amsop
          proj_amson(ip,ikmj) = proj_amson(ip,ikmj) + fac6(0) * rot_amson
          proj_istz(ip,ikmj)  = proj_istz(ip,ikmj)  + fac6(0) * rot_istz
          proj_istz2(ip,ikmj) = proj_istz2(ip,ikmj) + fac6(0) * rot_istz2
          proj_ist2(ip,ikmj)  = proj_ist2(ip,ikmj)  + fac6(0) * rot_ist2
          proj_ra2p(ip,ikmj)  = proj_ra2p(ip,ikmj)  + fac6(0) * rot_ra2p
          proj_ra2n(ip,ikmj)  = proj_ra2n(ip,ikmj)  + fac6(0) * rot_ra2n
          proj_ra2m(ip,ikmj)  = proj_ra2m(ip,ikmj)  + fac6(0) * rot_ra2m
        
          proj_Q1mp(:,ip,ikmj) = proj_Q1mp(:,ip,ikmj) + fac6(-1:1) * rot_Q1mp(:)
          proj_Q1mn(:,ip,ikmj) = proj_Q1mn(:,ip,ikmj) + fac6(-1:1) * rot_Q1mn(:)
          proj_Q2mp(:,ip,ikmj) = proj_Q2mp(:,ip,ikmj) + fac6(-2:2) * rot_Q2mp(:)
          proj_Q2mn(:,ip,ikmj) = proj_Q2mn(:,ip,ikmj) + fac6(-2:2) * rot_Q2mn(:)
          proj_Q3mp(:,ip,ikmj) = proj_Q3mp(:,ip,ikmj) + fac6(-3:3) * rot_Q3mp(:)
          proj_Q3mn(:,ip,ikmj) = proj_Q3mn(:,ip,ikmj) + fac6(-3:3) * rot_Q3mn(:)
          proj_M1ma(:,ip,ikmj) = proj_M1ma(:,ip,ikmj) + fac6(-1:1) * rot_M1ma(:)
          proj_M2ma(:,ip,ikmj) = proj_M2ma(:,ip,ikmj) + fac6(-2:2) * rot_M2ma(:)
        
          proj_occn(:,:,ip,ikmj) = proj_occn(:,:,ip,ikmj) + fac6(0) * & 
                                                            rot_occn(:,:)
        endif
      enddo
    enddo
  enddo
enddo

end subroutine integrate_rotmatelem

!------------------------------------------------------------------------------!
! subroutine open_file_rotmatelem                                              !
!                                                                              !
! Opens the file containing the rotated matrix elements calculated in a prev-  !
! ious calculation and reads the parameters used to perform said calculation.  !
!------------------------------------------------------------------------------!
subroutine open_file_rotmatelem

integer :: i=0, ialloc=0
character(14) :: prefix='rotmatelem.bin'
character(100) :: tmp_hamil
character(len=:), allocatable :: file_rotated
logical :: is_exist 
!cmpi character(10) :: nrank

!!! Writes the basic parameters 
if ( misc_frot == 1 ) then

  !cmpi write(nrank,'(i10)') paral_myrank
  !cmpi i = len_trim(adjustl(nrank))
  allocate (character(14+i) :: file_rotated)
  file_rotated = prefix 
  !cmpi file_rotated = prefix // trim(adjustl(nrank))

  open(utro, file=file_rotated, status='replace', action='write', &
       access='stream', form='unformatted', convert='little_endian') 

  !cmpi if ( paral_myrank == 0 ) then        
  tmp_hamil = adjustl(hamil_file)

  write(utro) tmp_hamil, hamil_type, misc_cutrot, & 
              HOsh_dim, HOsh_na, &
              bogo_labelL, bogo_labelR, bogo_npaL, bogo_npaR, &
              phiZ_dim, phiN_dim, phiA_dim, pnp_facpi, &
              alpJ_dim, betJ_dim, gamJ_dim, &
              parP_dim, &
              is_pnpZ, is_pnpN, is_pnpA, is_phiZ, is_phiN, is_phiA, &
              is_ampJ, is_ampMJ, is_ampKJ, is_alpJ, is_betJ, is_gamJ, &
              is_papP, is_parP, &  
              good_Z, good_N, good_A, &
              good_J, good_MKJ, &
              good_P, & 
              good_T, good_MKT, &
              is_good_Z, is_good_N, is_good_A, is_good_ZN, & 
              is_good_J, is_good_MKJ, &
              is_good_P, & 
              is_good_T, is_good_MKT
  !cmpi endif

!!! Reads the basic parameters 
elseif ( misc_frot == 2 ) then

  allocate (character(14) :: file_rotated)
  file_rotated = prefix

  inquire(file=file_rotated, exist=is_exist)
  if ( .not. is_exist ) then
    print "(a,a,a)", "The file containing the rotated matrix elements &
           &(file_rotated) = ", file_rotated, " can not be found!"
    stop
  endif 

  open(utro, file=file_rotated, status='old', action='read', &
       access='stream', form='unformatted', convert='little_endian') 

  read(utro) tmp_hamil, hamil_type, frot_cutrot, & 
             HOsh_dim

  allocate( HOsh_na(HOsh_dim), stat=ialloc )
  if ( ialloc /= 0 ) stop 'Error during allocation of shells (frot)'

  read(utro) HOsh_na, & 
             bogo_labelL, bogo_labelR, bogo_npaL, bogo_npaR, &
             phiZ_dim, phiN_dim, phiA_dim, pnp_facpi, &
             alpJ_dim, betJ_dim, gamJ_dim, &
             parP_dim, &
             is_pnpZ, is_pnpN, is_pnpA, is_phiZ, is_phiN, is_phiA, &
             is_ampJ, is_ampMJ, is_ampKJ, is_alpJ, is_betJ, is_gamJ, &
             is_papP, is_parP, &  
             good_Z, good_N, good_A, &
             good_J, good_MKJ, &
             good_P, & 
             good_T, good_MKT, &
             is_good_Z, is_good_N, is_good_A, is_good_ZN, & 
             is_good_J, is_good_MKJ, &
             is_good_P, & 
             is_good_T, is_good_MKT

  i = len_trim(adjustl(tmp_hamil))
  allocate (character(i) :: frot_hamil)
  frot_hamil = trim(adjustl(tmp_hamil))

  if ( .not. is_phiZ ) phiZ_dim = 0
  if ( .not. is_phiN ) phiN_dim = 0
  if ( .not. is_phiA ) phiA_dim = 0
  if ( .not. is_alpJ ) alpJ_dim = 0
  if ( .not. is_betJ ) betJ_dim = 0
  if ( .not. is_gamJ ) gamJ_dim = 0
  if ( .not. is_parP ) parP_dim = 0

endif

end subroutine open_file_rotmatelem

!------------------------------------------------------------------------------!
! subroutine readwrite_rotmatelem                                              !
!                                                                              ! 
! Reads or writes the rotated expectation values at each angle from file.      ! 
!                                                                              !  
! Input: iang = position in the loop over all angles                           !
!------------------------------------------------------------------------------!
subroutine readwrite_rotmatelem(iang)

integer(i64), intent(in) :: iang
integer(i64) :: head, offset=0, mypos
integer(i64), dimension(3) :: base(0:2)

!!! Number of octets written per angle for each physics case
head = 364 + 4 * HOsh_dim
base(0) = 956 + 32 * HOsh_dim
!base(1) = 140
!base(2) = 92 

!!! Determines the position before reading
if ( iang == loop_ini ) then
  !cmpi offset = paral_myoffset
  mypos = 1 + head + offset * base(misc_phys)
else
  inquire(utro,pos=mypos)
endif 

!!! Writes/reads the angles and rotated matrix elements
if ( misc_frot == 1 ) then
  write(utro) phiZ, phiN, phiA, &
              alpJ, betJ, gamJ, &
              parP

  if ( misc_phys == 0 ) then
    write(utro) rot_over,  rot_pari,  rot_ener,  &
                rot_prot,  rot_neut,  rot_prot2, rot_neut2, & 
                rot_amjz,  rot_amjx2, rot_amjy2, rot_amjz2, &
                rot_amsop, rot_amson, & 
                rot_istz,  rot_istx2, rot_isty2, rot_istz2, &
                rot_ra2p,  rot_ra2n, rot_ra2m, &
                rot_Q1mp, rot_Q1mn, rot_Q2mp, rot_Q2mn, rot_Q3mp, rot_Q3mn, &                 
                rot_M1ma, rot_M2ma, rot_occn
  endif
elseif ( misc_frot == 2 ) then
  read(utro, pos=mypos) phiZ, phiN, phiA, & 
                        alpJ, betJ, gamJ, &
                        parP

  if ( misc_phys == 0 ) then
    read(utro) rot_over,  rot_pari,  rot_ener, & 
               rot_prot,  rot_neut,  rot_prot2, rot_neut2, & 
               rot_amjz,  rot_amjx2, rot_amjy2, rot_amjz2, & 
               rot_amsop, rot_amson, &
               rot_istz,  rot_istx2, rot_isty2, rot_istz2, &
               rot_ra2p,  rot_ra2n, rot_ra2m, &
               rot_Q1mp, rot_Q1mn, rot_Q2mp, rot_Q2mn, rot_Q3mp, rot_Q3mn, &
               rot_M1ma, rot_M2ma, rot_occn  
  endif
endif

!!! Closes the file if this is the last angle
if ( iang == loop_fin ) then
  close(utro, status='keep')
endif

end subroutine readwrite_rotmatelem

!------------------------------------------------------------------------------!
! subroutine print_rotmatelem                                                  !
!                                                                              ! 
! Prints information obtained from the file containing the rotated matrix      ! 
! elements about the parameters used for the original calculation.             ! 
!------------------------------------------------------------------------------!
subroutine print_rotmatelem

character(3), dimension(2) :: ch_good_Z, ch_good_N, ch_good_A, ch_good_ZN,  &
                              ch_good_J, ch_good_MKJ, ch_good_P, ch_good_T, &
                              ch_good_MKT
character(len=*), parameter :: format1 = "(1a20,2x,1a)" , &
                               format2 = "(1a20,1x,1es10.3)" , &
                               format3 = "(1a20,1x,1i19,1x,1i19)" , &
                               format4 = "(1a20,9x,1i3,17x,1i3)", &
                               format5 = "(1a20,6x,1a3,1x,1i5,11x,1a3,1x,1i5)",&
                               format6 = "(1a20,6x,1a3,17x,1a3)" 

ch_good_Z   = ' no'
ch_good_N   = ' no'
ch_good_A   = ' no'
ch_good_ZN  = ' no'
ch_good_J   = ' no'
ch_good_P   = ' no'
ch_good_MKJ = ' no'
ch_good_T   = ' no'
ch_good_MKT = ' no'

where( is_good_Z   ) ch_good_Z   = "yes"
where( is_good_N   ) ch_good_N   = "yes"
where( is_good_A   ) ch_good_A   = "yes"
where( is_good_ZN  ) ch_good_ZN  = "yes"
where( is_good_J   ) ch_good_J   = "yes"
where( is_good_P   ) ch_good_P   = "yes"
where( is_good_MKJ ) ch_good_MKJ = "yes"
where( is_good_T   ) ch_good_T   = "yes"
where( is_good_MKT ) ch_good_MKT = "yes"

print '(60("%"),/,19x,"ROTATED MATRIX ELEMENTS",18x,/,60("%"),//, &
      & 3x,"Description",15x,"Values",/,40("-"))'
print format1, 'Hamiltonian name    ', frot_hamil
print format2, 'Cutoff rot. overlaps', frot_cutrot
print format3, 'Label of states L/R ', bogo_labelL, bogo_labelR
print format4, 'Total number parity ', bogo_npaL, bogo_npaR
print format5, 'Good Z              ', ch_good_Z(1), good_Z(1), &
                                       ch_good_Z(2), good_Z(2)
print format5, 'Good N              ', ch_good_N(1), good_N(1), & 
                                       ch_good_N(2), good_N(2)
print format5, 'Good A              ', ch_good_A(1), good_A(1), &
                                       ch_good_A(2), good_A(2)
print format6, 'Separate Z/N        ', ch_good_ZN(1), ch_good_ZN(2)
print format5, 'Good J              ', ch_good_J(1), good_J(1), &
                                       ch_good_J(2), good_J(2)
print format5, 'Good MKJ            ', ch_good_MKJ(1), good_MKJ(1), &
                                       ch_good_MKJ(2), good_MKJ(2)
print format5, 'Good P              ', ch_good_P(1), good_P(1), & 
                                       ch_good_P(2), good_P(2)
print format5, 'Good T              ', ch_good_T(1), good_T(1), &
                                       ch_good_T(2), good_T(2)
print format5, 'Good MKT            ', ch_good_MKT(1), good_MKT(1), &
                                       ch_good_MKT(2), good_MKT(2)
print*,' '

end subroutine print_rotmatelem

!------------------------------------------------------------------------------!
! subroutine generate_rotation_spatial                                         !
!                                                                              ! 
! Build the transformation matrix corresponding to the rotation over Euler     ! 
! angles. In the HO basis, the Euler rotation reads                            ! 
! | c > = | n_c l_c j_c m_c mt_c >                                             !
! | d > = | n_d l_d j_d m_d mt_d >                                             !
! < c | R(a,b,g) | d > = e^{i a m_c} * d^{j_c}_{m_c m_d}(b) * e^{i g m_d} *    ! 
!                        delta_{n_c n_d} * delta_{j_c j_d} * delta_{l_c l_d} * ! 
!                        delta_{mt_c mt_d}                                     ! 
!                                                                              !  
! Input: ndim = dimension of the sp basis                                      !
!        xaplha,xbeta,xgamma = Euler angles                                    !  
!        RI = initial rotation matrix                                          !  
!                                                                              !  
! Output: RF = final rotation matrix                                           !  
!------------------------------------------------------------------------------!
subroutine generate_rotation_spatial(xalpha,xbeta,xgamma,RI,RF,ndim)

integer, intent(in) :: ndim
real(r64), intent(in) :: xalpha, xbeta, xgamma
complex(r64), dimension(:,:), intent(in) :: RI
complex(r64), dimension(:,:), intent(out) :: RF
integer :: i, j, jm, jt, jshe, ij, im, it, ishe
real(r64) :: wigbet
complex(r64) :: wigalp, wiggam
complex(r64), dimension(ndim,ndim) :: D 

D = zzero

do j = 1, ndim
  jm = HOsp_2mj(j)     
  jt = HOsp_2mt(j) 
  jshe = jt * HOsp_sh(j)

  do i = 1, ndim
    ij = HOsp_2j(i)
    im = HOsp_2mj(i)     
    it = HOsp_2mt(i) 
    ishe = it * HOsp_sh(i)

    if ( jshe == ishe ) then  
      wigalp = exp(-zimag * xalpha * im / 2.0d0)
      wiggam = exp(-zimag * xgamma * jm / 2.0d0) 
      call dWigner(ij,im,jm,xbeta,wigbet)
      D(i,j)  = wigalp * wigbet * wiggam
    endif

  enddo
enddo

call zgemm('n','n',ndim,ndim,ndim,zone,RI,ndim,D,ndim,zzero,RF,ndim)

end subroutine generate_rotation_spatial

!------------------------------------------------------------------------------!
! subroutine generate_rotation_parity                                          !
!                                                                              ! 
! Multiply the matrix RI with the parity transformation, the result is RF.     ! 
! In the HO basis, the parity operation reads                                  ! 
! | a > = | n_a l_a j_a m_a mt_a >                                             !
! | b > = | n_b l_b j_b m_b mt_b >                                             !
! < a | P | b > = (-1)**l_a delta_ab                                           ! 
!                                                                              !  
! Input: ndim = dimension of the sp basis                                      !
!        ip = parity angle                                                     !  
!        RI = initial rotation matrix                                          !  
!                                                                              !  
! Output: RF = final rotation matrix                                           !  
!------------------------------------------------------------------------------!
subroutine generate_rotation_parity(ip,RI,RF,ndim)

integer, intent(in) :: ip, ndim
complex(r64), dimension(ndim,ndim), intent(in) :: RI
complex(r64), dimension(ndim,ndim), intent(out) :: RF
integer :: i, j

if ( ip == 1 ) then !!! no parity operation
  RF = RI
else                !!! parity operation
  do j = 1, ndim
    do i = 1, ndim
      RF(i,j) = (-1)**(HOsp_l(i)) * RI(i,j)
    enddo
  enddo
endif

end subroutine generate_rotation_parity

!------------------------------------------------------------------------------!
! subroutine generate_transfo_rotation_gauge                                    !
!                                                                              ! 
! Multiply the matrix RI with the protons+neutrons gauge transformation, the   ! 
! result is RF.                                                                ! 
! In the HO basis, the gauge rotation reads                                    ! 
! | a > = | n_a l_a j_a m_a mt_a >                                             !
! | b > = | n_b l_b j_b m_b mt_b >                                             !
! < a | R(phi_p, phi_n) | b > =  e^(i phi_{mt_a}) delta_ab                     ! 
!------------------------------------------------------------------------------!
subroutine generate_rotation_gauge(xphiZ,xphiN,xphiA,RI,RF,ndim)

integer, intent(in) :: ndim
real(r64), intent(in) :: xphiA, xphiZ, xphiN
complex(r64), dimension(ndim,ndim), intent(in) :: RI
complex(r64), dimension(ndim,ndim), intent(out) :: RF
integer :: i, j
complex(r64) :: phaseA, phaseZ, phaseN 

phaseZ = exp(-zimag * xphiZ) 
phaseN = exp(-zimag * xphiN) 
phaseA = exp(-zimag * xphiA) 

do j = 1, ndim
  do i = 1, ndim/2
    RF(i,j) = phaseA * phaseZ * RI(i,j)
  enddo

  do i = 1+ndim/2, ndim
    RF(i,j) = phaseA * phaseN * RI(i,j)
  enddo
enddo

end subroutine generate_rotation_gauge

!------------------------------------------------------------------------------!
! subroutine rotate_wavefunction                                               !
!                                                                              !
! Computes the Bogoliubov matrices of the rotated state                        !
!   UF = R UI                                                                  !
!   VF = R^* VI                                                                !
! Note that here R is an arbitrary rotation matrix (e.g. Euler+Parity+Gauge)   !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        UI,VI = initial Bogoliubov matrices                                   !
!        R = rotation matrix                                                   !
!                                                                              !
! Output: UF,VF = final Bogoliubov matrices                                    ! 
!------------------------------------------------------------------------------!
subroutine rotate_wavefunction(R,UI,VI,UF,VF,ndim)        

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: R, UI, VI
complex(r64), dimension(ndim,ndim), intent(out) :: UF, VF
complex(r64), dimension(ndim,ndim) :: Rc

Rc = conjg(R)

call zgemm('n','n',ndim,ndim,ndim,zone, R,ndim,UI,ndim,zzero,UF,ndim)
call zgemm('n','n',ndim,ndim,ndim,zone,Rc,ndim,VI,ndim,zzero,VF,ndim)

end subroutine rotate_wavefunction 

!------------------------------------------------------------------------------!
! subroutine calculate_thouless                                                !
!                                                                              ! 
! According to Ring&Schuck (equation E54, p619) the transition densities can be!  
! calculated as                                                                !  
!        rho_ij = <0|a^+_j a_i|1>   = V_1^* (U^T)^-1 V_0^T                     !  
!    kappa01_ij = <0|a_j a_i|1>     = V_1^* (U^T)^-1 U_0^T                     !  
!  kappa10_ij^* = <0|a^+_i a^+_j|1> = - U_1^* (U^T)^-1 V_0^T                   !
! where                                                                        ! 
!  U = U_1^dagger U_0 + V_1^dagger V_0                                         ! 
!                                                                              ! 
! In our case, U_1 and V_1 are to be replaced by the rotated transformations   ! 
! Ubar and Vbar. We then define                                                ! 
!  Utilde = Ubar (U^dagger)^-1                                                 ! 
!  Vtilde = Vbar (U^dagger)^-1                                                 ! 
! to be injected in the routine CALC_DENSITIES                                 ! 
!                                                                              !  
! Input: ndim = dimension of the sp basis                                      !
!        UL,VL = left  Bogoliubov matrices                                     !  
!        UR,VR = right Bogoliubov matrices                                     !  
!------------------------------------------------------------------------------!
subroutine calculate_thouless(UL,VL,UR,VR,ndim)                 

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: UL, VL, UR, VR
integer :: info1, info2
integer, dimension(ndim) :: ipiv
complex(r64), dimension(ndim) :: zwork
complex(r64), dimension(ndim,ndim) :: UD 

!!! U^dagger
call zgemm('c','n',ndim,ndim,ndim,zone,VL,ndim,VR,ndim,zzero,UD,ndim)
call zgemm('c','n',ndim,ndim,ndim,zone,UL,ndim,UR,ndim, zone,UD,ndim)

!!! (U^dagger)^-1  (then UD will contain its inverse)
call zgetrf(ndim,ndim,UD,ndim,ipiv,info1)
call zgetri(ndim,UD,ndim,ipiv,zwork,ndim,info2)

!!! Utilde and Vtilde
call zgemm('n','n',ndim,ndim,ndim,zone,UR,ndim,UD,ndim,zzero,Util,ndim)
call zgemm('n','n',ndim,ndim,ndim,zone,VR,ndim,UD,ndim,zzero,Vtil,ndim)
   
end subroutine calculate_thouless

!------------------------------------------------------------------------------!
! subroutine print_projmatelem_states                                          !
!                                                                              !
!------------------------------------------------------------------------------!
subroutine print_projmatelem_states  

integer :: j, mj, kj, ikmj, p, ip, nj, nmj, nkj, np, &
           protv, neutv, nuclv, amjzv, istzv
real(r64) :: eps=1.d-12, xj, xt
complex(r64) :: over, ener, pari, prot, neut, nucl, prot2, neut2, nucl2, & 
                amjz, amjz2, amj2, amsop, amson, istz, istz2, ist2, &
                sumtot1, sumtotE, sumcomp1, sumcompE
character(4) :: chJ, chMJ, chKJ
character(len=*), parameter :: format1 = "(2x,1a3,2x,1a3,2x,1a3,2x,1a3,1x, & 
                                          &1f12.8,1x,1f12.5,1f12.7,1i4,1f12.7, &
                                          &1i4,1f12.7,1i4,2f10.5,1i4,3f10.5, &
                                          &1i4)", &
                               format2 = "(2x,1a3,2x,1a3,1x,2f12.8,2f12.5)"

protv = 0
neutv = 0
nuclv = 0
amjzv = 0
istzv = 0

if ( (-1)**amp_2jmin == 1 ) then
  chJ  = '   J'
  chMJ = '  MJ'
  chKJ = '  KJ'
else
  chJ  = '   J'
  chMJ = '2*MJ'
  chKJ = '2*KJ'
endif

print '(/,60("%"),/,18x,"PROJECTED MATRIX ELEMENTS",17x,/,60("%"),/)'

!!! Print individually each component        
print "('All non-vanishing projected components',/,38('='),//, &
       &'v(x) = floor(log10(variance(x)))',//, &
       &1x,1a4,1x,1a4,1x,1a4,4x,'P',1x,'|',5x,'1',5x,'|',6x,'E',5x, &
       &'|',5x,'Z',7x,'v |',5x,'N',7x,'v |',5x,'A',7x,'v |',4x,'J',4x, &
       &'|',4x,'Jz',5x,'v |',4x,'P',4x,'|',4x,'T',4x,'|',4x,'Tz',5x,'v',&
       &/,152('-'))", &
       chJ,chMJ,chKJ

!!! Print individually each component        
do j = amp_2jmin, amp_2jmax, 2
  do p = pap_pmin, pap_pmax, 2
    ip = tabid_P(p) 
    do mj = -j, j, 2
      do kj = -j, j, 2
        ikmj = tabid_J(kj,mj,j) 

        over = proj_over(ip,ikmj)
        if ( abs(over) > eps ) then 
          ener  = proj_ener(ip,ikmj)  / over
          pari  = proj_pari(ip,ikmj)  / over
          prot  = proj_prot(ip,ikmj)  / over
          neut  = proj_neut(ip,ikmj)  / over
          prot2 = proj_prot2(ip,ikmj) / over
          neut2 = proj_neut2(ip,ikmj) / over
          amjz  = proj_amjz(ip,ikmj)  / over
          amjz2 = proj_amjz2(ip,ikmj) / over
          amj2  = proj_amj2(ip,ikmj)  / over
          amsop = proj_amsop(ip,ikmj) / over
          amson = proj_amson(ip,ikmj) / over
          istz  = proj_istz(ip,ikmj)  / over
          istz2 = proj_istz2(ip,ikmj) / over
          ist2  = proj_ist2(ip,ikmj)  / over

          nucl  = prot + neut
          nucl2 = 2.0d0 * (prot2 + neut2 - 2.0d0*istz2)

          !Remark: could add the parity: var(P) = 1 - <P>**2/<1>
          if ( misc_part /= 1 ) then
            protv = floor(log10(abs(prot2-prot**2)))
            neutv = floor(log10(abs(neut2-neut**2)))
            nuclv = floor(log10(abs(nucl2-nucl**2)))
            amjzv = floor(log10(abs(amjz2-amjz**2)))
            istzv = floor(log10(abs(istz2-istz**2)))
          endif
        else
          cycle
        endif
        xj = 0.5d0 * (-1 + sqrt(1+4*abs(amj2)))
        xt = 0.5d0 * (-1 + sqrt(1+4*abs(ist2)))

        if ( is_ampJ ) then 
          nj = j
          if ( (-1)**nj == 1 ) nj = j/2
        else
          nj = 1000
        endif
        
        if ( is_ampMJ ) then 
          nmj = mj
          if ( (-1)**nmj == 1 ) nmj = mj/2
        else
          nmj = 1000
        endif
        
        if ( is_ampKJ ) then 
          nkj = kj
          if ( (-1)**nkj == 1 ) nkj = kj/2
        else
          nkj = 1000
        endif
        
        if ( is_papP ) then
          np = p
        else
          np = 1000
        endif

        print format1, string3(nj), string3(nmj), string3(nkj), string3(np), & 
              real(over), real(ener), real(prot), protv, real(neut), neutv, &
              real(nucl), nuclv, &
              xj, real(amjz), amjzv, real(pari), xt, real(istz), istzv
      enddo
    enddo
  enddo
enddo

!!! Print sum of components for a given J+P
print "(/,'Sum of projected components for J/P',/,35('='),//, &        
       &1x,1a4,4x,'P',1x,'|',11x,'1',11x,'|',11x,'E',10x,/,59('-'))",chJ

sumtot1 = zzero
sumtotE = zzero

do j = amp_2jmin, amp_2jmax, 2
  do p = pap_pmin, pap_pmax, 2
    ip = tabid_P(p) 
    do kj = -j, j, 2
      ikmj = tabid_J(kj,kj,j) 
      sumtot1 = sumtot1 + proj_over(ip,ikmj) 
      sumtotE = sumtotE + proj_ener(ip,ikmj) 
    enddo
  enddo
enddo
if ( abs(sumtot1) > eps ) sumtotE = sumtotE / sumtot1

do p = pap_pmin, pap_pmax, 2
  ip = tabid_P(p) 
  do j = amp_2jmin, amp_2jmax, 2
    sumcomp1 = zzero 
    sumcompE = zzero 
    do kj = -j, j, 2
      ikmj = tabid_J(kj,kj,j) 
      sumcomp1 = sumcomp1 + proj_over(ip,ikmj) 
      sumcompE = sumcompE + proj_ener(ip,ikmj)
    enddo
    
    if ( abs(sumcomp1) < eps ) cycle 
    sumcompE = sumcompE / sumtot1
    
    if ( is_ampJ ) then 
      nj = j
      if ( (-1)**nj == 1 ) nj = j/2
    else
      nj = 1000
    endif
    
    if ( is_papP ) then
      np = p
    else
      np = 1000
    endif

    print format2, string3(nj), string3(np), real(sumcomp1), aimag(sumcomp1), & 
                                             real(sumcompE), aimag(sumcompE)
  enddo
enddo
print "(59('-'),/,4x,'Total',2x,2f12.8,2f12.5)", &
       real(sumtot1), aimag(sumtot1), real(sumtotE), aimag(sumtotE)

!!! Print sum of components for a given K
print "(/,'Sum of projected components for KJ/P',/,36('='),//, &        
       &1x,1a4,4x,'P',1x,'|',11x,'1',11x,'|',11x,'E',10x,/,59('-'))",chKJ

do p = pap_pmin, pap_pmax, 2
  ip = tabid_P(p) 
  do kj = -amp_2jmax, amp_2jmax, 2 
    sumcomp1 = zzero 
    sumcompE = zzero 
    do j = amp_2jmin, amp_2jmax, 2
      if ( j < abs(kj) ) cycle
      ikmj = tabid_J(kj,kj,j) 
      sumcomp1 = sumcomp1 + proj_over(ip,ikmj) 
      sumcompE = sumcompE + proj_ener(ip,ikmj)
    enddo

    if ( abs(sumcomp1) < eps ) cycle 
    sumcompE = sumcompE / sumtot1

    if ( is_ampKJ ) then 
      nkj = kj
      if ( (-1)**nkj == 1 ) nkj = kj/2
    else
      nkj = 1000
    endif
    
    if ( is_papP ) then
      np = p
    else
      np = 1000
    endif
    
    print format2, string3(nkj), string3(np), real(sumcomp1), aimag(sumcomp1), & 
                                              real(sumcompE), aimag(sumcompE)
  enddo
enddo
print "(59('-'),/,4x,'Total',2x,2f12.8,2f12.5)", &
       real(sumtot1), aimag(sumtot1), real(sumtotE), aimag(sumtotE)

end subroutine print_projmatelem_states 

!------------------------------------------------------------------------------!
! subroutine write_projmatelem_states                                          !
!                                                                              ! 
! Writes in a file the expectation values of a selection of operators for all  !
! projected states.                                                            ! 
! BB: I may add N^2, Z^2, A, A^2, Jz, Tz                                       ! 
!------------------------------------------------------------------------------!
subroutine write_projmatelem_states

integer :: j, mj, kj, ikmj, p, ip 

open(utst, file='projmatelem_states.bin', status='replace', action='write', &
           form='unformatted')

do j = amp_2jmin, amp_2jmax, 2           
  do mj = -j, j, 2            
    do kj = -j, j, 2         
      ikmj = tabid_J(kj,mj,j) 
      do p = pap_pmin, pap_pmax, 2
        ip = tabid_P(p)
        write(utst) bogo_labelL, bogo_labelR, j, mj, kj, p, &
                    real(proj_over(ip,ikmj)),  real(proj_ener(ip,ikmj)), &
                    real(proj_pari(ip,ikmj)),  real(proj_prot(ip,ikmj)), &
                    real(proj_neut(ip,ikmj)),  real(proj_amj2(ip,ikmj)), &
                    real(proj_amsop(ip,ikmj)), real(proj_amson(ip,ikmj)), &
                    real(proj_ist2(ip,ikmj)),  real(proj_ra2p(ip,ikmj)), &
                    real(proj_ra2n(ip,ikmj)),  real(proj_ra2m(ip,ikmj))
      enddo
    enddo
  enddo
enddo

close(utst, status='keep')

end subroutine write_projmatelem_states

!------------------------------------------------------------------------------!
! subroutine write_projmatelem_occnumb                                         !
!                                                                              ! 
! Writes in a file the occupation numbers for all projected states.            !
!                                                                              ! 
! Remark: I could change the other writing routines to be like this one.       ! 
! Also, this could be merged with write_projmatelem_states, but might result   ! 
! in a large file.                                                             ! 
!------------------------------------------------------------------------------!
subroutine write_projmatelem_occnumb   

integer :: j, mj, kj, ikmj, p, ip, projdim

open(utoc, file='projmatelem_occnumb.bin', status='replace', action='write', &
           form='unformatted')

!!! Counts the number of states
projdim = 0

do j = amp_2jmin, amp_2jmax, 2           
  projdim = projdim + (j + 1)**2
enddo

projdim = projdim * (pap_pmax - pap_pmin + 2) / 2

!!! Writes the occupation numbers
write(utoc) projdim, hamil_type, HOsh_dim, HOsh_na

do j = amp_2jmin, amp_2jmax, 2           
  do mj = -j, j, 2            
    do kj = -j, j, 2         
      ikmj = tabid_J(kj,mj,j) 
      do p = pap_pmin, pap_pmax, 2
        ip = tabid_P(p)
        write(utoc) bogo_labelL, bogo_labelR, j, mj, kj, p, &
                    real(proj_occn(:,:,ip,ikmj))
      enddo
    enddo
  enddo
enddo

close(utoc, status='keep')

end subroutine write_projmatelem_occnumb   

!------------------------------------------------------------------------------!
! subroutine write_projmatelem_electromagnetic                                 !
!                                                                              ! 
! Writes in a file the matrix elements for electromagnetic transitions between !
! the projected states.                                                        ! 
!------------------------------------------------------------------------------!
subroutine write_projmatelem_electromagnetic

integer :: lambda, mu, ji, mji, kji, jf, kjf, pi, ikmj, ip, uteX, utmX, ialloc=0
real(r64) :: cg
complex(r64) :: redmat_p, redmat_n, redmat_a
complex(r64), dimension(:,:,:), allocatable :: AP, AN, AA

!!!
!!! Electric multipoles
!!!

open(ute1, file='projmatelem_E1.bin', status='replace', action='write', &
           form='unformatted')
open(ute2, file='projmatelem_E2.bin', status='replace', action='write', &
           form='unformatted')
open(ute3, file='projmatelem_E3.bin', status='replace', action='write', &
           form='unformatted')

allocate( AP(-3:3,tabid_P_dim,tabid_J_dim), AN(-3:3,tabid_P_dim,tabid_J_dim), &
          stat=ialloc, source=zzero )
if ( ialloc /= 0 ) stop 'Error during allocation of ELM matrices'

do lambda = 2, 6, 2
  if ( lambda == 2 ) then
    uteX = ute1
    AP(-lambda/2:lambda/2,:,:) = proj_Q1mp(-lambda/2:lambda/2,:,:)
    AN(-lambda/2:lambda/2,:,:) = proj_Q1mn(-lambda/2:lambda/2,:,:)
  elseif ( lambda == 4 ) then
    uteX = ute2
    AP(-lambda/2:lambda/2,:,:) = proj_Q2mp(-lambda/2:lambda/2,:,:)
    AN(-lambda/2:lambda/2,:,:) = proj_Q2mn(-lambda/2:lambda/2,:,:)
  else
    uteX = ute3
    AP(-lambda/2:lambda/2,:,:) = proj_Q3mp(-lambda/2:lambda/2,:,:)
    AN(-lambda/2:lambda/2,:,:) = proj_Q3mn(-lambda/2:lambda/2,:,:)
  endif

  do pi = pap_pmin, pap_pmax, 2
    ip = tabid_P(pi)                                                       
    do ji = amp_2jmin, amp_2jmax, 2                                                           
      do kji = -ji, ji, 2                                                            
        do jf = max(amp_2jmin,abs(ji-lambda)), min(amp_2jmax,ji+lambda), 2 
          do kjf = -jf, jf, 2                                                        
            redmat_p = zzero                                                 
            redmat_n = zzero                                                 
            do mji = -ji, ji, 2                                                        
              ikmj = tabid_J(kji,mji,ji)                                                       
              do mu = -lambda, lambda, 2                                                        
                call ClebschGordan(ji,lambda,jf,mji,mu,kjf,cg)                                        
                redmat_p = redmat_p + cg * AP(mu/2,ip,ikmj)
                redmat_n = redmat_n + cg * AN(mu/2,ip,ikmj)
              enddo ! mu                                                         
            enddo ! mji                                                         
            redmat_p = redmat_p * sqrt(jf + one)
            redmat_n = redmat_n * sqrt(jf + one)
            write(uteX) bogo_labelL, bogo_labelR, &
                        jf, kjf, pi*(-1)**(lambda/2), ji, kji, pi, &
                        real(redmat_p), real(redmat_n)
          enddo ! kjf                                                          
        enddo ! jf                                                            
      enddo ! kji                                                              
    enddo ! ji                                                              
  enddo ! pi                                                              
enddo ! lambda                                                           


close(ute1, status='keep')
close(ute2, status='keep')
close(ute3, status='keep')

!!!
!!! Magnetic multipoles
!!!

open(utm1, file='projmatelem_M1.bin', status='replace', action='write', &
           form='unformatted')
open(utm2, file='projmatelem_M2.bin', status='replace', action='write', &
           form='unformatted')

allocate( AA(-2:2,tabid_P_dim,tabid_J_dim), stat=ialloc, source=zzero )
if ( ialloc /= 0 ) stop 'Error during allocation of ELM matrices'

do lambda = 2, 4, 2
  if ( lambda == 2 ) then
    utmX = utm1
    AA(-lambda/2:lambda/2,:,:) = proj_M1ma(-lambda/2:lambda/2,:,:)
  else                                                       
    utmX = utm2
    AA(-lambda/2:lambda/2,:,:) = proj_M2ma(-lambda/2:lambda/2,:,:)
  endif

  do pi = pap_pmin, pap_pmax, 2
    ip = tabid_P(pi)                                                       
    do ji = amp_2jmin, amp_2jmax, 2                                                           
      do kji = -ji, ji, 2                                                            
        do jf = max(amp_2jmin,abs(ji-lambda)), min(amp_2jmax,ji+lambda), 2
          do kjf = -jf, jf, 2                                                        
            redmat_a = zzero                                                 
            do mji = -ji, ji, 2                                                        
              ikmj = tabid_J(kji,mji,ji)
              do mu = -lambda, lambda, 2                                                        
                call ClebschGordan(ji,lambda,jf,mji,mu,kjf,cg)
                redmat_a = redmat_a + cg * AA(mu/2,ip,ikmj)
              enddo ! mu                                                         
            enddo ! mji                                                         
            redmat_a = redmat_a * sqrt(jf + one)
            write(utmX) bogo_labelL, bogo_labelR, &
                        jf, kjf, pi*(-1)**(1+lambda/2), ji, kji, pi, &
                        real(redmat_a)
          enddo ! kjf                                                          
        enddo ! jf                                                            
      enddo ! kji                                                              
    enddo ! ji                                                              
  enddo ! pi                                                              
enddo ! lambda                                                           

close(utm1, status='keep')
close(utm2, status='keep')

end subroutine write_projmatelem_electromagnetic

!------------------------------------------------------------------------------!
! subroutine print_files_complementary                                         !
!                                                                              ! 
! Prints the names of the complementary files created during the run.          ! 
!------------------------------------------------------------------------------!
subroutine print_files_complementary

integer :: i=0

print '(/,60("%"),/,20x,"COMPLEMENTARY FILES",21x,/,60("%"),/)'
print '(6x,"Description",16x,"File",/,48("-"))'

if ( hamil_read == 0 ) then
  if ( misc_phys == 0 ) then
    print*,"Hamiltonian reduced    : ", hamil_fred
  endif
endif

if ( misc_frot == 1 ) then 
!cmpi   i = 1
  if ( i == 0 ) then
    print*,"Rotated mat. elem.     : rotmatelem.bin"
  else
    print*,"Rotated mat. elem.     : rotmatelem.bin*" 
  endif
endif

if ( misc_phys == 0 ) then  
  print*,"Proj. mat. elem. states    : projmatelem_states.bin"
  print*,"Proj. mat. elem. occ. numb.: projmatelem_occnumb.bin"
  print*,"Proj. mat. elem. E1        : projmatelem_E1.bin"
  print*,"Proj. mat. elem. E2        : projmatelem_E2.bin"
  print*,"Proj. mat. elem. E3        : projmatelem_E3.bin"
  print*,"Proj. mat. elem. M1        : projmatelem_M1.bin"
  print*,"Proj. mat. elem. M2        : projmatelem_M2.bin"
endif

end subroutine print_files_complementary

!------------------------------------------------------------------------------!
! subroutine reduce_projmatelem                                                !
!                                                                              !
! Sum all arrays containing projected rotated matrix elements from all         !
! processes into process 0.                                                    !
!                                                                              !
! Input: ndimP = dimension projection on parity                                !
!        ndimJ = dimension projection on angular momementum                    !
!        sdim  = dimension of HO shells                                        !
!------------------------------------------------------------------------------!
!cmpi subroutine reduce_projmatelem(ndimP,ndimJ,sdim)
!cmpi 
!cmpi integer, intent(in) :: ndimP, ndimJ, sdim
!cmpi integer :: ndim1, ndim3, ndim5, ndim7, sdim2, ierr=0
!cmpi complex(r64), dimension(ndimP,ndimJ) :: tmp_over,  tmp_pari,  tmp_ener,  &
!cmpi                                         tmp_prot,  tmp_neut,  tmp_prot2, &
!cmpi                                         tmp_neut2, tmp_amjz,  tmp_amjz2, &
!cmpi                                         tmp_amj2,  tmp_amsop, tmp_amson, &
!cmpi                                         tmp_istz,  tmp_istz2, tmp_ist2, &
!cmpi                                         tmp_ra2p,  tmp_ra2n, tmp_ra2m
!cmpi complex(r64), dimension(-1:1,ndimP,ndimJ) :: tmp_Q1mp, tmp_Q1mn, tmp_M1ma
!cmpi complex(r64), dimension(-2:2,ndimP,ndimJ) :: tmp_Q2mp, tmp_Q2mn, tmp_M2ma
!cmpi complex(r64), dimension(-3:3,ndimP,ndimJ) :: tmp_Q3mp, tmp_Q3mn
!cmpi complex(r64), dimension(sdim,2,ndimP,ndimJ) :: tmp_occn
!cmpi 
!cmpi ndim1 = ndimP * ndimJ
!cmpi ndim3 = 3 * ndim1             
!cmpi ndim5 = 5 * ndim1             
!cmpi ndim7 = 7 * ndim1             
!cmpi sdim2 = 2 * sdim * ndim1      
!cmpi 
!cmpi if ( paral_myrank == 0 ) then
!cmpi   if ( misc_phys == 0 ) then
!cmpi     tmp_over  = zzero
!cmpi     tmp_pari  = zzero
!cmpi     tmp_ener  = zzero
!cmpi     tmp_prot  = zzero
!cmpi     tmp_neut  = zzero
!cmpi     tmp_prot2 = zzero
!cmpi     tmp_neut2 = zzero
!cmpi     tmp_amjz  = zzero
!cmpi     tmp_amjz2 = zzero
!cmpi     tmp_amj2  = zzero
!cmpi     tmp_amsop = zzero
!cmpi     tmp_amson = zzero
!cmpi     tmp_istz  = zzero
!cmpi     tmp_istz2 = zzero
!cmpi     tmp_ist2  = zzero
!cmpi     tmp_ra2p  = zzero
!cmpi     tmp_ra2n  = zzero
!cmpi     tmp_ra2m  = zzero
!cmpi     tmp_Q1mp  = zzero
!cmpi     tmp_Q1mn  = zzero
!cmpi     tmp_Q2mp  = zzero
!cmpi     tmp_Q2mn  = zzero
!cmpi     tmp_Q3mp  = zzero
!cmpi     tmp_Q3mn  = zzero
!cmpi     tmp_M1ma  = zzero
!cmpi     tmp_M2ma  = zzero
!cmpi   endif
!cmpi endif
!cmpi 
!cmpi if ( misc_phys == 0 ) then
!cmpi   call mpi_reduce(proj_over,tmp_over,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_pari,tmp_pari,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_ener,tmp_ener,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_prot,tmp_prot,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_neut,tmp_neut,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_prot2,tmp_prot2,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_neut2,tmp_neut2,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_amjz,tmp_amjz,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_amjz2,tmp_amjz2,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_amj2,tmp_amj2,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_amsop,tmp_amsop,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_amson,tmp_amson,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_istz,tmp_istz,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_istz2,tmp_istz2,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_ist2,tmp_ist2,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_ra2p,tmp_ra2p,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_ra2n,tmp_ra2n,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_ra2m,tmp_ra2m,ndim1,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_Q1mp,tmp_Q1mp,ndim3,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_Q1mn,tmp_Q1mn,ndim3,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_Q2mp,tmp_Q2mp,ndim5,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_Q2mn,tmp_Q2mn,ndim5,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_Q3mp,tmp_Q3mp,ndim7,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_Q3mn,tmp_Q3mn,ndim7,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_M1ma,tmp_M1ma,ndim3,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_M2ma,tmp_M2ma,ndim5,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi   call mpi_reduce(proj_occn,tmp_occn,sdim2,mpi_double_complex,mpi_sum, &
!cmpi                   0,mpi_comm_peers,ierr)
!cmpi endif
!cmpi 
!cmpi if ( paral_myrank == 0 ) then
!cmpi   if ( misc_phys == 0 ) then
!cmpi     proj_over  = tmp_over
!cmpi     proj_pari  = tmp_pari
!cmpi     proj_ener  = tmp_ener
!cmpi     proj_prot  = tmp_prot
!cmpi     proj_neut  = tmp_neut
!cmpi     proj_prot2 = tmp_prot2
!cmpi     proj_neut2 = tmp_neut2
!cmpi     proj_amjz  = tmp_amjz
!cmpi     proj_amjz2 = tmp_amjz2
!cmpi     proj_amj2  = tmp_amj2
!cmpi     proj_amsop = tmp_amsop
!cmpi     proj_amson = tmp_amson
!cmpi     proj_istz  = tmp_istz
!cmpi     proj_istz2 = tmp_istz2
!cmpi     proj_ist2  = tmp_ist2
!cmpi     proj_ra2p  = tmp_ra2p
!cmpi     proj_ra2n  = tmp_ra2n
!cmpi     proj_ra2m  = tmp_ra2m
!cmpi     proj_Q1mp  = tmp_Q1mp
!cmpi     proj_Q1mn  = tmp_Q1mn
!cmpi     proj_Q2mp  = tmp_Q2mp
!cmpi     proj_Q2mn  = tmp_Q2mn
!cmpi     proj_Q3mp  = tmp_Q3mp
!cmpi     proj_Q3mn  = tmp_Q3mn
!cmpi     proj_M1ma  = tmp_M1ma
!cmpi     proj_M2ma  = tmp_M2ma
!cmpi     proj_occn  = tmp_occn
!cmpi   endif
!cmpi endif
!cmpi 
!cmpi end subroutine reduce_projmatelem

END MODULE Projection 
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
