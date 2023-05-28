!==============================================================================!
! PROGRAM TAURUS_pav                                                           !
!                                                                              !
! This program performs the symmetry projection (N,Z,A,J,MJ,P) of real general !
! quasi-particle states represented in a spherical harmonic ocillator basis.   !
!                                                                              !
! Licence: GNU General Public License version 3 or later                       !
! DOI: 10.5281/zenodo.XXXXXXX                                                  !
! Repository: github.com/project-taurus/taurus_pav                             !
!                                                                              !
! Article describing the code:                                                 !
!==============================================================================!
PROGRAM TAURUS_pav

!cmpi use MPI            
!cmpi use Parallelization
use Nucleus  
use Projection    
use Initialization

implicit none

integer(i64) :: iang
!cmpi integer :: ierr=0

!!!
!!! INITIALIZATION
!!!

!cmpi call mpi_init(ierr)
!cmpi call mpi_comm_size(mpi_comm_world,paral_worldsize,ierr)
!cmpi call mpi_comm_rank(mpi_comm_world,paral_myrank,ierr)

call print_version        
call read_input
!cmpi call set_parallel_teams

if ( misc_frot < 2 ) then
  !!! Sets the basis and basic operators
  call set_nucleus
  call set_basis   
  call set_angularmomentum    
  call set_isospin         
  if ( misc_phys == 0 ) then 
    if ( misc_part == 0 ) then 
      call set_hamiltonian
    endif
    call set_radius(misc_oper)
    call set_multipoles(misc_oper)
  endif

  !!! Sets the wave functions and checks the symmetries
  call set_wavefunctions  
  call set_fields(misc_phys)
  call read_wavefunctions(HOsp_dim)
  !cmpi if ( paral_myrank == 0 ) then        
  call print_wavefunctions
  !cmpi endif
  call check_symmetries(HOsp_dim) 
else 
  !cmpi if ( paral_myrank == 0 ) then        
  call print_rotmatelem 
  !cmpi endif
endif
call check_compatibilities

!!! Sets the integrals and projected matrix elements
if ( misc_frot == 1 ) then
  call open_file_rotmatelem
endif
call set_projection_angles
call set_projection_matelem

!!!
!!! PROJECTING ON N,Z,A,J,MJ,P
!!!
do iang = loop_ini, loop_fin

  !!! Reads the rotated matrix elements from file and integrates them
  if ( misc_frot == 2 ) then
    call readwrite_rotmatelem(iang)
    if ( (.not. is_papP) .and. (parP == 2) ) cycle
    if ( abs(rot_over) < misc_cutrot ) cycle 
    call integrate_rotmatelem
    cycle
  endif

  !!! Performs the rotation and computes non-diagonal densities
  call reset_rotmatelem
  call determine_angle(iang)
  call generate_rotation_spatial(angle_alpJ(alpJ),angle_betJ(betJ), &
                                 angle_gamJ(gamJ),Rid,ROT,HOsp_dim)
  call generate_rotation_parity(parP,zone*ROT,ROT,HOsp_dim)
  call generate_rotation_gauge(angle_phiZ(phiZ),angle_phiN(phiN), &
                               angle_phiA(phiA),zone*ROT,ROT,HOsp_dim)
  call calculate_overlap(bogo_occL,bogo_empL,bogo_ovvaL,bogo_vouL,bogo_zDL, &
                         bogo_occR,bogo_empR,bogo_ovvaR,bogo_vouR,bogo_zDR, &
                         ROT,rot_over,HOsp_dim)
  if ( abs(rot_over) < misc_cutrot ) cycle 

  ! Cycles if one computes only the overlap
  if ( misc_part == 1 ) then
    call integrate_rotmatelem
    cycle
  endif
 
  call rotate_wavefunction(ROT,bogo_zUR,bogo_zVR,Ubar,Vbar,HOsp_dim)
  call calculate_thouless(bogo_zUL,bogo_zVL,Ubar,Vbar,HOsp_dim)     
  call calculate_densities(bogo_zUL,bogo_zVL,Util,Vtil,dens_rhoLR, &
                           dens_kappaLR,dens_kappaRL,HOsp_dim)
  if ( misc_phys == 0 ) then
    if ( misc_part == 0 ) then 
      call calculate_fields(dens_rhoLR,dens_kappaLR,field_gammaLR,field_hspLR, &
                            field_deltaLR,HOsp_dim)
      call calculate_expectval_energy(dens_rhoLR,dens_kappaRL,field_gammaLR, &
                                      field_deltaLR,rot_ener,HOsp_dim)
    endif
    call calculate_radius(dens_rhoLR,dens_kappaLR,dens_kappaRL,rot_ra2p, & 
                          rot_ra2n,rot_ra2m,HOsp_dim)
    call calculate_particlenumber(0,dens_rhoLR,dens_kappaLR,dens_kappaRL, &
                                 rot_prot,rot_neut,rot_prot2,rot_neut2,HOsp_dim)
    call calculate_angularmomentum(0,dens_rhoLR,dens_kappaLR,dens_kappaRL, &
                                   rot_amjx,rot_amjy,rot_amjz,rot_amjx2, &
                                   rot_amjy2,rot_amjz2,rot_amsop,rot_amson, &
                                   HOsp_dim)
    call calculate_isospin(0,dens_rhoLR,dens_kappaLR,dens_kappaRL,rot_istx, &
                           rot_isty,rot_istz,rot_istx2,rot_isty2,rot_istz2, &
                           HOsp_dim)
    call calculate_multipoles(dens_rhoLR,rot_Q1mp,rot_Q1mn,rot_Q2mp,rot_Q2mn, &
                              rot_Q3mp,rot_Q3mn,rot_M1ma,rot_M2ma,HOsp_dim)
    call calculate_occupation_number(dens_rhoLR,rot_occn,HOsp_dim,HOsh_dim)
    call generate_rotation_parity(2,zone*ROT,ROT,HOsp_dim)
    call calculate_overlap(bogo_occL,bogo_empL,bogo_ovvaL,bogo_vouL,bogo_zDL, &
                           bogo_occR,bogo_empR,bogo_ovvaR,bogo_vouR,bogo_zDR, &
                           ROT,rot_pari,HOsp_dim)
  endif

  !!! Integrate the rotated matrix elements to obtain the projected ones and
  !!! writes them to file depending on option
  call integrate_rotmatelem
  if ( misc_frot == 1 ) then
    call readwrite_rotmatelem(iang)
  endif 

enddo

!!! Reduces the projected matrix elemepts (MPI)
!cmpi if ( (paral_teams > 1) .and. (paral_myteamrank == 0) ) then        
!cmpi   call reduce_projmatelem(tabid_P_dim,tabid_J_dim,HOsh_dim)
!cmpi endif

!!!
!!! WRITING/PRINTING THE PROJECTED MATRIX ELEMENTS
!!!

!cmpi if ( paral_myrank == 0 ) then       
if ( misc_phys == 0 ) then 
  call write_projmatelem_states 
  call write_projmatelem_occnumb
  call write_projmatelem_electromagnetic
  call print_projmatelem_states 
endif

call print_files_complementary

print '(/,"This is the end, my only friend, the end.")'  
!cmpi endif

!!! In principle, this barrier is useless. In practice, however, it proved
!!! necessary with certain mpi implementations
!cmpi call mpi_barrier(mpi_comm_world,ierr)

!cmpi call mpi_finalize(ierr)

END PROGRAM TAURUS_pav   
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
