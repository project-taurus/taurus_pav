!==============================================================================!
! MODULE Initialization                                                        !
!                                                                              !
! This module contains the variables and routines related to the reading of    !
! the input parameters and files.                                              !
!                                                                              !
! List of routines and functions:                                              !
! subroutine print_version                                                     !
! subroutine read_input                                                        !
! subroutine print_input                                                       !
! subroutine broadcast_inputs (MPI)                                            !
! subroutine check_input                                                       !
! subroutine open_files_hamiltonian                                            !
!==============================================================================!
MODULE Initialization

!cmpi use MPI            
!cmpi use Parallelization
use Projection

implicit none

integer :: init_parP_dim
character(30), dimension(40) :: input_names ! Name of inputs

CONTAINS 

!------------------------------------------------------------------------------!
! subroutine print_version                                                     !
!                                                                              !
! Prints the panel with the version, the description and the logo of the code. !
! This is the first action of the code.                                        !
!------------------------------------------------------------------------------!
subroutine print_version

!cmpi if ( paral_myrank == 0 ) then        
print '(" __________________________________________________________ ",/, &
      & "|                                                          |",/, &
      & "|  (______)                                                |",/, &
      & "|  <(0  0)>   TAURUS_pav, version 2023.05.28               |",/, &
      & "|    (°°)                                                  |",/, &
      & "|                                                          |",/, &
      & "| This code performs the symmetry projections (N,Z,J,MJ,P) |",/, &
      & "| of real general Bogoliubov quasi-particle states         |",/, &
      & "| states represented in a spherical harmonic oscillator    |",/, &
      & "| basis.                                                   |",/, &
      & "|                                                          |",/, &
      & "| Licence: GNU General Public License version 3 or later   |",/, &
      & "| DOI: https://doi.org/10.5281/zenodo.XXXXXXX              |",/, &
      & "| Git: https://github.com/project-taurus/taurus_pav.git    |",/, &
      & "|                                                          |",/, &
      & "| Contributors: B. Bally, T. Rodríguez                     |",/, &
      & "|__________________________________________________________|",/)' 
!cmpi endif

end subroutine print_version

!------------------------------------------------------------------------------!
! subroutine read_input                                                        !
!                                                                              !
! Reads the input parameters from the input file and perform somes checks to   !
! see if they have approriate values. Also, opens the files containing the     !
! information on the Hamiltonian and model space.                              !
!------------------------------------------------------------------------------!
subroutine read_input

integer :: i, dummy_teamssize
character(3) :: app2b='.2b'
character(4) :: appsho='.sho', app01b='.01b', appred='.red', appcom='.com'
character(100) :: hamil_dummy
character(len=*), parameter :: format1 = "(1a)", &
                               format2 = "(1a30, 1a100)", &
                               format3 = "(1a30, 1i1)", &
                               format4 = "(1a30, 1i5)", &
                               format5 = "(1a30, 1es10.3)"

!!! Initialization
is_pnpZ  = .false.
is_pnpN  = .false.
is_pnpA  = .false.
is_ampJ  = .false.
is_ampMJ = .false.
is_ampKJ = .false.
is_papP  = .false.

is_phiZ = .false.
is_phiN = .false.
is_phiA = .false.
is_alpJ = .false.
is_betJ = .false.
is_gamJ = .false.
is_parP = .false.

is_good_1   = .false.
is_good_Z   = .false.
is_good_N   = .false.
is_good_A   = .false.
is_good_ZN  = .false.
is_good_J   = .false.
is_good_MKJ = .false.
is_good_T   = .false.
is_good_MKT = .false.

good_Z   = 0
good_N   = 0
good_A   = 0
good_J   = 0
good_MKJ = 0
good_P   = 0
good_T   = 0
good_MKT = 0

!!! Reads the input parameters
!cmpi if ( paral_myrank == 0 ) then        
read(uti,format1) input_names(1)
read(uti,format1) input_names(2)
read(uti,format2) input_names(3),  hamil_dummy
read(uti,format3) input_names(4),  hamil_com  
read(uti,format3) input_names(5),  hamil_read
read(uti,format4) input_names(6),  dummy_teamssize
read(uti,format1) input_names(7) 
read(uti,format1) input_names(8) 
read(uti,format1) input_names(9) 
read(uti,format3) input_names(10), misc_phys
read(uti,format3) input_names(11), misc_part
read(uti,format3) input_names(12), misc_oper
read(uti,format3) input_names(13), misc_frot
read(uti,format5) input_names(14), misc_cutrot
read(uti,format3) input_names(15), seed_text   
read(uti,format5) input_names(16), seed_occeps
read(uti,format3) input_names(17), seed_allemp
read(uti,format1) input_names(18)
read(uti,format1) input_names(19)
read(uti,format1) input_names(20)
read(uti,format4) input_names(21), valence_Z
read(uti,format4) input_names(22), valence_N
read(uti,format4) input_names(23), phiZ_dim
read(uti,format4) input_names(24), phiN_dim
read(uti,format4) input_names(25), phiA_dim
read(uti,format3) input_names(26), pnp_nosimp
read(uti,format1) input_names(27)
read(uti,format1) input_names(28)
read(uti,format1) input_names(29)
read(uti,format4) input_names(30), amp_2jmin
read(uti,format4) input_names(31), amp_2jmax
read(uti,format4) input_names(32), alpJ_dim
read(uti,format4) input_names(33), betJ_dim
read(uti,format4) input_names(34), gamJ_dim
read(uti,format3) input_names(35), amp_nosimp
read(uti,format1) input_names(36)
read(uti,format1) input_names(37)
read(uti,format1) input_names(38)
read(uti,format3) input_names(39), parP_dim   
read(uti,format3) input_names(40), pap_nosimp
!cmpi endif

!!! MPI adjustment and broadcast
!cmpi paral_teamssize = max(dummy_teamssize,1) 
!cmpi call broadcast_inputs(hamil_dummy)

!!! Reads file with rotated matrix elements to get the proper dimensions
if ( misc_frot == 2 ) then
  init_parP_dim = parP_dim
  call open_file_rotmatelem
  valence_A = valence_Z + valence_N
endif

!!! Assigns the names for the hamiltonian files and opens them
i = len_trim(adjustl(hamil_dummy))

allocate (character(i  ) :: hamil_file)
allocate (character(i+4) :: hamil_fsho)
allocate (character(i+4) :: hamil_f01b)
allocate (character(i+3) :: hamil_f2b )
allocate (character(i+4) :: hamil_fred)
allocate (character(i+4) :: hamil_fcom)

hamil_file = trim(adjustl(hamil_dummy))
hamil_fsho = hamil_file // appsho
hamil_f01b = hamil_file // app01b
hamil_f2b  = hamil_file // app2b
hamil_fred = hamil_file // appred
hamil_fcom = hamil_file // appcom

!!! Defines the real indices for the loop/arrays. This is needed in case of no
!!! projections. The values may be changed by the routine check_symmetries.

! Angular momentum
if ( alpJ_dim + betJ_dim + gamJ_dim == 0 ) then
  amp_2jmin = 0 
  amp_2jmax = 0
endif

if ( alpJ_dim > 0 ) then 
  is_alpJ  = .true.
  is_ampMJ = .true.
endif

if ( betJ_dim > 0 ) then
  is_betJ = .true.
  is_ampJ = .true.
endif

if ( gamJ_dim > 0 ) then 
  is_gamJ  = .true.
  is_ampKJ = .true.
endif

alpJ_dim = max(1,alpJ_dim)
betJ_dim = max(1,betJ_dim)
gamJ_dim = max(1,gamJ_dim)

! Particle number
if ( phiZ_dim > 0 ) then
  is_phiZ = .true.
  is_pnpZ = .true.
endif

if ( phiN_dim > 0 ) then
  is_phiN = .true.
  is_pnpN = .true.
endif

if ( phiA_dim > 0 ) then
  is_phiA = .true.
  is_pnpA = .true.
endif

phiZ_dim = max(1,phiZ_dim)
phiN_dim = max(1,phiN_dim)
phiA_dim = max(1,phiA_dim)

! Parity
if ( parP_dim > 0 ) then
  is_papP = .true.
  is_parP = .true. 
  pap_pmin = -1
  pap_pmax =  1
else
  pap_pmin = 1
  pap_pmax = 1
endif
parP_dim = min(2,parP_dim + 1)

if ( (misc_frot == 2) .and. (init_parP_dim == 0) ) then
  is_papP = .false.
  pap_pmin = 1
  pap_pmax = 1
endif

! Neutrinoless double beta decay: from 0+ to 0+
if ( misc_phys == 2 ) then
  amp_2jmin = 0
  amp_2jmax = 0
  pap_pmin = 1
  pap_pmax = 1
endif

!!! Cutoffs 
if ( seed_occeps < tiny(zero) ) seed_occeps = 1.d-8
if ( misc_cutrot < tiny(zero) ) misc_cutrot = 1.d-16

!!! Performs some tests on the value of the inputs and link the units for
!!! the Hamiltonian files
!cmpi if ( paral_myrank == 0 ) then        
call print_input
!cmpi endif
call check_input
call open_file_sho
if ( misc_phys == 0 ) then
  call open_files_hamiltonian
endif

end subroutine read_input

!------------------------------------------------------------------------------!
! subroutine print_input                                                       !
!                                                                              !
! Prints the input parameters at the beginning of the calculation in the same  !
! format such that it can be copied and reused in an input file.               !
!------------------------------------------------------------------------------!
subroutine print_input

integer :: dummy_teamssize=0, tmp_alpJ_dim, tmp_betJ_dim, tmp_gamJ_dim, &
           tmp_phiZ_dim, tmp_phiN_dim, tmp_phiA_dim, tmp_parP_dim
character(5) :: paral_teamssize_ch, valence_Z_ch, valence_N_ch, phiZ_dim_ch, &
                phiN_dim_ch, phiA_dim_ch, amp_2jmin_ch, amp_2jmax_ch, &
                alpJ_dim_ch, betJ_dim_ch, gamJ_dim_ch 
character(10) :: misc_cutrot_ch, seed_occeps_ch
character(len=*), parameter :: format1 = "(1a)", &
                               format2 = "(1a30, 1a)", &
                               format3 = "(1a30, 1i1)", &
                               format4 = "(1a30, 1a5)", &
                               format5 = "(1a30, 1a10)"

!!! Sets the (printed) dimensions to 0 if no projection 
tmp_phiZ_dim = 0 
tmp_phiN_dim = 0 
tmp_phiA_dim = 0 
tmp_alpJ_dim = 0 
tmp_betJ_dim = 0 
tmp_gamJ_dim = 0 
tmp_parP_dim = 0 

if ( is_phiZ ) tmp_phiZ_dim = phiZ_dim
if ( is_phiN ) tmp_phiN_dim = phiN_dim
if ( is_phiA ) tmp_phiA_dim = phiA_dim
if ( is_alpJ ) tmp_alpJ_dim = alpJ_dim
if ( is_betJ ) tmp_betJ_dim = betJ_dim
if ( is_gamJ ) tmp_gamJ_dim = gamJ_dim
if ( is_parP ) tmp_parP_dim = parP_dim-1

!!! Formats the variable to eliminate the unpleasant blanck spaces
!cmpi dummy_teamssize = paral_teamssize 
write(paral_teamssize_ch,'(1i5)') dummy_teamssize
paral_teamssize_ch = adjustl(paral_teamssize_ch)

write(misc_cutrot_ch,'(1es10.3)') misc_cutrot
write(seed_occeps_ch,'(1es10.3)') seed_occeps
misc_cutrot_ch = adjustl(misc_cutrot_ch)
seed_occeps_ch = adjustl(seed_occeps_ch)

write(valence_Z_ch,'(1i5)') valence_Z
write(valence_N_ch,'(1i5)') valence_N
write(phiZ_dim_ch,'(1i5)') tmp_phiZ_dim
write(phiN_dim_ch,'(1i5)') tmp_phiN_dim
write(phiA_dim_ch,'(1i5)') tmp_phiA_dim
valence_Z_ch = adjustl(valence_Z_ch)
valence_N_ch = adjustl(valence_N_ch)
phiZ_dim_ch = adjustl(phiZ_dim_ch)
phiN_dim_ch = adjustl(phiN_dim_ch)
phiA_dim_ch = adjustl(phiA_dim_ch)

write(amp_2jmin_ch, '(1i5)') amp_2jmin
write(amp_2jmax_ch, '(1i5)') amp_2jmax
write(alpJ_dim_ch,'(1i5)') tmp_alpJ_dim
write(betJ_dim_ch, '(1i5)') tmp_betJ_dim 
write(gamJ_dim_ch,'(1i5)') tmp_gamJ_dim
amp_2jmin_ch  = adjustl(amp_2jmin_ch )
amp_2jmax_ch  = adjustl(amp_2jmax_ch )
alpJ_dim_ch = adjustl(alpJ_dim_ch)
betJ_dim_ch  = adjustl(betJ_dim_ch )
gamJ_dim_ch = adjustl(gamJ_dim_ch)

!!! Prints the input parameters
print '(60("%"),/,22x,"INPUT PARAMETERS",22x,/,60("%"),/)'
write(uto,format1) input_names(1)
write(uto,format1) input_names(2)
write(uto,format2) input_names(3),  trim(adjustl(hamil_file))
write(uto,format3) input_names(4),  hamil_com  
write(uto,format3) input_names(5),  hamil_read 
write(uto,format4) input_names(6),  paral_teamssize_ch
write(uto,format1) input_names(7) 
write(uto,format1) input_names(8) 
write(uto,format1) input_names(9) 
write(uto,format3) input_names(10), misc_phys
write(uto,format3) input_names(11), misc_part
write(uto,format3) input_names(12), misc_oper
write(uto,format3) input_names(13), misc_frot
write(uto,format5) input_names(14), misc_cutrot_ch
write(uto,format3) input_names(15), seed_text   
write(uto,format5) input_names(16), seed_occeps_ch
write(uto,format3) input_names(17), seed_allemp
write(uto,format1) input_names(18)
write(uto,format1) input_names(19)
write(uto,format1) input_names(20)
write(uto,format4) input_names(21), valence_Z_ch
write(uto,format4) input_names(22), valence_N_ch
write(uto,format4) input_names(23), phiZ_dim_ch
write(uto,format4) input_names(24), phiN_dim_ch
write(uto,format4) input_names(25), phiA_dim_ch
write(uto,format3) input_names(26), pnp_nosimp
write(uto,format1) input_names(27)
write(uto,format1) input_names(28)
write(uto,format1) input_names(29)
write(uto,format4) input_names(30), amp_2jmin_ch
write(uto,format4) input_names(31), amp_2jmax_ch
write(uto,format4) input_names(32), alpJ_dim_ch
write(uto,format4) input_names(33), betJ_dim_ch
write(uto,format4) input_names(34), gamJ_dim_ch
write(uto,format3) input_names(35), amp_nosimp
write(uto,format1) input_names(36)
write(uto,format1) input_names(37)
write(uto,format1) input_names(38)
write(uto,format3) input_names(39), tmp_parP_dim   
write(uto,format3) input_names(40), pap_nosimp
print*,' '

end subroutine print_input

!------------------------------------------------------------------------------!
! subroutine braodcast_inputs                                                  !
!                                                                              !
! Broadcasts the inputs read by process rank = 0 to all others. This is safer  !
! I think than making all the processes read the input file.                   !
!------------------------------------------------------------------------------!
!cmpi subroutine broadcast_inputs(hamil_dummy)

!cmpi character(100), intent(in) :: hamil_dummy
!cmpi integer :: ierr=0

!!! Hamiltonian
!cmpi call mpi_bcast(hamil_dummy,100,mpi_character,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(hamil_com,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(hamil_read,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(paral_teamssize,1,mpi_integer,0,mpi_comm_world,ierr)

!!! Miscellaneous   
!cmpi call mpi_bcast(misc_phys,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(misc_part,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(misc_frot,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(misc_cutrot,1,mpi_double_precision,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(seed_text,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(seed_occeps,1,mpi_double_precision,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(seed_allemp,1,mpi_integer,0,mpi_comm_world,ierr)

!!! Particle number
!cmpi call mpi_bcast(valence_Z,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(valence_N,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(phiZ_dim,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(phiN_dim,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(phiA_dim,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(pnp_nosimp,1,mpi_integer,0,mpi_comm_world,ierr)

!!! Angular momentum
!cmpi call mpi_bcast(amp_2jmin,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(amp_2jmax,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(alpJ_dim,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(betJ_dim,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(gamJ_dim,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(amp_nosimp,1,mpi_integer,0,mpi_comm_world,ierr)

!!! Parity            
!cmpi call mpi_bcast(parP_dim,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(pap_nosimp,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(init_parP_dim,1,mpi_integer,0,mpi_comm_world,ierr)

!cmpi end subroutine broadcast_inputs

!------------------------------------------------------------------------------!
! subroutine check_input                                                       !
!                                                                              !
! Checks the values of the input parameters to see if they are appropriate.    !
! If not, the code will stop.                                                  !
!------------------------------------------------------------------------------!
subroutine check_input

integer :: ierror
!cmpi integer :: ierr=0

!!! Counter for the number of errors (should be 0 at the end)
ierror = 0

!!!
!!! Hamiltonian    
!!!

!cmpi if ( paral_myrank == 0 ) then
if ( hamil_com > 1 ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option for the COM correction (hamil_com) = ", & 
         hamil_com," should be 0 or 1."
endif 

if ( hamil_read > 1 ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to read an uncoupled file (hamil_read) = ", & 
         hamil_read," should be 0 or 1."
endif 

!cmpi if ( hamil_read /= 1 ) then
!cmpi   ierror = ierror + 1
!cmpi   print "(a,1i1,a)","The option to read an uncoupled file (hamil_read) &  
!cmpi         &= ",hamil_read," should be 1 when doing MPI calculations."
!cmpi endif

!cmpi if ( (paral_teamssize < 0) .or. (paral_teamssize > paral_worldsize) ) then
!cmpi   ierror = ierror + 1
!cmpi   print "(a,1i1,a)","The numer of processes per team (MPI) = ", & 
!cmpi          paral_teamssize," should be positive and smaller than the &
!cmpi         &total number of processes."
!cmpi endif 

!!!
!!! Miscellaneous
!!!

if ( misc_phys > 0 ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The physical case studied (misc_phys) = ", & 
         misc_phys," should be 0 in this version of the code."
endif 

if ( misc_part > 2 ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The part of the computation (misc_part) = ", & 
         misc_part," should be 0, 1 or 2."
endif 

if ( misc_oper > 1 ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to read the matrix elements of the operators &
         &(misc_oper) = ", misc_oper," should be 0 or 1."
endif 

if ( misc_frot > 2 ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to write/read the rotated matrix elements &
        &(misc_frot) = ",misc_frot," should be 0, 1 or 2."
endif 

if ( misc_cutrot < 0.0d0 ) then
  ierror = ierror + 1
  print "(a,1es10.3,a)","The cutoff for the rotated norm overlaps &
        &(misc_cutrot) = ", misc_cutrot," should be positive."
endif 

if ( seed_text > 1 ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to read binary/text wave functions &
        &(seed_text) = ",seed_text," should be 0 or 1."
endif 

if ( seed_occeps < 0.0d0 ) then
  ierror = ierror + 1
  print "(a,1es10.3,a)","The cutoff for the occupied single-particle states &
        &(seed_occeps) = ", seed_occeps," should be positive."
endif 

if ( seed_allemp > 1 ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option for including all the empty states &
        &(seed_allemp) = ", seed_allemp," should be 0 or 1."
endif 

!!!
!!! Particle number
!!!

if ( valence_Z < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The number of active protons (valence_Z) = ", & 
         valence_Z," should be positive."
endif 

if ( valence_N < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The number of active neutrons (valence_N) = ", & 
         valence_N," should be positive."
endif 

if ( (abs(valence_Z) + abs(valence_N)) == 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The number of active nucleons (valence_Z + valence_N) &
        &= ",valence_Z+valence_N," should be strictly positive."
endif 

if ( phiZ_dim < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The number of gauge angles for protons (phiZ_dim) = ", & 
         phiZ_dim," should be positive."
endif 

if ( phiN_dim < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The number of gauge angles for neutrons (phiN_dim) = ", & 
         phiN_dim," should be positive."
endif 

if ( phiA_dim < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The number of gauge angles for nucleons (phiA_dim) = ", & 
         phiA_dim," should be positive."
endif 

if ( pnp_nosimp > 1 ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to disable the simplifications for Z/N/A &
        & (pnp_nosimp) = ",pnp_nosimp," should be 0 or 1."
endif 

!!!
!!! Angular momentum
!!!

if ( amp_2jmin < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The minimum angular momentum 2J (amp_2jmin) = ", & 
         amp_2jmin," should be positive."
endif 

if ( amp_2jmax < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The maximum angular momentum 2J (amp_2jmax) = ", & 
         amp_2jmax," should be positive."
endif 

if ( amp_2jmax < amp_2jmin ) then
  ierror = ierror + 1
  print "(a,1i5,a,1i5,a)","The maximum angular momentum 2J (amp_2jmax) = ", & 
         amp_2jmax," should be greater than or equal to the minimum angular &
         &momentum 2J (amp_2Jmin) =",amp_2jmin,"."
endif 

if ( (-1)**amp_2jmin /= (-1)**amp_2jmax ) then
  ierror = ierror + 1
  print "(a,1i5,a,1i5,a)","The minimum angular momentum 2J (amp_2jmin) = ", & 
         amp_2jmin," should have the same number parity as the maximum angular &
         &momentum 2J (amp_2jmax) = ",amp_2jmax,"."
endif 

if ( (is_ampMJ .or. is_ampJ .or. is_ampKJ) .and. &
     ((-1)**amp_2jmin /= (-1)**(valence_N + valence_Z)) ) then
  ierror = ierror + 1
  print "(a,1i5,a,1i5,a)","The minimum angular momentum 2J (amp_2jmin) = ", & 
         amp_2jmin," should have the same number parity as the nucleon number &
         &(valence_Z + valence_N) = ",valence_Z + valence_N,"."
endif 

if ( alpJ_dim < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The number of Euler angles for alpha (alpJ_dim) = ", & 
         alpJ_dim," should be positive."
endif 

if ( betJ_dim < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The number of Euler angles for beta (betJ_dim) = ", & 
         betJ_dim," should be positive."
endif 

if ( gamJ_dim < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The number of Euler angles for gamma (gamJ_dim) = ", & 
         gamJ_dim," should be positive."
endif 

if ( amp_nosimp > 1 ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to disable the simplifications for J/MJ/KJ &
        &(amp_nosimp) = ",amp_nosimp," should be 0 or 1."
endif 

!!!
!!! Parity
!!!

if ( parP_dim-1 > 1 ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to perform parity projection (parP_dim) = ", & 
         parP_dim," should be 0 or 1."
endif 

if ( pap_nosimp > 1 ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to disable the simplifications for P &
        & (pap_nosimp) = ",pap_nosimp," should be 0 or 1."
endif 
!cmpi endif 

!!!
!!! Stops the code if an error has been found in the input file
!!!

!cmpi call mpi_bcast(ierror,1,mpi_integer,0,mpi_comm_world,ierr)

if ( ierror /= 0 ) then
!cmpi   if ( paral_myrank == 0 ) then
  print "(a,1i2,a)", "The code has dectected ",ierror," problem(s) with the &
        &input parameters and will stop. Please check the manual."
!cmpi   endif 
  stop 
endif

end subroutine check_input

!------------------------------------------------------------------------------!
! subroutine open_file_sho                                                     !
!                                                                              !
! Opens the SHO file to link the unit if it exists. If not the code will stop. !
! BB: contrarily to subroutine check_input, when performing MPI calculations   !
! all the processes will print the error messages. I might want to change that !
! latter but it might be of advantage to keep it in case only part of the      !
! processes cannot access the files.                                           !
!------------------------------------------------------------------------------!
subroutine open_file_sho             

integer :: ierror
logical :: is_exist            

!!! Counter for the number of errors (should be 0 at the end)
ierror = 0

!!!
!!! SHO file 
!!!
inquire (file=hamil_fsho, exist=is_exist)

if ( is_exist ) then
  open(uth, file=hamil_fsho, status='old', action='read', form='formatted')
  read(uth,'(a)') hamil_name
  read(uth,*) hamil_type
  rewind(uth)
else
  ierror = ierror + 1
  print "(a,a,a)", "The main file (hamil_fsho) = ", hamil_fsho, &
        " can not be found."
endif 

!!!
!!! Stops the code if an error has been found with the hamiltonian files
!!!
if ( ierror /= 0 ) then
  print "(a,1i1,a)", "The code has dectected ",ierror," problem(s) with the &
        &hamiltonian files and will stop. Please check the files."
  stop 
endif

end subroutine open_file_sho

!------------------------------------------------------------------------------!
! subroutine open_files_hamiltonian                                            !
!                                                                              !
! Opens the hamiltonian files to link the units if they exist. If not, the     !
! code will stop.                                                              !
! BB: same comment as in open_file_sho.                                        !
!------------------------------------------------------------------------------!
subroutine open_files_hamiltonian

integer :: ierror, iwarn
character(100) :: hname1, hname2, hnamer
logical :: is_exist            

!!! Counter for the number of errors (should be 0 at the end)
ierror = 0
iwarn  = 0

!!! Reduced file
if ( hamil_read == 1 ) then
  inquire (file=hamil_fred, exist=is_exist)
  if ( is_exist ) then
    open(uthr, file=hamil_fred, status='old', action='read', access='stream',&
         form='unformatted')              
    read(uthr) hnamer
    if ( hnamer /= hamil_name ) then
      iwarn = iwarn + 1
      print "(a,a,a,a)", "Warning: the name in the reduced file = ", &
             trim(adjustl(hnamer)), &
            " does not correspond with the one of the main file = ", &
             trim(adjustl(hamil_name))     
    endif
  else  
    ierror = ierror + 1
    print "(a,a,a)", "The uncoupled binary file (hamil_fred) = ",hamil_fred, &
          " can not be found."
  endif
endif

!!!
!!! 0+1 and 2-body files
!!!
select case ( hamil_type*(1-hamil_read) )
  case (0:2)
    continue

  case (3:4)
    !!! 1-body
    if ( hamil_type == 3 ) then
      inquire (file=hamil_f01b, exist=is_exist)
      if ( is_exist ) then
        open(uth1, file=hamil_f01b, status='old', action='read', & 
             form='formatted')
        read(uth1,'(a)') hname1
        rewind(uth1)

        if ( hname1 /= hamil_name ) then
          iwarn = iwarn + 1
          print "(a,a,a,a)", "Warning: the name in the 0+1-body file = ", &
                 trim(adjustl(hname1)), &
                " does not correspond with the one of the main file = ",  &
                 trim(adjustl(hamil_name))     
        endif
      else  
        ierror = ierror + 1
        print "(a,a,a)", "The 0+1-body file (hamil_f01b) = ",hamil_f01b, &
              " can not be found."
      endif
    endif

    !!! 2-body
    inquire (file=hamil_f2b, exist=is_exist)
    if ( is_exist ) then
      open(uth2, file=hamil_f2b, status='old', action='read', form='formatted')
      read(uth2,'(a)') hname2
      rewind(uth2)
      if ( hname2 /= hamil_name ) then
        iwarn = iwarn + 1
        print "(a,a,a,a)", "Warning: the name in the 2-body file = ",  &  
               trim(adjustl(hname2)), &
              " does not correspond with the one of the main file = ", &
               trim(adjustl(hamil_name))     
      endif
    else 
      ierror = ierror + 1
      print "(a,a,a)", "The 2-body file (hamil_f2b) = ",hamil_f2b, &
            " can not be found."
    endif

    !!! Center of mass (2-body)
    if ( hamil_com == 1 ) then
      inquire (file=hamil_fcom, exist=is_exist)
      if ( is_exist ) then
        open(uthc, file=hamil_fcom, status='old', action='read', & 
             form='formatted')
      else  
        ierror = ierror + 1
        print "(a,a,a)", "The center-of-mass file (hamil_fcom) = ",hamil_fcom, &
              " can not be found."
      endif
    endif

  case default
    ierror = ierror + 1
    print "(a,1i1,a)", "The hamiltonian format (hamil_type) = ",hamil_type, & 
          " should be between 1 and 4."
end select

!!!
!!! Stops the code if an error has been found with the hamiltonian files
!!!
if ( ierror /= 0 ) then
  print "(a,1i1,a)", "The code has dectected ",ierror," problem(s) with the &
        &hamiltonian files and will stop. Please check the files."
  stop 
endif

if ( iwarn /= 0 ) print*,' '

end subroutine open_files_hamiltonian

END MODULE Initialization
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
