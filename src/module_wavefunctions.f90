!==============================================================================!
! MODULE Wavefunctions                                                         !
!                                                                              !
! This module contains the variables and routines related to the wave funct-   !
! ions (including their densities).                                            !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_wavefunctions                                               !
! - subroutine read_wavefunctions                                              !
! - subroutine construct_canonical_basis                                       !
! - subroutine calculate_densities                                             !
! - subroutine calculate_overlap                                               !
! - subroutine calculate_norm                                                  !
!==============================================================================!
MODULE Wavefunctions 

use Basis

implicit none

!!! Parameters that determines the seed wave function
integer :: seed_text, & ! format of seed                 
           seed_allemp  ! option to include the empty states (overlap)
real(r64) :: seed_occeps ! cutoff for occupied single-particle states 

!!! Bogoliubov matrices
integer :: bogo_occL, & ! number of fully-occupied states in left  state
           bogo_occR, & !   "    "    "      "       "    "  right   "
           bogo_empL, & !   "    "    "    empty     "    "  left    "
           bogo_empR, & !   "    "    "      "       "    "  right   "
           bogo_npaL, & ! number parity of left  state 
           bogo_npaR    !   "       "   "  right   "
integer(i64) :: bogo_labelL, & ! label for the left  state
                bogo_labelR    !   "    "   "  right   "
real(r64) :: bogo_ovvaL, & ! overlap paired part of left  state with vacuum
             bogo_ovvaR    !    "      "     "   "  right   "    "     "
real(r64), dimension (:), allocatable :: bogo_vouL, & ! v/u for left  state
                                         bogo_vouR    !  "   "  right   "
complex(r64), dimension(:,:), allocatable :: bogo_zUL,  & ! left state U  
                                             bogo_zVL,  & !  "     "   V
                                             bogo_zDL,  & !  "     "   can. bas.
                                             bogo_zUR,  & ! right state U  
                                             bogo_zVR,  & !  "      "   V
                                             bogo_zDR     !  "      "  can. bas.
                
!!! Density matrices
complex(r64), dimension(:,:), allocatable :: dens_rhoLR,   & ! <L|a+a|R> 
                                             dens_kappaLR, & ! <L|aa|R>
                                             dens_kappaRL    ! <L|a+a+|R>^*

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_wavefunctions                                                 !
!                                                                              !
! Initializes the arrays related to the wave functions                         !
!------------------------------------------------------------------------------!
subroutine set_wavefunctions

integer :: ialloc=0

!!! Wavefunctions
allocate( bogo_zUL(HOsp_dim,HOsp_dim),  bogo_zUR(HOsp_dim,HOsp_dim), &
          bogo_zVL(HOsp_dim,HOsp_dim),  bogo_zVR(HOsp_dim,HOsp_dim), &
          bogo_zDL(HOsp_dim,HOsp_dim),  bogo_zDR(HOsp_dim,HOsp_dim), &
          stat=ialloc ) 
if ( ialloc /= 0 ) stop 'Error during allocation of wave functions'

bogo_zUL  = zzero
bogo_zVL  = zzero
bogo_zDL  = zzero
bogo_zUR  = zzero
bogo_zVR  = zzero
bogo_zDR  = zzero

!!! Densities
allocate( dens_rhoLR(HOsp_dim,HOsp_dim), dens_kappaLR(HOsp_dim,HOsp_dim), &
          dens_kappaRL(HOsp_dim,HOsp_dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of densities'

dens_rhoLR   = zzero
dens_kappaLR = zzero
dens_kappaRL = zzero

end subroutine set_wavefunctions

!------------------------------------------------------------------------------!
! subroutine read_wavefunctions                                                !
!                                                                              ! 
! This subroutine allocates the arrays related to the wave functions, reads    !
! the initial wave functions from left.wf and right.wf, and determine their    !
! canonical bases.                                                             !
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      !
!------------------------------------------------------------------------------!
subroutine read_wavefunctions(ndim)

integer, intent(in) :: ndim
integer :: i, j, icheck, HOsh_dimL, HOsh_dimR
integer, dimension(:), allocatable :: HOsh_naL, HOsh_naR   
real(r64), dimension (ndim,ndim) :: UL, VL, UR, VR
complex(r64), dimension (ndim,ndim) :: ULc, VLc, URc, VRc 
character(4) :: filetype  
character(11) :: fileform, filenameL
character(12) :: filenameR 
logical :: is_exist, is_binary

!!! Determines the name of the file (binary or text file)
select case (seed_text)
  case(0)
    is_binary = .true.
  case(1)
    is_binary = .false.
end select

if ( is_binary ) then
  fileform = 'unformatted'
  filetype = '.bin'
else
  fileform = 'formatted'
  filetype = '.txt'
endif

filenameL = 'left_wf' // filetype
filenameR = 'right_wf' // filetype

!!! Checks if the wave functions exist, otherwise stops the run 
inquire (file=filenameL, exist=is_exist)
if ( is_exist .eqv. .false. ) then
  print '(/,"The file containing the left wave function is missing.")'
  stop 
endif 

inquire (file=filenameR, exist=is_exist)
if ( is_exist .eqv. .false. ) then
  print '(/,"The file containing the right wave function is missing.")'
  stop 
endif 

!!! Opens the files
open(utlw, file=filenameL, status='old', action='read', form=fileform)
open(utrw, file=filenameR, status='old', action='read', form=fileform)

!!! Reads the model space 
if ( is_binary ) then
  read(utlw) HOsh_dimL
  allocate(HOsh_naL(HOsh_dimL))
  read(utlw) (HOsh_naL(i), i=1,HOsh_dimL)

  read(utrw) HOsh_dimR
  allocate(HOsh_naR(HOsh_dimR))
  read(utrw) (HOsh_naR(i), i=1,HOsh_dimR)
else
  read(utlw,*) HOsh_dimL
  allocate(HOsh_naL(HOsh_dimL))
  do i = 1, HOsh_dimL
    read(utlw,*) HOsh_naL(i)
  enddo

  read(utrw,*) HOsh_dimR
  allocate(HOsh_naR(HOsh_dimR))
  do i = 1, HOsh_dimR
    read(utrw,*) HOsh_naR(i)
  enddo
endif

!!! Stops the run if the model spaces of the wave func. and interaction differ    
icheck = 0

if ( HOsh_dimL /= HOsh_dim ) icheck = icheck + 1
do i = 1, min(HOsh_dim,HOsh_dimL)
  if ( HOsh_naL(i) /= HOsh_na(i) ) icheck = icheck + 1
enddo

if ( icheck /= 0 ) then
  print '(/,"The model space of the left wave function is not consistent", & 
        & " with the one of the Hamiltonian.")'
  print*, 'Hamil:', HOsh_dim, (HOsh_na(i), i=1,HOsh_dim)
  print*, 'State:', HOsh_dimL, (HOsh_naL(i), i=1,HOsh_dimL)
  stop
endif

if ( HOsh_dimR /= HOsh_dim ) icheck = icheck + 1
do i = 1, min(HOsh_dim,HOsh_dimR)
  if ( HOsh_naR(i) /= HOsh_na(i) ) icheck = icheck + 1
enddo

if ( icheck /= 0 ) then
  print '(/,"The model space of the right wave function is not consistent", & 
        & " with the one of the Hamiltonian.")'
  print*, 'Hamil:', HOsh_dim, (HOsh_na(i), i=1,HOsh_dim)
  print*, 'State:', HOsh_dimR, (HOsh_naR(i), i=1,HOsh_dimR)
  stop
endif

!!! Reads the wave function                                                  
if ( is_binary ) then
  read(utlw) bogo_labelL   
  read(utlw) ((UL(j,i), j=1,HOsp_dim), i=1,HOsp_dim)
  read(utlw) ((VL(j,i), j=1,HOsp_dim), i=1,HOsp_dim)

  read(utrw) bogo_labelR   
  read(utrw) ((UR(j,i), j=1,HOsp_dim), i=1,HOsp_dim)
  read(utrw) ((VR(j,i), j=1,HOsp_dim), i=1,HOsp_dim)
else
  read(utlw,*) bogo_labelL  
  do i = 1, HOsp_dim
    do j = 1, HOsp_dim
      read(utlw,*) UL(j,i)
    enddo
  enddo
  do i = 1, HOsp_dim
    do j = 1, HOsp_dim
      read(utlw,*) VL(j,i)
    enddo
  enddo

  read(utrw,*) bogo_labelR  
  do i = 1, HOsp_dim
    do j = 1, HOsp_dim
      read(utrw,*) UR(j,i)
    enddo
  enddo
  do i = 1, HOsp_dim
    do j = 1, HOsp_dim
      read(utrw,*) VR(j,i)
    enddo
  enddo
endif

bogo_zUL = zone * UL
bogo_zVL = zone * VL

bogo_zUR = zone * UR
bogo_zVR = zone * VR

close (utlw, status='keep')
close (utrw, status='keep')

deallocate(HOsh_naL,HOsh_naR)

!!! Determine the canonical basis of the original state and  
!!! count the number of fully occupied states (used for the overlap
!!! calculation later on)
call construct_canonical_basis(UL,VL,ULc,VLc,bogo_zDL,bogo_ovvaL,bogo_occL, &
                               bogo_empL,HOsp_dim)
call construct_canonical_basis(UR,VR,URc,VRc,bogo_zDR,bogo_ovvaR,bogo_occR, &
                               bogo_empR,HOsp_dim)

allocate(bogo_vouL(HOsp_dim-bogo_occL))
allocate(bogo_vouR(HOsp_dim-bogo_occR))

bogo_vouL = zero
bogo_vouR = zero

j = 1
do i = 1+bogo_empL, HOsp_dim-bogo_occL
  if ( (-1)**j == -1 ) then 
    bogo_vouL(i) =  real(VLc(i,i+1) / ULc(i,i))
  else
    bogo_vouL(i) =  real(VLc(i,i-1) / ULc(i,i))
  endif
  j = j + 1
enddo

j = 1
do i = 1+bogo_empR, HOsp_dim-bogo_occR
  if ( (-1)**j == -1 ) then 
    bogo_vouR(i) =  real(VRc(i,i+1) / URc(i,i))
  else
    bogo_vouR(i) =  real(VRc(i,i-1) / URc(i,i))
  endif
  j = j + 1
enddo

!! Slater determinants
if ( HOsp_dim-bogo_occL-bogo_empL == 0 ) bogo_ovvaL = 1.0d0
if ( HOsp_dim-bogo_occR-bogo_empR == 0 ) bogo_ovvaR = 1.0d0

!!! Computes the number parity of the input wave functions
bogo_npaL = (-1)**bogo_occL
bogo_npaR = (-1)**bogo_occR

end subroutine read_wavefunctions 

!------------------------------------------------------------------------------!
! subroutine print_wavefunctions                                               !
!                                                                              !
! Prints the informations about the left and right wavefunctions.              !
!------------------------------------------------------------------------------!
subroutine print_wavefunctions

character(len=*), parameter :: format1 = "(1a20,1x,1i19,1x,1i19)" , &
                               format2 = "(1a20,8x,1i3,17x,1i3)" 

print '(/,60("%"),/,23x,"WAVE FUNCTIONS",23x,/,60("%"),//, &
      & 3x,"Description",15x,"Left",15x,"Right",/,60("-"))'
print format1, 'Label of state      ', bogo_labelL, bogo_labelR
print format2, 'Total number parity ', bogo_npaL, bogo_npaR
print format2, 'No. of fully occ. sp', bogo_occL, bogo_occR

end subroutine print_wavefunctions 

!------------------------------------------------------------------------------!
! subroutine construct_canonical_basis                                         !
!                                                                              !
! Constructs the canonical basis of the Bogoliubov quasiparticle states with   !
! matrices U,V. Also, computes the number of fully occupied and empty states   !
! in this basis, which is needed to compute the overlap later.                 !
! The algorithm used is the following:                                         !
! 1) We first diagonalize rho                                                  !
! 2) We transforms kappa in the basis that diagonalizes rho                    !
! 3) If they are non-canonical block in kappa (due to degenerate eigenvalues   !
!    of rho), we perform a Schur decomposition on these blocks Given kappa is  !
!    real, it will automatically put it in its canonical form.                 !
!                                                                              !
! Remark: the determination of non-canoncial block can only be numerical and   !
!         therefore may fail in some cases.                                    !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        U,V  Bovoliubov matrices                                              !
!                                                                              !
! Output: zUc,zVc = U and V in the canonical basis                             !
!         zDc = unitary transformation to change the basis                     !
!         ovacc = overlap with the bare vacuum of the fully-paired part of wf  !
!         nocc,nemp = number of fully occupied/empty states in canon. basis    !
!------------------------------------------------------------------------------!
subroutine construct_canonical_basis (U,V,zUc,zVc,zDc,ovacc,nocc,nemp,ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim,ndim), intent(in) :: U, V
integer, intent(out) :: nocc, nemp
real(r64), intent(out) :: ovacc       
complex(r64), dimension(ndim,ndim), intent(out) :: zUc, zVc, zDc
integer :: i, j, k, m, info, sdim, ialloc=0
integer, dimension(ndim) :: submax, subdim
real(r64) :: occu_u, occu_v, occu_v2, eps
real(r64), dimension(ndim) :: eigen_rho
real(r64), dimension(3*ndim-1) :: work1   
real(r64), dimension(:), allocatable :: wr, wi, work2
real(r64), dimension(ndim,ndim) :: rho, kappa, Dc, rhoc, kappac, A1, A2 
real(r64), dimension(:,:), allocatable :: B1, vs 
logical, dimension(:), allocatable :: bwork

!!! Cutoff for occupied single-particle states
eps = seed_occeps

!!! Builds the density rho, then diagonalizes it, then transforms rho in the  
!!! basis where it is diagonal.
call dgemm('n','t',ndim,ndim,ndim,one,V,ndim,V,ndim,zero,rho,ndim)

Dc = rho
call dsyev('V','U',ndim,Dc,ndim,eigen_rho,work1,3*ndim-1,info)

call dgemm('t','n',ndim,ndim,ndim,one,Dc,ndim,rho,ndim,zero,A1,ndim)
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,Dc,ndim,zero,rhoc,ndim)

!!! Counts the number of occupied/empty single-particle states. To determine if 
!!! a state is occupied/empty, we use a cutoff: eps.
!!! Then computes the overlap of the fully paired part with the bare vacuum. As
!!! the quantity occu_u can be quite small, we build \sum log(occu_u) rather 
!!! than \product occu_u.
nocc = 0
nemp = 0
ovacc = 0.0d0

do i = 1, ndim
  occu_v2 = abs(rhoc(i,i))
  occu_u  = sqrt(abs(1.d0 - occu_v2))  
  if ( occu_v2 >= 1.0d0 - eps ) then
    nocc = nocc + 1    
  elseif ( (occu_v2 <= eps) .and. (seed_allemp == 0) ) then
    nemp = nemp + 1
  else
    ovacc = ovacc + log(occu_u) 
  endif
enddo

!!! Builds the density kappa, then builds it in the basis that diagonlizes rho.
!!! At this point, there is no guaranty that kappa will be in its canonical form
!!! yet. There could be degeneracies in the v2,u2 that would prevent a full     
!!! canonical structure.
call dgemm('n','t',ndim,ndim,ndim,one,V,ndim,U,ndim,zero,kappa,ndim) 
call dgemm('t','n',ndim,ndim,ndim,one,Dc,ndim,kappa,ndim,zero,A1,ndim)
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,Dc,ndim,zero,kappac,ndim)

!!! Checks if the kappa is already in its canonical form by counting the number
!!! of non-zero marix in a column of kappac. If it is greater than one (i.e. 
!!! the dimension of the subspace subdim > 2), we need to further reduce kappac.
submax = 0

do i = 1+nemp, ndim-nocc
  submax(i) = i 
  do j = i+1, ndim-nocc
    if ( abs(kappac(j,i)) > sqrt(eps)*1.0d-2 ) submax(i) = j
  enddo
enddo

!!! Tries to determine the dimension of each subspace
k = 1 + nemp

do while ( k < ndim - nocc )
  j = submax(k)
  subdim = 0
  subdim(k:j) = submax(k:j)
  m = maxval(subdim)
  submax(k:j) = m
  if ( m == j ) k = j + 1
enddo

k = 1 + nemp

do while ( k < ndim - nocc )
  subdim(k) = submax(k) - k + 1
  if ( (subdim(k) > 2) .and. (mod(subdim(k),2) /= 0) ) then
    !print*, k, subdim(k)
    subdim(k) = 2
    print*,"Warning: subspace of odd dimension when building kappa canonical. &
           &Intenting the calculation assuming dim = 2."
  endif
  k = k + subdim(k)
enddo

!!! Transforms the degenerates subspaces in their Schur form, which automat-
!!! ically puts them in their canonical form (kappa being real)
A1 = zero
do i = 1, ndim
  A1(i,i) = one
enddo

k = 1 + nemp

do while ( k < ndim - nocc )
  m = subdim(k)
  !!! Schur factorization of this block 
  if ( m > 2 ) then
    !print*,"Diangostic purpose (Schur), dim = :", m
    allocate(B1(m,m), wr(m), wi(m), work2(3*m), bwork(m), vs(m,m), &
             stat=ialloc)
    if ( ialloc /= 0 ) stop 'Error during allocation for Schur factorization'
    
    B1(1:m,1:m) = kappac(k:k+m-1,k:k+m-1)
    call dgees('V','N',lapack_sel,m,B1,m,sdim,wr,wi,vs,m,work2,3*m,bwork,info)
  
    ! store the transformation for this block   
    A1(k:k+m-1,k:k+m-1) = vs(1:m,1:m)  
    kappac(k:k+m-1,k:k+m-1) = B1(1:m,1:m)  
    deallocate(B1,wr,wi,work2,bwork,vs)
  endif
  
  k = k + m
enddo    

!!! Update the transformation D and recomptues kappa in the canonical basis
!!! if necessary
if ( maxval(subdim) > 2 ) then
  call dgemm('n','n',ndim,ndim,ndim,one,Dc,ndim,A1,ndim,zero,A2,ndim)
  Dc = A2
endif

zDc = zone * Dc

!!! Construct U and V in the canonical basis. Note that as we know rho/kappa in
!!! the canonical basis, we do not need to use D or construct C (à la BMZ).
!!! Nevertheless, we have to be careful about the sign of the matrix elements
!!! to reproduce kappac.
zUc = zzero
zVc = zzero

k = 0

do i = 1, ndim
  occu_v2 = abs(rhoc(i,i))
  occu_v = sqrt(occu_v2)

  zUc(i,i) = sqrt(abs(1.d0 - occu_v2))

  if ( occu_v2 > 1.0d0 - eps ) then   ! occupied states
    zVc(i,i) = occu_v * zone          
  elseif ( occu_v2 <= eps ) then      ! empty states
    zVc(i,i) = occu_v * zone          
  else                                ! paired states
    if ( mod(k,2) == 0 ) then
      zVc(i,i+1) = sign(occu_v,kappac(i,i+1)) * zone 
    else
      zVc(i,i-1) = sign(occu_v,kappac(i,i-1)) * zone
    endif
    k = k + 1
  endif
enddo

end subroutine construct_canonical_basis

!------------------------------------------------------------------------------!
! subroutine calculate_densities                                               !
!                                                                              !
! Calculates the (possibly non-diagonal) densities according to                !
!     rhoLR =   V_R^* V_L^T                                                    !
!   kappaLR =   V_R^* U_L^T                                                    !
! kappaRL^* = - U_R^* V_L^T = V_L U_R^dagger                                   !
! Be careful that what is called kappaRL in the code is actually kappaRL^*     !
!                                                                              !
! Actually, the projection routine will inject the complex conjugate of Utilde !
! and Vtilde because we use in that case the formulae (that can be found in    !
! T.R. Rodriguez master thesis)                                                !
!     rhoLR =   Vtilde V_L^T                                                   !
!   kappaLR =   Vtilde U_L^T                                                   !
! kappaRL^* = - Utilde V_L^T = V_L Utilde^T                                    !
! where                                                                        !
! Utilde =  U_L^* + V_L Aphi                                                   !
! Vtilde =  V_L^* + U_L Aphi                                                   !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        UL,VL = Left Bogoliubov matrices                                      !
!        UR,VR = Right Bogoliubov matrices                                     !
!                                                                              !
! Output: rhoLR,kappaLR,kappaRL = densities                                    !
!------------------------------------------------------------------------------!
subroutine calculate_densities(UL,VL,UR,VR,rhoLR,kappaLR,kappaRL,ndim)
                                                                                 
integer, intent(in) :: ndim                                                                  
complex(r64), dimension(ndim,ndim), intent(in) :: UL, VL, UR, VR
complex(r64), dimension(ndim,ndim), intent(out) :: rhoLR, kappaLR, kappaRL
complex(r64), dimension(ndim,ndim) :: URc, VRc

URc = conjg(UR)                                                                 
VRc = conjg(VR)                                                                 
                                                                                 
call zgemm('n','t',ndim,ndim,ndim,zone,VRc,ndim,VL,ndim,zzero,rhoLR,ndim)             
call zgemm('n','t',ndim,ndim,ndim,zone,VRc,ndim,UL,ndim,zzero,kappaLR,ndim)             
call zgemm('n','t',ndim,ndim,ndim,zone,VL,ndim,URc,ndim,zzero,kappaRL,ndim)             
                                                                                 
end subroutine calculate_densities

!------------------------------------------------------------------------------!
! subroutine calculate_overlap                                                 !
!                                                                              !
! Computes the overlap <L|ROT|R> for general Bogoliubov quasiparticle states   !
! |L> and |R> and where ROT is a general rotation matrix (e.g. spatial + gauge !
! + parity).                                                                   !
! The overlap is computed through the pfaffian formula found in the reference  !
! Avez.2012.PhysRevC.85.034325.                                                !
! The routine to computes the pfaffian of a skew-symmetric matrix is taken     !
! from the reference Wimmer.2012.ACM.TransMathSoftware.38.30.                  !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        occL = number of occupied states in the canonical basis of |L>        !
!        empL = number of empty    states in the canonical basis of |L>        !
!        ovvaL = overlap of |0> with the single-particule vacuum               !
!        vouL = values of v/u in the canonical basis of |L>                    !
!        DL = matrix D of the BMZ transformation for |L> (canon. basis)        !
!        occR, empR, ovvaR, vouR, DR = same but for the state |R>              !
!        ROT = rotation matrix                                                 !
!                                                                              !
! Ouput: overlap = complex overlap <L|ROT|R>                                   !
!------------------------------------------------------------------------------!
subroutine calculate_overlap(occL,empL,ovvaL,vouL,DL, &
                             occR,empR,ovvaR,vouR,DR,ROT,overlap,ndim)

integer, intent(in) :: ndim, occL, empL, occR, empR 
real(r64), intent(in) :: ovvaL, ovvaR
real(r64), dimension(ndim-occL), intent(in) :: vouL                 
real(r64), dimension(ndim-occR), intent(in) :: vouR                 
complex(r64), dimension(ndim,ndim), intent(in) :: DL, DR, ROT
complex(r64), intent(out) :: overlap
integer :: i, j, nreg, nempm, nblocL, nblocR, nL, nR, nE, info, islaL, islaR, &
           ialloc=0
complex(r64) :: detR, normvac, over_module, over_phase
real(r64) :: sgn, fac1, fac2, fac3, signfac
complex(r64), dimension(ndim,ndim) :: DRc, Rn
complex(r64), dimension(:,:), allocatable :: Mreg, RnE, Rinv
integer   , dimension(:), allocatable :: ipiv, iwork        
real(r64), dimension(:), allocatable :: rwork        
complex(r64), dimension(:), allocatable :: zwork, zwork2

!!! Computes the dimensions for the regularized matrix M, and computes the 
!!! phase factor accordingly 
nempm = min(empR,empL) 
nblocL = ndim - occL - nempm
nblocR = ndim - occR - nempm
nreg = nblocR + nblocL ! = 2*(n-nempm) - (occL+occR)

nL = ndim - occL
nR = ndim - occR
nE = ndim - nempm
fac1 = (-1)**(nE*(nE+1)/2)
fac2 = (-1)**(occL*(occL-1)/2)
fac3 = (-1)**(occL*(occL+1)/2 + occR*(occR+1)/2 + occL*(nL+nE) + occR*nR)
signfac = fac1 * fac2 * fac3

!!! When constructing the canonical basis, we store ovac = log(prod u) such 
!!! that we use sqrt(a*b) = exp(1/2 [log(a) + log(b)])
normvac = zone * (0.5d0 * (ovvaL + ovvaR))

!!! Check if we are dealing with Slater determinants
islaL = 0
if ( nL-empL == 0) islaL = 1
islaR = 0
if ( nR-empR == 0) islaR = 2

select case ( islaL+islaR )
  case (1)
    if ( occL > ndim-empR ) then
      overlap = zzero
      return
    endif
  case (2)
    if ( occR > ndim-empL ) then
      overlap = zzero
      return
    endif
  case (3)
    if ( occL == occR ) then
      overlap = zone  
    else
      overlap = zzero
      return
    endif
  case default
    continue
end select

!!! Determine the overlap matrix R between the canonical bases of the original
!!! left state and the rotated right state. The inversion and determinant 
!!! calculations of R are done using LU factorization
allocate( RnE(nE,nE), Rinv(nE,nE), ipiv(nE), zwork(nE), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation (1) in calculate overlap'

call zgemm('n','n',ndim,ndim,ndim,zone,ROT,ndim,DR,ndim,zzero,DRc,ndim)
call zgemm('c','n',ndim,ndim,ndim,zone,DL,ndim,DRc,ndim,zzero,Rn,ndim)

RnE(1:nE,1:nE) = Rn(1+nempm:ndim,1+nempm:ndim) 
Rinv = RnE

call zgetrf(nE,nE,Rinv,nE,ipiv,info)
if (info /= 0) then
  print*,'In calculate_overlap got info = ',info,' from zgetrf'
  stop 
endif

detR = zone
do i = 1, nE  
  detR = detR * Rinv(i,i)
enddo
 
sgn = one   
do i= 1, nE
  if (ipiv(i) /= i) then
    sgn = -sgn
  endif
enddo
detR = sgn * detR
         
call zgetri(nE,Rinv,nE,ipiv,zwork,nE,info)
if (info /= 0) then
  print*,'In calculate_overlap got info = ',info,' from zgetri'
  stop 
endif

!!! Construct the M matrix 
if ( islaL+islaR /= 3 ) then
  allocate( Mreg(nreg,nreg), zwork2(nreg**2), iwork(nreg), rwork(nreg), &
           stat=ialloc )
  if ( ialloc /= 0 ) stop 'Error during allocation (2) in calculate overlap'

  Mreg = zzero

  ! Upper left corner
  do i = 1, nblocR
    if ( (-1)**i == -1 ) then 
      Mreg(i,i+1) = vouR(i+nempm) 
    else
      Mreg(i,i-1) = vouR(i+nempm) 
    endif
  enddo

  ! Lower right corner
  do j = 1, nblocL
    i = j + nblocR
    if ( (-1)**i == -1 ) then 
      Mreg(i,i+1) = -vouL(j+nempm) 
    else
      Mreg(i,i-1) = -vouL(j+nempm) 
    endif
  enddo

  ! Upper right and bottom left corners
  do j = 1, nblocL
    do i = 1, nblocR
      Mreg(i,j+nblocR) = -Rinv(i,j) 
      Mreg(j+nblocR,i) =  Rinv(i,j) 
    enddo
  enddo

  ! Computes the pfaffian using routines from M. Wimmer
  call zskpfa('u','p',nreg,Mreg,nreg,overlap,iwork,zwork2,nreg**2,rwork,info)

  deallocate(Mreg,zwork2,iwork,rwork)
endif 


!!! Final value of the overlap
if ( abs(overlap) > 0.0d0 ) then
  over_phase = zimag * atan2(aimag(overlap),real(overlap))
  over_module = zone * log(abs(overlap)) 

  overlap = detR * signfac * exp(over_phase + over_module + normvac)
endif

deallocate(Rinv,RnE,ipiv,zwork)

end subroutine calculate_overlap

!------------------------------------------------------------------------------!
! subroutine calculate_norm                                                    !
!                                                                              ! 
! Computes the module of the overlap using the Onishi formula                  ! 
!   |<L|R>| = sqrt(|det(U)|) = sqrt(|det(UR^\dagger UL + VR^\dagger VL)|)      ! 
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      !
!        UL,VL = Bogoliubov matrices of the left  state                        ! 
!        UR,VR =      "         "    "   "  right   "                          ! 
!                                                                              ! 
! Output: norm = norm between the two states                                   ! 
!------------------------------------------------------------------------------!
subroutine calculate_norm(UL,VL,UR,VR,norm,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: UL, VL, UR, VR
real(r64), intent(out) :: norm 
integer :: i, info
integer, dimension(ndim) :: ipiv
complex(r64) :: detU 
complex(r64), dimension(ndim,ndim) :: UI       

!!! U
call zgemm('c','n',ndim,ndim,ndim,zone,UR,ndim,UL,ndim,zzero,UI,ndim)
call zgemm('c','n',ndim,ndim,ndim,zone,VR,ndim,VL,ndim, zone,UI,ndim)

call zgetrf(ndim,ndim,UI,ndim,ipiv,info)
if (info /= 0) then
  print*,' In calculate_norm got info = ',info,' from zgetrf'
  !stop 
endif

detU = zone
do i = 1, ndim
  detU = detU * UI(i,i)
enddo

norm = sqrt(abs(detU))

end subroutine calculate_norm

END MODULE Wavefunctions 
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
