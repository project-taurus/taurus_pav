!==============================================================================!
! MODULE Radius                                                                !
!                                                                              !
! This module contains the variables and routines related to radius squared    !
! one-body operator.                                                           !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_radius                                                      !
! - subroutine calculate_radius                                                !
! - subroutine calculate_fields_r1r2                                           !
!==============================================================================!
MODULE Radius 

use Basis
!cmpi use MPI

implicit none 
public

!!! Factors for proton/neutron/matter depending on com corrections or not
real(r64), dimension(3,-1:1) :: radius_fac1b, & ! one-body
                                radius_fac2b    ! two-body

!!! Matrix elements in m-scheme (HO basis)  
complex(r64), dimension(:,:), allocatable :: radius_r2p, & ! r^2 protons  (1b)
                                             radius_r2n    ! r^2 neutrons (1b)
real(rH2), dimension(:), allocatable :: radius_H2 ! 2-body part 
integer(i64) :: radius_H2dim, radius_H2dim_all    ! number of 2BME stored 
integer(i16), dimension(:), allocatable :: radius_abcd  ! indices of 2BME   
integer(i8), dimension(:), allocatable :: radius_trperm ! time reversal permut.

!!! Fields in the single-particle basis
complex(r64), dimension(:,:,:), allocatable :: radius_gammaLR, & ! Gamma^LR
                                               radius_deltaLR    ! Delta^LR

!!! Other quantities
logical :: is_r1r2 ! switch to compute 2-body operator r1r2

!!! Private routines
private :: calculate_fields_r1r2, read_radius_reduced

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_radius                                                        !
!                                                                              ! 
! Defines the HO matrix elements of the radius square one-body operator.       ! 
!                                                                              ! 
! Input: opt_oper = option to read the matrix elements from a file             !
!------------------------------------------------------------------------------!
subroutine set_radius(opt_oper)

integer, intent(in) :: opt_oper
integer :: ia, ib, ja, jb, mja, mjb, la, lb, utn, ialloc=0
complex(r64) :: xrad
character(24) :: namefile
character(26) :: namefiler1r2
logical :: is_exist

is_exist = .false.

allocate( radius_r2p(HOsp_dim/2,HOsp_dim/2), &
          radius_r2n(HOsp_dim/2,HOsp_dim/2), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of square radius operator'

radius_r2p = zzero
radius_r2n = zzero

!!! Factors 
radius_fac1b(1,-1) = one / valence_Z 
radius_fac1b(1,+1) = zero
radius_fac1b(2,-1) = zero
radius_fac1b(2,+1) = one / valence_N
radius_fac1b(3,-1) = one / valence_A
radius_fac1b(3,+1) = one / valence_A

!!! Reads the one-body part from a file
if ( opt_oper == 1 ) then 
  namefile = "matelem_radius_r2_1b.bin"
  inquire(file=namefile, exist=is_exist)
 
  if ( is_exist ) then
    open(newunit=utn, file=namefile, status='old', action='read', &
         form='unformatted')
 
    do ib = 1, HOsp_dim/2
      do ia = 1, HOsp_dim/2
        read(utn) radius_r2p(ia,ib), radius_r2n(ia,ib)
      enddo
    enddo
 
    close(utn, status="keep")
    return
  endif
endif

!!! Standard definition
do ib = 1, HOsp_dim/2
  lb = HOsp_l(ib)
  jb = HOsp_2j(ib)
  mjb = HOsp_2mj(ib)
  do ia = 1, HOsp_dim/2
    la = HOsp_l(ia)
    ja = HOsp_2j(ia)
    mja = HOsp_2mj(ia)
    if ( (ja /= jb) .or. (mja /= mjb) .or. (la /= lb) ) cycle

    xrad = radial_integral(ia,ib,2) * zone
    radius_r2p(ia,ib) = xrad
    radius_r2n(ia,ib) = xrad
  enddo
enddo

!!! Reads the two-body part from a file
is_r1r2 = .false.

if ( opt_oper == 1 ) then 
  namefiler1r2 = "matelem_radius_r1r2_2b.bin"
  inquire(file=namefiler1r2, exist=is_exist)

  if ( is_exist ) then
    is_r1r2 = .true.

    !!! Sets the factors: 1-body
    radius_fac1b(1,-1) = ( valence_A * (valence_A - 2) + valence_Z ) * one / &
                      ( valence_Z ) 
    radius_fac1b(1,+1) = one 
    radius_fac1b(2,-1) = one 
    radius_fac1b(2,+1) = ( valence_A * (valence_A - 2) + valence_N ) * one / &
                      ( valence_N ) 
    radius_fac1b(3,-1) = ( valence_A - 1 ) * one
    radius_fac1b(3,+1) = ( valence_A - 1 ) * one

    radius_fac1b = radius_fac1b / ( valence_A**2)

    !!! Sets the factors: 2-body
    radius_fac2b(1,-1) = +( 2*valence_A - valence_Z ) * one / ( valence_Z )
    radius_fac2b(2,-1) = - one
    radius_fac2b(3,-1) = + one
                        
    radius_fac2b(1, 0) = +( valence_A - valence_Z ) * one / ( valence_Z )
    radius_fac2b(2, 0) = +( valence_A - valence_N ) * one / ( valence_N )
    radius_fac2b(3, 0) = + one
                        
    radius_fac2b(1,+1) = - one 
    radius_fac2b(2,+1) = +( 2*valence_A - valence_N ) * one / ( valence_N )
    radius_fac2b(3,+1) = + one 

    ! Remark: minus sign compated to the paper of Cipollone
    radius_fac2b = radius_fac2b * 2 * hbarmass / ( HO_hw * valence_A**2 )

    !!! Sets the fields
    allocate( radius_gammaLR(3,HOsp_dim,HOsp_dim), &
              radius_deltaLR(3,HOsp_dim,HOsp_dim), &
              stat=ialloc ) 
    if ( ialloc /= 0 ) stop 'Error during allocation of fields'
    
    radius_gammaLR = zzero
    radius_deltaLR = zzero
    
    !!! Reads the matrix elements from the file 
    call read_radius_reduced
    
  endif
endif

end subroutine set_radius     

!------------------------------------------------------------------------------!
! subroutine calculate_radius                                                  !
!                                                                              ! 
! Expectation value for the one-body squared radius operator                   ! 
! < r^2 > = Tr(r^2 rho) + 1/2 Tr(gammaLR * rhoLR) + 1/2 Tr(deltaLR * kappaRL)  ! 
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      !
!        rhoLR,kappaLR,kappRL = transition densities                           !
!                                                                              ! 
! Output: ra2_x  = < r^2 >_x (x: proton, neutron, matter)                      !
!------------------------------------------------------------------------------!
subroutine calculate_radius(rhoLR,kappaLR,kappaRL,ra2p,ra2n,ra2m,ndim)                                       

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL
complex(r64), intent(out) :: ra2p, ra2n, ra2m
integer :: hdim, i, j
complex(r64) :: ra2p_bare, ra2n_bare, ra2_sum
complex(r64), dimension(ndim/2,ndim/2) :: rho_p, rho_n, A1, A2
complex(r64), dimension(ndim,ndim) :: A3, A4

!!! One-body part
hdim = ndim/2

do j = 1, hdim
  do i = 1, hdim 
    rho_p(i,j) = rhoLR(i,j)
    rho_n(i,j) = rhoLR(i+hdim,j+hdim)
  enddo
enddo

call zgemm('n','n',hdim,hdim,hdim,zone,radius_r2p,hdim,rho_p,hdim,zzero,A1,hdim)
call zgemm('n','n',hdim,hdim,hdim,zone,radius_r2n,hdim,rho_n,hdim,zzero,A2,hdim)

ra2p_bare = zzero
ra2n_bare = zzero

do i = 1, hdim
  ra2p_bare = ra2p_bare + A1(i,i)
  ra2n_bare = ra2n_bare + A2(i,i)
enddo             

ra2p = radius_fac1b(1,-1) * ra2p_bare + radius_fac1b(1,+1) * ra2n_bare
ra2n = radius_fac1b(2,-1) * ra2p_bare + radius_fac1b(2,+1) * ra2n_bare
ra2m = radius_fac1b(3,-1) * ra2p_bare + radius_fac1b(3,+1) * ra2n_bare

!!! Two-body part
if ( is_r1r2 ) then
  call calculate_fields_r1r2(rhoLR,kappaLR,radius_gammaLR,radius_deltaLR,ndim)

  do i = 1, 3
    call zgemm('n','n',ndim,ndim,ndim,zone,radius_gammaLR(i,:,:),ndim,rhoLR, &
               ndim,zzero,A3(:,:),ndim)
    call zgemm('n','t',ndim,ndim,ndim,zone,radius_deltaLR(i,:,:),ndim,kappaRL, &
               ndim,zzero,A4(:,:),ndim)

    ra2_sum = zzero
    do j = 1, ndim
      ra2_sum = ra2_sum + A3(j,j) + A4(j,j)
    enddo

    if ( i == 1 ) then
      ra2p = ra2p + 0.5d0 * ra2_sum
    elseif ( i == 2 ) then
      ra2n = ra2n + 0.5d0 * ra2_sum
    else
      ra2m = ra2m + 0.5d0 * ra2_sum
    endif
  enddo
endif

end subroutine calculate_radius

!------------------------------------------------------------------------------!
! subroutine calculate_fields_r1r2                                             !
!                                                                              !
! Calculates the HFB fields h, Gamma an Delta for the two-body operator r1r2   !
! which are used to compute the center-of-mass correction.                     !
!                                                                              !
! The structure is similar to the one of the energy fields, but with the       !
! differeence that we compute different fields for proton, neutron and matter  !
! radii using the factors defined in: Cipollone.2015.PRC.92.014306.            !
!                                                                              !
! Order for the indices:                                                       !
! 1: proton                                                                    !
! 2: neutron                                                                   !
! 3: matter                                                                    !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        rhoLR,kappaLR = transition densities                                  !
!                                                                              !
! Output: gammaLR,deltaLR = transisition fields                                !
!------------------------------------------------------------------------------!
subroutine calculate_fields_r1r2(rhoLR,kappaLR,gammaLR,deltaLR,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR
complex(r64), dimension(3,ndim,ndim), intent(out) :: gammaLR, deltaLR 
integer :: i, j, ia, ib, ic, id, it, perm
integer(i64) :: kk
real(r64), dimension(3) :: h2b, f2b
!cmpi integer :: ierr=0
!cmpi complex(r64), dimension(3,ndim,ndim) :: gammaLR_red, deltaLR_red

gammaLR = zzero
deltaLR = zzero

!$OMP PARALLEL DO FIRSTPRIVATE(rhoLR,kappaLR) &
!$OMP             PRIVATE(ia,ib,ic,id,h2b,f2b,perm,it) &
!$OMP             REDUCTION(+:gammaLR,deltaLR)
do kk = 1, radius_H2dim  
  ia = radius_abcd(1+4*(kk-1))
  ib = radius_abcd(2+4*(kk-1))
  ic = radius_abcd(3+4*(kk-1))
  id = radius_abcd(4+4*(kk-1))
  it = ( HOsp_2mt(ia) + HOsp_2mt(ib) ) / 2
  h2b(:) = radius_fac2b(:,it) * radius_H2(kk)
  perm = radius_trperm(kk)

  !!! Loop on time reversal
  do it = 1, 2 
    if ( it == 2 ) then
      if ( HOsp_2mj(ia) + HOsp_2mj(ib) == 0 ) cycle
      call find_timerev(perm,ia,ib,ic,id)
      h2b = sign(one,perm*one) * h2b
    endif

    !!! Faster than using if ((a /= c).or.(b /= d))
    f2b = h2b * (1 - kdelta(ia,ic) * kdelta(ib,id))

    !!! Calculation of Gamma^LR
    gammaLR(:,ia,ic) = gammaLR(:,ia,ic) + h2b(:) * rhoLR(id,ib)
    gammaLR(:,ia,id) = gammaLR(:,ia,id) - h2b(:) * rhoLR(ic,ib)
    gammaLR(:,ic,ia) = gammaLR(:,ic,ia) + f2b(:) * rhoLR(ib,id)
    gammaLR(:,ic,ib) = gammaLR(:,ic,ib) - f2b(:) * rhoLR(ia,id)
   
    !!! Calculation of Delta^LR
    deltaLR(:,ib,ia) = deltaLR(:,ib,ia) + h2b(:) * kappaLR(id,ic)
    deltaLR(:,id,ic) = deltaLR(:,id,ic) + f2b(:) * kappaLR(ib,ia)

  enddo
enddo
!$OMP END PARALLEL DO

!!! Reduces the values for the processes in the same team                        
!cmpi if ( paral_myteamsize > 1 ) then
!cmpi   call mpi_reduce(gammaLR,gammaLR_red,3*ndim**2,mpi_double_complex, &        
!cmpi                   mpi_sum,0,mpi_comm_team,ierr)                            
!cmpi   call mpi_reduce(deltaLR,deltaLR_red,3*ndim**2,mpi_double_complex, &        
!cmpi                   mpi_sum,0,mpi_comm_team,ierr)                            
!cmpi   gammaLR = gammaLR_red                                                    
!cmpi   deltaLR = deltaLR_red                                                    
!cmpi endif          

!!! Symmetry from sums (included here)       
gammaLR = 2.0d0 * gammaLR
  
!!! Skew symmetry 
do j = 1, HOsp_dim
  do i = 1, j-1      
    deltaLR(:,i,j) = -1.0d0 * deltaLR(:,j,i) 
  enddo
enddo

end subroutine calculate_fields_r1r2  

!------------------------------------------------------------------------------!
! subroutine read_radius_reduced                                               !
!                                                                              !
! Reads the reduced radius two-body matrix elements. Very similar to the       !
! equivalent routine for the energy.                                           !
!------------------------------------------------------------------------------!
subroutine read_radius_reduced
  
integer(i32) :: bytes_H2, utn, ialloc=0
integer(i64) :: kk, pos_ini, pos_abcd, pos_H2, pos_perm
character(26) :: namefiler1r2
!cmpi integer(i64) :: divide, rest
  
namefiler1r2 = "matelem_radius_r1r2_2b.bin"
  
open(newunit=utn, file=namefiler1r2, status='old', action='read', &
     access='stream', form='unformatted')
    
read(utn) radius_H2dim_all

radius_H2dim = radius_H2dim_all

!!! Determines the position of the arrays in the file. This is not needed in a 
!!! sequential run but I prefer to have the same code structure as in the 
!!! MPI version.
inquire(utn, pos=pos_ini)

if ( rH2 == r32 ) then 
  bytes_H2 = 4
else
  bytes_H2 = 8
endif

pos_abcd = pos_ini
pos_H2 = pos_ini + radius_H2dim_all * 8 
pos_perm = pos_ini + radius_H2dim_all * (8 + bytes_H2)

!!! MPI: distribution of matrix elements among the membrs of the teams. 
!!! First computes the number of matrix elements and then their position in the
!!! files. 
!cmpi divide = radius_H2dim_all / paral_myteamsize
!cmpi rest = modulo(radius_H2dim_all,paral_myteamsize)

!cmpi radius_H2dim = divide
!cmpi if ( paral_myteamrank < rest ) then
!cmpi   radius_H2dim = radius_H2dim + 1
!cmpi endif

!cmpi pos_abcd = pos_ini + paral_myteamrank * radius_H2dim * 8 
!cmpi pos_H2 = pos_ini + radius_H2dim_all * 8 + paral_myteamrank * &
!cmpi                                           radius_H2dim * bytes_H2
!cmpi pos_perm = pos_ini + radius_H2dim_all * (8 + bytes_H2)  & 
!cmpi            + paral_myteamrank * radius_H2dim * 1

!cmpi if ( paral_myteamrank >= rest ) then
!cmpi   pos_abcd = pos_abcd + rest * 8
!cmpi   pos_H2 = pos_H2 + rest * bytes_H2
!cmpi   pos_perm = pos_perm + rest 
!cmpi endif

!!! Reads the arrays related to the two-body matrix elements
allocate ( radius_H2(radius_H2dim), radius_abcd(4*radius_H2dim), & 
           radius_trperm(radius_H2dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of array of indices'

read(utn, pos=pos_abcd) (radius_abcd(kk), kk=1, 4*radius_H2dim)
read(utn, pos=pos_H2) (radius_H2(kk), kk=1, radius_H2dim)
read(utn, pos=pos_perm) (radius_trperm(kk), kk=1, radius_H2dim)
  
close(utn, status='keep')

end subroutine read_radius_reduced

END MODULE Radius    
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
