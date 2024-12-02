!==============================================================================!
! MODULE Fields                                                                !
!                                                                              !
! This module contains the variables and routines related to the fields and    !
! the matrix elements of H in the QP basis.                                    !
!                                                                              !
! Remark: at some point, this module should be merged with the Hamiltonian     !
! one.                                                                         !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_fields                                                      !
! - subroutine calculate_fields                                                !
!==============================================================================!
MODULE Fields           

use Wavefunctions  
use Hamiltonian    
!cmpi use MPI

implicit none
public

!!! Fields in the single-particle basis
complex(r64), dimension(:,:), allocatable :: field_hspLR,   & ! h^LR         
                                             field_gammaLR, & ! Gamma^LR
                                             field_deltaLR    ! Delta^LR 

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_fields                                                        !
!                                                                              !
! Allocates the arrays for the fields.                                         !
!                                                                              !
! Input: opt_phys = option controlling the physics case                        !
!------------------------------------------------------------------------------!
subroutine set_fields(opt_phys)

integer, intent(in) :: opt_phys
integer :: ialloc=0

if ( opt_phys == 0 ) then
  allocate( field_hspLR(HOsp_dim,HOsp_dim), field_gammaLR(HOsp_dim,HOsp_dim), &
            field_deltaLR(HOsp_dim,HOsp_dim), stat=ialloc ) 
  if ( ialloc /= 0 ) stop 'Error during allocation of fields'

  field_hspLR   = zzero
  field_gammaLR = zzero
  field_deltaLR = zzero
endif

end subroutine set_fields

!------------------------------------------------------------------------------!
! subroutine calculate_fields                                                  !
!                                                                              !
! Calculates the HFB fields h, Gamma an Delta which are then used to compute   !
! other quantities of interest (in particular the energy).                     !
!                                                                              !
! What is calculated:                                                          !
!   h^LR_{ac} = Gamma^LR_{ac} + t_{ac}                                         !
!   Gamma^LR_{ac} = sum_{bd} V_{abcd} rho^LR_{db}                              !
!   Delta^LR_{ab} = 1/2 sum_{cd}  V_{abcd} kappa^LR_{cd}                       !
!                 =     sum_{c<d} V_{abcd} kappa^LR_{cd}                       !
!                                                                              !
! The fields have the symmetry properties:                                     !
!   Delta^LR_{ab} = - Delta^LR_{ba}    (Skew symmetry)                         !
!                                                                              !
! The actual calculation below uses these symmetries and the fact that only    !
! a subset of the matrix elements of the interaction are stored.               !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        rhoLR,kappaLR = transition densities                                  !
!                                                                              !
! Output: gammaLR,deltaLR = transisition fields                                !
!------------------------------------------------------------------------------!
subroutine calculate_fields(rhoLR,kappaLR,gammaLR,hspLR,deltaLR,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR
complex(r64), dimension(ndim,ndim), intent(out) :: gammaLR, hspLR, deltaLR 
integer :: i, j, ia, ib, ic, id, it, perm
integer(i64) :: kk
real(r64) :: h2b, f2b
!cmpi integer :: ierr=0
!cmpi complex(r64), dimension(ndim,ndim) :: gammaLR_red, deltaLR_red

gammaLR = zzero
deltaLR = zzero

!$OMP PARALLEL DO FIRSTPRIVATE(rhoLR,kappaLR) &
!$OMP             PRIVATE(ia,ib,ic,id,h2b,f2b,perm,it) &
!$OMP             REDUCTION(+:gammaLR,deltaLR)
do kk = 1, hamil_H2dim  
  ia = hamil_abcd(1+4*(kk-1))
  ib = hamil_abcd(2+4*(kk-1))
  ic = hamil_abcd(3+4*(kk-1))
  id = hamil_abcd(4+4*(kk-1))
  h2b = hamil_H2(kk)
  perm = hamil_trperm(kk)



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
    gammaLR(ia,ic) = gammaLR(ia,ic) + h2b * rhoLR(id,ib)
    gammaLR(ia,id) = gammaLR(ia,id) - h2b * rhoLR(ic,ib)
    gammaLR(ic,ia) = gammaLR(ic,ia) + f2b * rhoLR(ib,id)
    gammaLR(ic,ib) = gammaLR(ic,ib) - f2b * rhoLR(ia,id)
   
    !!! Calculation of Delta^LR
    deltaLR(ib,ia) = deltaLR(ib,ia) + h2b * kappaLR(id,ic)
    deltaLR(id,ic) = deltaLR(id,ic) + f2b * kappaLR(ib,ia)

  enddo
enddo
!$OMP END PARALLEL DO

stop

!!! Reduces the values for the processes in the same team                        
!cmpi if ( paral_myteamsize > 1 ) then
!cmpi   call mpi_reduce(gammaLR,gammaLR_red,ndim**2,mpi_double_complex, &        
!cmpi                   mpi_sum,0,mpi_comm_team,ierr)                            
!cmpi   call mpi_reduce(deltaLR,deltaLR_red,ndim**2,mpi_double_complex, &        
!cmpi                   mpi_sum,0,mpi_comm_team,ierr)                            
!cmpi   gammaLR = gammaLR_red                                                    
!cmpi   deltaLR = deltaLR_red                                                    
!cmpi endif          
       
gammaLR = 2.0d0 * gammaLR
  
!!! h = Gamma + 1body (BB: not used right now)
do j = 1, HOsp_dim
  do i = 1, HOsp_dim
    hspLR(i,j) = gammaLR(i,j) + hamil_H1(i,j)
  enddo
enddo

!!! Skew symmetry 
do j = 1, HOsp_dim
  do i = 1, j-1      
    deltaLR(i,j) = -1.0d0 * deltaLR(j,i) 
  enddo
enddo

end subroutine calculate_fields

END MODULE Fields        
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
