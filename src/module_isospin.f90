!==============================================================================!
! MODULE Isospin                                                               !
!                                                                              !
! This module contains the variables and routines related to the isospin.      !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_isospin                                                     !
! - subroutine calculate_isospin                                               !
!==============================================================================!
MODULE Isospin

use Basis

implicit none 
public

complex(r64), dimension(:,:), allocatable :: isospin_Tx,  & ! Tx
                                             isospin_Ty,  & ! Ty
                                             isospin_Tz,  & ! Tz
                                             isospin_Tx2, & ! Tx^2
                                             isospin_Ty2, & ! Ty^2
                                             isospin_Tz2    ! Tz^2

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_isospin                                                       !
!                                                                              ! 
! Defines the HO matrix elements of the isospin operators in cartesian coord-  ! 
! inates.                                                                      ! 
!------------------------------------------------------------------------------!
subroutine set_isospin       

integer :: i, j, ialloc=0

allocate( isospin_Tx(HOsp_dim,HOsp_dim), isospin_Tx2(HOsp_dim,HOsp_dim), &                     
          isospin_Ty(HOsp_dim,HOsp_dim), isospin_Ty2(HOsp_dim,HOsp_dim), &                     
          isospin_Tz(HOsp_dim,HOsp_dim), isospin_Tz2(HOsp_dim,HOsp_dim), &                     
!         isospin_Tp(HOsp_dim,HOsp_dim), isospin_Tm(HOsp_dim,HOsp_dim),  &
          stat=ialloc )                                                          
if ( ialloc /= 0 ) stop 'Error during allocation of isospin operators'

isospin_Tx  = zzero
isospin_Tx2 = zzero
isospin_Ty  = zzero
isospin_Ty2 = zzero
isospin_Tz  = zzero
isospin_Tz2 = zzero
!isospin_Tp = zzero
!isospin_Tm = zzero

!!! Tx, Ty, Tz
do i = 1, HOsp_dim/2
  j = i + HOsp_dim/2
  isospin_Tx(i,j) = 0.5d0 
  isospin_Tx(j,i) = 0.5d0 
  isospin_Ty(i,j) =  zimag * 0.5d0
  isospin_Ty(j,i) = -zimag * 0.5d0
  isospin_Tz(i,i) = HOsp_2mt(i) * 0.5d0
  isospin_Tz(j,j) = HOsp_2mt(j) * 0.5d0
enddo 

!!! Tx^2, Ty^2, Tz^2
do i = 1, HOsp_dim
  isospin_Tx2(i,i)= 0.25d0
  isospin_Ty2(i,i)= 0.25d0
  isospin_Tz2(i,i)= 0.25d0
enddo

end subroutine set_isospin    

!------------------------------------------------------------------------------!
! subroutine calculate_isospin                                                 !
!                                                                              ! 
! Expectation values of the angular momentum operator and its square in        ! 
! cartesian coordinates.                                                       ! 
! < T_i > = Tr(t_i rho)                                                        ! 
! < T_i^2 > = Tr(t_i^2 rho) + Tr(t_i^T kappa10^dagger t_i kappa01)             ! 
!             - Tr([t_i rho]^2) + Tr(t_i rho)^2                                ! 
! < T^2 > = < T_x^2 > + < T_y^2 > + < T_z^2 >                                  ! 
!                                                                              ! 
! To speed-up the calculation, we could avoid the calculations of T_z using    ! 
! T_z = (N-Z)/2                                                                ! 
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      !
!        io = 0 computes < T_i >, < T_i^2 > with i = x,y,z                     !
!           = 1 computes < T_i >                                               !
!        rhoLR,kappaLR,kappaRL = transition densities                          !
!                                                                              !
! Output: isoti, isoti2 = < T_i >, < T_i^2 > with i = x,y,z                    !
!------------------------------------------------------------------------------!
subroutine calculate_isospin(io,rhoLR,kappaLR,kappaRL, & 
                             isotx,isoty,isotz,isotx2,isoty2,isotz2,ndim)

integer, intent(in) :: io, ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL
complex(r64), intent(out) :: isotx, isoty, isotz, isotx2, isoty2, isotz2  
integer :: i, j, imin, imax
complex(r64) :: tr1, tr2, tr5, tr6, ttmp, t2tmp
complex(r64), dimension(ndim,ndim) :: Ti, Ti2, A1, A2, A3, A4, A5, A6

!!! Determines what components have to be calculated
if ( io == 0 ) then 
  imin = 1
  imax = 3
elseif ( io == 1 ) then
  imin = 3
  imax = 3
elseif ( io == 2 ) then
  imin = 1
  imax = 2
else
  write(uto,*) 'In calculate_isospin: wrong argument io = ', io 
  stop 
endif

!!! Computes the expectation values
isotx  = zzero
isoty  = zzero
isotz  = zzero
isotx2 = zzero
isoty2 = zzero
isotz2 = zzero

do i = imin, imax

  if ( i == 1 ) then       ! X
    Ti  = isospin_Tx
    Ti2 = isospin_Tx2
  elseif ( i == 2 ) then   ! Y
    Ti  = isospin_Ty
    Ti2 = isospin_Ty2
  else                     ! Z
    Ti  = isospin_Tz
    Ti2 = isospin_Tz2
  endif

  call zgemm('n','n',ndim,ndim,ndim,zone, Ti,ndim,  rhoLR,ndim,zzero,A1,ndim)
  call zgemm('n','n',ndim,ndim,ndim,zone,Ti2,ndim,  rhoLR,ndim,zzero,A2,ndim)
  call zgemm('t','t',ndim,ndim,ndim,zone, Ti,ndim,kappaRL,ndim,zzero,A3,ndim)   
  call zgemm('n','n',ndim,ndim,ndim,zone, Ti,ndim,kappaLR,ndim,zzero,A4,ndim)   
  call zgemm('n','n',ndim,ndim,ndim,zone, A1,ndim,     A1,ndim,zzero,A5,ndim)   
  call zgemm('n','n',ndim,ndim,ndim,zone, A3,ndim,     A4,ndim,zzero,A6,ndim)

  tr1 = zzero
  tr2 = zzero
  tr5 = zzero
  tr6 = zzero

  do j = 1, ndim
    tr1 = tr1 + A1(j,j)
    tr2 = tr2 + A2(j,j)
    tr5 = tr5 + A5(j,j)
    tr6 = tr6 + A6(j,j)
  enddo

  ttmp  = tr1 
  t2tmp = tr1**2 + tr2 - tr5 + tr6

  if ( i == 1 ) then
    isotx  = ttmp
    isotx2 = t2tmp
  elseif ( i == 2 ) then
    isoty  = ttmp
    isoty2 = t2tmp
  else
    isotz  = ttmp
    isotz2 = t2tmp
  endif

enddo

end subroutine calculate_isospin

END MODULE Isospin   
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
