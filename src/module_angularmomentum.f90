!==============================================================================!
! MODULE AngularMomentum                                                       !
!                                                                              !
! This module contains the variables and routines related to angular momentum  !
! operators.                                                                   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_angularmomentum                                             !
! - subroutine calculate_angularmomentum                                       !
!==============================================================================!
MODULE AngularMomentum

use Basis

implicit none 
public

complex(r64), dimension(:,:), allocatable :: angumome_Jx,  & ! Jx
                                             angumome_Jy,  & ! Jy
                                             angumome_Jz,  & ! Jz
                                             angumome_Jx2, & ! Jx^2
                                             angumome_Jy2, & ! Jy^2
                                             angumome_Jz2, & ! Jz^2
                                             angumome_so     ! spin-orbit   

complex(r64), dimension(:,:,:), allocatable :: angumome_S,   & ! S_mu (=-1,0,1)
                                               angumome_L,   & ! L_mu 
                                               angumome_Sgr, & ! S_mu grad(rY)
                                               angumome_Lgr    ! L_mu grad(rY)

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_angularmomentum                                               !
!                                                                              ! 
! Defines the HO matrix elements of the angular-momentum operators.            ! 
!------------------------------------------------------------------------------!
subroutine set_angularmomentum

integer :: i, j, ni, nj, ji, jj, li, lj, mji, mjj, mli, mlj, ij, im, &
           ms1, ms2, mu, msp, mp, msi, msj, ialloc=0
real(r64) :: xj, xm, cb1, cb2, cb3, cb4, cb5, cb6, cb7, cb8, xl, xml, xms, &
           fac1, fac2, fac3, fac4
complex(r64), dimension(:,:), allocatable :: Jp,  Jm

!!! TOTAL ANGULAR MOMENTUM
allocate( angumome_Jx(HOsp_dim,HOsp_dim), angumome_Jx2(HOsp_dim,HOsp_dim), &
          angumome_Jy(HOsp_dim,HOsp_dim), angumome_Jy2(HOsp_dim,HOsp_dim), &
          angumome_Jz(HOsp_dim,HOsp_dim), angumome_Jz2(HOsp_dim,HOsp_dim), &
          Jp(HOsp_dim,HOsp_dim), Jm(HOsp_dim,HOsp_dim), &
          stat=ialloc )                                                          
if ( ialloc /= 0 ) stop 'Error during allocation of angular momentum operators'

angumome_Jx  = zzero
angumome_Jy  = zzero
angumome_Jz  = zzero
angumome_Jx2 = zzero
angumome_Jy2 = zzero
angumome_Jz2 = zzero
Jp = zzero
Jm = zzero

do i = 1, HOsp_dim
  ij = HOsp_2j(i)
  im = HOsp_2mj(i)
  xj = ij / 2.0d0
  xm = im / 2.0d0
  angumome_Jz(i,i) = xm 
  angumome_Jz2(i,i)= xm**2
  if ( im /= -ij ) Jp(i, i+1) = sqrt( xj*(xj+1) - xm*(xm-1) )
  if ( im /=  ij ) Jm(i, i-1) = sqrt( xj*(xj+1) - xm*(xm+1) )
enddo

angumome_Jx = zone  * ( Jm + Jp) / 2.0d0 
angumome_Jy = zimag * ( Jm - Jp) / 2.0d0 

call zgemm('n','n',HOsp_dim,HOsp_dim,HOsp_dim,zone,angumome_Jx,HOsp_dim, &
           angumome_Jx,HOsp_dim,zzero,angumome_Jx2,HOsp_dim)
call zgemm('n','n',HOsp_dim,HOsp_dim,HOsp_dim,zone,angumome_Jy,HOsp_dim, &
           angumome_Jy,HOsp_dim,zzero,angumome_Jy2,HOsp_dim)

deallocate(Jp, Jm)

!!! SPIN
allocate(angumome_S(HOsp_dim/2,HOsp_dim/2,-1:+1), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of angular momentum operators'

angumome_S = zzero

do i = 1, HOsp_dim/2
  ni = HOsp_n(i)                                                                
  ji = HOsp_2j(i)                                                                   
  li = 2*HOsp_l(i)                                                                   
  mji = HOsp_2mj(i)                                                                 
  do j = 1, HOsp_dim/2
    nj = HOsp_n(j)                                                              
    jj = HOsp_2j(j)                                                                 
    lj = 2*HOsp_l(j)                                                                 
    mjj = HOsp_2mj(j)                                                               
    if ( (ni /= nj) .or. (li /= lj) ) cycle                    
                                                                                 
    do mli = -li, li, 2                                                          
      do ms1 = -1, 1, 2                                                          
        do ms2 = -1, 1, 2                                                        
          xms = ms2 / 2.0d0                           
          call ClebschGordan(li,1,ji,mli,ms1,mji,cb1)                            
          call ClebschGordan(li,1,jj,mli,ms2,mjj,cb2)  
          cb3 = sqrt(3.d0/4.d0 - xms*(xms+1.0)) * kdelta(ms1,ms2+2)
          angumome_S(i,j,+1) = angumome_S(i,j,+1) + cb1 * cb2 * cb3 
          cb3 = xms * kdelta(ms1,ms2)
          angumome_S(i,j, 0) = angumome_S(i,j, 0) + cb1 * cb2 * cb3 
          cb3 = sqrt(3.d0/4.d0 - xms*(xms-1.0)) * kdelta(ms1,ms2-2)
          angumome_S(i,j,-1) = angumome_S(i,j,-1) + cb1 * cb2 * cb3 
        enddo                                                                    
      enddo                                                                      
    enddo                                                                        
                                                                                 
  enddo                                                                          
enddo            

!!! ORBITAL
allocate(angumome_L(HOsp_dim/2,HOsp_dim/2,-1:+1), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of angular momentum operators'

angumome_L = zzero

do i = 1, HOsp_dim/2
  ni = HOsp_n(i)                                                                
  ji = HOsp_2j(i)                                                                   
  li = 2*HOsp_l(i)                                                                   
  mji = HOsp_2mj(i)                                                                 
  do j = 1, HOsp_dim/2
    nj = HOsp_n(j)                                                              
    jj = HOsp_2j(j)                                                                 
    lj = 2*HOsp_l(j)                                                                 
    xl = HOsp_l(j) * one                                                           
    mjj = HOsp_2mj(j)                                                               
    if ( (ni /= nj) .or. (li /= lj) ) cycle                    
                                                                                 
    do ms1 = -1, 1, 2                                                          
      do mli = -li, li, 2                                                          
        do mlj = -li, li, 2
          xml = mlj / 2.0d0                           
          call ClebschGordan(li,1,ji,mli,ms1,mji,cb1)                            
          call ClebschGordan(li,1,jj,mlj,ms1,mjj,cb2)                            
          cb3 = sqrt(xl*(xl+1.0) - xml*(xml+1.0)) * kdelta(mli,mlj+2)
          angumome_L(i,j,+1) = angumome_L(i,j,+1) + cb1 * cb2 * cb3 
          cb3 = xml * kdelta(mli,mlj)
          angumome_L(i,j, 0) = angumome_L(i,j, 0) + cb1 * cb2 * cb3 
          cb3 = sqrt(xl*(xl+1.0) - xml*(xml-1.0)) * kdelta(mli,mlj-2)
          angumome_L(i,j,-1) = angumome_L(i,j,-1) + cb1 * cb2 * cb3 
        enddo                                                                    
      enddo                                                                      
    enddo                                                                        
                                                                                 
  enddo                                                                          
enddo            

!!! Conversion from ladder to spherical operators
do mu = -1, 1, 2
  angumome_S(:,:,mu) = (-mu/sqrt(2.0d0)) *  angumome_S(:,:,mu) 
  angumome_L(:,:,mu) = (-mu/sqrt(2.0d0)) *  angumome_L(:,:,mu) 
enddo

!!! L grad(r^1 Y_1m) and S grad(r^1 Y_1m)
! (this part can be easily generalized to handle any lambda) 
allocate(angumome_Lgr(HOsp_dim/2,HOsp_dim/2,-2:+2), &
         angumome_Sgr(HOsp_dim/2,HOsp_dim/2,-2:+2), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of angular momentum operators'

angumome_Lgr = zzero
angumome_Sgr = zzero

do i = 1, HOsp_dim/2
  ji = HOsp_2j(i)                                                                   
  li = 2*HOsp_l(i)                                                                   
  mji = HOsp_2mj(i)                                                                 
  do j = 1, HOsp_dim/2
    jj = HOsp_2j(j)                                                                 
    lj = 2*HOsp_l(j)                                                                 
    xl = HOsp_l(j) * one                                                           
    mjj = HOsp_2mj(j)                                                               
    
    call ClebschGordan(lj,2,li,0,0,0,cb1)   
    fac1 = cb1 * sqrt(3*(lj+one)/((li+one)*4*pi)) * radial_integral_even(i,j,1)  
    if ( cb1 == zero ) cycle
                                                                             
    do mu = -2, 2, 2                                                          
      do mp = -2, 2, 2                                                          
        do mlj = -lj, lj, 2
          xml = mlj / 2.d0
          fac2 = -(mu/2) * sqrt((xl*(xl+one)-xml*(xml+mu/2))/2.d0) &
                 + xml * kdelta(mu,0) 
          do mli = -li, li, 2
            do msi = -1, 1, 2 
              call ClebschGordan(li,1,ji,mli,msi,mji,cb2)    
              call ClebschGordan(lj,1,jj,mlj,msi,mjj,cb3)    
              call ClebschGordan(lj,2,li,mlj+mu,mp,mli,cb4)
              call ClebschGordan(lj,2,li,mlj,mp,mli,cb5)
              fac3 = fac2 * cb2 * cb3 * cb4    
              fac4 = cb2 * cb5 * sqrt(3.d0/4.d0)
              do msp = -4, 4, 2
                call ClebschGordan(2,2,4,mp,mu,msp,cb6)    
                angumome_Lgr(i,j,msp/2) = angumome_Lgr(i,j,msp/2) + fac3 * cb6
                do msj = -1, 1, 2
                  call ClebschGordan(lj,1,jj,mlj,msj,mjj,cb7)    
                  call ClebschGordan(1,2,1,msj,mu,msi,cb8)    
                  angumome_Sgr(i,j,msp/2) = angumome_Sgr(i,j,msp/2) + & 
                                             fac4 * cb6 * cb7 * cb8
                enddo !msj
              enddo !msp
            enddo !msi
          enddo !mli
        enddo !mlj
      enddo !mp
    enddo !mu
    angumome_Lgr(i,j,:) = angumome_Lgr(i,j,:) * fac1
    angumome_Sgr(i,j,:) = angumome_Sgr(i,j,:) * fac1
  enddo !j
enddo !i

!!! SPIN-ORBIT
allocate( angumome_so(HOsp_dim,HOsp_dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of angular momentum operators'

angumome_so = zzero

do i = 1, HOsp_dim
  ji = HOsp_2j(i)                                                                   
  li = HOsp_l(i) 
  xj = ji / 2.d0                                                 
  angumome_so(i,i) = xj*(xj+1) - li*(li+1) - 0.75d0
enddo

end subroutine set_angularmomentum

!------------------------------------------------------------------------------!
! subroutine calculate_angularmomentum                                         !
!                                                                              ! 
! Expectation values of the angular momentum operator and its square in        ! 
! cartesian coordinates.                                                       ! 
! < J_i > = Tr(j_i rho)                                                        ! 
! < J_i^2 > = Tr(j_i^2 rho) + Tr(j_i^T kappa10^dagger j_i kappa01)             ! 
!             - Tr([j_i rho]^2) + Tr(j_i rho)^2                                ! 
! < J^2 > = < J_x^2 > + < J_y^2 > + < J_z^2 >                                  ! 
! < ls > = Tr(ls rho)                                                          ! 
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      !
!        io = 0 computes < J_i >, < J_i^2 > with i = x,y,z                     !
!           = 1 computes < J_z >, < J_z^2 >                                    !
!        rhoLR,kappaLR,kappaRL = transition densities                          !
!                                                                              !
! Output: amji, amji2 = < J_i >, < J_i^2 > with i = x,y,z                      !
!         amsoi = < ls >  with i = p,n                                         !
!------------------------------------------------------------------------------!
subroutine calculate_angularmomentum(io,rhoLR,kappaLR,kappaRL,amjx,amjy,amjz, &
                                     amjx2,amjy2,amjz2,amsop,amson,ndim)

integer, intent(in) :: io, ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL
complex(r64), intent(out) :: amjx, amjy, amjz, amjx2, amjy2, amjz2, amsop, amson
integer :: i, j, imin
complex(r64) :: tr1, tr2, tr5, tr6, jtmp, j2tmp
complex(r64), dimension(ndim,ndim) :: Ji, Ji2, A1, A2, A3, A4, A5, A6

!!! Determines what components have to be calculated
if ( io == 0 ) then 
  imin = 1
elseif ( io == 1 ) then
  imin = 3
else
  write(uto,*) 'In calculate_angularmomentum: wrong argument io = ', io 
  stop 
endif

!!! Computes the expectation values of J_i and J_i^2
amjx  = zzero
amjy  = zzero
amjz  = zzero
amjx2 = zzero
amjy2 = zzero
amjz2 = zzero

do i = imin, 3

  if ( i == 1 ) then       ! X
    Ji  = angumome_Jx
    Ji2 = angumome_Jx2
  elseif ( i == 2 ) then   ! Y
    Ji  = angumome_Jy
    Ji2 = angumome_Jy2
  else                     ! Z
    Ji  = angumome_Jz
    Ji2 = angumome_Jz2
  endif

  call zgemm('n','n',ndim,ndim,ndim,zone, Ji,ndim,  rhoLR,ndim,zzero,A1,ndim)
  call zgemm('n','n',ndim,ndim,ndim,zone,Ji2,ndim,  rhoLR,ndim,zzero,A2,ndim)
  call zgemm('t','t',ndim,ndim,ndim,zone, Ji,ndim,kappaRL,ndim,zzero,A3,ndim)   
  call zgemm('n','n',ndim,ndim,ndim,zone, Ji,ndim,kappaLR,ndim,zzero,A4,ndim)   
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

  jtmp  = tr1 
  j2tmp = tr1**2 + tr2 - tr5 + tr6

  if ( i == 1 ) then
    amjx  = jtmp
    amjx2 = j2tmp
  elseif ( i == 2 ) then
    amjy  = jtmp
    amjy2 = j2tmp
  else
    amjz  = jtmp
    amjz2 = j2tmp
  endif

enddo 

!!! Computes the expectation value of the spin-orbit 
amsop = zzero
amson = zzero

call zgemm('n','n',ndim,ndim,ndim,zone,angumome_so,ndim,rhoLR,ndim,zzero,A1, &
           ndim)

do i = 1, ndim/2
  amsop = amsop + A1(i,i)
  amson = amson + A1(i+ndim/2,i+ndim/2)
enddo

end subroutine calculate_angularmomentum

END MODULE AngularMomentum
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
