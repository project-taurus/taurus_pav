!==============================================================================!
! MODULE Multipoles                                                            !
!                                                                              !
! This module contains the variables and routines related to multipole operat- !
! ors (up to Q44 for now).                                                     !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_multipoles                                                  !
! - subroutine calculate_multipoles                                            !
!==============================================================================!
MODULE Multipoles

use Basis
use AngularMomentum

implicit none 
public

complex(r64), dimension(:,:,:), allocatable :: multipole_Q1p, & ! Q1 protons 
                                               multipole_Q1n, & ! Q1 neutrons
                                               multipole_Q2p, & ! Q2 protons 
                                               multipole_Q2n, & ! Q2 neutrons
                                               multipole_Q3p, & ! Q3 protons 
                                               multipole_Q3n, & ! Q3 neutrons
                                               multipole_M1p, & ! M1 protons
                                               multipole_M1n, & ! M1 neutrons
                                               multipole_M2p, & ! M2 protons
                                               multipole_M2n    ! M2 neutrons

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_multipoles                                                    !
!                                                                              !
! Defines the HO matrix elements of the multipoles operators. The following    !
! formula is used:                                                             !
! <a|Q_lambda,m|b> = b^lambda * R^lambda_ab * sqrt((2lb+1)*(2lambda+1)/(4pi *  !
!                    (2la+1))) * CG(lb,0,lambda,0,la,0) *                      !
!                    sum_ms,mla,mlb CG(la,mla,1/2,ms,ja,ma) *                  !
!                       CG(lb,mlb,1/2,ms,jb,mb) * CG(lb,mlb,lambda,m,la,mla)   !
! where b is the oscillator length, R is the radial part, and CGs are Clebsch- !
! Gordan coefficients.                                                         !
!                                                                              !
! Input: opt_oper = option to read the matrix elements from a file             !
!------------------------------------------------------------------------------!
subroutine set_multipoles(opt_oper)

integer, intent(in) :: opt_oper
integer :: a, b, ja, jb, la, lb, mja, mjb, mla, mlb, ms, lambda, mu, ialloc=0, &
           utn
real(r64) :: fac, cb1, cb2, cb3
complex(r64) :: zsum
complex(r64), dimension(:,:,:), allocatable :: Qp, Qn
character(27) :: namefile
logical, dimension(3) :: is_exist

allocate( multipole_Q1p(HOsp_dim/2,HOsp_dim/2,-1:1), & 
          multipole_Q1n(HOsp_dim/2,HOsp_dim/2,-1:1), & 
          multipole_Q2p(HOsp_dim/2,HOsp_dim/2,-2:2), & 
          multipole_Q2n(HOsp_dim/2,HOsp_dim/2,-2:2), & 
          multipole_Q3p(HOsp_dim/2,HOsp_dim/2,-3:3), & 
          multipole_Q3n(HOsp_dim/2,HOsp_dim/2,-3:3), & 
          multipole_M1p(HOsp_dim/2,HOsp_dim/2,-1:1), & 
          multipole_M1n(HOsp_dim/2,HOsp_dim/2,-1:1), & 
          multipole_M2p(HOsp_dim/2,HOsp_dim/2,-2:2), & 
          multipole_M2n(HOsp_dim/2,HOsp_dim/2,-2:2), & 
          stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of mulitpole operators'

multipole_Q1p = zzero
multipole_Q1n = zzero
multipole_Q2p = zzero
multipole_Q2n = zzero
multipole_Q3p = zzero
multipole_Q3n = zzero
multipole_M1p = zzero
multipole_M1n = zzero
multipole_M2p = zzero
multipole_M2n = zzero

!!! 
!!! Electric moments
!!! 

!!! Reads from file if opt_oper = 1
is_exist = .false.

if ( opt_oper == 1 ) then 
  allocate( Qp(HOsp_dim/2,HOsp_dim/2,-3:3), Qn(HOsp_dim/2,HOsp_dim/2,-3:3),  & 
            stat=ialloc )
  if ( ialloc /= 0 ) stop 'Error during allocation of mulitpole operators'

  do lambda = 1, 3
    if ( lambda == 1 ) then
      namefile = "matelem_multipole_Q1_1b.bin"
    elseif ( lambda == 2 ) then
      namefile = "matelem_multipole_Q2_1b.bin"
    else
      namefile = "matelem_multipole_Q3_1b.bin"
    endif
  
    inquire(file=namefile, exist=is_exist(lambda))
 
    if ( is_exist(lambda) ) then
      open(newunit=utn, file=namefile, status='old', action='read', &
           form='unformatted')
 
      do mu = -lambda, lambda
        do b = 1, HOsp_dim/2
          do a = 1, HOsp_dim/2
            read(utn) Qp(a,b,mu), Qn(a,b,mu) 
          enddo
        enddo
      enddo

      if ( lambda == 1 ) then
        multipole_Q1p(:,:,-lambda:lambda) = Qp(:,:,-lambda:lambda)
        multipole_Q1n(:,:,-lambda:lambda) = Qn(:,:,-lambda:lambda)
      elseif ( lambda == 2 ) then
        multipole_Q2p(:,:,-lambda:lambda) = Qp(:,:,-lambda:lambda)
        multipole_Q2n(:,:,-lambda:lambda) = Qn(:,:,-lambda:lambda)
      else
        multipole_Q3p(:,:,-lambda:lambda) = Qp(:,:,-lambda:lambda)
        multipole_Q3n(:,:,-lambda:lambda) = Qn(:,:,-lambda:lambda)
      endif

      close(utn, status="keep")
    endif ! is_exist
  enddo ! lambda
  deallocate( Qp, Qn )
endif ! opt_oper


!!! Standard definition
do a = 1, HOsp_dim/2
  ja = HOsp_2j(a)
  la = 2*HOsp_l(a)
  mja = HOsp_2mj(a)
  do b = 1, HOsp_dim/2
    jb = HOsp_2j(b)
    lb = 2*HOsp_l(b)
    mjb = HOsp_2mj(b)
    do lambda = 2, 6, 2
      if ( is_exist(lambda/2) ) cycle
      if ( mod(la+lb+lambda,4) /= 0 ) cycle
      call ClebschGordan(lb,lambda,la,0,0,0,cb1)
      fac = cb1 * radial_integral_even(a,b,lambda/2) * &
            sqrt( ((lb+1)*(lambda+1)) / (4*pi*(la+1)) )
      do mu = -lambda, lambda, 2
        zsum = zzero
        do ms = -1, 1, 2    
          do mla = -la, la, 2    
            do mlb = -lb, lb, 2    
              call ClebschGordan(la,1,ja,mla,ms,mja,cb1)
              call ClebschGordan(lb,1,jb,mlb,ms,mjb,cb2)
              call ClebschGordan(lb,lambda,la,mlb,mu,mla,cb3)
              zsum = zsum + cb1 * cb2 * cb3
            enddo
          enddo
        enddo
        zsum = zsum * fac
        if ( lambda == 2 ) then
          multipole_Q1p(a,b,mu/2) = zsum
          multipole_Q1n(a,b,mu/2) = zsum
        elseif ( lambda == 4 ) then
          multipole_Q2p(a,b,mu/2) = zsum
          multipole_Q2n(a,b,mu/2) = zsum
        elseif ( lambda == 6 ) then
          multipole_Q3p(a,b,mu/2) = zsum
          multipole_Q3n(a,b,mu/2) = zsum
        endif
      enddo
    enddo 
  enddo
enddo

!!! 
!!! Magnetic moments
!!! 

!!! Reads from file if opt_oper = 1
is_exist = .false.

if ( opt_oper == 1 ) then 
  allocate( Qp(HOsp_dim/2,HOsp_dim/2,-3:3), Qn(HOsp_dim/2,HOsp_dim/2,-3:3),  & 
            stat=ialloc )
  if ( ialloc /= 0 ) stop 'Error during allocation of mulitpole operators'

  do lambda = 1, 2
    if ( lambda == 1 ) then
      namefile = "matelem_multipole_M1_1b.bin"
    else
      namefile = "matelem_multipole_M2_1b.bin"
    endif
  
    inquire(file=namefile, exist=is_exist(lambda))
 
    if ( is_exist(lambda) ) then
      open(newunit=utn, file=namefile, status='old', action='read', &
           form='unformatted')
 
      do mu = -lambda, lambda
        do b = 1, HOsp_dim/2
          do a = 1, HOsp_dim/2
            read(utn) Qp(a,b,mu), Qn(a,b,mu) 
          enddo
        enddo
      enddo

      if ( lambda == 1 ) then
        multipole_M1p(:,:,-lambda:lambda) = Qp(:,:,-lambda:lambda)
        multipole_M1n(:,:,-lambda:lambda) = Qn(:,:,-lambda:lambda)
      else
        multipole_M2p(:,:,-lambda:lambda) = Qp(:,:,-lambda:lambda)
        multipole_M2n(:,:,-lambda:lambda) = Qn(:,:,-lambda:lambda)
      endif

      close(utn, status="keep")
    endif ! is_exist
  enddo ! lambda
  deallocate( Qp, Qn )
endif ! opt_oper

!!! M1    
if ( .not. is_exist(1) ) then                     
  do mu = -1, 1
    do b = 1, HOsp_dim/2
      do a = 1, HOsp_dim/2
        multipole_M1p(a,b,mu) = sqrt(3/(4*pi)) * &
                                ( gyro_glp * angumome_L(a,b,mu) + &
                                  gyro_gsp * angumome_S(a,b,mu) )
        multipole_M1n(a,b,mu) = sqrt(3/(4*pi)) * &
                                ( gyro_gln * angumome_L(a,b,mu) + &
                                  gyro_gsn * angumome_S(a,b,mu) )
      enddo
    enddo
  enddo
endif

!!! M2                         
if ( .not. is_exist(2) ) then                     
  do mu = -2, 2
    do b = 1, HOsp_dim/2
      do a = 1, HOsp_dim/2
        multipole_M2p(a,b,mu) = sqrt(10.d0) * &
                                ( gyro_glp * angumome_Lgr(a,b,mu) * 2.d0/3 + &
                                  gyro_gsp * angumome_Sgr(a,b,mu) )
        multipole_M2n(a,b,mu) = sqrt(10.d0) * & 
                                ( gyro_gln * angumome_Lgr(a,b,mu) * 2.d0/3 + &
                                  gyro_gsn * angumome_Sgr(a,b,mu))
      enddo
    enddo
  enddo
endif

end subroutine set_multipoles 

!------------------------------------------------------------------------------!
! subroutine calculate_multipoles                                              !
!                                                                              ! 
! Computes the expectation values of the electromagnetic multipole operators.  ! 
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      !
!        rhoLR = transition densities                                          !
!                                                                              ! 
! Output: QX(mu) = < Q_X,mu > (protons and neutrons separetely)                !
!         MX(mu) = < M_X,mu > (nucleons)                                       !
!------------------------------------------------------------------------------!
subroutine calculate_multipoles(rhoLR,Q1p,Q1n,Q2p,Q2n,Q3p,Q3n,M1a,M2a,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR
complex(r64), dimension(-1:1), intent(out) :: Q1p, Q1n, M1a
complex(r64), dimension(-2:2), intent(out) :: Q2p, Q2n, M2a
complex(r64), dimension(-3:3), intent(out) :: Q3p, Q3n
integer :: hdim, i, j, lambda, mu
complex(r64) :: tr1, tr2
complex(r64), dimension(ndim/2,ndim/2) :: A1, A2, A3, A4, rho_p, rho_n

hdim = ndim/2

do j = 1, hdim
  do i = 1, hdim 
    rho_p(i,j)  = rhoLR(i,j)
    rho_n(i,j)  = rhoLR(i+hdim,j+hdim)
  enddo
enddo

!!! Electric moments               
do lambda = 1, 3 
  do mu = -lambda, lambda

    if ( lambda == 1 ) then
      A1(:,:) = multipole_Q1p(:,:,mu)
      A2(:,:) = multipole_Q1n(:,:,mu)
    elseif ( lambda == 2 ) then
      A1(:,:) = multipole_Q2p(:,:,mu)
      A2(:,:) = multipole_Q2n(:,:,mu)
    else
      A1(:,:) = multipole_Q3p(:,:,mu)
      A2(:,:) = multipole_Q3n(:,:,mu)
    endif

    call zgemm('n','n',hdim,hdim,hdim,zone,A1,hdim,rho_p,hdim,zzero,A3,hdim)
    call zgemm('n','n',hdim,hdim,hdim,zone,A2,hdim,rho_n,hdim,zzero,A4,hdim)

    tr1 = zzero
    tr2 = zzero
    do i = 1, hdim   
      tr1 = tr1 + A3(i,i)
      tr2 = tr2 + A4(i,i)
    enddo 

    if ( lambda == 1 ) then
      Q1p(mu) = tr1
      Q1n(mu) = tr2
    elseif ( lambda == 2 ) then
      Q2p(mu) = tr1
      Q2n(mu) = tr2
    else
      Q3p(mu) = tr1
      Q3n(mu) = tr2
    endif
  enddo
enddo

!!! Magnetic moments
! I have defined the magnetic moments using the spherical operators rather than
! the ladder ones (e.g. L+).
! Here are the relations between the two
! L(+) = -1/sqrt(2) L+
! L(-) =  1/sqrt(2) L-

do lambda = 1, 2
  do mu = -lambda, lambda
    if ( lambda == 1 ) then
      A1(:,:) = multipole_M1p(:,:,mu)
      A2(:,:) = multipole_M1n(:,:,mu)
    else
      A1(:,:) = multipole_M2p(:,:,mu)
      A2(:,:) = multipole_M2n(:,:,mu)
    endif

    call zgemm('n','n',hdim,hdim,hdim,zone,A1,hdim,rho_p,hdim,zzero,A3,hdim)
    call zgemm('n','n',hdim,hdim,hdim,zone,A2,hdim,rho_n,hdim,zzero,A4,hdim)
    
    tr1 = zzero
    tr2 = zzero
    do i = 1, hdim   
      tr1 = tr1 + A3(i,i)
      tr2 = tr2 + A4(i,i)
    enddo 
    
    if ( lambda == 1 ) then
      M1a(mu) = tr1 + tr2
    else
      M2a(mu) = tr1 + tr2
    endif
  enddo 
enddo 

end subroutine calculate_multipoles

END MODULE Multipoles
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
