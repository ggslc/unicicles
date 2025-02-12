subroutine add_snowpacks(nsmax                         &
                        ,dzsnow                        &
                        ,ds                            &
                        ,sice                          & 
                        ,sliq                          &
                        ,tsnow                         &
                        ,rgrainl                       & 
                        ,snowdepth,nsnow               &
                        ,ds_in                         &
                        ,sice_in                       &
                        ,sliq_in                       &
                        ,tsnow_in                      &
                        ,rgrainl_in                    &
                        )
USE water_constants_mod, ONLY :                                              &
!  imported scalar parameters
 tm              ! Temperature at which fresh water freezes
!

USE c_perma, ONLY :                                               &
!  imported scalar parameters
 hcapi                                                            &
                  !  Specific heat capacity of ice (J/kg/K)
,hcapw           !  Specific heat capacity of water (J/kg/K)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nsmax

REAL, INTENT(IN) :: dzsnow(nsmax)

INTEGER, INTENT(INOUT) :: nsnow

REAL, INTENT(INOUT) ::                                              &
 ds(nsmax),ds_in(nsmax)                                          &
,sice(nsmax),sice_in(nsmax)                                          &
,sliq(nsmax),sliq_in(nsmax)                                          &
,tsnow(nsmax),tsnow_in(nsmax)                                        &
,rgrainl(nsmax),rgrainl_in(nsmax)                                    &
,snowdepth




! obviously mostly this is relayersnow. Most of its actually needed for
! adding mass in consistently, it seems, not taking it out!
REAL, PARAMETER :: thin_snow_limit = 1.0e-12
                  ! Maximum snow thickness (m) that is neglected
                  ! during relayering. All contributions
                  ! (mass, energy etc) from that snow are
                  ! neglected.

! Local scalars
INTEGER ::                                                        &
 i                                                                &
                       ! Land point index
,iznew                                                            &
                       ! layer index
,izz                                                              &
                       ! layer index
,k                                                                &
                       ! Tile point index
,n                                                                &
                       ! Snow layer index
,new                                                              &
                       ! layer index
,old                   ! layer index

REAL ::                                                           &
 csnow                                                    &
                       ! Areal heat capacity of layer (J/K/m2)
,oldremains                                                       &
,remains                                                       &
                       ! remaining depth in an old layer (m)
,wt                                                            &
                       ! weight given to a layer value
,denom

! Local arrays
INTEGER ::                                                        &
 nold

REAL ::                                                           &
 d0(1:nsmax)                                             &
                              ! Layer thicknesses before adjustment (m)
!                             ! D0(:,0) represents new snow if NSNOW>0,
!                             ! otherwise it is all snow.
,e(1:nsmax)                                                       &
                     ! Internal energy before adjustment (J/m2)
,newremains(nsmax)                                                &
                     ! available (unfilled) depth in new layer (m)
,r(1:nsmax)                                                       &
                     ! Grain size before adjustment (kg/m2)
,s(1:nsmax)                                                       &
                     ! Ice content before adjustment (kg/m2)
,w(1:nsmax)                                                       &
                     ! Liquid content before adjustment (kg/m2)
,u(nsmax)            ! Layer energy contents


nold = nsmax
d0 = ds+ds_in
!-----------------------------------------------------------------------
! Store previous snow layer energy and mass contents
!-----------------------------------------------------------------------
r(:)=0.
s(:)=0.
w(:)=0.
DO n=1,nold
!    csnow = sice(n)*hcapi + sliq(n)*hcapw
!    e(n) = csnow * ( tsnow(n) - tm )
  e(n) = (sice(n)*hcapi + sliq(n)*hcapw) * ( tsnow(n) - tm )    &
       + (sice_in(n)*hcapi + sliq_in(n)*hcapw) * ( tsnow_in(n) - tm )
  denom=1.
  if (ds(n) .gt. 0) r(n)=r(n)+rgrainl(n)
  if (ds_in(n) .gt. 0) r(n)=r(n)+rgrainl_in(n)
  if (ds(n) .gt. 0 .AND. ds_in(n) .gt. 0) denom=2.
  r(n)=r(n)/denom
  s(n) = sice(n)+sice_in(n)
  w(n) = sliq(n)+sliq_in(n)
END DO

!-----------------------------------------------------------------------
! Calculate snowdepth
!-----------------------------------------------------------------------
snowdepth = 0.
DO n=1,nsmax
  snowdepth = snowdepth + d0(n)
END DO


!-----------------------------------------------------------------------
! Divide snowpack into new layers
!-----------------------------------------------------------------------
nsnow = 0
ds(:) = 0.0

!   Only divide into layers if depth is >= a threshold.
! layersnow 1D
  remains = snowdepth

  DO n=1,nsmax
    ds(n) = dzsnow(n)
    remains = remains - dzsnow(n)
    IF ( remains <= dzsnow(n).OR. n == nsmax ) THEN
      ds(n) = ds(n) + remains
      remains=0
      EXIT
    END IF
  END DO

  nsnow = n

!!-----------------------------------------------------------------------
!! Initialise accumulations for new layer values.
!!-----------------------------------------------------------------------
  u(:) = 0.
  sice(:) = 0.
  sliq(:) = 0.
  rgrainl(:) = 0.

!!-----------------------------------------------------------------------
!! Set the state of the new layers.
!!-----------------------------------------------------------------------
!!    Initialise with all new layers empty.
  newremains(1:nsnow) = ds(1:nsnow)
!!     Start by filling top new layer.
  iznew = 1
!
!!     Loop over the old layers.
    DO old=1,nold

!!       All of this old layer remains to be reassigned to new layer(s).
      oldremains = d0(old)
!
!!       Point to first new layer with remaining space.
      izz = iznew
!
!!       Loop over new layers with remaining space.
      DO new=izz,nsnow
!
        IF ( oldremains > newremains(new) ) THEN
!!-----------------------------------------------------------------------
!! The remaining depth in the new layer will be exhausted by some or
!! all of the remaining depth from the old layer.
!!-----------------------------------------------------------------------
!
!!           Decrement old layer by the remaining space in new layer.
          oldremains = oldremains - newremains(new)
!
!!           Add properties from old layer to accumulation for new layer.
!!           Note that wt is <= 1 since here we have oldRemains>newRemain
!!           and oldRemains <= d0.
          IF ( d0(old) > thin_snow_limit ) THEN
            wt =  newremains(new) / d0(old)
            u(new) = u(new) + e(old) * wt
            sice(new) = sice(new) + s(old) * wt
            sliq(new) = sliq(new) + w(old) * wt
            rgrainl(new) = rgrainl(new) +     &
                             r(old) * newremains(new)
          END IF

!!           Update the pointer to the next new layer with space.
          izz = new + 1

        ELSE

!!-----------------------------------------------------------------------
!! The old layer will be exhausted by this increment.
!!-----------------------------------------------------------------------
!!           Decrement available space in the new layer.
          newremains(new) = newremains(new) - oldremains
!!           Add properties from old layer to accumulation for new layer.
          IF ( d0(old) > thin_snow_limit ) THEN
            wt = oldremains /  d0(old)
            u(new) = u(new) + e(old) * wt
            sice(new) = sice(new) + s(old) * wt
            sliq(new) = sliq(new) + w(old) * wt
            rgrainl(new) = rgrainl(new) + &
                            r(old) * oldremains
          END IF
!!           Proceed to the next old layer by exiting from the new layer 
          EXIT
        END IF
      END DO  !  new layers
!!       Update pointer to the next new layer with space.
      iznew = izz
    END DO  !  old layers

!!-----------------------------------------------------------------------
!! Diagnose layer temperatures and densities.
!-----------------------------------------------------------------------
    DO n=1,nsnow
      csnow = sice(n)*hcapi + sliq(n)*hcapw
      tsnow(n) = tm + u(n) / csnow
      rgrainl(n) = rgrainl(n) / ds(n)
    END DO

END
