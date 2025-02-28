subroutine subtract_snowpacks(nsmax,dzsnow                        &
                             ,ds, sice, sliq                      &
                             ,tsnow,rgrainl                       &
                             ,snowdepth,nsnow                     &
                             ,sice0,glimmer_ice_density           &
                             )

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                            &
 nsmax

INTEGER, INTENT(INOUT) ::                                         &
 nsnow

REAL, INTENT(IN) :: dzsnow(nsmax)

REAL, INTENT(IN) ::                                               &
 sice0                                                            & 
,glimmer_ice_density

REAL, INTENT(INOUT) ::                                            &
 snowdepth

REAL, INTENT(INOUT) ::                                            &
 ds(nsmax)                                                        &
,sice(nsmax)                                                      &
,sliq(nsmax)                                                      &
,tsnow(nsmax)                                                     &
,rgrainl(nsmax)

! Local scalars
INTEGER ::                                                        &
 n                                                                &
                     ! Snow layer index
,nold 

REAL ::                                                           &
 d0(1:nsmax)                                                      &
                   ! Layer thicknesses before adjustment (m)
,s(1:nsmax)                                                       &
                   ! Ice content before adjustment (kg/m2)
,w(1:nsmax)                                                       &
                   ! Liquid content before adjustment (kg/m2)
,remains                                                          &
,mass_n

nold = nsnow
d0(1:nsnow) = ds(1:nsnow)
s(1:nsnow) = sice(1:nsnow)
w(1:nsnow) = sliq(1:nsnow)

d0(nsnow) = d0(nsnow) + sice0 / glimmer_ice_density
!much easier to treat mass loss basically ice/water agnostic in this framework
!some energy violation from implied phase changes
s(nsnow)=s(nsnow)+sice0


!-----------------------------------------------------------------------
! Calculate snowdepth
!-----------------------------------------------------------------------
snowdepth = 0.
DO n=1,nsnow
  snowdepth = snowdepth + d0(n)
END DO

!-----------------------------------------------------------------------
! Divide snowpack into new layers
!-----------------------------------------------------------------------
nsnow = 0
ds(:) = 0.0

! Only divide into layers if depth is >= a threshold.
! layersnow 1D
remains = snowdepth

DO n=1,nsmax
  ds(n) = dzsnow(n)
  remains = remains - dzsnow(n)
  IF ( remains <= dzsnow(n) .OR. n == nsmax) THEN
    ds(n) = ds(n) + remains
    remains=0
    EXIT
  END IF
END DO

nsnow = n

! subtract mass from layers until you have something left
remains=0.
do n=nold,1,-1
  mass_n=s(n)+w(n)
  if (mass_n+remains .gt. 0 .and. ds(n).gt.0) then
      s(n)=s(n)+s(n)*remains/(mass_n)
      w(n)=w(n)+w(n)*remains/(mass_n)
      EXIT
  else
    s(n)=0.
    w(n)=0.
    remains=remains+mass_n
  endif
end do
sice(1:) = s(1:)
sliq(1:) = w(1:)

END
