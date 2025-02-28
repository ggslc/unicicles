SUBROUTINE  adjust_stdalone(                &
                         nx                 &
                        ,ny                 &
                        ,ntiles             &
                        ,nsmax              &
                        ,dzsnow             &
                        ,elev_start         &
                        ,elev_levels        &
                        ,tile_frac          &
                        ,tile_frac_old      &
                        ,nsnow              &
                        ,nonice_snowdepth   &
                        ,snow_tile          &
                        ,snowdepth          &
                        ,rho_snow_grnd      &
                        ,rgrain             &
                        ,sice               &
                        ,sliq               &
                        ,rgrainl            &
                        ,tsnow              &
                        ,rho_snow           &
                        ,ds                 &
                        ,smb_from_nisnow    &
                         )

USE water_constants_mod, ONLY :             &
 tm              ! Temperature at which fresh water freezes

IMPLICIT NONE

REAL, PARAMETER  :: icetile_tol=1e-4
REAL, PARAMETER  :: r0 = 50
REAL, PARAMETER  :: rho_snow_const=250.
REAL, PARAMETER  :: glimmer_ice_density=910.
!REAL, PARAMETER  :: ice_reset_snowmass=glimmer_ice_density*100.

INTEGER,INTENT(IN) :: nx,ny,ntiles,nsmax
INTEGER,INTENT(IN) :: elev_start,elev_levels

INTEGER,INTENT(INOUT) :: nsnow(nx,ny,ntiles)

REAL,INTENT(IN) :: dzsnow(nsmax)

REAL,INTENT(INOUT) ::                                                           &
 tile_frac(nx,ny,ntiles)                                       &
,tile_frac_old(nx,ny,ntiles)                                   &
,nonice_snowdepth(nx,ny,ntiles)                                &
,snow_tile(nx,ny,ntiles)                                       &
,snowdepth(nx,ny,ntiles)                                       &
,rho_snow_grnd(nx,ny,ntiles)                                   &
,rgrain(nx,ny,ntiles)                                          &
,sice(nx,ny,ntiles,nsmax)                                      &
,sliq(nx,ny,ntiles,nsmax)                                      &
,rgrainl(nx,ny,ntiles,nsmax)                                   &
,tsnow(nx,ny,ntiles,nsmax)                                     &
,rho_snow(nx,ny,ntiles,nsmax)                                  &
,ds(nx,ny,ntiles,nsmax)                                        &
,smb_from_nisnow(nx,ny,ntiles)
!

!LOCAL
INTEGER ::                                                        &
 i,j,l,e                                                          &
,n,n_ni,n_i                                                    &
,ns                                                               &
,elev_ni_start, elev_ni_end                                       &
,elev_i_start, elev_i_end                                         &
,init_i, init_ni                                                  &
,grow_i, grow_ni                                                  &
,nsnow_work

REAL ::                                                           &
 diff                                                             &
,new_elev_temp                                                    &
,new_elev_q                                                       &
,area_factor                                                      &
,telev_nrg_old,telev_nrg

REAL ::                                                           &
 ds_work(nsmax)                                                   &
,sice_work(nsmax)                                                 &
,sliq_work(nsmax)                                                 &
,tsnow_work(nsmax)                                                &
,rgrainl_work(nsmax)                                              &
,snowdepth_work                                                   &
,snow_tile_work                                                   &
,rho_snow_grnd_work


elev_i_start = elev_start
elev_i_end   = elev_i_start+elev_levels-1

elev_ni_start    = elev_start+elev_levels
elev_ni_end      = elev_ni_start+elev_levels-1

do j=1,ny
do i=1,nx

  diff=0.
  smb_from_nisnow(i,j,:)=0.

  do n=elev_i_start,elev_ni_end
    diff=diff+abs(tile_frac_old(i,j,n)-tile_frac(i,j,n))
  enddo 

!      tile_fractions have changed in this gridpoint
!      loop over elevations, do the two types within by hand for each elev
    do n_i=elev_i_start,elev_i_end
      n_ni=n_i+elev_levels

      !write statements only really useful for a 1x1 grid!
      !write(6,*)"change in fracs, tile:",n_ni, n_i
      !write(6,*)"old               ",tile_frac_old(i,j,n_ni),tile_frac_old(i,j,n_i)
      !write(6,*)"new               ",tile_frac(i,j,n_ni),tile_frac(i,j,n_i)
      !write(6,*)"ni_snowdepth      ",nonice_snowdepth(i,j,n_ni)

      !write(6,*)"b4: ni_snow",tile_frac_old(i,j,n_ni) &
      !                       ,snow_tile(i,j,n_ni)     &
      !                       ,nsnow(i,j,n_ni) &
      !                       ,snowdepth(i,j,n_ni) &
      !                       ,rho_snow_grnd(i,j,n_ni) &
      !                       ,tile_frac_old(i,j,n_ni)*snow_tile(i,j,n_ni)

      !write(6,*)"    i_snow",tile_frac_old(i,j,n_i) &
      !                       ,snow_tile(i,j,n_i)-ice_reset_snowmass &
      ! ,tile_frac_old(i,j,n_i)*(snow_tile(i,j,n_i)-ice_reset_snowmass)

  if (abs(diff) .gt. icetile_tol) then
      grow_ni=0
      init_ni=0
      grow_i=0
      init_i=0
      if (tile_frac(i,j,n_ni) .gt. tile_frac_old(i,j,n_ni)) grow_ni= 1
      if (tile_frac(i,j,n_ni) .lt. tile_frac_old(i,j,n_ni)) grow_ni=-1
      if (tile_frac(i,j,n_i)  .gt. tile_frac_old(i,j,n_i))  grow_i=  1
      if (tile_frac(i,j,n_i)  .lt. tile_frac_old(i,j,n_i))  grow_i= -1

      if( grow_ni .eq. 1 .AND. tile_frac_old(i,j,n_ni).lt. icetile_tol) init_ni=1
      if( grow_i  .eq. 1 .AND. tile_frac_old(i,j,n_i) .lt. icetile_tol) init_i=1

      !write(6,*)"         grow_i,     grow_ni,    init_i,    init_ni"
      !write(6,*)grow_i,grow_ni,init_i,init_ni

      if (init_ni .gt. 0) then
        snow_tile(i,j,n_ni)=0.
        nsnow(i,j,n_ni)    =0.
        snowdepth(i,j,n_ni)=0.
        rho_snow_grnd(i,j,n_ni)=rho_snow_const
        rgrain(i,j,n_ni)=r0
      
        ds(i,j,n_ni,:)=0.
        sice(i,j,n_ni,:)=0.
        sliq(i,j,n_ni,:)=0.
        tsnow(i,j,n_ni,:)=tm
        rgrainl(i,j,n_ni,:)=r0
        rho_snow(i,j,n_ni,:)=rho_snow_const
      end if


      !changing fractions. Start by rearranging the ni_snowpack to the new fractions
      !newly initialised points need no rearranging

        !if ni_frac growing, spread it out. Quantities /m^2 need to change
        !Area-integrated quantities are unchanged as the area has increased
        !factor for change to new area=tile_frac/tile_frac_old
        if (grow_ni .eq. 1 .AND. init_ni .lt. 1) then
          area_factor=tile_frac(i,j,n_ni)/tile_frac_old(i,j,n_ni)
          if (nsnow(i,j,n_ni) .gt. 0) then
            call scale_snowpack(nsmax,dzsnow                                  &
                               ,ds(i,j,n_ni,:),sice(i,j,n_ni,:),sliq(i,j,n_ni,:)    &
                               ,tsnow(i,j,n_ni,:),rgrainl(i,j,n_ni,:)       &
                               ,snowdepth(i,j,n_ni),nsnow(i,j,n_ni)         &
                               ,area_factor )

        !  Can the above reduce a pack to 0?
          else
        ! only tile mass/m^2 needs changing
            snow_tile(i,j,n_ni)=snow_tile(i,j,n_ni)/area_factor
          end if !what sort of snowpack are we changing?

        end if !pre-existing non_ice tile growing


        !mass loss from ni_snowpacks
        if (grow_ni .eq. -1 .AND. snow_tile(i,j,n_ni) .gt. dzsnow(1)  &
                            .AND. tile_frac(i,j,n_i)  .gt. icetile_tol) then
        !if ni_frac shrinking, keep ni_snowpack/m^2 the same, but add the 
        !lost mass implied by the change of area to the ice_snowpack if there is
        !one
          area_factor=tile_frac(i,j,n_i)/(tile_frac_old(i,j,n_ni)-tile_frac(i,j,n_ni))

          smb_from_nisnow(i,j,n_i)=smb_from_nisnow(i,j,n_i)+snow_tile(i,j,n_ni)/area_factor

        end if !redistribution of ni_snow mass

        if (tile_frac(i,j,n_ni) .lt. icetile_tol) then ! if the new ni frac has -> 0
          snow_tile(i,j,n_ni)=0.                       ! we shouldn't have any snow
          !other terms handled by general regularisation below
        end if


        !Now, add any ni_snow changes from Glimmer's internal SMB-creation or
        !destruction of ice points. At this point, nonice_snowdepth is an array of
        !ice-density depth anomalies applicable to a whole UM gridbox (or just
        !the nonice fraction thereof?)

  endif !is there a tile change in this gridbox?

        if (tile_frac(i,j,n_ni) .gt. icetile_tol) then !it should be if there's a
                                                       !Glimmer anomaly to
                                                       !apply, but not
                                                       !guaranteed under
                                                       !BISICLES
          !this assumes Glimmer has derived a depth for its anomaly by spreading mass across the
          !whole area of UM gridbox as it interpolated/binned. Not sure this is
          !actally true - investigate.
          if (nonice_snowdepth(i,j,n_ni) .gt. dzsnow(1)) then
  
            !the ice-sheet is adding mass. Use the ice_snowpack for the properties. If we
            !didn't have multilayer snow before, we probably will now
            if (nsnow(i,j,n_ni) .eq. 0) then
             !need to init a multilayer pack
             ds(i,j,n_ni,:)=0.
             ds(i,j,n_ni,1)=dzsnow(1)
             sice(i,j,n_ni,:)=0.
             sice(i,j,n_ni,1)=snow_tile(i,j,n_ni)
             sliq(i,j,n_ni,:)=0.
             tsnow(i,j,n_ni,:)=tm
             rgrainl(i,j,n_ni,:)=r0
             nsnow(i,j,n_ni)=1
             snowdepth(i,j,n_ni)=dzsnow(1)
            endif

            !there really *should* be valid icepack properties at this elev to use...
            ds_work(:)=ds(i,j,n_i,:)
            sice_work(:)=sice(i,j,n_i,:)
            sliq_work(:)=sliq(i,j,n_i,:)
            tsnow_work(:)=tsnow(i,j,n_i,:)
            rgrainl_work(:)=rgrainl(i,j,n_i,:)
            nsnow_work=nsnow(i,j,n_i)
            snowdepth_work=snowdepth(i,j,n_i)

            !Make a snowpack of ni_snowdepth*gbm_area mass with the ice_tile
            !snowpack/m^2 properties on the ni_tile area so we can add it
            !to what's already on the ni_tile - note the area_factor assumption 
            !that this was mass in this elevation (used to be whole gridbox)

            !scale icetile snowpack to match the mass we need to add, on the
            !nonice fraction
            area_factor=sum(sice_work(:)+sliq_work(:)) &
                       /(nonice_snowdepth(i,j,n_ni)*glimmer_ice_density)  &
                       *tile_frac(i,j,n_ni)/(tile_frac(i,j,n_ni)+tile_frac(i,j,n_i))

            call scale_snowpack(nsmax,dzsnow,ds_work(:),sice_work(:),sliq_work(:)    &
                                     ,tsnow_work(:),rgrainl_work(:)             &
                                     ,snowdepth_work,nsnow_work                 &
                                     ,area_factor )

            call add_snowpacks(nsmax                            &
                          ,dzsnow &
                          ,ds(i,j,n_ni,:)                         &
                          ,sice(i,j,n_ni,:)                       &
                          ,sliq(i,j,n_ni,:)                       &
                          ,tsnow(i,j,n_ni,:)                      &
                          ,rgrainl(i,j,n_ni,:)                    &
                          ,snowdepth(i,j,n_ni)                    &
                          ,nsnow(i,j,n_ni)                        &
                          ,ds_work(:),sice_work(:),sliq_work(:) &
                          ,tsnow_work(:),rgrainl_work(:)        &
                          )
          endif
          if (nonice_snowdepth(i,j,n_ni) .lt. -1.*dzsnow(1)) then

            area_factor= tile_frac(i,j,n_ni)/(tile_frac(i,j,n_ni)+tile_frac(i,j,n_i))
            !we're subtracting mass. This is difficult in a layer-by-layer sense
            !so just pull (ice-density) mass out of the base. Energy is hard too,
            !but we can try to work out the non-con implications at least
            !
            !initial check. Does this exhaust our snowpack, or at least change
            !it to 0-layer?
            snowdepth(i,j,n_ni)=snowdepth(i,j,n_ni) + (nonice_snowdepth(i,j,n_ni)/area_factor)
            snow_tile(i,j,n_ni)=snow_tile(i,j,n_ni) + (nonice_snowdepth(i,j,n_ni)/area_factor)*glimmer_ice_density

            if (snowdepth(i,j,n_ni) .lt. dzsnow(1) .OR. (snow_tile(i,j,n_ni).lt.dzsnow(1)*glimmer_ice_density) ) then
            !nsnow=0 now if it wasn't before
              sice(i,j,n_ni,:)=0.
              sliq(i,j,n_ni,:)=0.
              ds(i,j,n_ni,:)=0
              rgrainl(i,j,n_ni,:)=r0
              tsnow(i,j,n_ni,:)=tm
              snow_tile(i,j,n_ni)=max(                                     &
                snow_tile(i,j,n_ni)+(nonice_snowdepth(i,j,n_ni)/area_factor) &
                *GLIMMER_ICE_DENSITY                                     &
                ,0.)
              snowdepth(i,j,n_ni)=0.
              nsnow(i,j,n_ni)=0.
            else
            !we've still got multi-layer snow even after the subtraction

               call subtract_snowpacks(nsmax,dzsnow                       &
                              ,ds(i,j,n_ni,:),sice(i,j,n_ni,:),sliq(i,j,n_ni,:) &
                              ,tsnow(i,j,n_ni,:),rgrainl(i,j,n_ni,:)          &
                              ,snowdepth(i,j,n_ni),nsnow(i,j,n_ni)            &
                              ,nonice_snowdepth(i,j,n_ni)/area_factor       &
                               *glimmer_ice_density,glimmer_ice_density   &
                              )
            endif !has the subtraction left us with multi-snow, or not?

          endif !is there an icesheet ni_snowpack anomaly to deal with
        endif !is there actually tile_frac to apply it to?!

        !have a simple routine to regularise all the snowpack vars from 
        !ds, sice, sliq, tsnow and rgrainl alone, to make sure diagnostics
        !etc are happy? snowdepth and nsnow should already be OK from routines

        if (nsnow(i,j,n_ni) .gt. 0 .and. ds(i,j,n_ni,1).ge.dzsnow(1)) then
          !multilayer. Make sure everything is consistent
          snow_tile(i,j,n_ni)=0.
          snowdepth(i,j,n_ni)=0.
          do ns=1,nsmax
            if (ns .le. nsnow(i,j,n_ni)) then
              rho_snow(i,j,n_ni,ns)=(sice(i,j,n_ni,ns)+sliq(i,j,n_ni,ns))/ds(i,j,n_ni,ns)
              snow_tile(i,j,n_ni)=snow_tile(i,j,n_ni)+sice(i,j,n_ni,ns)+sliq(i,j,n_ni,ns)
              snowdepth(i,j,n_ni)=snowdepth(i,j,n_ni)+ds(i,j,n_ni,ns)
            else
              rgrainl(i,j,n_ni,ns) = r0
              rho_snow(i,j,n_ni,ns)=0.
              sice(i,j,n_ni,ns)=0.
              sliq(i,j,n_ni,ns)=0.
              tsnow(i,j,n_ni,ns)=tm
            endif
          end do
          rgrain(i,j,n_ni) = rgrainl(i,j,n_ni,1)
          rho_snow_grnd(i,j,n_ni) = snow_tile(i,j,n_ni)/snowdepth(i,j,n_ni)
        else
          !0layer - only snow_tile is valid really
          nsnow(i,j,n_ni)=0.
          rgrain(i,j,n_ni)=r0
          rho_snow_grnd(i,j,n_ni)=rho_snow_const
          snowdepth(i,j,n_ni)=0
          rgrainl(i,j,n_ni,:) = r0
          rho_snow(i,j,n_ni,:)=0.
          sice(i,j,n_ni,:)=0.
          sliq(i,j,n_ni,:)=0.
          tsnow(i,j,n_ni,:)=tm
        end if

        !REPEAT FOR ICE TILE
        if (nsnow(i,j,n_i) .gt. 0 .and. ds(i,j,n_i,1).ge.dzsnow(1)) then
          !multilayer. Make sure everything is consistent
          snow_tile(i,j,n_i)=0.
          snowdepth(i,j,n_i)=0.
          do ns=1,nsmax
            if (ns .le. nsnow(i,j,n_i)) then
              rho_snow(i,j,n_i,ns)=(sice(i,j,n_i,ns)+sliq(i,j,n_i,ns))/ds(i,j,n_i,ns)
              snow_tile(i,j,n_i)=snow_tile(i,j,n_i)+sice(i,j,n_i,ns)+sliq(i,j,n_i,ns)
              snowdepth(i,j,n_i)=snowdepth(i,j,n_i)+ds(i,j,n_i,ns)
            else
              rgrainl(i,j,n_i,ns) = r0
              rho_snow(i,j,n_i,ns)=0.
              sice(i,j,n_i,ns)=0.
              sliq(i,j,n_i,ns)=0.
              tsnow(i,j,n_i,ns)=tm
            endif
          end do
          rgrain(i,j,n_i) = rgrainl(i,j,n_i,1)
          rho_snow_grnd(i,j,n_i) = snow_tile(i,j,n_i)/snowdepth(i,j,n_i)
        end if

    end do !loop over elevations


enddo !gridpoints
enddo !gridpoints

END SUBROUTINE
