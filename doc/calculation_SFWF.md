from https://bb.cgd.ucar.edu/cesm/threads/which-variables-to-include-in-ocean-surface-salt-budget.2385/:

##  FW formulation
SFWF - QFLUX/latent_heat_fusion/1.e4 = PREC_F + EVAP_F + ROFF_F + IOFF_F + MELT_F 
                                       + SALT_F*(sflux_factor/salinity_factor) 
                                       - QFLUX/latent_heat_fusion/1.e4    
                                       + (any weak or strong salinity restoring)

MELT_F represents any FW flux computed by the ice model and sent to the ocean (associated with ice melt/growth)

SALT_F represents any SALT flux computed by the ice model and sent to the ocean (associated with the fact that ice melt/growth requires that salt be added/removed from the ocean because ice has a constant salinity of ~4 psu!). This salt flux must of course be converted to appropriate units; in the above equation, it is converted from (kg SALT/m^2/s) to (kg FW/m^2/s)

The QFLUX term isn't included in SFWF because it's not considered a flux exchanged between model components -- it is the FW flux associated with frazil ice formation which occurs within the ocean model!

POP code in `forcing_coupled.F90`; specifically the line

```fortran
  if (sfc_layer_type == sfc_layer_varthick .and.   &
    .not. lfw_as_salt_flx) then
    !*** compute fresh water flux (cm/s)
    (...)
  else  ! convert fresh water to virtual salinity flux

  STF,             &! surface tracer fluxes at current timestep
! l. 801ff
  STF(:,:,2,iblock) = RCALCT(:,:,iblock)*(  &
                      (PREC_F(:,:,iblock) + EVAP_F(:,:,iblock) +  &
                      MELT_F(:,:,iblock) + ROFF_F(:,:,iblock) + IOFF_F(:,:,iblock))*salinity_factor + &
                      SALT_F(:,:,iblock)*sflux_factor)
! CALCT          ,&! flag=true if point is an ocean point
! RCALCT         ,&! real equiv of CALCT,U to use as more (./code/source/grid.F90)
```
in `forcing.F90`
```fortran
SFWF = STF/salinity_factor
```
salinity_factor = -0.00347
sflux_factor = 0.1 -- Convert Salt Flux to Salt Flux (in model units)

SFWF = (PREC_F + EVAP_F + MELT_F + ROFF_F + IOFF_F)*salinity_factor
     +  SALT_F(:,:,iblock)*sflux_factor

IOFF_F: Ice Runoff Flux from Coupler due to Land-Model Snow Capping, units='kg/m^2/s', grid_loc='2110'

##  virtual salt flux formulation:  see `POP_ConstantsMod.F90`