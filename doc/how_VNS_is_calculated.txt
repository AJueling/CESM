# ----------------------------------------------------------------------
# How is VNS exactly calculated
# code from `advection.F90`
# ----------------------------------------------------------------------

   real (POP_r8), dimension(nx_block,ny_block) :: & 
     UTE,UTW,VTN,VTS,  &! tracer flux velocities across E,W,N,S faces

   integer (POP_i4), parameter :: &
      tadvect_centered = 1,  &! centered leap-frog tracer advection


   real (POP_r8), dimension(:,:,:,:), allocatable :: &
      VTN_to_VVEL_N        ! converts VTN to VVEL_E


   call define_tavg_field(tavg_VN_TRACER(2),'VNS',3,                   &
                          long_name='Salt Flux in grid-y direction',   &
                          units='gram/gram/s', grid_loc='3121')

# ----------------------------------------------------------------------
# VTN
# ----------------------------------------------------------------------

   if (k == 1 .or. .not. luse_lw_lim) then
      call comp_flux_vel(k,UUU,VVV,WTK, &
                         UTE,UTW,VTN,VTS,WTKB,this_block)
   else
      VTN  = FLUX_VEL_prev(:,:,3,bid)
   end if


 subroutine comp_flux_vel(k,UUU,VVV,WTK, &
                          UTE,UTW,VTN,VTS,WTKB,this_block)

! !DESCRIPTION:
!  This routine computes the vertical velocity at bottom of a layer.
!  It also returns the tracer flux velocities across the lateral faces.
      UUU, VVV            ! U,V for this block at current time

         VTN(i,j) = p5*(VVV(i  ,j  ,k)*DXU(i  ,j    ,bid)*  &
                                       DZU(i  ,j  ,k,bid) + &
                        VVV(i-1,j  ,k)*DXU(i-1,j    ,bid)*  &
                                       DZU(i-1,j  ,k,bid))

# ----------------------------------------------------------------------
# FVN
# ----------------------------------------------------------------------

      if (partial_bottom_cells) then
         FVN =  p5*VTN*TAREA_R(:,:,bid)/DZT(:,:,k,bid)
      else
         FVN =  p5*VTN*TAREA_R(:,:,bid)
      endif

# ----------------------------------------------------------------------
# VNS
# ----------------------------------------------------------------------

      do n=1,nt
         if (tadvect_itype(n) == tadvect_centered) then
            if (tavg_requested(tavg_VN_TRACER(n))) then
               WORK = FVN*(        TRCR(:,:,k,n) + &
                           eoshift(TRCR(:,:,k,n),dim=2,shift=1))
               call accumulate_tavg_field(WORK,tavg_VN_TRACER(n),bid,k)
            endif

         else
            if (tavg_requested(tavg_VN_TRACER(n))) then
               WORK = c2*FVN*TRACER_N(:,:,n)
               call accumulate_tavg_field(WORK,tavg_VN_TRACER(n),bid,k)
            endif


# ----------------------------------------------------------------------
# my summary
# ----------------------------------------------------------------------

# volume flux between at (nlon_t,nlat_u) with DXZU = DXU*DZU  [cm^3/s]
VTN = [DXZU(i-1)*VVEL(i-1) + DXZU(i)*VVEL(i) ]/2

# half the northward rate at (nlon_t,nlat_u) with TVOL = DZT*TAREA  [1/s]
# FVN = VTN/TVOL/2

# northward advection of Salt at (nlon_t,nlat_u)  [g/kg/s]
# VNS = FVN*(SALT(j)+SALT(j+1))