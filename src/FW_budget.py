""" Atlantic Freswater Budget under Climate Change

this script generates the freshwater budget files for the CESM HIGH and LOW simulations
figures are then generated in jupyter notebooks

SALT bias figure creates in `SALT_bias.ipynb`
AMOC trends in `AMOC.ipynb`

necessary fields:
1. SFWF: PREC, EVAP, ROFF, salt terms,
2. meridional transport: VVEL, SALT, VNS
"""

import os
import sys
import numpy as np
import xarray as xr

from paths import CESM_filename, path_prace, path_results, file_RMASK_ocn_low, file_ex_ocn_ctrl
from timeseries import IterateOutputCESM
from xr_DataArrays import xr_AREA
from xr_regression import ocn_field_regression
from aa_derivation_fields import DeriveField

""" functions to create certain files """

lat_bands = [(45,60), (10,45), (-10,10), (-34,10), (-34,60), (60,90)]

def trex(fn, fct, kwargs={}):
    """  try: existence / except: create
    checks whether file exists already and if not, creates it
    input:
    fn      .. file name of file that function is supposed to create
    fct     .. function that will create data stored as `fn`
    kwargs  .. keyword arguments for function `fct`
    """
    assert type(fn) in [str, list]
    assert hasattr(fct, '__call__')
    assert type(kwargs)==dict
    try:
        if type(fn)==list:
            for f in fn: assert os.path.exists(f), f' now creating file:  {fn[46:]}'
            print(f' files exists:  {fn[0][46:]} ...')
        elif type(fn)==str:
            assert os.path.exists(fn), f' now creating file:  {fn[46:]}'
            print(f' file exists:  {fn[46:]}')
    except:
        fct(**kwargs)
    return

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def make_lat_dict():
    """ creates dictionary of nlat and nlon of E-W/N-S sections """
    sections = ['Drake', 'Agulhas', 'Med', '30W', 'Bering', '34S', '10S', 'Eq', '10N', '26N', '45N', '60N', '65N']
    section_dict = {}
    for r, run in enumerate(['lr1', 'rcp']):
        ds = [dsl, dsh][r]
        for s in sections:
            if   s=='Drake'  : lat, lon = (-71, -55  ), 293
            elif s=='Agulhas': lat, lon = (-75, -34.5),  20
            elif s=='Med'    : lat, lon = ( 30,  40  ), 364.5
            elif s=='30W'    : lat, lon = (-79,  68  ), 330
            elif s=='Bering' : lat, lon =  67         , (180  , 200)
            elif s=='34S'    : lat, lon = -34.5       , (306  ,  20)
            elif s=='10S'    : lat, lon = -10         , (322  ,  15)
            elif s=='Eq'     : lat, lon =   0         , (309  ,  10)
            elif s=='10N'    : lat, lon =  10         , (298  , 346)
            elif s=='26N'    : lat, lon =  26.5       , (279.5, 347)
            elif s=='45N'    : lat, lon =  45         , (298  , 359)
            elif s=='60N'    : lat, lon =  60         , (280  ,  [22,19][r])
            elif s=='65N'    : lat, lon =  65         , (315  ,  [22,19][r])
                
            if type(lat)==tuple:  # N-S section
                nlon1, nlat1 = find_indices_near_coordinate(ds, lat[0], lon)
                nlon2, nlat2 = find_indices_near_coordinate(ds, lat[1], lon)
                nlon = int((nlon1+nlon2)/2)
                section_dict[s] = {'nlon':int(nlon), 'nlat':slice(int(nlat1),int(nlat2))}
            elif type(lon)==tuple:  # E-W section
                nlon1, nlat1 = find_indices_near_coordinate(ds, lat, lon[0])
                nlon2, nlat2 = find_indices_near_coordinate(ds, lat, lon[1])
                nlat = int((nlat1+nlat2)/2)
                section_dict[s] = {'nlon':slice(int(nlon1), int(nlon2)), 'nlat':int(nlat)}
        save_obj(section_dict, f'{path_results}/sections/section_dict_{["low","high"][r]}')
    return

def make_SFWF_means(run, fn):
    """ calculates mean of SFWF fields"""
    # ctrl: 2min 14s, lpd: 3 sec
    if run=='ctrl':  yy = 200
    elif run=='lpd':  yy = 500
    for i, qs in enumerate([['EVAP_F', 'PREC_F', 'ROFF_F'],['SALT_F']]):
        name = '_'.join(qs)
        ds_list = [xr.open_dataset(f'{path_prace}/{run}/ocn_yrly_{name}_0{y}.nc') for y in np.arange(yy,yy+30)]
        fn = f'{path_prace}/{run}/{name}_{run}_mean_{years[0]}-{years[-1]+1}.nc'
        xr.concat(ds_list, dim='time').mean(dim='time').to_netcdf()
    return

def make_SFWF_trends(run):
    """ calculate linear trends of SFWF """
    for qs in [['EVAP_F', 'PREC_F', 'ROFF_F'],['SALT_F']]:
        name = '_'.join(qs)
        for i, (y,m,f) in enumerate(IterateOutputCESM(domain='ocn', run=run, tavg='yrly', name=name)):
            ds_ = xr.open_dataset(f, decode_times=False)
            if i==0:  ds = []
            ds.append(ds_)
        ds = xr.concat(ds, dim='time')
        for q in qs:
            ds[q].to_netcdf(f'{path_prace}/{run}/{q}_yrly_{run}.nc')
            trend = ocn_field_regression(xa=ds[q], run=run)
            trend[0].to_netcdf(f'{path_prace}/{run}/{q}_yrly_trend_{run}.nc')
    return

def make_SFWF_surface_int_dict():
    """ calculate SFWF surface integrals for specific regions """
    RMASK_ocn = xr.open_dataset(file_ex_ocn_ctrl, decode_times=False).REGION_MASK
    RMASK_low = xr.open_dataset(file_RMASK_ocn_low, decode_times=False).REGION_MASK
    # Atlantic + Labrador + GIN, no Med
    Atl_MASK_ocn = xr.DataArray(np.in1d(RMASK_ocn, [6,8,9]).reshape(RMASK_ocn.shape),
                                dims=RMASK_ocn.dims, coords=RMASK_ocn.coords)
    Atl_MASK_low = xr.DataArray(np.in1d(RMASK_low, [6,8,9]).reshape(RMASK_low.shape),
                                dims=RMASK_low.dims, coords=RMASK_low.coords)
    AREA_ocn = xr_AREA(domain='ocn')
    AREA_low = xr_AREA(domain='ocn_low')

    ds_ctrl = xr.open_dataset(f'{path_prace}/ctrl/EVAP_F_PREC_F_ROFF_F_ctrl_mean_200-230.nc')
    ds_lpd  = xr.open_dataset(f'{path_prace}/lpd/EVAP_F_PREC_F_ROFF_F_lpd_mean_500-530.nc')
    Sm_ctrl = xr.open_dataset(f'{path_prace}/lpd/SALT_F_ctrl_mean_200-230.nc').SALT_F
    Sm_lpd  = xr.open_dataset(f'{path_prace}/lpd/SALT_F_lpd_mean_500-530.nc').SALT_F
    Et_rcp = xr.open_dataarray(f'{path_prace}/rcp/EVAP_F_yrly_trend.nc')
    Pt_rcp = xr.open_dataarray(f'{path_prace}/rcp/PREC_F_yrly_trend.nc')
    Rt_rcp = xr.open_dataarray(f'{path_prace}/rcp/ROFF_F_yrly_trend.nc')
    St_rcp = xr.open_dataarray(f'{path_prace}/rcp/SALT_F_yrly_trend.nc')
    Et_lr1 = xr.open_dataarray(f'{path_prace}/lr1/EVAP_F_yrly_trend.nc')
    Pt_lr1 = xr.open_dataarray(f'{path_prace}/lr1/PREC_F_yrly_trend.nc')
    Rt_lr1 = xr.open_dataarray(f'{path_prace}/lr1/ROFF_F_yrly_trend.nc')
    St_lr1 = xr.open_dataarray(f'{path_prace}/lr1/SALT_F_yrly_trend.nc')

    for i, sim in enumerate(['HIGH', 'LOW']):
        d = {}
        (Em, Pm, Rm) = [(ds_ctrl.EVAP_F, ds_ctrl.PREC_F, ds_ctrl.ROFF_F),
                        (ds_lpd.EVAP_F, ds_lpd.PREC_F, ds_lpd.ROFF_F)][i]
        (Et, Pt, Rt) = [(Et_rcp, Pt_rcp, Rt_rcp),
                        (Et_lr1, Pt_lr1, Rt_lr1)][i]
        AREA = [AREA_ocn, AREA_low][i]
        RMASK = [RMASK_ocn, RMASK_low][i]
        MASK = [Atl_MASK_ocn, Atl_MASK_low][i]

        for (latS, latN) in lat_bands:
            if latN>60:  # + Arctic & Hudson Bay
                MASK = MASK + RMASK.where(RMASK_ocn==10) + RMASK.where(RMASK_ocn==11)  
            MASK_ = MASK.where(Pm.ULAT<latN).where(Pm.ULAT>latS)
            AREA_total = AREA.where(MASK_).sum()
            
            # integrals of mean
            Pmi     = (Pm.where(MASK_)*AREA).sum().values
            Emi     = (Em.where(MASK_)*AREA).sum().values
            Rmi     = (Rm.where(MASK_)*AREA).sum().values
            d['Pmi_mmd'] = Pmi/AREA_total.values*24*3600       # [kg/s] -> [mm/d]
            d['Emi_mmd'] = Emi/AREA_total.values*24*3600
            d['Rmi_mmd'] = Rmi/AREA_total.values*24*3600
            d['Pmi_Sv']  = Pmi/1e9                             # [kg/m^2/s] -> [Sv]
            d['Emi_Sv']  = Emi/1e9
            d['Rmi_Sv']  = Rmi/1e9
            d['Tmi_mmd'] = d['Pmi_mmd'] + d['Emi_mmd'] + d['Rmi_mmd']
            d['Tmi_Sv']  = d['Pmi_Sv']  + d['Emi_Sv']  + d['Rmi_Sv']
            
            # integrals of trends
            Pti     = (Pt.where(MASK_)*AREA).sum().values
            Eti     = (Et.where(MASK_)*AREA).sum().values
            Rti     = (Rt.where(MASK_)*AREA).sum().values
            d['Pti_mmd'] = Pti/AREA_total.values*24*3600*365*100  # [mm/d/100yr]
            d['Eti_mmd'] = Eti/AREA_total.values*24*3600*365*100
            d['Rti_mmd'] = Rti/AREA_total.values*24*3600*365*100
            d['Pti_Sv']  = Pti/1e9*365*100                        # [Sv/100yr]
            d['Eti_Sv']  = Eti/1e9*365*100
            d['Rti_Sv']  = Rti/1e9*365*100
            d['Tti_mmd'] = d['Pti_mmd'] + d['Eti_mmd'] + d['Rti_mmd']
            d['Tti_Sv']  = d['Pti_Sv']  + d['Eti_Sv']  + d['Rti_Sv']
            
            print(f'\n{latS}N to {latN}N,   {AREA_total.values:4.2E} m^2\n')
            print('             PREC  EVAP  ROFF    SUM')

            print(f'[mm/d]      {d["Pmi_mmd"]:5.2f} {d["Emi_mmd"]:5.2f} {d["Rmi_mmd"]:5.2f} {d["Tmi_mmd"]:5.4f}')
            print(f'[mm/d/100y] {d["Pti_mmd"]:5.2f} {d["Eti_mmd"]:5.2f} {d["Rti_mmd"]:5.2f} {d["Tti_mmd"]:5.4f}')
            print(f'[Sv]        {d["Pmi_Sv"]:5.2f} {d["Emi_Sv"]:5.2f} {d["Rmi_Sv"]:5.2f} {d["Tmi_Sv"]:5.4f}')
            print(f'[Sv/100y]   {d["Pti_Sv"]:5.2f} {d["Eti_Sv"]:5.2f} {d["Rti_Sv"]:5.2f} {d["Tti_Sv"]:5.4f}')
            print(f'[%/100y]    {d["Pti_Sv"]/d["Pmi_Sv"]*100:5.1f} {d["Eti_Sv"]/d["Emi_Sv"]*100:5.1f} {d["Rti_Sv"]/d["Rmi_Sv"]*100:5.1f} {d["Tti_Sv"]/d["Tmi_Sv"]*100:5.1f}')
            fn = f'{path_results}/SFWF/Atlantic_SFWF_integrals_{sim}_{latS}N_{latN}N'
            save_obj(d, fn)
    return

def make_ddt_SALT_int():

    return


if __name__=='__main__':
    # 0.   geometry:  create dicts with nlat of specific latitudes
    fn = [f'{path_results}/sections/section_dict_low.pkl',
          f'{path_results}/sections/section_dict_high.pkl']
    fct = make_lat_dict
    trex(fn=fn, fct=fct)

    # 1.   annually averaged files
    fields_ = [['SALT','VNS','UES'],          # d/dt, transport
               ['EVAP_F','PREC_F','ROFF_F'],  # SFWF (FW)
               ['SALT_F'],                    # SFWF (SALT)
              ]
    for fields in fields_:
        name = '_'.join(fields)
        for run in ['ctrl', 'rcp', 'lpd', 'lr1']:
            fct = DeriveField(run).yrly_avg_nc
            if run=='ctrl':   years = np.arange(200,230)
            elif run=='lpd':  years = np.arange(500,530)
            else:             years = np.arange(2000,2101)

            fn = []
            for y in years:
                fn.append(CESM_filename(domain='ocn',run=run,y=y,m=0, name=name))

            kwargs = dict(domain='ocn', fields=fields, years=years)
            trex(fn=fn, fct=fct, kwargs=kwargs)

    # 2.   SFWF
    # 2.1  CTRL means and RCP trends for surface fields
    for run in ['ctrl', 'rcp', 'lpd', 'lr1']:
        if run=='ctrl':   years = np.arange(200,230)
        elif run=='lpd':  years = np.arange(500,530)
        else:             years = np.arange(2000,2101)

        if run in ['ctrl', 'lpd']:  # means
            fct = make_SFWF_means
            fn = [f'{path_prace}/{run}/EVAP_F_PREC_F_ROFF_F_{run}_mean_{years[0]}-{years[-1]+1}.nc',
                  f'{path_prace}/{run}/SALT_F_{run}_mean_{years[0]}-{years[-1]+1}.nc']
            kwargs = dict(run=run, fn=fn)
            trex(fn=fn, fct=fct, kwargs=kwargs)
        
        elif run in ['rcp', 'lr1']:  # trends
            for q in ['EVAP_F','PREC_F','ROFF_F', 'SALT_F']:
                fn = [f'{path_prace}/{run}/{q}_yrly_{run}.nc',
                      f'{path_prace}/{run}/{q}_yrly_trend_{run}.nc']
                fct = make_SFWF_trends
                kwargs = dict(run=run)
                trex(fn=fn, fct=fct, kwargs=kwargs)


    # 2.2  surface area integrals: means and trends
    fct = make_SFWF_surface_int_dict
    fn = []
    for i, sim in ['HIGH', 'LOW']:
        for (latS, latN) in lat_bands:
            fn.append(f'{path_prace}/SFWF/Atlantic_SFWF_integrals_{sim}_{latS}N_{latN}N')
    trex(fn=fn, fct=fct)


    # 3.   meridional transport terms
    # 3.1  as a function of latitude


    # 3.2  select latitudes, transport divergence: means and trends


    # 4.1  d/dt (F/S)


    # 4.2  F/S_{diff} = residual


""" Figures in jupyter notebooks """
# 4:  mean/trend maps:  `SFWF.ipynb`

# 6:  regional budget:  `Atl_FW_budget.ipynb`

# 6.1:  surface integrals
#     f'{path_prace}/{ctrl/lpd}/EVAP_PREC_ROFF_{ctrl/lpd}_mean_200-230.nc' 
#     f'{path_prace}/{rcp/lr1}/{EVAP/PREC/ROFF}_yrly_trend.nc'
#     calculation of surface integrals in .ipynb 
#         -> f'{path_results}/SFWF/Atlantic_SFWF_integrals_{sim}_{latS}N_{latN}N'

# 6.2:  d/dt

# 6.3:  meridional transport
#     f'{path_prace}/{run}/FW_SALT_fluxes_{run}.nc'

# 7:  latitudinal dependence of transport terms:  `M_ov.ipynb`



