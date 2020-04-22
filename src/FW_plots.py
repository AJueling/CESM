""" plotting freshwater/salt fluxes """

import numpy as np
import pickle
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt

from paths import path_results, path_prace
from FW_budget import load_obj, lat_bands
from mimic_alpha import make_color_transparent as mct
from xr_regression import ocn_field_regression, xr_linear_trend

# term                 where is the data stored                         to improve
# d/dt(SALT)           f'{path_results}/SALT/SALT_integrals_ctrl.nc'    use proper DZT and nlats as boundaries
# SFWF                                                                  add other terms, use proper nlats as boundaries
# fluxes/convergences                                                   include Baltic Sea

# color scheme: 
# fluxes/convergences   C0 (ov); C1 (az); C2 (eddy); C3 (total), C8 (Med/BS)
# SFWF                  C4 (P+E+R); C6 (P+E); C7 (R)
# mixing                C5
# d/dt(SALT)            C9

S0 = 35.
lat_bounds = ['90N', '60N', '45N', '10N', '10S', '34S']
nabla, Delta, mix, deg = r'$\nabla$', r'$\Delta$', r'$_{mix}$', r'$\!^\circ\!$'
ov, az, ed, to = r'$_{ov}$', r'$_{az}$', r'$_{eddy}$', r'$_{total}$'
xrkw = {'decode_times':False}                 # xr.open_* keywords
vbkw = dict(y=0, width=.1, align='edge')      # vertical bar keywords
hbkw = dict(left=0, height=.1, align='edge')  # horizontal bar keywords

def nlat_at(sim, lat):
    """ returns nlat of high or low res grid for certain latitudes """
    assert sim in ['HIGH', 'LOW']
    if type(lat)==str:
        assert lat in lat_bounds
        lat_str = lat
        if lat_str[-1]=='N':  lat = float(lat_str[:-1])
        if lat_str[-1]=='S':  lat = -float(lat_str[:-1])
    elif type(lat)==float or type(lat)==int:
        if lat<0: lat_str = str(-lat)+'S'
        else:     lat_str = str(lat)+'N'
        assert lat_str in lat_bounds
    else:
        print(f'lat {lat}, type is {type(lat)}')
    
    sdict = load_obj(f'{path_results}/sections/section_dict_{sim.lower()}')
    nlat = sdict[lat_str]['nlat']
    if sim=='LOW' and lat==-34:  nlat = 85
    if sim=='LOW' and lat==60:   nlat = 348
    if sim=='HIGH' and lat==-34:   nlat = 811
    if sim=='HIGH' and lat==60:   nlat = 1867
    return nlat

def get_ddt_SALT(sim, latS, latN):
    """ calculate d(S)/dt term in [kg/s] 
    as the difference between the last and first years
    """
    if sim=='HIGH':   run, rcp = 'ctrl', 'rcp'
    elif sim=='LOW':  run, rcp = 'lpd' , 'lr1'
    # files created in `FW_budget.py` `make_SALT_vol)_int`
    salt_ctl = xr.open_dataset(f'{path_results}/SALT/SALT_integrals_{run}.nc', **xrkw)
    salt_rcp = xr.open_dataset(f'{path_results}/SALT/SALT_integrals_{rcp}.nc', **xrkw)
    assert len(salt_ctl.time)==101  # dataset contains only integrals for 101 years concurrent to rcp
    if latS==-34 and latN==60:
        for i, (latS_, latN_) in enumerate([(-34,-10),(-10,10),(10,45),(45,60)]):
            n1 = f'SALT_0-1000m_timeseries_{latS_}N_{latN_}N'
            n2 = f'SALT_below_1000m_timeseries_{latS_}N_{latN_}N'
            if i==0:
                salt_ctl_ = (salt_ctl[n1] + salt_ctl[n2])
                salt_rcp_ = (salt_rcp[n1] + salt_rcp[n2])
            else:
                salt_ctl_ += (salt_ctl[n1] + salt_ctl[n2])
                salt_rcp_ += (salt_rcp[n1] + salt_rcp[n2])
        ddtS_ctl = (salt_ctl_.isel(time=30)-salt_ctl_.isel(time=0)).values/30/365/24/3600
        ddtS_rcp = (salt_rcp_.isel(time=-1)-salt_rcp_.isel(time=0)).values/101/365/24/3600
    else:

        n1 = f'SALT_0-1000m_timeseries_{latS}N_{latN}N'
        n2 = f'SALT_below_1000m_timeseries_{latS}N_{latN}N'
        # multiply by 100 because saved as g_S/kg_W * m^3 = 1e-2 kg_S
        salt_ctl_ = (salt_ctl[n1] + salt_ctl[n2])
        salt_rcp_ = (salt_rcp[n1] + salt_rcp[n2])
        ddtS_ctl = (salt_ctl_.isel(time=30)-salt_ctl_.isel(time=0)).values/30/365/24/3600
        ddtS_rcp = (salt_rcp_.isel(time=-1)-salt_rcp_.isel(time=0)).values/101/365/24/3600
    return ddtS_ctl, ddtS_rcp

def get_SFWF(sim, quant, latS, latN):
    """ reads the integrated surface freshwater fluxes between latS and latN from pickled object """
    assert sim in ['HIGH', 'LOW'] and  quant in ['SALT', 'FW']
    assert latS in [-34, -10, 10, 45, 60]
    assert latN in [-10, 10, 45, 60, 90]

    d = load_obj(f'{path_results}/SFWF/Atlantic_SFWF_integrals_{sim}_{latS}N_{latN}N')
    # if quant=='SALT':  fac = -.00347e9  # 1e9 kg/s/Sv * salinity_factor of POP2 \approx -S0/1e6*.1?
    if quant=='SALT':  fac = -35e6  # 1e9 kg/s/Sv * salinity_factor of POP2 \approx -S0/1e6*.1?
    elif quant=='FW':  fac = 1.
    mean, trend = {}, {}
    
    mean['PE']   = (d[f'Pmi_Sv']+d[f'Emi_Sv'])*fac
    mean['R']    = d[f'Rmi_Sv']*fac
    mean['SFWF'] = d[f'Tmi_Sv']*fac

    trend['PE']  = (d[f'Pti_Sv']+d[f'Eti_Sv'])*fac
    trend['P']   = d[f'Pti_Sv']*fac
    trend['E']   = d[f'Eti_Sv']*fac
    trend['R']   = d[f'Rti_Sv']*fac
    trend['SFWF'] = d[f'Tti_Sv']*fac
    return mean, trend

def get_fluxes(sim, quant, lat=None, latS=None, latN=None):
    """ fluxes at `lat` and neg. flux convergences between latS and latN"""
    assert sim in ['HIGH', 'LOW'] and  quant in ['SALT', 'FW']
    assert latS in [-34, -10, 10, 45, None]
    assert latN in [-10, 10, 45, 60, None]
    if sim=='HIGH':   run, rcp = 'ctrl', 'rcp'
    elif sim=='LOW':  run, rcp = 'lpd', 'lr1'

    dso = xr.open_dataset(f'{path_prace}/Mov/FW_SALT_fluxes_{run}.nc', **xrkw)
    dst = xr.open_dataset(f'{path_prace}/Mov/FW_SALT_fluxes_{rcp}.nc', **xrkw)
    if lat is not None:   nlat = nlat_at(sim, lat)
    if latS is not None:  nlat_S, nlat_N =  nlat_at(sim, latS), nlat_at(sim, latN)
    # select 30 year mean of fluxes fro branch-off year 200/500
    if sim=='HIGH':   dso = dso.isel(time=slice(200,230)).mean('time')
    elif sim=='LOW':  dso = dso.isel(time=slice(346,376)).mean('time')

    fluxes, fac = {}, 365*100
    if quant=='SALT' and lat is not None:
        fluxes['ov']  = dso.Sov.isel(nlat_u=nlat)
        fluxes['az']  = dso.Saz.isel(nlat_u=nlat)
        fluxes['ed']  = dso.Se.isel(nlat_u=nlat)
        fluxes['to']  = dso.St.isel(nlat_u=nlat)
        
        fluxes['tov'] = fac*xr_linear_trend(dst.Sov.isel(nlat_u=nlat))
        fluxes['taz'] = fac*xr_linear_trend(dst.Saz.isel(nlat_u=nlat))
        fluxes['ted'] = fac*xr_linear_trend(dst.Se.isel(nlat_u=nlat))
        fluxes['tto'] = fac*xr_linear_trend(dst.St.isel(nlat_u=nlat))

    elif quant=='FW' and lat is not None:
        fluxes['ov']  = dso.Fov.isel(nlat_u=nlat)
        fluxes['az']  = dso.Faz.isel(nlat_u=nlat)
        fluxes['ed']  = -dso.Se.isel(nlat_u=nlat)/35e6  # [kg/s] -> [Sv]
        fluxes['to']  = (dso.Fov+dso.Faz+fluxes['ed']).isel(nlat_u=nlat)

        fluxes['tov'] = fac*xr_linear_trend(dst.Fov.isel(nlat_u=nlat))
        fluxes['taz'] = fac*xr_linear_trend(dst.Faz.isel(nlat_u=nlat))
        fluxes['ted'] = fac*xr_linear_trend(dst.Se.isel(nlat_u=nlat)/35e6)
        fluxes['tto'] = fac*xr_linear_trend((dst.Fov+dst.Faz-dst.Se/35e6).isel(nlat_u=nlat))

    conv = {}
    if quant=='SALT' and latS is not None:
        conv['ov']  = dso.Sov.isel(nlat_u=nlat_S) - dso.Sov.isel(nlat_u=nlat_N)
        conv['az']  = dso.Saz.isel(nlat_u=nlat_S) - dso.Saz.isel(nlat_u=nlat_N)
        conv['ed']  = dso.Se.isel(nlat_u=nlat_S)  - dso.Se.isel(nlat_u=nlat_N)
        conv['to']  = dso.St.isel(nlat_u=nlat_S)  - dso.St.isel(nlat_u=nlat_N)
        
        conv['tov'] = fac*xr_linear_trend(dst.Sov.isel(nlat_u=nlat_S) - dst.Sov.isel(nlat_u=nlat_N))
        conv['taz'] = fac*xr_linear_trend(dst.Saz.isel(nlat_u=nlat_S) - dst.Saz.isel(nlat_u=nlat_N))
        conv['ted'] = fac*xr_linear_trend(dst.Se.isel(nlat_u=nlat_S)  - dst.Se.isel(nlat_u=nlat_N))
        conv['tto'] = fac*xr_linear_trend(dst.St.isel(nlat_u=nlat_S)  - dst.St.isel(nlat_u=nlat_N))
        
    elif quant=='FW' and latS is not None:
        conv['ov']  = dso.Fov.isel(nlat_u=nlat_S) - dso.Fov.isel(nlat_u=nlat_N)
        conv['az']  = dso.Faz.isel(nlat_u=nlat_S) - dso.Faz.isel(nlat_u=nlat_N)
        conv['ed']  = (dso.Se.isel(nlat_u=nlat_S) - dso.Se.isel(nlat_u=nlat_N))/35e6
        conv['to']  = (dso.Fov+dso.Faz-dso.Se/35e6).isel(nlat_u=nlat_S) - \
                      (dso.Fov+dso.Faz-dso.Se/35e6).isel(nlat_u=nlat_N)

        conv['tov'] = fac*xr_linear_trend(dst.Fov.isel(nlat_u=nlat_S) - dst.Fov.isel(nlat_u=nlat_N))
        conv['taz'] = fac*xr_linear_trend(dst.Faz.isel(nlat_u=nlat_S) - dst.Faz.isel(nlat_u=nlat_N))
        conv['ted'] = fac*xr_linear_trend(dst.Se.isel(nlat_u=nlat_S)  - dst.Se.isel(nlat_u=nlat_N))/35e6
        conv['tto'] = fac*xr_linear_trend((dst.Fov+dst.Faz-dst.Se/35e6).isel(nlat_u=nlat_S) - \
                                          (dst.Fov+dst.Faz-dst.Se/35e6).isel(nlat_u=nlat_N))
    return fluxes, conv

def get_BS_Med(sim):
    """ Bering Strait and Striat of Gibraltat transports:
    mean of control and trend of rcp simulations
    """
    if sim=='HIGH':  ctl, rcp = 'ctrl', 'rcp'
    elif sim=='LOW':  ctl, rcp = 'lpd', 'lr1'

    ds_ctl = xr.open_mfdataset(f'{path_prace}/Mov/FW_Med_Bering_{ctl}*.nc', concat_dim='time', decode_times=False)
    ds_rcp = xr.open_mfdataset(f'{path_prace}/Mov/FW_Med_Bering_{rcp}*.nc', concat_dim='time', decode_times=False)

    return ds_ctl.mean('time')

def draw_region_labels(fig):
    for i, lb in enumerate([f'34{deg}S',f'10{deg}S',f'10{deg}N',f'45{deg}N',f'60{deg}N','BS',f'34{deg}S',f'60{deg}N']):
        fig.text(3/32+i*1/8,.95, lb, va='center', ha='center', fontsize=8, color='grey')
    return


def FW_summary_plot(quant):
    """ plots main budget terms for certain regions
    d/dt = conv + SFWF + diffusion
    """
    assert quant in ['FW', 'SALT']
    if quant=='FW':      q, Q = 0, 'F'
    elif quant=='SALT':  q, Q = 1, 'S'

    f = plt.figure(figsize=(8,2))
    draw_region_labels(f)  

    for i, (latS, latN) in enumerate(lat_bands):

        #region: define axes
        if i==0:  # Atlantic 34S-60N
            ax = f.add_axes([13.5/16,.05,1/8,.8])
            for s in ['right', 'top', 'bottom']:  ax.spines[s].set_visible(False)
            ax.xaxis.set_ticks([])
            ax.axhline(c='k', lw=.5)
            # ax.set_ylabel(f'{quant} fluxes in [{["Sv","kg/s"][q]}]')
            vert_lim = [(-1.2,1.2),(-3.9e7,3.9e7)][q]
            ax.set_xlim((-.25,1.2))
            ax.set_ylim(vert_lim)

        else:
            ax = f.add_axes([(2*i+1)/16-1.5/16,.05,1/8,.8])
            for s in ['right', 'top', 'bottom']:  ax.spines[s].set_visible(False)
            ax.xaxis.set_ticks([])
            ax.axhline(c='k', lw=.5)
            if i==1:  ax.set_ylabel(f'{quant} fluxes in [{["Sv","kg/s"][q]}]')
            else:     ax.yaxis.set_ticklabels([])
            vert_lim = [(-.9,.9),(-3.2e7,3.2e7)][q]
            ax.set_xlim((-.25,1.2))
            ax.set_ylim(vert_lim)
        #endregion
        
        for s, sim in enumerate(['HIGH', 'LOW']):  

            #region: d/dt
            if quant=='SALT':  fac = 1
            elif quant=='FW':  fac = -1e-6/S0
            ddtS_ctl, ddtS_rcp = get_ddt_SALT(sim=sim, latS=latS, latN=latN)
            bdt = ax.bar(x=s/10, height=ddtS_ctl*fac, **vbkw, color=mct('C9', alpha=1-.3*s), label=r'$\Delta_t$S')
            ax.arrow(x=.05+s/10, y=0, dy=ddtS_rcp*fac, dx=0)
            #endregion
            
            #region: surface fluxes
            sfwf_mean, sfwf_trend = get_SFWF(sim=sim, quant=quant, latS=latS, latN=latN)
            bsf = ax.bar(x=.25+s/10 , height=sfwf_mean['SFWF'], **vbkw, color=mct('C4', 1-.3*s), label='SFWF')
            ax.arrow(x=.3+s/10, y=sfwf_mean['SFWF'], dy=sfwf_trend['SFWF'], dx=0)
            #endregion

            #region: meridional flux convergences incl. BS/Med
            inflow = get_BS_Med(sim)
            if latS==60:  # Arctic
                fluxes, conv = get_fluxes(sim=sim, quant=quant, lat=latS, latS=None, latN=None)
                conv['to'] = fluxes['to'] + inflow[f'{Q}_BS'].values/[1e6,1][q]  # [m^3/s] -> [Sv]
            else:
                fluxes, conv = get_fluxes(sim=sim, quant=quant, lat=latS, latS=latS, latN=latN)
                if latS<35 and latN>35:  # Med
                    conv['to'] -= inflow[f'{Q}_Med'].values/[1e6,1][q]  # [m^3/s] -> [Sv]
                bct = ax.bar(x=.5+s/10, height=conv['to'], **vbkw, color=mct('C3', 1-.3*s), label=f'-{nabla}{Q}')
                ax.arrow(x=.55+s/10, y=conv['to'], dx=0, dy=conv['tto'])
            #endregion
            
            #region: mixing term
            diffusion = ddtS_ctl*fac - conv['to'] - sfwf_mean['SFWF']
            bdi = ax.bar(x=.75+s/10, height=diffusion, **vbkw, color=mct('C5', 1-.3*s), label=f'{Q}{mix}')
            #endregion
            
    #region: legend, numbering
    ax = f.add_axes([11.1/16,.2,1/8,.55])
    ax.axis('off')
    ax.legend(handles=[bct, bsf, bdt, bdi], fontsize=7, loc='center')
    f.text(.01,.9,['(a)','(c)'][q])
    #endregion

    plt.savefig(f'{path_results}/Mov/{["FW","SALT"][q]}_region_budget_total.eps')
    return


def FW_region_plot(quant):
    """ plots SFWF, merid. transport and its convergence for specific regions
    quant .. Fw or SALT
    """
    assert quant in ['FW', 'SALT']
    if quant=='FW':      q, Q, labels = 0, 'F', [r'$-\nabla F_{oc}$', r'$F_{oc}$', r'$\int F_{atm}$']
    elif quant=='SALT':  q, Q, labels = 1, 'S', [r'$-\nabla S_{oc}$', r'$S_{oc}$', r'$\int S_{atm}$']
    lat_bounds = ['90N', '60N', '45N', '10N', '10S', '34S']

    f = plt.figure(figsize=(8,3))
    draw_region_labels(f)

    dd = {}
    for i, (latS, latN) in enumerate(lat_bands):
        #region: define axes
        if i==0:  # Atlantic 34S-60N
            ax = f.add_axes([13.5/16,.2,1/8,.55])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])

            axb = f.add_axes([13.5/16,-.05,1/8,.5])
            axt = f.add_axes([13.5/16,.5,1/8,.5])
            axl = f.add_axes([12/16,.2,3/16,.55])
            axr = f.add_axes([14/16,.2,3/16,.55])
            vert_lim = [(-2,2),(-8e7,8e7)][q]
            hor_lim = [(-1,1),(-3e8,3e8)][q]
        else:
            ax = f.add_axes([(2*i+1)/16-1.5/16,.2,1/8,.55])
            if i==5:  ax.spines['right'].set_visible(False)
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])

            axb = f.add_axes([(2*i+1)/16-1.5/16,-.05,1/8,.5])
            axt = f.add_axes([(2*i+1)/16-1.5/16,.5,1/8,.5])
            axl = f.add_axes([(2*i+1)/16-3/16,.2,3/16,.55])
            axr = f.add_axes([(2*i+1)/16,.2,3/16,.55])

            vert_lim = [(-1.5,1.5),(-4e7,4e7)][q]
            hor_lim = [(-.5,.5),(-1.2e8,1.2e8)][q]

        for ax in [axl, axr]:
            ax.set_xlim(hor_lim)
            ax.set_ylim((-.35,1.3))

        for ax in [axb, axt, axl, axr]:
            # ax.patch.set_alpha(0)
            ax.axis('off')
            ax.xaxis.set_ticks([])

        for ax in [axt, axb]:
            ax.set_ylim(vert_lim)
            ax.set_xlim((-.25,1.45))
        #endregion

        for s, sim in enumerate(['HIGH', 'LOW']):
            nlat_ = nlat_at(sim, latS)
            run, rcp = ['ctrl','lpd'][s], ['rcp', 'lr1'][s]

            #region: d/dt  [axb]
            if quant=='SALT':  fac = 1
            elif quant=='FW':  fac = -1e-6/S0
            ddtS_ctl, ddtS_rcp = get_ddt_SALT(sim=sim, latS=latS, latN=latN)
            bdt = axb.bar(x=s/10, height=ddtS_ctl*fac, **vbkw, color=mct('C9', alpha=1-.3*s), label=r'$\Delta_t$S')
            axb.arrow(x=s/10, y=0, dx=ddtS_rcp*fac, dy=0)
            #endregion

            #region: meridional flux (convergences) means/trends, incl. BS/Med  [axl, axr, axb]
            inflow = get_BS_Med(sim)
            if latS==60:
                fluxes, conv = get_fluxes(sim=sim, quant=quant, lat=latS, latS=None, latN=None)
                conv['ov'], conv['az'], conv['ed'] = 0, 0, 0
                conv['to'] = fluxes['to'] + inflow[f'{Q}_BS'].values/[1e6,1][q]  # [m^3/s] -> [Sv]
                conv['tov'], conv['taz'], conv['ted'], conv['tto'] = 0, 0, 0, fluxes['tto']
                axt.bar(x=.8+s/10, height=-inflow[f'{Q}_BS'].values/[1e6,1][q], **vbkw, color=mct('C8', 1-.3*s))
            else:
                fluxes, conv = get_fluxes(sim=sim, quant=quant, lat=latS, latS=latS, latN=latN)
                if latS<35 and latN>35:  # Med
                    conv['to'] -= inflow[f'{Q}_Med'].values/[1e6,1][q]  # [m^3/s] -> [Sv] 
                    axt.bar(x=.8+s/10, height=inflow[f'{Q}_Med'].values/[1e6,1][q], **vbkw, color=mct('C8', 1-.3*s))

            bov = axl.barh(y=.5 +s/10, width=fluxes['ov'], **hbkw, color=mct('C0', 1-.3*s), label=f'(-{nabla}){Q}{ov}')
            baz = axl.barh(y=.25+s/10, width=fluxes['az'], **hbkw, color=mct('C1', 1-.3*s), label=f'(-{nabla}){Q}{az}')
            bed = axl.barh(y=0  +s/10, width=fluxes['ed'], **hbkw, color=mct('C2', 1-.3*s), label=f'(-{nabla}){Q}{ed}')
            bto = axl.barh(y=.75+s/10, width=fluxes['to'], **hbkw, color=mct('C3', 1-.3*s), label=f'(-{nabla}){Q}{to}')

            axl.arrow(x=fluxes['ov'], y=.55+s/10, dx=fluxes['tov'], dy=0)
            axl.arrow(x=fluxes['az'], y=.3 +s/10, dx=fluxes['taz'], dy=0)
            axl.arrow(x=fluxes['ed'], y=.05+s/10, dx=fluxes['ted'], dy=0)
            axl.arrow(x=fluxes['to'], y=.8 +s/10, dx=fluxes['tto'], dy=0)
            
            if i==0:  # draw 60 North values in Atlantic box
                fluxes_, conv_ = get_fluxes(sim=sim, quant=quant, lat=latN, latS=latS, latN=latN)

                axr.barh(y=.5 +s/10, width=fluxes_['ov'], **hbkw, color=mct('C0', 1-.3*s))
                axr.barh(y=.25+s/10, width=fluxes_['az'], **hbkw, color=mct('C1', 1-.3*s))
                axr.barh(y=0  +s/10, width=fluxes_['ed'], **hbkw, color=mct('C2', 1-.3*s))
                axr.barh(y=.75+s/10, width=fluxes_['to'], **hbkw, color=mct('C3', 1-.3*s))

                axr.arrow(x=fluxes_['ov'], y=.55+s/10, dx=fluxes_['tov'], dy=0)
                axr.arrow(x=fluxes_['az'], y=.3 +s/10, dx=fluxes_['taz'], dy=0)
                axr.arrow(x=fluxes_['ed'], y=.05+s/10, dx=fluxes_['ted'], dy=0)
                axr.arrow(x=fluxes_['to'], y=.8 +s/10, dx=fluxes_['tto'], dy=0)

            axb.bar(x=.5 +s/10, height=conv['ov'], **vbkw, color=mct('C0', 1-.3*s))
            axb.bar(x=.75+s/10, height=conv['az'], **vbkw, color=mct('C1', 1-.3*s))
            axb.bar(x=1. +s/10, height=conv['ed'], **vbkw, color=mct('C2', 1-.3*s))
            axb.bar(x=.25+s/10, height=conv['to'], **vbkw, color=mct('C3', 1-.3*s))
        
            axb.arrow(x= .55+s/10, y=conv['ov'], dx=0, dy=conv['tov'])
            axb.arrow(x= .8 +s/10, y=conv['az'], dx=0, dy=conv['taz'])
            axb.arrow(x=1.05+s/10, y=conv['ed'], dx=0, dy=conv['ted'])
            axb.arrow(x= .3 +s/10, y=conv['to'], dx=0, dy=conv['tto'])


            #endregion

            #region: surface fluxes [axt]
            sfwf_mean, sfwf_trend = get_SFWF(sim=sim, quant=quant, latS=latS, latN=latN)

            bpe = axt.bar(x=.55+s/10, height=-sfwf_mean['PE']  , **vbkw, color=mct('C6', 1-.3*s), label='P-E')
            br  = axt.bar(x=.3+s/10 , height=-sfwf_mean['R']   , **vbkw, color=mct('C7', 1-.3*s), label='R')
            bsf = axt.bar(x=.0+s/10 , height=-sfwf_mean['SFWF'], **vbkw, color=mct('C4', 1-.3*s), label='SFWF')

            axt.arrow(x=.6+s/10 , y=-sfwf_mean['PE']  , dx=0, dy=-sfwf_trend['PE'] )
            axt.arrow(x=.62+s/10, y=0                 , dx=0, dy=-sfwf_trend['P']  )
            axt.arrow(x=.58+s/10, y=0                 , dx=0, dy=-sfwf_trend['E']  )
            axt.arrow(x=.35+s/10, y=-sfwf_mean['R']   , dx=0, dy=-sfwf_trend['R']  )
            axt.arrow(x=.05+s/10, y=-sfwf_mean['SFWF'], dx=0, dy=-sfwf_trend['SFWF'])
            #endregion

            #region: diffusion [axt]
            diffusion = ddtS_ctl*fac - sfwf_mean['SFWF'] - conv['to']
            bdi = axt.bar(x=1+s/10, height=diffusion, **vbkw, color=mct('C5', alpha=1-.3*s), label=f'{Q}{mix}')
            #endregion

    #region: legend, numbering, scales
    ax = f.add_axes([11.1/16,.2,1/8,.55])
    ax.axis('off')
    ax.legend(handles=[bov, baz, bed, bto, bpe, br, bsf, bdt, bdi], fontsize=7, loc='center')

    # axa = f.add_axes([0.04,.08,.125,.125])
    # axa.set_xlim((0,[.5,2e2][q]))    
    # axa.set_ylim((0,[3,1.5e2][q]/2))
    # axa.text(.1,.1,['Sv\nSv/100yr','kt/s\nkt/s/100yr'][q], transform=axa.transAxes, fontsize=8)
    # axa.patch.set_alpha(0)
    # axa.spines['top'].set_visible(False)
    # axa.spines['right'].set_visible(False)
    f.text(.01,.9,['(b)','(d)'][q])
    #endregion

    plt.savefig(f'{path_results}/Mov/{["FW","SALT"][q]}_region_budget.eps')
    return

