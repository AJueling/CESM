import numpy as np
import pickle
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt

from paths import path_results, path_prace
from FW_budget import load_obj, lat_bands
from mimic_alpha import make_color_transparent as mct
from xr_regression import ocn_field_regression, xr_linear_trend


lat_bounds = ['90N', '60N', '45N', '10N', '10S', '34S']
nabla, Delta, diff = r'$\nabla$', r'$\Delta$', r'$_{diff}$'
xrkw = {'decode_times':False}              # xr.open_* keywords
vbkw = dict(y=0, width=.1, align='edge')   # vertical bar keywords
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

# term                 where is the data stored                         to improve
# d/dt(SALT)           f'{path_results}/SALT/SALT_integrals_ctrl.nc'    use proper DZT and nlats as boundaries
# SFWF                                                                  add other terms, use proper nlats as boundaries
# fluxes          
# flux divergences

# colors: 
# d/dt(SALT)           C9
# SFWF                 C6 (P+E); C7 (R); C4 (P+E+R)
# fluxes/divergences   C0 (ov); C1 (az); C2 (eddy); C3 (total)
# diffusion            C5

def get_ddt_SALT(sim, latS, latN):
    """ calculate d(S)/dt term in """
    if sim=='HIGH':   run, rcp = 'ctrl', 'rcp'
    elif sim=='LOW':  run, rcp = 'lpd' , 'lr1'
    salt_ctl = xr.open_dataset(f'{path_results}/SALT/SALT_integrals_{run}.nc', **xrkw)
    salt_rcp = xr.open_dataset(f'{path_results}/SALT/SALT_integrals_{rcp}.nc', **xrkw)
    assert len(salt_ctl.time)==101  # dataset contains only integrals for 101 years concurrent to rcp
    n1 = f'SALT_0-1000m_timeseries_{latS}N_{latN}N'
    n2 = f'SALT_below_1000m_timeseries_{latS}N_{latN}N'
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
    if quant=='SALT':  fac = -.00347e9  # 1e9 kg/s/Sv * salinity_factor of POP2
    elif quant=='FW':  fac = 1.
    mean, trend = {}, {}
    
    mean['PE']   = (d[f'Pmi_Sv']+d[f'Emi_Sv'])*fac
    mean['R']    = d[f'Rmi_Sv']*fac
    mean['PER']  = d[f'Tmi_Sv']*fac

    trend['PE']  = (d[f'Pti_Sv']+d[f'Eti_Sv'])*fac
    trend['P']   = d[f'Pti_Sv']*fac
    trend['E']   = d[f'Eti_Sv']*fac
    trend['R']   = d[f'Rti_Sv']*fac
    trend['PER'] = d[f'Tti_Sv']*fac
    return mean, trend

def get_fluxes(sim, quant,lat=None, latS=None, latN=None):
    """ fluxes at `lat` and neg. flux divergences between latS and latN"""
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
        fluxes['ed']  = 0

        fluxes['tov'] = fac*xr_linear_trend(dst.Fov.isel(nlat_u=nlat))
        fluxes['taz'] = fac*xr_linear_trend(dst.Faz.isel(nlat_u=nlat))
        fluxes['ted'] = 0

        fluxes['to']  = (dso.Fov+dso.Faz).isel(nlat_u=nlat)
        fluxes['tto'] = fac*xr_linear_trend((dst.Fov+dst.Faz).isel(nlat_u=nlat))

    div = {}
    if quant=='SALT' and latS is not None:
        div['ov']  = dso.Sov.isel(nlat_u=nlat_S) - dso.Sov.isel(nlat_u=nlat_N)
        div['az']  = dso.Saz.isel(nlat_u=nlat_S) - dso.Saz.isel(nlat_u=nlat_N)
        div['ed']  = dso.Se.isel(nlat_u=nlat_S)  - dso.Se.isel(nlat_u=nlat_N)
        div['to']  = dso.St.isel(nlat_u=nlat_S)  - dso.St.isel(nlat_u=nlat_N)
        
        div['tov'] = fac*xr_linear_trend(dst.Sov.isel(nlat_u=nlat_S) - dst.Sov.isel(nlat_u=nlat_N))
        div['taz'] = fac*xr_linear_trend(dst.Saz.isel(nlat_u=nlat_S) - dst.Saz.isel(nlat_u=nlat_N))
        div['ted'] = fac*xr_linear_trend(dst.Se.isel(nlat_u=nlat_S)  - dst.Se.isel(nlat_u=nlat_N))
        div['tto'] = fac*xr_linear_trend(dst.St.isel(nlat_u=nlat_S)  - dst.St.isel(nlat_u=nlat_N))
        
    elif quant=='FW' and latS is not None:
        div['ov']  = dso.Fov.isel(nlat_u=nlat_S) - dso.Fov.isel(nlat_u=nlat_N)
        div['az']  = dso.Faz.isel(nlat_u=nlat_S) - dso.Faz.isel(nlat_u=nlat_N)
        div['ed']  = 0

        div['tov'] = fac*xr_linear_trend(dst.Fov.isel(nlat_u=nlat_S) - dst.Fov.isel(nlat_u=nlat_N))
        div['taz'] = fac*xr_linear_trend(dst.Faz.isel(nlat_u=nlat_S) - dst.Faz.isel(nlat_u=nlat_N))
        div['ted'] = 0

        div['to']  = (dso.Fov+dso.Faz).isel(nlat_u=nlat_S) - (dso.Fov+dso.Faz).isel(nlat_u=nlat_N)
        div['tto'] = fac*xr_linear_trend((dst.Fov+dst.Faz).isel(nlat_u=nlat_S) - (dst.Fov+dst.Faz).isel(nlat_u=nlat_N))
    return fluxes, div


def FW_summary_plot(quant):
    """ plots main budget terms for certain regions
    d/dt = -div + SFWF + diffusion
    """
    assert quant in ['FW', 'SALT']
    if quant=='FW':      q=0
    elif quant=='SALT':  q=1

    f = plt.figure(figsize=(6.4,2))
    for i in range(6):  f.text(.98-.9*i/5,.95, lat_bounds[i], va='center', ha='center', fontsize=8, color='grey')
    for i, (latS, latN) in enumerate(lat_bands):

        #region: define axes
        ax = f.add_axes([.1+.9*i/5,.05,.8/5,.8])
        for s in ['right', 'top', 'bottom']:  ax.spines[s].set_visible(False)
        ax.xaxis.set_ticks([])
        ax.axhline(c='k', lw=.5)
        if i==0:  ax.set_ylabel(f'{quant} fluxes in [{["Sv","kt/s"][q]}]')
        else:     ax.yaxis.set_ticklabels([])
        vert_lim = [(-2,2),(-2.5e7,2.5e7)][q]
        ax.set_xlim((-.25,1.2))
        ax.set_ylim(vert_lim)
        #endregion
        
        for s, sim in enumerate(['HIGH', 'LOW']):  

            # d/dt
            if quant=='SALT':
                ddtS_ctl, ddtS_rcp = get_ddt_SALT(sim=sim, latS=latS, latN=latN)
                ax.bar(x=s/10, height=ddtS_ctl, **vbkw, color=mct('C9', alpha=1-.3*s), label=r'$\Delta_t$S')
                ax.arrow(x=.05+s/10, y=0, dy=ddtS_rcp, dx=0)
                print(ddtS_ctl)
            
            # surface fluxes
            sfwf_mean, sfwf_trend = get_SFWF(sim=sim, quant=quant, latS=latS, latN=latN)
            ax.bar(x=.25+s/10 , height=sfwf_mean['PER'], **vbkw, color=mct('C4', 1-.3*s), label='P+E+R')
            ax.arrow(x=.3+s/10, y=sfwf_mean['PER'], dy=sfwf_trend['PER'], dx=0)

            # meridional flux divergences
            if latS==60:
                fluxes, div = get_fluxes(sim=sim, quant=quant, lat=latS, latS=None, latN=None)
                div['to'] = fluxes['to']
            else:
                fluxes, div = get_fluxes(sim=sim, quant=quant, lat=latS, latS=latS, latN=latN)
                ax.bar(x=.5+s/10, height=div['to'], **vbkw, color=mct('C3', 1-.3*s), label=f'-{nabla}{quant}')
                ax.arrow(x=.55+s/10, y=div['to'], dx=0, dy=div['tto'])
            
            # diff term
            if quant=='SALT':  diffusion = ddtS_ctl - div['to'] - sfwf_mean['PER']
            elif quant=='FW':  diffusion = div['to'] - sfwf_mean['PER']
            ax.bar(x=.75+s/10, height=diffusion, **vbkw, color=mct('C5', 1-.3*s), label=f'{quant}{diff}')
            
            # legend
            if i==2 and s==0:  ax.legend(fontsize=6)

    # plt.savefig(f'{path_results}/Mov/{["FW","SALT"][q]}_region_budget_total.eps')
    return


def FW_region_plot(quant):
    """ plots SFWF, merid. transport and its divergence for specific regions
    quant .. Fw or SALT
    """
    assert quant in ['FW', 'SALT']
    if quant=='FW':      q, labels =0,  [r'$-\nabla F_{oc}$', r'$F_{oc}$', r'$\int F_{atm}$']
    elif quant=='SALT':  q, labels = 1, [r'$-\nabla S_{oc}$', r'$S_{oc}$', r'$\int S_{atm}$']
    lat_bounds = ['90N', '60N', '45N', '10N', '10S', '34S']

    f = plt.figure(figsize=(6.4,3))

    # latitude labels
    for i in range(3):  f.text(.92, [.1,.5,.9][i], labels[i], va='center', fontsize=8, color='grey')
    for i in range(6):  f.text(.9-.8*i/5,.77, lat_bounds[i], va='center', ha='center', fontsize=8, color='grey')

    dd = {}
    for i, (latS, latN) in enumerate(lat_bands):
        print(latS, latN)
        #region: define axes
        ax = f.add_axes([.1+.8*i/5,.2,.8/5,.55])
        if i==4:  ax.spines['right'].set_visible(False)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])

        vert_lim = [(-3.5,3.5),(-3e7,3e7)][q]
        axb = f.add_axes([.1+.8*i/5,-.05,.8/5,.5])
        axb.set_xlim((-.25,1.2))
        axb.set_ylim(vert_lim)
        
        axt = f.add_axes([.1+.8*i/5,.5,.8/5,.5])
        axt.set_xlim((-.25,1))
        axt.set_ylim(vert_lim)
        
        axl = f.add_axes([-.025+.8*i/5,.2,.25,.55])
        axl.set_xlim([(-.5,.5),(-2e8,2e8)][q])
        axl.set_ylim((-.35,1.3))
        for ax in [axb, axt, axl]:
            ax.patch.set_alpha(0)
            ax.axis('off')
            ax.xaxis.set_ticks([])
        #endregion

        for s, sim in enumerate(['HIGH', 'LOW']):
            nlat_ = nlat_at(sim, latS)
            run, rcp = ['ctrl','lpd'][s], ['rcp', 'lr1'][s]

            # d/dt  [axb]
            if quant=='SALT':
                ddt = get_ddt_SALT(sim=sim, latS=latS, latN=latN)
                axb.bar(x=s/10, height=ddt, **vbkw, color=mct('C9', alpha=1-.3*s), label=r'$\Delta_t$S')

            #region: meridional flux (divergences) means/trends  [axb]
            if latS==60:
                fluxes, div = get_fluxes(sim=sim, quant=quant, lat=latS, latS=None, latN=None)
                div['ov'], div['az'], div['ed'], div['to'] = fluxes['ov'], fluxes['az'], fluxes['ed'], fluxes['to']
                div['tov'], div['taz'], div['ted'], div['tto'] = fluxes['tov'], fluxes['taz'], fluxes['ted'], fluxes['tto']
                
            else:
                fluxes, div = get_fluxes(sim=sim, quant=quant, lat=latS, latS=latS, latN=latN)

            axl.barh(y=0  +s/10, width=fluxes['ov'], **hbkw, color=mct('C0', 1-.3*s))
            axl.barh(y=.25+s/10, width=fluxes['az'], **hbkw, color=mct('C1', 1-.3*s))
            axl.barh(y=.5 +s/10, width=fluxes['ed'], **hbkw, color=mct('C2', 1-.3*s))
            axl.barh(y=.75+s/10, width=fluxes['to'], **hbkw, color=mct('C3', 1-.3*s))

            axl.arrow(x=fluxes['ov'], y=.05+s/10, dx=fluxes['tov'], dy=0)
            axl.arrow(x=fluxes['az'], y=.3 +s/10, dx=fluxes['taz'], dy=0)
            axl.arrow(x=fluxes['ed'], y=.55+s/10, dx=fluxes['ted'], dy=0)
            axl.arrow(x=fluxes['to'], y=.8 +s/10, dx=fluxes['tto'], dy=0)
            
            ax.bar(x=0  +s/10, height=div['ov'], **vbkw, color=mct('C0', 1-.3*s))
            ax.bar(x=.25+s/10, height=div['az'], **vbkw, color=mct('C1', 1-.3*s))
            ax.bar(x=.5 +s/10, height=div['ed'], **vbkw, color=mct('C2', 1-.3*s))              
            ax.bar(x=.75+s/10, height=div['to'], **vbkw, color=mct('C3', 1-.3*s))
        
            ax.arrow(x=.05+s/10, y=div['ov'], dx=0, dy=div['tov'])
            ax.arrow(x=.3 +s/10, y=div['az'], dx=0, dy=div['taz'])
            ax.arrow(x=.55+s/10, y=div['ed'], dx=0, dy=div['ted'])
            ax.arrow(x=.8 +s/10, y=div['to'], dx=0, dy=div['tto'])
            #endregion

            #region: surface fluxes [axt]
            sfwf_mean, sfwf_trend = get_SFWF(sim=sim, quant=quant, latS=latS, latN=latN)

            axt.bar(x=.55+s/10, height=-sfwf_mean['PE'] , **vbkw, color=mct('C6', 1-.3*s))
            axt.bar(x=.3+s/10 , height=-sfwf_mean['R']  , **vbkw, color=mct('C7', 1-.3*s))
            axt.bar(x=.0+s/10 , height=-sfwf_mean['PER'], **vbkw, color=mct('C4', 1-.3*s), label='P+E+R')

            axt.arrow(x=.6+s/10 , y=-sfwf_mean['PE'] , dx=0, dy=-sfwf_trend['PE'] )
            axt.arrow(x=.62+s/10, y=0                , dx=0, dy=-sfwf_trend['P']  )
            axt.arrow(x=.58+s/10, y=0                , dx=0, dy=-sfwf_trend['E']  )
            axt.arrow(x=.35+s/10, y=-sfwf_mean['R']  , dx=0, dy=-sfwf_trend['R']  )
            axt.arrow(x=.05+s/10, y=-sfwf_mean['PER'], dx=0, dy=-sfwf_trend['PER'])
            #endregion

            #  labels: PER, ov/az/eddy/total
            if i==0 and s==0:  
                axt.text(0.25,.95, 'P+E+R', transform=axt.transAxes, ha='center', va='top', rotation=90, fontsize=8)                
                axt.text(0.5 ,.95, 'R'    , transform=axt.transAxes, ha='center', va='top', rotation=90, fontsize=8)                
                axt.text(0.75,.95, 'P+E'  , transform=axt.transAxes, ha='center', va='top', rotation=90, fontsize=8)                
                axl.text([-.4,-1.2e8][q],.1 , ['F','S'][q]+r'$_{ov}$'   , va='center', fontsize=8)
                axl.text([-.4,-1.2e8][q],.35, ['F','S'][q]+r'$_{az}$'   , va='center', fontsize=8)
                axl.text([-.4,-1.2e8][q],.6 , ['',r'S$_{eddy}$'][q]     , va='center', fontsize=8)                
                axl.text([-.4,-1.2e8][q],.85, ['F','S'][q]+r'$_{total}$', va='center', fontsize=8)
    
    #region: axes legend
    axa = f.add_axes([0.04,.08,.125,.125])
    axa.set_xlim((0,[.5,2e2][q]))    
    axa.set_ylim((0,[3,1.5e2][q]/2))
    axa.text(.1,.1,['Sv\nSv/100yr','kt/s\nkt/s/100yr'][q], transform=axa.transAxes, fontsize=8)
    axa.patch.set_alpha(0)
    axa.spines['top'].set_visible(False)
    axa.spines['right'].set_visible(False)
    #endregion
    plt.savefig(f'{path_results}/Mov/{["FW","SALT"][q]}_region_budget.eps')
    return

