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
nabla = r'$\nabla$'
xrkw = {'decode_times':False} 
diff = r'$_{diff}$'

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


def FW_summary_plot(quant):
    """ plots main budget terms for certain regions
    d/dt = -div + SFWF + diffusion
    """
    assert quant in ['FW', 'SALT']
    if quant=='FW':      q=0
    elif quant=='SALT':  q=1

    salt_ctrl = xr.open_dataset(f'{path_results}/SALT/SALT_integrals_ctrl.nc', decode_times=False)
    salt_rcp  = xr.open_dataset(f'{path_results}/SALT/SALT_integrals_rcp.nc' , decode_times=False)
    salt_lpd  = xr.open_dataset(f'{path_results}/SALT/SALT_integrals_lpd.nc' , decode_times=False)
    salt_lr1  = xr.open_dataset(f'{path_results}/SALT/SALT_integrals_lr1.nc' , decode_times=False)

    f = plt.figure(figsize=(6.4,2))
    for i in range(6):  f.text(.98-.9*i/5,.95, lat_bounds[i], va='center', ha='center', fontsize=8, color='grey')
    for i, (latS, latN) in enumerate(lat_bands):

        # ss
        ax = f.add_axes([.1+.9*i/5,.05,.8/5,.8])
        for s in ['right', 'top', 'bottom']:  ax.spines[s].set_visible(False)
        ax.xaxis.set_ticks([])
        ax.axhline(c='k', lw=.5)
        if i==0:
            ax.set_ylabel(f'{quant} fluxes in [{["Sv","kt/s"][q]}]')
        else:
            ax.yaxis.set_ticklabels([])
        vert_lim = [(-2,2),(-2.5e7,2.5e7)][q]
        ax.set_xlim((-.25,1.2))
        ax.set_ylim(vert_lim)
        
        for s, sim in enumerate(['HIGH', 'LOW']):
            run, rcp = ['ctrl','lpd'][s], ['rcp', 'lr1'][s]

            # d/dt
            if quant=='SALT':
                n1 = f'SALT_0-1000m_timeseries_{latS}N_{latN}N'
                n2 = f'SALT_below_1000m_timeseries_{latS}N_{latN}N'
                salt = [salt_ctrl[n1]+salt_ctrl[n2], salt_lpd[n1]+salt_lpd[n2]][s]
                ddt = (salt.isel(time=30)-salt.isel(time=0)).values/30/365/24/3600
                ax.bar(x=s/10, y=0, width=.1, height=ddt, label=r'$\Delta_t$S')
            
            # surface fluxes
            d = load_obj(f'{path_results}/SFWF/Atlantic_SFWF_integrals_{sim}_{latS}N_{latN}N')
            if quant=='SALT':  fac = -.00347e9  # 1e9 kg/s/Sv * salinity_factor of POP2
            elif quant=='FW':  fac = 1.
            sfwf = (d[f'Tmi_Sv'])*fac
            ax.bar(x=.25+s/10 , height=sfwf, width=.1,
                   align='edge', color=mct(color='C3', alpha=1-.3*s),
                   label='P+E+R')
            ax.arrow(x=.3+s/10, y=sfwf, dy=d[f'Tti_Sv']*fac, dx=0)

            # meridional flux divergences
            dso = xr.open_dataset(f'{path_prace}/Mov/FW_SALT_fluxes_{run}.nc', **xrkw)
            dst = xr.open_dataset(f'{path_prace}/Mov/FW_SALT_fluxes_{rcp}.nc', **xrkw)
            if latN!=90:  
                nlat_S, nlat_N = nlat_at(sim, latS), nlat_at(sim, latN)
                if quant=='SALT':
                    to  = dso.St.isel(nlat_u=nlat_S) - dso.St.isel(nlat_u=nlat_N)
                    tto = dst.St.isel(nlat_u=nlat_S) - dst.St.isel(nlat_u=nlat_N)
                elif quant=='FW':
                    to  = (dso.Fov+dso.Faz).isel(nlat_u=nlat_S) - (dso.Fov+dso.Faz).isel(nlat_u=nlat_N)
                    tto = (dst.Fov+dst.Faz).isel(nlat_u=nlat_S) - (dst.Fov+dst.Faz).isel(nlat_u=nlat_N)
                ax.bar(x=.5+s/10, y=0, height=to, width=.1, align='edge',\
                       color=mct(color='C6', alpha=.9-.3*s), label=f'-{nabla}{quant}')
#                 ax.arrow(x=.55+s/10, y=to, dx=0, dy=xr_linear_trend(tto).values*365*100)
            
            
            # diff term
            if quant=='SALT':  diffusion = ddt - to - sfwf
            elif quant=='FW':  diffusion = to - sfwf
            ax.bar(x=.75+s/10, y=0, height=diffusion, width=.1, align='edge',\
                   color=mct(color='C7', alpha=.9-.3*s), label=f'{quant}{diff}')
            
            if i==2 and s==0:  ax.legend(fontsize=8)

    plt.savefig(f'{path_results}/Mov/{["FW","SALT"][q]}_region_budget_total.eps')
    return


def FW_region_plot(quant):
    """ plots SFWF, merid. transport and its divergence for specific regions
    quant .. Fw or SALT
    """
    assert quant in ['FW', 'SALT']
    if quant=='FW':
        q=0
        labels = [r'$-\nabla F_{oc}$', r'$F_{oc}$', r'$\int F_{atm}$']
    elif quant=='SALT':
        q=1
        labels = [r'$-\nabla S_{oc}$', r'$S_{oc}$', r'$\int S_{atm}$']
    lat_bounds = ['90N', '60N', '45N', '10N', '10S', '34S']

    f = plt.figure(figsize=(6.4,3))
    for i in range(3):  f.text(.92, [.1,.5,.9][i], labels[i], va='center', fontsize=8, color='grey')
    for i in range(6):  f.text(.9-.8*i/5,.77, lat_bounds[i], va='center', ha='center', fontsize=8, color='grey')

    dd = {}
    for i, (latS, latN) in enumerate(lat_bands):
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

        for s, sim in enumerate(['HIGH', 'LOW']):
            nlat_ = nlat_at(sim, latS)
            run, rcp = ['ctrl','lpd'][s], ['rcp', 'lr1'][s]
            dso = xr.open_dataset(f'{path_prace}/Mov/FW_SALT_fluxes_{run}.nc',\
                                  decode_times=False).mean(dim='time').isel(nlat_u=nlat_)
            dst = xr.open_dataset(f'{path_prace}/Mov/FW_SALT_fluxes_{rcp}.nc',\
                                  decode_times=False).isel(nlat_u=nlat_)

            # meridional flux means
            if quant=='SALT':  ov, az, ed, to = dso.Sov, dso.Saz, dso.Se, dso.Sov+dso.Saz+dso.Se
            elif quant=='FW':  ov, az, to = dso.Fov, dso.Faz, dso.Fov+dso.Faz
            axl.barh(y=0  +s/10, width=ov, height=.1, align='edge', color=mct(color='C4', alpha=1-.3*s))
            axl.barh(y=.25+s/10, width=az, height=.1, align='edge', color=mct(color='C5', alpha=1-.3*s))
            axl.barh(y=.75+s/10, width=to, height=.1, align='edge', color=mct(color='C6', alpha=1-.3*s))
            if quant=='SALT':  axl.barh(y=.5 +s/10, width=ed, height=.1, align='edge', color=mct(color='C7', alpha=.9-.3*s))

            # meridional flux trends
            if quant=='SALT':  tov, taz, ted, tto = dst.Sov, dst.Saz, dst.Se, dst.Sov+dst.Saz+dst.Se
            elif quant=='FW':  tov, taz, tto = dst.Fov, dst.Faz, dst.Fov+dst.Faz
            axl.arrow(x=ov, y=.05+s/10, dx=xr_linear_trend(tov)*365*100, dy=0)
            axl.arrow(x=az, y=.3 +s/10, dx=xr_linear_trend(taz)*365*100, dy=0)
            axl.arrow(x=to, y=.8 +s/10, dx=xr_linear_trend(tto)*365*100, dy=0)
            if quant=='SALT':  axl.arrow(x=ed, y=.55 +s/10, dx=xr_linear_trend(ted)*365*100, dy=0)
            
            # divergences
            if i==0:  # temporarily store southern boundary values
                dd[f'ov_{sim}'], dd[f'az_{sim}'], dd[f'to_{sim}'], ax_ = ov.values, az.values, to.values, axb
                dd[f'tov_{sim}'], dd[f'taz_{sim}'], dd[f'tto_{sim}'] = tov.copy(), taz.copy(), tto.copy()
                if quant=='SALT':  dd[f'ed_{sim}'], dd[f'ted_{sim}'] = ed.values, ted.copy()
            else:
                ax_.bar(x=0  +s/10, y=0, height=dd[f'ov_{sim}']-ov.values,\
                        width=.1, align='edge', color=mct(color='C4', alpha=.9-.3*s))
                ax_.bar(x=.25+s/10, y=0, height=dd[f'az_{sim}']-az.values,\
                        width=.1, align='edge', color=mct(color='C5', alpha=.9-.3*s))
                ax_.bar(x=.75+s/10, y=0, height=dd[f'to_{sim}']-to.values,\
                        width=.1, align='edge', color=mct(color='C6', alpha=.9-.3*s))
                if quant=='SALT':  ax_.bar(x=.5 +s/10, y=0, height=dd[f'ed_{sim}']-ed.values,\
                                           width=.1, align='edge', color='C7', alpha=.9-.3*s)
            
                ax_.arrow(x=.05+s/10, y=dd[f'ov_{sim}']-ov.values,\
                          dx=0, dy=xr_linear_trend(dd[f'tov_{sim}']-tov).values*365*100)
                ax_.arrow(x=.3 +s/10, y=dd[f'az_{sim}']-az.values,\
                          dx=0, dy=xr_linear_trend(dd[f'taz_{sim}']-taz).values*365*100)
                ax_.arrow(x=.8 +s/10, y=dd[f'to_{sim}']-to.values,\
                          dx=0, dy=xr_linear_trend(dd[f'tto_{sim}']-tto).values*365*100)
                if quant=='SALT':  ax_.arrow(x=.55+s/10, y=dd[f'ed_{sim}']-ed.values,\
                                             dx=0, dy=xr_linear_trend(dd[f'ted_{sim}']-ted).values*365*100)
                
                dd[f'ov_{sim}'], dd[f'az_{sim}'], dd[f'to_{sim}'], ax_ = ov.values, az.values, to.values, axb
                dd[f'tov_{sim}'], dd[f'taz_{sim}'], dd[f'tto_{sim}'] = tov.copy(), taz.copy(), tto.copy()
                if quant=='SALT':  dd[f'ed_{sim}'], dd[f'ted_{sim}'] = ed.values, ted.copy()

            # surface fluxes
            d = load_obj(f'{path_results}/SFWF/Atlantic_SFWF_integrals_{sim}_{latS}N_{latN}N')
            if quant=='SALT':  fac = -.00347e9  # 1e9 kg/s/Sv * salinity_factor of POP2

                # axt.bar(x=.05+s/10 , height=-d[f'Smi_kgs'], width=.1, align='edge', color='C2', alpha=1-.3*s)
                # axt.arrow(y=-d[f'Smi_kgs'], x=.05+s/10, dy=-d[f'Sti_kgs'], dx=0)
            elif quant=='FW':  fac = 1.
            axt.bar(x=.55+s/10, height=-(d[f'Pmi_Sv']+d[f'Emi_Sv'])*fac,\
                    width=.1, align='edge', color=mct(color='C1', alpha=1-.3*s))
            axt.bar(x=.3+s/10 , height=-(d[f'Rmi_Sv']             )*fac,\
                    width=.1, align='edge', color=mct(color='C2', alpha=1-.3*s))
            axt.bar(x=.0+s/10 , height=-(d[f'Tmi_Sv']             )*fac,\
                    width=.1, align='edge', color=mct(color='C3', alpha=1-.3*s))
            axt.arrow(y=-(d[f'Pmi_Sv']+d[f'Emi_Sv'])*fac, x=.6+s/10 , dy=-(d[f'Pti_Sv']+d[f'Eti_Sv'])*fac, dx=0)
            axt.arrow(y=0                               , x=.62+s/10, dy=-d[f'Pti_Sv']*fac               , dx=0)
            axt.arrow(y=0                               , x=.58+s/10, dy=-d[f'Eti_Sv']*fac               , dx=0)
            axt.arrow(y=-d[f'Rmi_Sv']*fac               , x=.35+s/10, dy=-d[f'Rti_Sv']*fac               , dx=0)
            axt.arrow(y=-d[f'Tmi_Sv']*fac               , x=.05+s/10, dy=-d[f'Tti_Sv']*fac               , dx=0)
         
            if i==0 and s==0:  #  labels: PER, ov/az/eddy/total
                axt.text(0.25,.95, 'P+E+R', transform=axt.transAxes, ha='center', va='top', rotation=90, fontsize=8)                
                axt.text(0.5 ,.95, 'R'    , transform=axt.transAxes, ha='center', va='top', rotation=90, fontsize=8)                
                axt.text(0.75,.95, 'P+E'  , transform=axt.transAxes, ha='center', va='top', rotation=90, fontsize=8)                
                axl.text([-.4,-1.2e8][q],.1 , ['F','S'][q]+r'$_{ov}$'   , va='center', fontsize=8)
                axl.text([-.4,-1.2e8][q],.35, ['F','S'][q]+r'$_{az}$'   , va='center', fontsize=8)
                axl.text([-.4,-1.2e8][q],.6 , ['',r'S$_{eddy}$'][q]     , va='center', fontsize=8)                
                axl.text([-.4,-1.2e8][q],.85, ['F','S'][q]+r'$_{total}$', va='center', fontsize=8)
    
    # legend
    axa = f.add_axes([0.04,.08,.125,.125])
    axa.set_xlim((0,[.5,2e2][q]))    
    axa.set_ylim((0,[3,1.5e2][q]/2))
    axa.text(.1,.1,['Sv\nSv/100yr','kt/s\nkt/s/100yr'][q], transform=axa.transAxes, fontsize=8)
    axa.patch.set_alpha(0)
    axa.spines['top'].set_visible(False)
    axa.spines['right'].set_visible(False)
    plt.savefig(f'{path_results}/Mov/{["FW","SALT"][q]}_region_budget.eps')
#     ax.arrow(x=0,dx=3)
    return

