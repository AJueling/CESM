import sys
sys.path.append("..")
from timeseries import yrly_avg_nc
opt = int(sys.argv[1])

if opt==1:
    domain = 'atm'
    run    = 'ctrl'
    fields = ['T', 'T850', 'U', 'V']
elif opt==2:
    domain = 'atm'
    run    = 'rcp'
    fields = ['T', 'T850', 'U', 'V']
    
elif opt==3:
    domain = 'ocn'
    run    = 'ctrl'
    fields = ['TEMP', 'PD']
elif opt==4:
    domain = 'ocn'
    run    = 'rcp'
    fields = ['TEMP', 'PD']
    
elif opt==5:
    domain = 'ocn_rect'
    run    = 'ctrl'
    fields = ['TEMP', 'PD']
elif opt==6:
    domain = 'ocn_rect'
    run    = 'rcp'
    fields = ['TEMP', 'PD']
elif opt==7:
    domain = 'ocn_rect'
    run    = 'ctrl'
    fields = ['SHF']
elif opt==8:
    domain = 'ocn_rect'
    run    = 'rcp'
    fields = ['SHF']
    
elif opt==9:
    domain = 'ocn'
    run    = 'ctrl'
    fields = ['UVEL', 'VVEL']
elif opt==10:
    domain = 'ocn'
    run    = 'rcp'
    fields = ['UVEL', 'VVEL']
    
    
yrly_avg_nc(domain=domain, run=run, fields=fields)