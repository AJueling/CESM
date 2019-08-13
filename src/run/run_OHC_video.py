import sys
sys.path.append("..")
import datetime

from cd_synthesis_OHC import SynthesizeOHC as SO

if __name__=="__main__":
    print(f'running OHC synthesis plots')
    print(f'{datetime.datetime.now()}\n\n')
    
    SO().plot_OHC_anomaly()
    
    print(f'\n\nfinished at\n{datetime.datetime.now()}')