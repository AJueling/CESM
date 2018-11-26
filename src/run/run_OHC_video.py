import sys
sys.path.append("..")
import datetime

run     = sys.argv[1]

from OHC_plots import OHC_vert_video

if __name__=="__main__":
    print(f'running CICE_XMXL plots run={run}')
    print(f'{datetime.datetime.now()}\n\n')
    
    OHC_vert_video(run=run)
    
    print(f'\n\nfinished at\n{datetime.datetime.now()}')