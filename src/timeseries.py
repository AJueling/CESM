import os
from analysis.paths import path_ocn_ctrl, path_ocn_rcp

class IterateOutputCESM:
    """ iterator over all CESM ctrl/rcp filenames
    automatically detects the last file
    
    example:
    >for filename in IterateOutputCESM('ctrl'):
    >    print(filename)
    """
    
    def __init__(self, run):
        self.run = run
        self.stop = False
        
        if run=='rcp':
            self.path  = path_ocn_rcp
            self.year  = 2000
            self.month = 1
        elif run=='ctrl':
            self.path  = path_ocn_ctrl
            self.year  = 200
            self.month = 1
        else:
            raise Exception('wrong run name, only rcp and ctrl are implemented')
            
    def file(self):
        time = f'{self.year:04}-{self.month:02}'
        if self.run=='rcp':
            filename = f'{self.path}/rcp8.5_co2_f05_t12.pop.h.{time}.nc'
        if self.run=='ctrl':
            filename = f'{self.path}/spinup_pd_maxcores_f05_t12.pop.h.{time}.nc'
            
        if os.path.exists(filename)==False:
            self.stop = True
        return filename
    
    def __iter__(self):
        return self

    def __next__(self):
        new_file = self.file()
        y = self.year
        m = self.month
        if self.stop:
            raise StopIteration
        else:
            if self.month==12:
                self.month = 1
                self.year +=1
            else:
                self.month += 1
            return y, m, new_file