import json
from importlib.resources import read_text

import pandas as pd

from utility import update_data_names
import data as pkg_static_data

class PpmsData:

    def __init__(self, file=None):
        
        self.acms_read_options = json.loads(read_text(pkg_static_data,
                                            'acms_read_options.json'))
        self.vsm_read_options = {}

        self.filename = file
        self.df = None
        self.header = None
        self.ftype = None
        self.read_into_df(self.filename)

    def read_into_df(self, file=None):
        
        try:
            assert not file is None
            header, df, ftype = self.read_ppms_file(file)
            assert not header is None
            res = self._validate_per_ftype(df, ftype)
            print(df)
            assert res
        except AssertionError:
            pass
        except FileNotFoundError:
            pass
        except UnicodeDecodeError:
            pass
        else:
            self.filename = file
            self.header = header
            self.df = df
            self.ftype = ftype
            
    def read_ppms_file(self, filename):
        
        header, df, ftype = None, None, None

        with open(filename, 'r') as f:
            d = f.readlines()

        head_flag, head_idx, data_flag, data_idx = self._locate_header_and_data(d)

        if (head_flag and data_flag):
            
            header = self._read_header(d, head_idx, data_idx)
            df = pd.read_csv(filename,
                             header=data_idx,
                             engine='python',
                             skip_blank_lines=False)
            ftype = self._read_ftype(d, head_idx)

        return header, df, ftype

    def _validate_per_ftype(self, df, ftype):
        
        res = False
        if ftype == 'VSM':
            res = self._validate_as_vsm(df)
        elif ftype == 'ACMS':
            res = self._validate_as_acms(df)
        return res

    def _validate_as_vsm(self, df):

        return True

    def _validate_as_acms(self, df):

        # 1 to make sure that none of the names in read_options were matched
        # more than once.
        # 2 to make sure that only Mp (and therefore Mpp) OR Xp (and therefore
        # Xpp) can appear at the same time. In the case that this is ever an
        # error, self.fill_df_data_values will have to be changed.
        try:
            summary = update_data_names(df, self.acms_read_options)
            counts = [val>1 for key, val in summary.items()]
            assert not any(counts) #1
            assert (summary['Mp (emu)']>0) != (summary['Xp (emu/Oe)']>0) #2
        except AssertionError:
            return False
        else:
            return True
    
    def _read_ftype(self, d, head_idx):

        ftype = d[head_idx].strip().split()[1]
        return ftype

    def _read_header(self, d, head_idx, data_idx):
        
        header = d[head_idx:data_idx-1]
        
        header = [h.strip().split(',') for h in header 
                  if not h.startswith(';') and h.startswith('INFO')]
        header_dict = {}
        for h in header:
            try:
                header_dict[h[2]] = h[1]
            except IndexError:
                continue
        header = header_dict
        
        return header      

    def _locate_header_and_data(self, d):
        
        head_flag, data_flag = False, False
        head_idx, data_idx = 0, 0

        for i, line in enumerate(d):
            if '[Header]' in line:
                head_idx = i+1
                head_flag = True
            elif '[Data]' in line:
                data_idx = i+1
                data_flag = True

        return head_flag, head_idx, data_flag, data_idx



#    def load_ppms_data(self):
#        
#            self.fill_df_data_values()
#            self.cleanup_loaded_ppms()
#            self.num_meas_freqs = len(set(self.raw_df['AC Frequency (Hz)']))
#            self.update_temp_subsets()
#            self.update_meas_temps()

if __name__ == '__main__':
    
    df = PpmsData(r'C:\Users\emilk\Documents\Uddannelse\PhD\python_library\process_ac_gui\ac_gui\prototyping\dy-dbm\20180209DyII_1000.dat')
