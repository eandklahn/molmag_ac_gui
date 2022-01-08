import json
from importlib.resources import read_text

import pandas as pd
import numpy as np

from utility import update_data_names
import data as pkg_static_data

class PpmsData:

    def __init__(self, file=None):
        
        self.acms_read_options = json.loads(read_text(pkg_static_data,
                                            'acms_read_options.json'))
        self.vsm_read_options = {}

        self.options = {'VSM':
                        {'validator': self._validate_as_vsm,
                         'attribs': self._set_vsm_attribs},
                        'ACMS':
                        {'validator': self._validate_as_acms,
                         'attribs': self._set_acms_attribs}
                        }

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
            self._clean_and_set_attribs_per_ftype()
            
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
        
        f = self.options[ftype]['validator']
        res = f(df)
        
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
    
    def _clean_and_set_attribs_per_ftype(self):

        self._cleanup_df()
        f = self.options[self.ftype]['attribs']
        f()

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

    def _fill_acms_data_values(self):
    
        if ('Xp (emu/Oe)' in self.df.columns and not ('Mp (emu)' in self.df.columns)):
            # Susceptibility exists in the data frame, but magnetisation does not
            Mp = self.df['Xp (emu/Oe)']*self.df['Magnetic Field (Oe)']
            Mpp = self.df['Xpp (emu/Oe)']*self.df['Magnetic Field (Oe)']
            Xp_idx = self.df.columns.get_loc('Xp (emu/Oe)')
            self.df.insert(Xp_idx, column='Mp (emu)', value=Mp)
            self.df.insert(Xp_idx+1, column='Mpp (emu)', value=Mpp)
            
        elif (not 'Xp (emu/Oe)' in self.df.columns and ('Mp (emu)' in self.df.columns)):
            # Magnetisation exists in the data frame, but susceptibility does not
            Xp = self.df['Mp (emu)']/self.df['Magnetic Field (Oe)']
            Xpp = self.df['Mpp (emu)']/self.df['Magnetic Field (Oe)']
            Mp_idx = self.df.columns.get_loc('Mp (emu)')
            self.df.insert(Mp_idx+2, column='Xp (emu/Oe)', value=Xp)
            self.df.insert(Mp_idx+3, column='Xpp (emu/Oe)', value=Xpp)

    def _set_vsm_attribs(self):

        pass

    def _set_acms_attribs(self):

        self._fill_acms_data_values()
        self.nfreq = len(set(self.df['AC Frequency (Hz)']))
        self.freq = np.array(sorted(list(set(self.df['AC Frequency (Hz)']))))
        self.update_temp_subsets()
        self.update_meas_temps()
    
    def _cleanup_df(self):
        
        # Drop columns where all values are NaN
        self.df.dropna(axis=1, how='all', inplace=True)
        # Removing instrument comment lines
        # Drop "Comment" column
        if 'Comment' in self.df.columns:
            self.df.drop(['Comment'], axis='columns', inplace=True)
        # Drop all rows where there is still a NaN value
        self.df.dropna(axis=0, inplace=True)
        
        # Make sure that the rows are named continuously
        old_indices = self.df.index.values
        new_indices = list(range(len(old_indices)))
        self.df.rename(index=dict(zip(old_indices, new_indices)),
                           inplace=True)

    def update_temp_subsets(self):
        
        self.T_subsets = []
        idx_list = [0]
        # Splits based on where the frequency "restarts"
        # Assumes that the frequency is always increasing within a measurement
        i=0
        old_val = 0
        while i<self.df.shape[0]:
            new_val = self.df['AC Frequency (Hz)'].iloc[i]
            if new_val<old_val:
                idx_list.append(i)
            else:
                pass
            old_val = new_val
            i+=1
        idx_list.append(self.df.shape[0])
        
        for n in range(len(idx_list)-1):
            self.T_subsets.append(self.df.iloc[idx_list[n]:idx_list[n+1]])

    def update_meas_temps(self):
        
        Ts = []
        for sub in self.T_subsets:
            Ts.append(sub['Temperature (K)'].mean())
        
        self.T = np.array(Ts)
        self.nT = len(self.T)
        self.Tmin = self.T.min()
        self.Tmax = self.T.max()

if __name__ == '__main__':
    
    ppmsdata = PpmsData(r'C:\Users\emilk\Documents\Uddannelse\PhD\python_library\process_ac_gui\ac_gui\prototyping\dy-dbm\20180209DyII_1000.dat')
    print(ppmsdata.freq)
