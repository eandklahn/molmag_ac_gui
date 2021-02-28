import pandas as pd

def read_ppms_file(filename):
    
    f = open(filename, 'r')
    d = f.readlines()
    f.close()
    
    found_header, found_data = False, False
    header_start, data_start = 0,0
    for i, line in enumerate(d):
        if '[Header]' in line:
            header_start = i+1
            found_header = True
        elif '[Data]' in line:
            data_start = i+1
            found_data = True
    
    if (found_header and found_data):
        header = d[header_start:data_start-1]
        header = [h.strip().split(',') for h in header if not h.startswith(';') and h.startswith('INFO')]
        header = {h[2]: h[1] for h in header}
        
        df = pd.read_csv(filename,
                        header=data_start,
                        engine='python')
    else:
        header, df = None, None
    
    return header, df

def get_ppms_column_name_matches(columns, options):

    matches = [x in columns for x in options]
    count = matches.count(True)
    if count>0:
        idx = matches.index(True)
        name = options[idx]
    else:
        name = None
    
    return count, name
   
if __name__ == '__main__':
    
    filename = input()
    h, df = read_ppms_file(filename)
    print(h)
    
def update_data_names(df, options):
    """
    # This function is supposed to update the names of the columns in raw_df, so that
    # the names conform to a standard to be used programwide.
    """
    
    summary = {}
    for key, val in options.items():
        
        count, name = get_ppms_column_name_matches(df, val)
        if count>0:
            df.rename(columns={name: key}, inplace=True)
        
        summary[key] = count
    
    return summary
    
    
        
    # First check that all data names can be found

    
    ## AC Frequency
    #freq_count, freq_name = get_ppms_column_name_matches(df, options['Frequency'])
    #if freq_count>0:
    #    self.raw_df.rename(columns={freq_name: "AC Frequency (Hz)"}, inplace=True)
    #
    ## AC Field
    #acfield_count, acfield_name = get_ppms_column_name_matches(self.raw_df.columns,
    #                                                           self.read_options['AC Field'])
    #if acfield_count>0:
    #    self.raw_df.rename(columns={acfield_name: "AC Amplitude (Oe)"}, inplace=True)
    #
    ## Magnetic field
    #magfield_count, magfield_name = get_ppms_column_name_matches(self.raw_df.columns,
    #                                                             self.read_options['Magfield'])
    #if magfield_count>0:
    #    self.raw_df.rename(columns={magfield_name: "Magnetic Field (Oe)"}, inplace=True)
    #    
    ## Mp
    #Mp_count, Mp_name = get_ppms_column_name_matches(self.raw_df.columns,
    #                                                 self.read_options['Mp'])
    #if Mp_count>0:
    #    self.raw_df.rename(columns={Mp_name: "Mp (emu)"}, inplace=True)
    #
    ## Mpp
    #Mpp_count, Mpp_name = get_ppms_column_name_matches(self.raw_df.columns,
    #                                                   self.read_options['Mpp'])
    #if Mpp_count>0:
    #    self.raw_df.rename(columns={Mpp_name: "Mpp (emu)"}, inplace=True)
    #
    ## Xp
    #Xp_count, Xp_name = get_ppms_column_name_matches(self.raw_df.columns,
    #                                                 self.read_options['Xp'])
    #if Xp_count>0:
    #    self.raw_df.rename(columns={Xp_name: "Xp (emu/Oe)"}, inplace=True)        
    #
    ## Mpp
    #Xpp_count, Xpp_name = get_ppms_column_name_matches(self.raw_df.columns,
    #                                                   self.read_options['Xpp'])
    #if Xpp_count>0:
    #    self.raw_df.rename(columns={Xpp_name: "Xpp (emu/Oe)"}, inplace=True)
    