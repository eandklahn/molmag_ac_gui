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
        header_dict = {}
        for h in header:
            try:
                header_dict[h[2]] = h[1]
            except IndexError:
                continue
        header = header_dict
        
        df = pd.read_csv(filename,
                         header=data_start,
                         engine='python',
                         skip_blank_lines=False)
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

def formatlabel(label): 
    #Lav dictonary: Mp fx som key, M' som value 
    
    if label == "Mp (emu)": 
        return "M' (emu)"
    elif label == "Mpp (emu)": 
        return "M'' (emu)"
    elif label == "Xp (emu/Oe)": 
        return r"$\chi'$ (emu/Oe)"
    elif label == "Xpp (emu/Oe)": 
        return r"$\chi''$ (emu/Oe)"
    elif label == "Mp_m (emu/mol)": 
        return r"$M_m'$ (emu/Oe)"
    elif label == "Mpp_m (emu/mol)":
        return r"$M_m''$ (emu/Oe)"
    elif label ==  "Xp_m (emu/(Oe*mol))":
        return r"$\chi_m'$ $(emu/(Oe⋅mol))$"
    elif label ==  "Xpp_m (emu/(Oe*mol))":
        return r"$\chi_m''$ $(emu/(Oe⋅mol))$"
    else: 
        return label


if __name__ == '__main__':
    
    filename = input()
    h, df = read_ppms_file(filename)
    print(h)
