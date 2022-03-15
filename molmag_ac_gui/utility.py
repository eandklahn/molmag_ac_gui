import pandas as pd
label_dict = {
        "Mp (emu)": "M' (emu)",
        "Mpp (emu)": "M'' (emu)", 
        "Mp_m (emu/mol)": r"$M_m'$ (emu/Oe)",
        "Mpp_m (emu/mol)": r"$M_m''$ (emu/Oe)",
        "Xp (emu/Oe)": r"$\chi'$ (emu/Oe)", 
        "Xpp (emu/Oe)": r"$\chi''$ (emu/Oe)", 
        "Xp_m (emu/(Oe*mol))": r"$\chi_m'$ $(emu/(Oe⋅mol))$", 
        "Xpp_m (emu/(Oe*mol))": r"$\chi_m''$ $(emu/(Oe⋅mol))$", 
        "XT_m (emu*K/Oe)": r"$\chi_m⋅T$ $(emu⋅K/Oe)$",
        "Moment (emu)":  "M (emu)", 
        "Moment_m (emu/mol)":  r"$M_m$ (emu/mol)", 
        "X_m (emu/(Oe*mol))": r"$\chi_m$ $(emu/(Oe⋅mol))$"
        }

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
        pre_header = d[:data_start]
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
        pre_header, header, df = None, None, None
        
    return pre_header, header, df

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

    try: 
        return label_dict[label]
    except KeyError: #If not in label dict
        return label



if __name__ == '__main__':
    
    filename = input()
    h, df = read_ppms_file(filename)
