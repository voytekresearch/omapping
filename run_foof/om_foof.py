


def foof_par(psd):
    """   """

    # Fit FOOF
    foof.model(freqs_ext, psd)
    
    # Store vals in tuple and return
    return (foof.chi_, foof.centers_, foof.powers_, foof.stdevs_)


def clean_file_list(files_in, string):
    """"   """
    
    files_out = []

    for i in range(0, len(files_in)):
        if(string in files_in[i]):
            files_out.append(files_in[i])
            
    return files_out


def extract_psd(psd, freqs, f_low, f_high):
    """   """

    # Drop frequencies below f_low
    f_low_mask = freqs > f_low
    freqs_ext = freqs[f_low_mask]
    psd_ext = psd[:, f_low_mask]

    # Drop frequencies above f_high
    f_high_mask = freqs_ext < f_high
    freqs_ext = freqs_ext[f_high_mask]
    psd_ext = psd_ext[:, f_high_mask]

    return psd_ext, freqs_ext
