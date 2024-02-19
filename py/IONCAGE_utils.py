__all__ = ["read_ioncage_output", "make_dlogD", "make_dlogv"]

import numpy as np
import sys

def find_occurrence_of_string(filename,string,listflag,first_idx):
# PURPOSE:
#     Looks for string in file 
# USAGE:
#     Given the path "filename" to a file, and a text "string" run
#     > find_occurrence_of_string(filename, string, listflag, first_idx)
#     where LISTFLAG controls if algorithm should search for 
#     first occurrence in file (listflag = "first") or list of all
#     occurrences in file, and first_idx controls the line number
#     of the first line (eg 1, 0 or something userdefined).
# OUTPUT:
#     Line number, or list of line numbers
# AUTHOR:
#     Jacob Svensmark, Feb. 2019

    import sys

    # Attempt to open file    
    try:
        fp = open(filename)
    except IOError:
        print("Could not read file:", filename)
        sys.exit()
    string_length = len(string)
    cnt = 0
    with fp:
        # Return line number of first occurrence in file
        if listflag == 'first':
            line = fp.readline()
            while not (string == line[0:string_length]):
                line = fp.readline()
                cnt += 1
            return cnt + first_idx
        # Return line numbers of all occurrences in file
        elif listflag == 'all':
            line_numbers = np.empty([1,0])
            for line in fp:
                if string == line[0:string_length]:
                    line_numbers = np.append(line_numbers,cnt)
                cnt +=1
            return line_numbers + first_idx
        # Raise error if listflag is wrong
        else:
            raise ValueError("Unknown listflag: " + listflag)
            sys.exit()


def read_ioncage_output(filename,**kwargs):

# PURPOSE:
#     Reads and returns output from the IONCAGE model
# USAGE:
#     Given the path of an IONCAGE outputfile, run
#     > t,d,v,N,n = read_ioncage_experiment(filename)
# OUTPUT:
#     t: Array of time values, one for each of the NS snapshots
#     d: Array of NTOT node diameters
#     v: Array of corresponding volume nodes
#     N: Array of size [NS,NTOT,3], where first index annotates time, second annotates 
#        aerosol node size and last index annotates neutral, positive and negative 
#        aerosols in that order.
#     n: Array of size [NS,3] where first index annotates time and  last index 
#        annotates neutral, positive and negative  aerosols. 
# AUTHOR:
#     Jacob Svensmark, Feb. 2019

    time_hack = True

    # Extract ntot:
    string = 'ntot'
    ntot_line = find_occurrence_of_string(filename,string,'first',0)
    string = ' --------------------------------------------' 
    snapshot_lines = find_occurrence_of_string(filename,string,'all',0)
    n_snapshots = len(snapshot_lines)

    string = ' //********************* Diameter/Volume nodes ********************//'
    d_line = find_occurrence_of_string(filename,string,'first',0)

    # Read NTOT value
    with open(filename) as fp:
        for i in range(0,20):
            line = next(fp).strip()
            if (i == ntot_line): # Read ntot value and allocate arrays
                ntot = int(line.rstrip().split()[1])

    # Allocate arrays
    d = np.zeros(ntot)
    v = np.zeros(ntot)
    t = np.zeros(n_snapshots)
    n = np.zeros([n_snapshots,3])
    N = np.zeros([ntot,n_snapshots,3])
    t_cnt = 0
    d_cnt = 0
    with open(filename) as fp:
        for i, line in enumerate(fp):
            if (i >= snapshot_lines[t_cnt]+12) and (i < snapshot_lines[t_cnt]+12+ntot): # Read Aerosols
                N[d_cnt,t_cnt,0],N[d_cnt,t_cnt,1],N[d_cnt,t_cnt,2] =  line.rstrip().split()[:3]
                d_cnt+=1
            elif (i == snapshot_lines[t_cnt]+9): # Read charged monomers
                n[t_cnt,1],n[t_cnt,2] = line.rstrip().split()[:2] 
            elif (i == snapshot_lines[t_cnt]+6): # Read neutral monomers
                n[t_cnt] = float(line)
            elif (i == snapshot_lines[t_cnt]+3): # Read time value 
                if time_hack==True:
                    if t_cnt < 2:
                        t[t_cnt] = float(line)
                    else:
                        t[t_cnt] = t[t_cnt-1]+t[1]
                else:
                    t[t_cnt] = float(line)
            elif (i >= d_line+2) and (i<d_line+2+ntot): # Read diameter and volume arrays
                d[i-(d_line+2)],v[i-(d_line+2)] = line.rstrip().split()[:2]
            elif (d_cnt==ntot) and (t_cnt<=n_snapshots): # Advance snapshot counter and reset d counter
                t_cnt+=1
                d_cnt=0
    return t,d,v,N,n    

def make_dlogD(v,N):

# PURPOSE:
#     Calculates the dN/dlogD derivative for IONCAGE output
# INPUT:
#     v: Array of Ntot volume values of length
#     N: Array of aerosol concentrations N[Ntot,NS,3], where NS is the number of snapshots (timesteps)
# RETURNS:
#    dN0/dlogD, dNp/dlogD, dNm/dlogD and the total dNT/dlogD, all og length Ntot
# AUTHOR:
#     Jacob Svensmark, Feb. 2019

    # Make Upper / lower bounds of bins
    ntot = len(v)
    aa  = np.zeros(ntot)
    bb  = np.zeros(ntot)
    fac = v[3]/v[2]
    for j in range(ntot):
        aa[j] = v[j]/np.sqrt(fac) # Logarithmic midpoints of bins
        bb[j] = v[j]*np.sqrt(fac)
    # Convert to diameters
    aa = np.log10(2.*(3./(4.*np.pi) * aa)**(1./3.))
    bb = np.log10(2.*(3./(4.*np.pi) * bb)**(1./3.))
    # Calculate dN / dlog(Dp) distribution
    nt = len(N[0,:,0])
    dN0_dlogD = np.zeros([ntot,nt])
    dNp_dlogD = np.zeros([ntot,nt])
    dNm_dlogD = np.zeros([ntot,nt])
    dNT_dlogD = np.zeros([ntot,nt])
    for i in range(nt):
        dN0_dlogD[:,i] = N[:,i,0] / (bb-aa)
        dNp_dlogD[:,i] = N[:,i,1] / (bb-aa)
        dNm_dlogD[:,i] = N[:,i,2] / (bb-aa)
        dNT_dlogD[:,i] = (N[:,i,0]+N[:,i,1]+N[:,i,2]) / (bb-aa)
    return  dN0_dlogD, dNp_dlogD, dNm_dlogD, dNT_dlogD


def make_dlogv(v,N):

# PURPOSE:
#     Calculates the d(N)/dlog(v) derivative for IONCAGE output
# INPUT:
#     v: Array of Ntot volume values of length
#     N: Array of aerosol concentrations N[Ntot,NS,3], where NS is the number of snapshots (timesteps)
#RETURNS:
#     dN0/dlogv, dNp/dlogv, dNm/dlogv and the total dNT/dlogv, all og length Ntot
# AUTHOR:
#     Jacob Svensmark, Feb. 2019

    # Make Upper / lower bounds of bins
    ntot = len(v)
    aa  = np.zeros(ntot)
    bb  = np.zeros(ntot)
    fac = v[3]/v[2]
    for j in range(ntot):
        aa[j] = v[j]/np.sqrt(fac) # Logarithmic midpoints of bins
        bb[j] = v[j]*np.sqrt(fac)
    # Calculate dN / dlog(Dp) distribution
    nt = len(N[0,:,0])
    dN0_dlogv = np.zeros([ntot,nt])
    dNp_dlogv = np.zeros([ntot,nt])
    dNm_dlogv = np.zeros([ntot,nt])
    dNT_dlogv = np.zeros([ntot,nt])
    for i in range(nt):
        dN0_dlogv[:,i] = N[:,i,0] / (bb-aa)
        dNp_dlogv[:,i] = N[:,i,1] / (bb-aa)
        dNm_dlogv[:,i] = N[:,i,2] / (bb-aa)
        dNT_dlogv[:,i] = (N[:,i,0]+N[:,i,1]+N[:,i,2]) / (bb-aa)
    return  dN0_dlogv, dNp_dlogv, dNm_dlogv, dNT_dlogv
