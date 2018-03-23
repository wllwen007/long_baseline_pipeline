#!/usr/bin/env python
import pyvo as vo
import os
import glob
import numpy as np
import argparse
import pyrap.tables as pt
from astropy.table import Column


## based on script by Leah Morabito
## updated by Kaspars and Atvars to work for lbcs
## helper functions taken from download_tgss_skymodel_target.py in prefactor pipeline

########################################################################
def grab_coo_MS(MS):
    """
    Read the coordinates of a field from one MS corresponding to the selection given in the parameters

    Parameters
    ----------
    MS : str
        Full name (with path) to one MS of the field

    Returns
    -------
    RA, Dec : "tuple"
        coordinates of the field (RA, Dec in deg , J2000)
    """

    # reading the coordinates ("position") from the MS
    # NB: they are given in rad,rad (J2000) 
    [[[ra,dec]]] = pt.table(MS+'/FIELD', readonly=True, ack=False).getcol('PHASE_DIR')

    # RA is stocked in the MS in [-pi;pi]
    # => shift for the negative angles before the conversion to deg (so that RA in [0;2pi])
    if ra<0:
        ra=ra+2*np.pi

    # convert radians to degrees
    ra_deg =  ra/np.pi*180.
    dec_deg = dec/np.pi*180.

    # and sending the coordinates in deg
    return ra_deg,dec_deg

########################################################################
def input2strlist_nomapfile(invar):
    """ from bin/download_IONEX.py
    give the list of MSs from the list provided as a string
    """

    str_list = None
    if type(invar) is str:
        if invar.startswith('[') and invar.endswith(']'):
            str_list = [f.strip(' \'\"') for f in invar.strip('[]').split(',')]
        else:
            str_list = [invar.strip(' \'\"')]
    elif type(invar) is list:
        str_list = [str(f).strip(' \'\"') for f in invar]
    else:
        raise TypeError('input2strlist: Type '+str(type(invar))+' unknown!')
    return str_list

def main(ms_input, ResultsFile, catalogue, Radius=1.5, DoDownload="True"):

    """
    Download the source list for the target field and 
        for lbcs find things with good P values
        for lotss (preliminary release) find things that are unresolved

    Parameters
    ----------
    ms_input : str
        String from the list (map) of the target MSs
    ResultsFile : str
        Full name (with path) to the skymodel; if YES is true, the skymodel will be downloaded here
    catalogue : str
        Name of the catalogue to use: supported lobos / lotss
    Radius : float (default = 1.5)
        Radius for the cone search in degrees
    DoDownload : str ("Force" or "True" or "False")
        Download or not the LBCS catalogue.
        "Force": download catalogue, delete existing if needed.
        "True" or "Yes": use existing file if it exists, download if it does not.
        "False" or "No": Do not download catalogue, raise an exception if file does not exist.
    
    """

    FileExists = os.path.isfile(ResultsFile)
    if (not FileExists and os.path.exists(ResultsFile)):
        raise ValueError("download_lbcs_catalogue: WTF! Path: \"%s\" exists but is not a file!"%(ResultsFile))
    download_flag = False
    if DoDownload.upper() == "FORCE":
        if FileExists:
            os.remove(ResultsFile)
        download_flag = True
    elif DoDownload.upper() == "TRUE" or DoDownload.upper() == "YES":
        if FileExists:
            print "USING the exising catalogue in "+ ResultsFile
            return
        else:
            download_flag = True
    elif DoDownload.upper() == "FALSE" or DoDownload.upper() == "NO":
         if FileExists:
            print "USING the exising catalogue in "+ ResultsFile
            return
         else:
            raise ValueError("download_lbcs_catalogue: Path: \"%s\" does not exist and LBCS download is disabled!"%(ResultsFile))


    # If we got here, then we are supposed to download the skymodel.
    assert download_flag == True # Jaja, belts and suspenders...
    print "DOWNLOADING LBCS catalogue to "+ ResultsFile

    # Reading a MS to find the coordinate (pyrap)
    RATar, DECTar=grab_coo_MS(input2strlist_nomapfile(ms_input)[0])
    mypos = [( float(RATar), float(DECTar) )]
    #mypos = grab_coo_MS(input2strlist_nomapfile(ms_input)[0])

    if catalogue == 'lobos':
        ## this is the lbcs database to query
        url = 'http://vo.astron.nl/lbcs/lobos/cone/scs.xml'
    
    elif catalogue == 'lotss':
        ## this is the tier 1 database to query
        # note - this is the preliminary release
        url = 'http://vo.astron.nl/lofartier1/q/cone/scs.xml'
    else:
        raise ValueError('unknkown catalogue specified:', catalogue)

    ## this works
    query = vo.dal.scs.SCSQuery( url )
    query['RA'] = float( RATar )
    query['DEC'] = float( DECTar )
    query.radius = float( Radius )
    t = query.execute()   
    ## this does not
    #t = vo.conesearch( url, pos=mypos, radius=float(Radius) )

    tb = t.votable.to_table()


    if catalogue == 'lobos':
        #### Workaround to sort and pick good calibrator info from tb array ###########

        counts=[]
        P_count  = 0                        
        for i in tb:
            b = i[5].count('P')      #### Count of 'P' - good baselines       
            counts.append(b)
            if b >=2:
                P_count = P_count + 1    #### To determine how many sources to take 
        print 'Good sources - ' + str(P_count)   
        inds = np.argsort(counts)
        tb_sorted =tb[inds[::-1]]
        len_array = len(tb_sorted)
        for i in range ((len_array-P_count)):
            len_array-=1
            tb_sorted.remove_row(len_array)

        ## remove duplicates
        tb_tmp = np.array( tb_sorted['raj2000','decj2000'] )
        result = [ idx for idx, item in enumerate( tb_tmp ) if item in tb_tmp[:idx] ]
        tb_sorted.remove_rows(result)
        
        ############ add dummy flux (spec ind columns) ################################
        
        tb_sorted.add_column(Column(np.ones(len(tb_sorted)), 'Total_flux'))
        tb_sorted.add_column(Column(np.nan*np.zeros(len(tb_sorted)), 'SpecInd'))
        
        ############ Print the new table with standard column names in it  ########################
        
        tb_out = tb_sorted['raj2000', 'decj2000', 'Total_flux', 'SpecInd', 'ObsID'] 
        tb_out.rename_column('raj2000', 'RA')
        tb_out.rename_column('decj2000', 'DEC')
        tb_out.rename_column('ObsID', 'ID')
        tb_out.write( ResultsFile, format='ascii.csv')    
   

    elif catalogue == 'lotss':

        ## find unresolved
        nsrcs = float( len( tb.columns['Resolved'] ) )
        unresolved_index = np.where( tb.columns['Resolved'] == 'U' )[0]
        perc_unres = len( unresolved_index ) / nsrcs * 100.
        print 'Percentage of sources which are unresolved: '+str( perc_unres )
        utb = tb[unresolved_index]
    
        ## sort by flux
        flux_sort = utb.argsort('Total_flux')
        utb_sorted = utb[ flux_sort[::-1] ]
        if AllFile is not None:
            utb_sorted.write( AllFile, format='ascii.csv' )
        
        ### add dummy spec ind
        utb_sorted.add_column(Column(np.nan*np.zeros(len(utb_sorted)), 'SpecInd'))
        
        ## keep only RA, DEC, Source_id and reorder the columns
        utb_sorted.keep_columns(['RA', 'DEC', 'Total_flux', 'SpecInd', 'Source_id'])

        utb_final = utb_sorted[ 'RA', 'DEC', 'Total_flux', 'SpecInd', 'Source_id' ]
        utb_final.rename_column('Source_id', 'ID')

        utb_final.write( ResultsFile, format='ascii.csv' )

    return




def main(ms_input, ResultsFile, Radius=1.5, DoDownload="True", AllFile=None):

    """
    Download the LoTSS skymodel for the target field

    Parameters
    ----------
    ms_input : str
        String from the list (map) of the target MSs
    ResultsFile : str
        Full name (with path) to the skymodel; if YES is true, the LOTSS skymodel will be downloaded here
    Radius : float (default = 1.5)
        Radius for the LOTSS cone search in degrees
    DoDownload : str ("Force" or "True" or "False")
        Download or not the LOTSS skymodel.
        "Force": download skymodel from LOTSS, delete existing skymodel if needed.
        "True" or "Yes": use existing skymodel file if it exists, download skymodel from 
                         LOTSS if it does not.
        "False" or "No": Do not download skymodel, raise an exception if skymodel
                         file does not exist.
    
    """

    FileExists = os.path.isfile(ResultsFile)
    if (not FileExists and os.path.exists(ResultsFile)):
        raise ValueError("download_lotss_skymodel: WTF! Path: \"%s\" exists but is not a file!"%(ResultsFile))
    download_flag = False
    if DoDownload.upper() == "FORCE":
        if FileExists:
            os.remove(ResultsFile)
        download_flag = True
    elif DoDownload.upper() == "TRUE" or DoDownload.upper() == "YES":
        if FileExists:
            print "USING the exising skymodel in "+ ResultsFile
            return
        else:
            download_flag = True
    elif DoDownload.upper() == "FALSE" or DoDownload.upper() == "NO":
         if FileExists:
            print "USING the exising skymodel in "+ ResultsFile
            return
         else:
            raise ValueError("download_lotss_skymodel: Path: \"%s\" does not exist and LOTSS download is disabled!"%(ResultsFile))

    # If we got here, then we are supposed to download the skymodel.
    assert download_flag == True # Jaja, belts and suspenders...
    print "DOWNLOADING LOTSS Skymodel for the target into "+ ResultsFile

    # Reading a MS to find the coordinate (pyrap)
    #[RATar,DECTar]=grab_coo_MS(input2strlist_nomapfile(ms_input)[0])
    #mypos = ( RATar, DECTar )
    RATar, DECTar = grab_coo_MS(input2strlist_nomapfile(ms_input)[0])




if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=' Download LBCS or LoTSS catalogue for a target field and find good calibrators (more than 2 P values)')

    parser.add_argument('MSfile', type=str, nargs='+', help='One (or more MSs) for which an LBCS catalogue will be download.')
    parser.add_argument('--Radius', type=float, help='Radius for the LBCS cone search in degrees')
    parser.add_argument('--Outfile', type=str, help='Filename to save the results to')

    args = parser.parse_args()
    radius=1.5
    if args.Radius:
        radius=args.Radius

    main(args.MSfile,args.Outfile,Radius=radius)
