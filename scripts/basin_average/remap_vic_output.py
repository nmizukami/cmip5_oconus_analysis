#!/usr/bin/env python
''' Process timeseries grid into mean areal timeseries for arbitrary polygons
    Depends on mapping file (i,j version) from poly2poly.py
    Rewritten from earlier script pieces to use vectorized operations -- massive speed up.
    Writes files that are in local time for the western US
    Converts met forecast raw units to SUMMA units
    Aug 2020, A. Wood.

    Modifications (this script version):
     * AWW Oct 2020 adapted more efficient/faster scripting to GMET datasets
     * divide input time var by 1e9 to get seconds since ... correct in output
     * from JS: reorders outputs:  uses a target attribute file to ensure that the HRUs are written to match
        the order of other inputs (eg attributes & params) for SUMMA (useful for specific basins);
        previous versions write out subsets of all HRUs in the order of spatial weights file (useful for large
        domains that must be split)
'''
# Notes:  on cheyenne: module load python/3.6.8; source /glade/u/apps/opt/ncar_pylib/ncar_pylib.csh
#                   (or ncar_pylib)
# =========================================================================
import sys, os, time
import argparse
import numpy as np
import xarray as xr
import pandas as pd

########################################################################
#                Subroutines / Functions                               #
########################################################################

TIME_DIM_NAME = 'time'
FILL_VALUE   = -999.

def process_command_line():
    '''Parse the commandline'''
    parser = argparse.ArgumentParser(description='remap WRFHydro 1km forcing data to HRU and write in netcdf')
    parser.add_argument('frc_nc',   help='path to WRFHydro 1km forcing netcdf.')
    parser.add_argument('wgt_nc',   help='path to mapping netcdf.')
    parser.add_argument('out_nc',   help='path to remapped forcing netcdf.')
    parser.add_argument('meta',     help='path to variable meta data.')
    parser.add_argument('hru_type', help='hru data type. e.g., str, int, int64, float, float32')
    # Optional arguments
    parser.add_argument('--idfile', action='store_true', default='allPolygon', help='path of file with list of reach ids.')

    return parser.parse_args()


class AutoVivification(dict):
    """Implementation of perl's autovivification feature to initialize structure."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
        return value


def getVar(inmeta):
    """ read metadata for hru and segment shp attributes """
    meta   = AutoVivification()
    with open(inmeta) as fp:
        for line in fp:
            cols=line.split(',')
            varname=cols[0].strip()
            meta[varname]['agg']       = cols[1].strip()
            meta[varname]['dtyp']      = cols[2].strip()
            meta[varname]['long_name'] = cols[3].strip()
            meta[varname]['units']     = cols[4].strip()
    return meta


def getMeta(forc_nc):
    """ get variable attributes and encodings in netcdf """
    attrs = {}; encodings={}
    with xr.open_dataset(forc_nc) as ds:
        for varname in ds.variables:
           attrs[varname] = ds[varname].attrs
           encodings[varname] = ds[varname].encoding
    return attrs, encodings


def writeNetCDFData(out_nc, hrus, dr_time, hru_type, remapped_data, var_meta, var_attrs, var_encodings, remap_idx):
    """ Write <vars>[time,hru] array in netCDF4 file,<fn> and variable of <varname> """

    dataset = xr.Dataset()

    for varname, meta in var_meta.items():
       foo = xr.DataArray(remapped_data[varname][:, remap_idx],
                          dims=['time', 'basinID'],
                          name=varname)

       foo.encoding = var_encodings[varname]
       foo.attrs    = var_attrs[varname]

       dataset[varname] = foo

    # HRU ID variables
    dataset['basinID'] = xr.DataArray(hrus[remap_idx], dims=['basinID'])
    dataset['basinID'].encoding = {'dtype': hru_type, '_FillValue': None}
    dataset['basinID'].attrs    = {'long_name': 'Basin ID'}

    dataset[TIME_DIM_NAME] = dr_time

    dataset.to_netcdf(out_nc, unlimited_dims='time')


def compAvgVal(matWgts, matIndex_i, matIndex_j, overlaps, xd, varname, default=FILL_VALUE):
    """Compute areal weighted avg value of <varname> in <nc_in> for all output polygons
       based on input/output overlap weights in <nc_wgt>"""

    # now read the input timeseries file
    print("-------------------")
    print("reading input timeseries data for %s " % varname)
    dataVals = xd[varname].values
    print("INFO: value at [0,1,1]: %f" % (dataVals[0,1,1]))

    array_shape = dataVals.shape
    nDims       = len(array_shape)
    nTimeSteps = array_shape[0]

    wgtedVals   = np.zeros((nTimeSteps, nOutPolys))
    matDataVals = np.zeros((nTimeSteps, nOutPolys, maxOverlaps))

    # reformat var data into regular matrix matching weights format (nOutPolygons, maxOverlaps)
    #   used advanced indexing to extract matching input grid indices
    for p in range(0, nOutPolys):
        if overlaps[p]>0:
            matDataVals[:, p, 0:overlaps[p]] = dataVals[:, matIndex_j[p, 0:overlaps[p]], matIndex_i[p, 0:overlaps[p]] ]
        else:
            matDataVals[:, p, 0] = default

    wgtedVals = np.sum(matDataVals * matWgts, axis=nDims-1)   # produces vector of weighted values

    print(" averaged var %s" % (varname))

    return wgtedVals

############################################
#                Main                      #
############################################
if __name__ == '__main__':

    # process command line
    args = process_command_line()

    # get variable meta data from ascii files
    var_meta = getVar(args.meta)

    # get variables attributes and encodings from input netCDF
    var_attrs, var_encodings = getMeta( args.frc_nc )

    # Get data from spatial weights file
    with xr.open_dataset(args.wgt_nc) as ds:
        wgtVals =       ds['weight'].values
        i_index =       ds['i_index'].values
        j_index =       ds['j_index'].values        # j_index --- LAT direction
        allOutPolyIDs = ds['polyid'].values
        overlaps =      ds['overlaps'].values
        #IDmask = xspw['IDmask'].values

        # Check data dimension size and sum of overlapping elements (expected to be the same)
        # grid2poly.py output skip data in data dimension, but poly2poly put missing values in data dimension
        dim_data_size = ds.dims['data']
        sum_overlaps  = overlaps.sum()
        flag = False
        if dim_data_size != sum_overlaps:
            flag = True
            print('WARNING: data dimension size (%d) is not equal to sum of overlapping elements (%d)'%(dim_data_size, sum_overlaps))

    remap_idx = np.arange(len(allOutPolyIDs))
    # OPTIONAL: subset mode
    if args.idfile != 'allPolygon':
        # open a netcdf file containing a vector of the target IDs (eg use the attributes file of target domain)
        #   these will be used to subset the spatial weights mappings and output data in desired id order
        dsTargIds = xr.open_dataset(idfile, decode_times=False)
        ids = dsTargIds.hruId.values

        # find matching spwts ID indices to target index file (this works but is slow; better method desired)
        #   remap_idx used to re-index data values (polys) to desired before writing
        remap_idx = [xp for (xi, x) in enumerate(ids) for (xp, y) in enumerate(allOutPolyIDs) if x==y]

    # start reading input data (dims of data variables: 3D [time, x, y])
    xd = xr.open_dataset(args.frc_nc, decode_cf=False)
    dr_time = xd[TIME_DIM_NAME]

    # set to zero based index and flip the N-S index (for read-in data array 0,0 is SW-corner)
    j_index = j_index - 1              # j_index starts at North, and at 1 ... make 0-based
    i_index = i_index - 1              # i_index starts at West, and a 1 ... make 0-based

    # number of target HRUs, and maximum number of overlapping grid boxes
    nOutPolys = len(allOutPolyIDs)
    maxOverlaps = overlaps.max()

    # assign weights and indices to a regular array (nOutPolys x maxOverlaps)
    matWgts    = np.zeros((nOutPolys, maxOverlaps), dtype='float32')
    matIndex_i = np.zeros((nOutPolys, maxOverlaps), dtype='int32')
    matIndex_j = np.zeros((nOutPolys, maxOverlaps), dtype='int32')
    ix2=0;
    for p in range(0, nOutPolys):
        if overlaps[p]>0:
            ix1 = ix2
            ix2 = ix1+overlaps[p]
        elif overlaps[p]==0 and flag: # WARNNG: grid2poly.py output skip data in data dimension, but poly2poly put missing values in data dimension
            ix1 = ix2
            ix2 = ix2+1
        elif overlaps[p]==0 and not flag:
            matWgts[p, 0:overlaps[p]] = 0
            matWgts[p, 0] = 1.0
            continue
        matWgts[p, 0:overlaps[p]]    = wgtVals[ix1:ix2]/wgtVals[ix1:ix2].sum()
        matIndex_i[p, 0:overlaps[p]] = i_index[ix1:ix2]
        matIndex_j[p, 0:overlaps[p]] = j_index[ix1:ix2]
    print("Pre-process weight and index arrays")

    # now use these weights in the avg val computation
    remapped_data = {}
    for var_name in var_meta.keys():
      remapped_data[var_name] = compAvgVal(matWgts, matIndex_i, matIndex_j, overlaps, xd, var_name)

    # write the output file
    writeNetCDFData(args.out_nc, allOutPolyIDs, dr_time, args.hru_type, remapped_data, var_meta, var_attrs, var_encodings, remap_idx)
    print("wrote output file %s " % args.out_nc)

# DONE
