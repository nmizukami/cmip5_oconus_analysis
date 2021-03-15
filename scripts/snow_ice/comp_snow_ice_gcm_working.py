#!/usr/bin/env python

import sys
import os
import time
import xarray as xr
import pandas as pd
import numpy as np

#from dask_mpi import initialize
#initialize()

downscales = ["BCSD"]
gcms       = ["ACCESS1-3","CanESM2","CCSM4","CSIRO-Mk3-6-0","GFDL-ESM2M","HadGEM2-ES","inmcm4","MIROC5","MPI-ESM-MR","MRI-CGCM3"]
rcps       = [45,85]
regions    = ["alaska"]

snow_variable = {'SWEmax': 'annual maximum SWE', 'DOY_SWEmax':'day of year with maximum SWE', 'DOY_MELT':'snow melt-out days'}

def subset(ds):
    variable = ['SWE','IWE']
    return ds[variable].isel(y=range(15,181), x=range(40,245))

if __name__ == "__main__":
    #client = Client()
    #print(client)

    yrs = range(1950,2100)
    for region in regions:
        main_dir = f'/glade/p/ral/hap/mizukami/oconus_hydro/{region}_run/output'

        ds1 = xr.open_dataset(os.path.join(main_dir, f'{region}_mask.nc'))
        mask = ds1['latitude'].notnull()

        for downscale in downscales:
            for gcm in gcms:
                for rcp in rcps:

                    sub_dir=os.path.join(main_dir,f'snow_ice/{downscale}/{gcm}/rcp{rcp}')
                    if not os.path.exists(sub_dir):
                        os.makedirs(sub_dir)

                    print(f'Processing {gcm} rcp{rcp} {downscale} data')

                    for iy, yr in enumerate(yrs):
                        innc  = os.path.join(main_dir, f'daily/{downscale}/{gcm}/rcp{rcp}/{gcm}_rcp{rcp}_BCSD_ws_{yr}.nc')
                        ds = xr.open_dataset(innc, chunks={"time": 365})
                        ds = subset(ds)
                        dr_max = ds['SWE'].resample(time='2QS').max('time').where(mask).isel(time=0)

                        dr_day_max = ds['SWE'].resample(time='2QS').argmax('time').where(mask)
                        print(dr_max)
                        sys.exit()

                    ds_snow = dr_max.to_dataset(name = 'SWEmax')
                    ds_snow['DOY_SWEmax'] = dr_day_max

                #    # clean up
                #    ds_snow = ds_snow.fillna(-999.0)

                    for var in ds_snow.variables:
                        ds_snow[var].attrs['long_name'] = snow_variable[var]
                        ds_snow[var].encoding['_FillValue'] = -999.0

                    ds_snow['time'].encoding['dtype']='int32'

                    outnc = os.path.join(sub_dir, f'{gcm}_rcp{rcp}_{downscale}_snow.nc')
                    ds_snow.load().to_netcdf(outnc, unlimited_dims=['time'])

                    sys.exit()
