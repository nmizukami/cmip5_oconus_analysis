#!/usr/bin/env python

import sys
import os
import time
import xarray as xr
import numpy as np

#from dask_mpi import initialize
#initialize()

np.seterr(divide='ignore')

downscales = ["BCSD"]
gcms       = ["ACCESS1-3","CanESM2","CCSM4","CSIRO-Mk3-6-0","GFDL-ESM2M","HadGEM2-ES","inmcm4","MIROC5","MPI-ESM-MR","MRI-CGCM3"]
rcps       = [85]
regions    = ["hawaii"] # hawaii or alaska

var_drop = {'alaska': ['height','yv','xv','tmax_anom_year','tmin_anom_year','pcp_anom_year'],
            'hawaii': ['height','tmax_anom_year','tmin_anom_year','pcp_anom_year']}

var_rename = {'alaska': {'yc':'latitude','xc':'longitude','ni':'x','nj':'y'},
              'hawaii': {}}


if __name__ == "__main__":
    #client = Client()
    #print(client)

    yrs = range(1950,2100)

    for region in regions:

        def preprocess(ds):
          ds = ds.drop_vars(var_drop[region])
          ds = ds.rename(var_rename[region])
          return ds

        main_dir = f'/glade/p/ral/hap/mizukami/oconus_hydro/{region}_run/met'

        ds1 = xr.open_dataset(os.path.join(main_dir, f'{region}_mask.nc'))
        mask = ds1['mask'].where(ds1['mask']==1).notnull()

        for downscale in downscales:
            for gcm in gcms:
                for rcp in rcps:

                    sub_dir=os.path.join(main_dir,f'monthly/{downscale}/{gcm}/rcp{rcp}')
                    if not os.path.exists(sub_dir):
                        os.makedirs(sub_dir)

                    print(f'Processing {gcm} rcp{rcp} {downscale} data')

                    downscale_lower = downscale.lower()
                    innc  = os.path.join(main_dir, f'daily/rcp{rcp}/{downscale_lower}_{gcm}_*_out.nc')
                    ds = xr.open_mfdataset(innc, chunks={"time": 365} ,preprocess=preprocess)

                    if region == 'alaska':
                        ds = ds.isel(y=range(15,181), x=range(40,245))

                    ds_mon = ds.resample(time="1MS").mean(dim="time").where(mask)

                    print(ds_mon)

                    # clean up
                    ds_mon = ds_mon.fillna(-999.0)

                    for var in ds.variables:
                        ds_mon[var].attrs = ds[var].attrs
                        ds_mon[var].encoding['_FillValue'] = -999.0

                    ds_mon['time'].encoding['dtype']='int32'

                    for yr in yrs:

                        outnc = os.path.join(sub_dir, f'{gcm}_rcp{rcp}_{downscale}_met_{yr}.nc')

                        ds_mon.sel(time=slice('%d-01-01'%yr, '%d-12-31'%(yr))).load().to_netcdf(outnc, unlimited_dims=['time'])
