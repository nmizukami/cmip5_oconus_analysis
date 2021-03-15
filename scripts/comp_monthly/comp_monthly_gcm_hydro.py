#!/usr/bin/env python

import sys
import os
import time
import xarray as xr
import numpy as np

#from dask_mpi import initialize
#initialize()

downscales = ["BCSD"]
#gcms       = ["ACCESS1-3","CanESM2","CCSM4","CSIRO-Mk3-6-0","GFDL-ESM2M","HadGEM2-ES","inmcm4","MIROC5","MPI-ESM-MR","MRI-CGCM3"]
gcms       = ["MRI-CGCM3"]
rcps       = [85]
varTypes   = ["wf"]
regions    = ["hawaii"]

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

                    sub_dir=os.path.join(main_dir,f'monthly/{downscale}/{gcm}/rcp{rcp}')
                    if not os.path.exists(sub_dir):
                        os.makedirs(sub_dir)
                    for varType in varTypes:

                        print(f'Processing {gcm} rcp{rcp} {downscale} {varType} data')
                        for yr in yrs:

                            innc  = os.path.join(main_dir, f'daily/{downscale}/{gcm}/rcp{rcp}/{gcm}_rcp{rcp}_{downscale}_{varType}_{yr}.nc')
                            outnc = os.path.join(sub_dir, f'{gcm}_rcp{rcp}_{downscale}_{varType}_{yr}.nc')

                            ds = xr.open_dataset(innc,chunks={"time": 365})
                            if region == 'alaska':
                                ds = ds.isel(y=range(15,181), x=range(40,245))

                            if varType == 'ws' or varType == 'eb':
                                ds_mon = ds.resample(time="1MS").mean(dim="time").where(mask)
                            elif varType == 'wf':
                                ds_mon = ds.resample(time="1MS").sum(dim="time").where(mask)

                            # clean up
                            ds_mon = ds_mon.fillna(-999.0)
                            ds_mon['latitude']  = ds_mon['latitude'].isel(time=0, drop=True)
                            ds_mon['longitude'] = ds_mon['longitude'].isel(time=0, drop=True)

                            for var in ds.variables:
                                ds_mon[var].attrs = ds[var].attrs
                                ds_mon[var].encoding['_FillValue'] = -999.0
                                if varType == 'wf' and var not in ['latitude','longitude','time']:
                                    ds_mon[var].attrs['units'] = 'mm/month'
                            ds_mon['time'].encoding['dtype']='int32'

                            ds_mon.load().to_netcdf(outnc, unlimited_dims=['time'])
