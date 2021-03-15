#!/usr/bin/env python

import sys
import os
import time
import xarray as xr
import numpy as np

#from dask_mpi import initialize
#initialize()

varTypes   = ["eb"]
regions    = ['alaska']

if __name__ == "__main__":

#    from distributed import Client
#    client = Client()
#    print(client)

    for region in regions:

        main_dir = f'/glade/p/ral/hap/mizukami/oconus_hydro/{region}_run/output'

        if region == 'alaska':
            yrs = range(1980,2017)
            obs = 'daymet'
        elif region == 'hawaii':
            yrs = range(1990,2015)
            obs = 'uh'

        ds1 = xr.open_dataset(os.path.join(main_dir, f'{region}_mask.nc'))
        mask = ds1['latitude'].notnull()

        sub_dir=os.path.join(main_dir, f'monthly/{obs}')
        if not os.path.exists(sub_dir):
            os.makedirs(sub_dir)
        for varType in varTypes:

            print(f'Processing {region} {obs} ata')
            for yr in yrs:

                innc  = os.path.join(main_dir, f'daily/{obs}/{varType}_{yr}-01.nc')
                outnc = os.path.join(sub_dir, f'{obs}_{varType}_{yr}.nc')

                ds = xr.open_dataset(innc,chunks={"time": 365})
                if region=='alaska':
                    ds = ds.isel(y=range(15,181), x=range(40,245))

                if varType=='ws' or varType=='eb':
                    ds_mon = ds.resample(time="1MS").mean(dim="time").where(mask)
                elif varType=='wf':
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
