#!/usr/bin/env python

import sys
import os
import time
import xarray as xr
import numpy as np

#from dask_mpi import initialize
#initialize()

regions    = ["hawaii"]

if __name__ == "__main__":
    #client = Client()
    #print(client)

    for region in regions:

        main_dir = f'/glade/p/ral/hap/mizukami/oconus_hydro/{region}_run/met'

        if region == 'alaska':
            yrs = range(1980,2017)
            obs = 'daymet'
        elif region == 'hawaii':
            yrs = range(1990,2015)
            obs = 'uh'

        ds1 = xr.open_dataset(os.path.join(main_dir, f'{region}_mask.nc'))
        mask = ds1['mask'].where(ds1['mask']==1).notnull()

        sub_dir=os.path.join(main_dir,f'monthly/{obs}')
        if not os.path.exists(sub_dir):
            os.makedirs(sub_dir)

        print(f'Processing {obs} data')
        for yr in yrs:

          innc  = os.path.join(main_dir, f'daily/{obs}/{obs}_met_{yr}.nc')
          outnc = os.path.join(sub_dir, f'{obs}_met_{yr}.nc')

          ds = xr.open_dataset(innc, chunks={"time": 365})
          if region == 'alaska':
              ds = ds.isel(y=range(15,181), x=range(40,245))

          ds_mon = ds.resample(time="1MS").mean(dim="time").where(mask)

          # clean up
          # need to overwrite latitude and longitude varialbes from daily file
          if region == 'alaska':
              ds_mon['latitude'] = ds['latitude']
              ds_mon['longitude'] = ds['longitude']

          ds_mon = ds_mon.fillna(-999.0)

          for var in ds.variables:
              ds_mon[var].attrs = ds[var].attrs
              ds_mon[var].encoding['_FillValue'] = -999.0

          ds_mon['time'].encoding['dtype']='int32'

          ds_mon.load().to_netcdf(outnc, unlimited_dims=['time'])
