#!/usr/bin/env python

import sys
import os
import xarray as xr
import numpy as np

if __name__ == "__main__":

    for region in ['alaska','hawaii']:
        main_dir = f'/glade/p/ral/hap/mizukami/oconus_hydro/{region}_run/output'

        ds1 = xr.open_dataset(os.path.join(main_dir, 'daily/BCSD/ACCESS1-3/rcp45/ACCESS1-3_rcp45_BCSD_wf_1950.nc'))

        if region=='alaska':
            ds2 = ds1[['latitude','longitude']].isel(y=range(15,181), x=range(40,245))
            ds2['mask'] = ds1['latitude'].isel(y=range(15,181), x=range(40,245)).notnull()
        elif region=='hawaii':
            ds2 = ds1[['latitude','longitude']]
            ds2['mask'] = ds1['latitude'].notnull()

        ds2.to_netcdf(os.path.join(main_dir,f'{region}_mask1.nc'))
