#!/usr/bin/env python

import sys
import os
import xarray as xr
import numpy as np

if __name__ == "__main__":

    main_dir = f'/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/met'

    #ds1 = xr.open_dataset(os.path.join(main_dir, 'daily/rcp45/bcsd_ACCESS1-3_historical_out.nc'))
    ds1 = xr.open_dataset(os.path.join(main_dir, 'daily/UH/UH_1990.nc'))
    ds2 = ds1[['lat','lon']]
    ds2['mask'] = ds1['tmax'].isel(time=0).notnull()
#    ds2 = ds2.drop_vars(['height','tmax_anom_year','tmin_anom_year','pcp_anom_year','time'])
    ds2 = ds2.drop_vars(['time'])

    ds2.to_netcdf(os.path.join(main_dir,'hawaii_met_mask.nc'))
