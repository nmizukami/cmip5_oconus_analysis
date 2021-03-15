#!/usr/bin/env python
'''
Only Alaska
'''

import sys, os, glob

downscales = ["BCSD"]
gcms       = ["ACCESS1-3","CanESM2","CCSM4","CSIRO-Mk3-6-0","GFDL-ESM2M","HadGEM2-ES","inmcm4","MIROC5","MPI-ESM-MR","MRI-CGCM3"]
rcps       = [45,85]
varTypes   = ["ws"]

region='alaska'

wgt_nc   = os.path.join('/glade/work/mizukami/py_mapping/alaska/mapping/spatialweights_grid_12km_to_MERIT_Hydro_alaska_basin.nc')

yrs = range(1950,2100)
main_dir = f'/glade/p/ral/hap/mizukami/oconus_hydro/{region}_run/output/monthly'

for downscale in downscales:
    for gcm in gcms:
        for rcp in rcps:

            sub_dir=os.path.join(main_dir,f'basin_average/{downscale}/{gcm}/rcp{rcp}')

            if not os.path.exists(sub_dir):
                os.makedirs(sub_dir)

            for varType in varTypes:

                print(f'Processing {gcm} rcp{rcp} {downscale} {varType} data')
                for yr in yrs:

                    innc  = os.path.join(main_dir, f'{downscale}/{gcm}/rcp{rcp}/{gcm}_rcp{rcp}_{downscale}_{varType}_{yr}.nc')
                    outnc = os.path.join(sub_dir, f'{gcm}_rcp{rcp}_{downscale}_{varType}_{yr}.nc')

                    if varType=='wf':
                        meta='flux.meta'
                    elif varType=='ws':
                        meta='state.meta'
                    elif varType=='eb':
                        meta='energy.meta'

                    cmd="./remap_vic_output.py %s %s %s ./%s int32"%(innc, wgt_nc, outnc, meta)
                    os.system(cmd)
