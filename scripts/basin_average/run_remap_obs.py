#!/usr/bin/env python
'''
Only Alaska
'''

import sys, os, glob

varTypes   = ["ws"]
regions    = ['alaska']

for region in regions:

    main_dir = f'/glade/p/ral/hap/mizukami/oconus_hydro/{region}_run/output/monthly'

    if region == 'alaska':
        wgt_nc = os.path.join('/glade/work/mizukami/py_mapping/alaska/mapping/spatialweights_grid_12km_to_MERIT_Hydro_alaska_basin.nc')
        yrs = range(1980,2017)
        obs = 'daymet'
    elif region == 'hawaii':
        yrs = range(1990,2015)
        obs = 'uh'

    sub_dir=os.path.join(main_dir,f'basin_average/{obs}')

    if not os.path.exists(sub_dir):
        os.makedirs(sub_dir)

    for varType in varTypes:

        print(f'Processing {obs} data')
        for yr in yrs:

            innc  = os.path.join(main_dir, f'{obs}/{obs}_{varType}_{yr}.nc')
            outnc = os.path.join(sub_dir, f'{obs}_{varType}_{yr}.nc')

            if varType=='wf':
                meta='flux.meta'
            elif varType=='ws':
                meta='state.meta'
            elif varType=='eb':
                meta='energy.meta'

            cmd="./remap_vic_output.py %s %s %s ./%s int32"%(innc, wgt_nc, outnc, meta)
            os.system(cmd)
