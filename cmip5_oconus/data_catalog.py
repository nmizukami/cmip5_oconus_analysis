
import os
import warnings
import glob

import xarray as xr

# TODO: make this more configurable
# GCMs
CMIP5_VIC_DAY_ROOT_DIR = {'grid': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/output/daily/BCSD',
                                   'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/output/daily/BCSD'},
                          'hrus': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/output/daily/basin_average/BCSD',
                                   'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/output/daily/basins_average/BCSD'}
                         }

CMIP5_VIC_MON_ROOT_DIR = {'grid': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/output/monthly/BCSD',
                                   'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/output/monthly/BCSD'},
                          'hrus': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/output/monthly/basin_average/BCSD',
                                   'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/output/monthly/basin_average/BCSD'}
                         }

CMIP5_MET_DAY_ROOT_DIR = {'grid': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/met',
                                   'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/met'},
                          'hrus': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/met/basin_average/BCSD',
                                   'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/met/basin_average/BCSD'}
                         }

CMIP5_MET_MON_ROOT_DIR = {'grid': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/met/monthly/BCSD',
                                   'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/met/monthly/BCSD'},
                          'hrus': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/met/monthly/basin_average/BCSD',
                                   'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/met/monthly/basin_average/BCSD'}
                         }

# observed
OBS_VIC_DAY_ROOT_DIR = {'grid': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/output/daily/daymet',
                                 'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/output/daily/uh'},
                        'hrus': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/output/daily/basin_average/daymet',
                                 'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/output/daily/basin_average/uh'}
                       }

OBS_VIC_MON_ROOT_DIR = {'grid': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/output/monthly/daymet',
                                 'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/output/monthly/uh'},
                        'hrus': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/output/monthly/basin_average/daymet',
                                 'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/output/monthly/basin_average/uh'}
                       }

OBS_MET_MON_ROOT_DIR = {'grid': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/met/monthly/daymet',
                                 'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/met/monthly/uh'},
                        'hrus': {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/met/monthly/basin_average/daymet',
                                 'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/met/monthly/basin_average/uh'}
                       }

FILE_TAG = {'NET_SHORT':'eb','NET_LONG':'eb','SENSIBLE':'eb','LATENT':'eb','GRND_FLUX':'eb','SOIL_TEMP1':'eb','SOIL_TEMP2':'eb','SOIL_TEMP3':'eb','ENERGY_ERROR':'eb',
            'SWE':'ws','IWE':'ws','SM1':'ws','SM2':'ws','SM3':'ws','WATER_ERROR':'ws',
            'total_runoff':'wf','RUNOFF':'wf','BASEFLOW':'wf','EVAP':'wf','PRCP':'wf', 'SNOW_MELT':'wf', 'GLACIER_MELT':'wf',
            'pcp':'met','tmean':'met', 'tmax':'met', 'tmin':'met', 'dtr':'met',
            }

DOMAIN = {'AK': '/glade/p/ral/hap/mizukami/oconus_hydro/alaska_run/output/alaska_mask.nc',
          'HI': '/glade/p/ral/hap/mizukami/oconus_hydro/hawaii_run/output/hawaii_mask.nc'}


DEFAULT_MON_HYDRO_VARS = ['PRCP', 'EVAP', 'total_runoff', 'SWE', 'SM1', 'SM2', 'SM3']
DEFAULT_DAY_HYDRO_VARS = ['total_runoff']

DEFAULT_MON_MET_VARS = ['pcp', 'tmean', 'tmin', 'tmax']
DEFAULT_DAY_MET_VARS = ['pcp', 'tmean', 'tmin', 'tmax']

KELVIN = 273.13
SEC_PER_DAY = 86400


class NoneError(Exception):
    pass

def progress(r):
    try:
        from tqdm import tqdm
        return tqdm(r)
    except ImportError:
        return r


def _calc_total_runoff(ds):
    if 'total_runoff' in ds:
        return ds['total_runoff']
    return ds['RUNOFF'] + ds['BASEFLOW']

# TODO --- make sure this is the common definition of t_mean
def _calc_t_mean(ds):
    if 'tmean' in ds:
        return ds['tmean']
    return (ds['tmin'] + ds['tmax']) / 2.

def _calc_dtr(ds):
    if 'dtr' in ds:
        return ds['dtr']
    return ds['tmax'] - ds['tmin']


def resample_data(ds, freq='MS', region=None, chunks=None):

    if region is None:
        raise NotImplementedError('No other regions are not provided')

    ds1 = xr.open_dataset(DOMAIN[region])
    mask = ds1['mask'].where(ds1['mask']==1).notnull()

    out = xr.Dataset()
    for name, da in ds.data_vars.items():
        if name in ['PRCP', 'EVAP', 'total_runoff', 'RUNOFF', 'BASEFLOW']:
            out[name] = da.resample(time=freq).sum(skipna=False)#.where(mask)
        else:
            # TODO: weight by days in month, or sum over year
            out[name] = da.resample(time=freq).mean(skipna=False)#.where(mask)

    if chunks is not None:
        out = out.persist().chunk(chunks)
    return out



# Wrappers - start
def load_monthly_historical_hydro_datasets(models=None,
                                           variables=DEFAULT_MON_HYDRO_VARS,
                                           region=None,
                                           dataType='grid',
                                           **kwargs):
    print('load_monthly_historical_hydro_datasets', flush=True)

    if region not in ['AK','HI']:
        raise NotImplementedError('No other regions are supported at this time')
    if dataType not in ['grid','hrus']:
        raise NotImplementedError('No other data types are supported at this time')

    data = {}
    data['gcm'] = load_monthly_cmip5_hydro_datasets('historical', models=models, variables=variables, region=region, dataType=dataType, **kwargs)
    data['obs'] = load_monthly_obs_hydro_datasets(variables=variables, region=region, dataType=dataType, **kwargs)

    return data


def load_daily_historical_hydro_datasets(models=None,
                                         variables=DEFAULT_DAY_HYDRO_VARS,
                                         region=None,
                                         dataType='grid',
                                         **kwargs):
    print('load_daily_historical_hydro_datasets', flush=True)

    if region not in ['AK','HI']:
        raise NotImplementedError('No other regions are supported at this time')
    if dataType not in ['grid','hrus']:
        raise NotImplementedError('No other data types are supported at this time')

    data = {}
    data['gcm'] = load_daily_cmip5_hydro_datasets('historical', models=models, variables=variables, region=region, dataType=dataType, **kwargs)
    data['obs'] = load_daily_obs_hydro_datasets(variables=variables, region=region, dataType=dataType, **kwargs)

    return data


def load_monthly_obs_hydro_datasets(variables=DEFAULT_MON_HYDRO_VARS,
                                    region=None,
                                    dataType='grid',
                                    **kwargs):
    if region not in ['AK','HI']:
        raise NotImplementedError('No other regions are supported at this time')
    if dataType not in ['grid','hrus']:
        raise NotImplementedError('No other data types are supported at this time')

    data = load_obs_dataset(OBS_VIC_MON_ROOT_DIR[dataType][region], variables=variables, **kwargs)

    # TODO: it would be better if we passed this info to the individual loaders
    if 'total_runoff' in variables:
        data['total_runoff'] = _calc_total_runoff(data)
    data = data[variables]
    return data


def load_daily_obs_hydro_datasets(variables=DEFAULT_DAY_HYDRO_VARS,
                                  region=None,
                                  dataType='grid',
                                  **kwargs):
    if region not in ['AK','HI']:
        raise NotImplementedError('No other regions are supported at this time')
    if dataType not in ['grid','hrus']:
        raise NotImplementedError('No other data types are supported at this time')

    data = load_obs_dataset(OBS_VIC_DAILY_ROOT_DIR[dataType][region], variables=variables, **kwargs)

    # TODO: it would be better if we passed this info to the individual loaders
    if 'total_runoff' in variables:
        data['total_runoff'] = _calc_total_runoff(data)
    data = data[variables]
    return data


def load_monthly_cmip5_hydro_datasets(scen, models=None,
                                      variables=DEFAULT_MON_HYDRO_VARS,
                                      region=None,
                                      dataType='grid',
                                      **kwargs):
    print('load_monthly_cmip5_hydro_datasets', flush=True)

    if region not in ['AK','HI']:
        raise NotImplementedError('No other regions are supported at this time')
    if dataType not in ['grid','hrus']:
        raise NotImplementedError('No other data types are supported at this time')

    data = load_cmip5_dataset(CMIP5_VIC_MON_ROOT_DIR[dataType][region], variables=variables, scen=scen, models=models, **kwargs)

    # TODO: it would be better if we passed this info to the individual loaders
    if 'total_runoff' in variables:
        data['total_runoff'] = _calc_total_runoff(data)

    data = data[variables]

    return data


def load_daily_cmip5_hydro_datasets(scen, models=None,
                                    variables=DEFAULT_DAY_HYDRO_VARS,
                                    region=None,
                                    dataType='grid',
                                    **kwargs):
    print('load_daily_cmip5_hydro_datasets', flush=True)

    if region not in ['AK','HI']:
        raise NotImplementedError('No other regions are supported at this time')
    if dataType not in ['grid','hrus']:
        raise NotImplementedError('No other data types are supported at this time')

    data = load_cmip5_dataset(CMIP5_VIC_DAY_ROOT_DIR[dataType][region], variables=variables, scen=scen, models=models, **kwargs)

    # TODO: it would be better if we passed this info to the individual loaders
    if 'total_runoff' in variables:
        data['total_runoff'] = _calc_total_runoff(data)
    data = dat[variables]
    return data


## Meteorology data
def load_monthly_historical_met_datasets(models=None,
                                         variables=DEFAULT_MON_MET_VARS,
                                         region=None,
                                         dataType='grid',
                                         **kwargs):
    print('load_monthly_historical_met_datasets', flush=True)

    if region not in ['AK','HI']:
        raise NotImplementedError('No other regions are supported at this time')
    if dataType not in ['grid','hrus']:
        raise NotImplementedError('No other data types are supported at this time')

    data = {}
    data['gcm'] = load_monthly_cmip5_met_datasets('historical', models=models, variables=variables, region=region, dataType=dataType, **kwargs)
    data['obs'] = load_monthly_obs_met_datasets(variables=variables, region=region, dataType=dataType, **kwargs)

    return data


def load_daily_historical_met_datasets(models=None,
                                       variables=DEFAULT_DAY_MET_VARS,
                                       region=None,
                                       dataType='grid',
                                       **kwargs):
    print('load_daily_historical_met_datasets', flush=True)

    if region not in ['AK','HI']:
        raise NotImplementedError('No other regions are supported at this time')
    if dataType not in ['grid','hrus']:
        raise NotImplementedError('No other data types are supported at this time')

    data = {}
    data['gcm'] = load_daily_cmip5_met_datasets('historical', models=models, variables=variables, region=region, dataType=dataType, **kwargs)
    data['obs'] = load_daily_obs_met_datasets(variables=variables, region=region, dataType=dataType, **kwargs)

    return data


def load_monthly_obs_met_datasets(variables=DEFAULT_MON_MET_VARS,
                                  region=None,
                                  dataType='grid',
                                  **kwargs):
    if region not in ['AK','HI']:
        raise NotImplementedError('No other regions are supported at this time')
    if dataType not in ['grid','hrus']:
        raise NotImplementedError('No other data types are supported at this time')

    data = load_obs_dataset(OBS_MET_MON_ROOT_DIR[dataType][region], variables=variables, **kwargs)

    if 'tmean' in variables:
        data['tmean'] = _calc_t_mean(data)

    if 'dtr' in variables:
        data['dtr'] = _calc_dtr(data)

    data = data[variables]

    return data


def load_monthly_cmip5_met_datasets(scen, models=None,
                                    variables=DEFAULT_MON_MET_VARS,
                                    region=None,
                                    dataType='grid',
                                    **kwargs):
    print('load_monthly_cmip5_met_datasets', flush=True)

    if region not in ['AK','HI']:
        raise NotImplementedError('No other regions are supported at this time')
    if dataType not in ['grid','hrus']:
        raise NotImplementedError('No other data types are supported at this time')

    data = load_cmip5_dataset(CMIP5_MET_MON_ROOT_DIR[dataType][region], variables=variables, scen=scen, models=models, **kwargs)

    if 'tmean' in variables:
        data['tmean'] = _calc_t_mean(data)

    if 'dtr' in variables:
        data['dtr'] = _calc_dtr(data)

    data = data[variables]

    return data


def load_daily_cmip5_met_datasets(scen, models=None,
                                  variables=DEFAULT_DAY_MET_VARS,
                                  region=None,
                                  dataType='grid',
                                  **kwargs):
    print('load_daily_cmip5_met_datasets', flush=True)

    if region not in ['AK','HI']:
        raise NotImplementedError('No other regions are supported at this time')
    if dataType not in ['grid','hrus']:
        raise NotImplementedError('No other data types are supported at this time')

    data = load_cmip5_dataset(CMIP5_MET_DAY_ROOT_DIR[dataType][region], variables=variables, scen=scen, models=models, **kwargs)

    if 'tmean' in variables:
        data['tmean'] = _calc_t_mean(data)

    if 'dtr' in variables:
        data['dtr'] = _calc_dtr(data)

    data = dat[variables]

    return data

# Wrappers - end


def drop_bound_varialbes(ds):
    drops = []
    for v in ['lon_bnds', 'lat_bnds', 'time_bnds']:
        if v in ds or v in ds.coords:
            drops.append(v)
    return ds.drop(drops)


def get_valid_years(scen):
    if 'hist' in scen:
        r = range(1950, 2006)
    else:
        r = range(2005, 2100)

    return list(map(str, list(r)))


def filter_files(files, valid_years):
    out = []
    for f in files:
        for y in valid_years:
            if y in f:
                out.append(f)
                break
    return list(set(out))


def load_cmip5_dataset(root, variables=None, scen='rcp85', models=None, **kwargs):
    print('load_cmip5_dataset', flush=True)

    if variables is None:
        raise NoneError("variables list not provided. Cannot access a 'None'.")

    valid_years = get_valid_years(scen)
    if 'hist' in scen:
        scen = 'rcp85'  # bcsd put historical in the rcp dataset
    elif 'hist_rcp45':
        scen = 'rcp45'  # bcsd put historical in the rcp dataset

    if models is None:
        models = os.listdir(root)
        models.sort()

    _file_tags = []
    for var in variables:
        _file_tags.append(FILE_TAG[var])
    file_tags = list(set(_file_tags))

    ds_list = []
    models_list = []
    for m in progress(models):

        for ix, tag in enumerate(file_tags):
            fpath = os.path.join(root, f'{m}', scen, f'*_{tag}_*nc')
            if ix==0:
                files = glob.glob(fpath)
            else:
                files += glob.glob(fpath)

        if not files:
            warnings.warn('no files to open: %s' % fpath)
            models.remove(m)
            continue

        files = filter_files(files, valid_years)
        try:
            ds_list.append(xr.open_mfdataset(files,
                                             preprocess=drop_bound_varialbes,
                                             **kwargs))
            models_list.append(m)
        except OSError:
            print('skipping %s' % m)

    ds = xr.concat(ds_list, dim=xr.Variable('gcm', models_list))

    for var in ['bounds_latitude', 'bounds_longitude']:
        if var in ds:
            ds = ds.drop(var)

    return ds

def load_obs_dataset(root, variables=None, **kwargs):
    print('load_obs_dataset', flush=True)

    if variables is None:
        raise NoneError("variables list not provided. Cannot access a 'None'.")

    _file_tags = []
    for var in variables:
        _file_tags.append(FILE_TAG[var])
    file_tags = list(set(_file_tags))

    files = []
    for tag in file_tags:
        fpath = os.path.join(root, f'*_{tag}_*nc')
        files.append(glob.glob(fpath))

    ds = xr.open_mfdataset(files, **kwargs)

    return ds


###. Not used
def load_daily_obs_meteorology(region=None, dataType=None, **kwargs):
    print('load_daily_obs_meteorology', flush=True)

    def preproc(ds):
        if 'latitude' in ds:
            # 1 or 2 files have different coordinate data so we fix that here
            ds = ds.rename({'latitude': 'lat', 'longitude': 'lon'})
        for var in ['bounds_latitude', 'bounds_longitude',
                    'longitude_bnds', 'latitude_bnds']:
            if var in ds:
                ds = ds.drop(var)
        ds['lon'] = ds['lon'].where(ds['lon'] <= 180, ds['lon'] - 360)
        return ds

    fpath = os.path.join(LIVNEH_MET_ROOT_DIR, '*nc')
    ds = xr.open_mfdataset(fpath, **kwargs)
    if 'longitude' in ds.coords:
        ds = ds.rename({'Prec': 'pcp', 'Tmin': 't_min', 'Tmax': 't_max'})
    else:
        ds = ds.rename({'Prec': 'pcp', 'Tmin': 't_min', 'Tmax': 't_max'})
    ds['t_mean'] = _calc_t_mean(ds)
    return ds[['t_mean', 'pcp']]
