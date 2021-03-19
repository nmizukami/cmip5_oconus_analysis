# CMIP5 OCONUS hydrologic projection Analysis
A collection of notebooks and tools for analyzing the CMIP5 OCONUS dataset


## Goals:

1.  Evaluate hydrologic changes over Alaska and Hawaii in early, mid, and later 21st century compared to historical period (1970-1999) 
2.	Assess how portrayals of hydrologic projections are different due to GCMs  

## Links

- NCAR Computational Hydrology: https://ncar.github.io/hydrology/

## Acknowledgements

This work is jointly supported by the US Army Corps of Engineers. NCAR is supported by the National Science Foundation.

## References:

- 

## additional info

To use cmip5_oconus modules:

Install the cmip5_oconus by

```bash
cd cmip5_oconus_analysis 
conda activate $ENVIRONMENT_NAME
pip install -e .
```
after that, cmip5_oconus module can be imported

```python
import cmip5_oconus 
```

or

Append to the system environment variable PYTHONPATH the pfafstetter directory 

```bash
PYTHONPATH=$PYTHONPATH:<path_to_cmip5_oconus directory>
export $PYTHONPATH
```
