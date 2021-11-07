[![Documentation Status](https://readthedocs.org/projects/xmca/badge/?version=latest)](https://xmca.readthedocs.io/en/latest/?badge=latest)

# About xMCA
xMCA, inspired from EOFs https://github.com/ajdawson/eofs,  is Maximum Covariance Analysis (sometimes also called SVD) in xarray. 

API Documentation: http://xmca.readthedocs.io/

# How to install
### via git
```
git clone https://github.com/Yefee/xMCA.git
cd xMCA
python setup.py install
```
# Future Plan
In next version, A Monte Carlo method will be added for statistical test.

# Example 1
MCA analysis for US surface air temperature and SST over the Pacific.
This example is taken from https://atmos.washington.edu/~breth/classes/AS552/matlab/lect/html/MCA_PSSTA_USTA.html


```python
from xMCA import xMCA
import xarray as xr
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
usta = xr.open_dataarray('xMCA/examples/data/USTA.nc').transpose(*['time', 'lat', 'lon'])
usta.name = 'USTA'
print(usta)
```

    <xarray.DataArray 'USTA' (time: 396, lat: 5, lon: 12)>
    array([[[-0.450303, -0.734848, ..., -4.270303, -2.69697 ],
            [ 1.066061,  2.691515, ..., -4.947273, -3.330303],
            ...,
            [      nan, -0.342424, ...,       nan,       nan],
            [      nan,       nan, ...,       nan,       nan]],
    
           [[ 1.524545,  1.370606, ..., -1.430303,  0.048485],
            [ 1.366364,  2.497273, ..., -0.593939, -0.079697],
            ...,
            [      nan,  0.695455, ...,       nan,       nan],
            [      nan,       nan, ...,       nan,       nan]],
    
           ...,
    
           [[ 1.077879,  0.630303, ..., -1.262727, -1.496364],
            [ 1.020606,  0.114848, ..., -0.786667, -0.573939],
            ...,
            [      nan,  1.65    , ...,       nan,       nan],
            [      nan,       nan, ...,       nan,       nan]],
    
           [[ 1.768182,  2.807879, ...,  0.885758,  0.618182],
            [ 1.555152,  3.435152, ..., -0.416667,  0.185152],
            ...,
            [      nan,  0.012121, ...,       nan,       nan],
            [      nan,       nan, ...,       nan,       nan]]])
    Coordinates:
      * lat      (lat) float64 47.5 42.5 37.5 32.5 27.5
      * lon      (lon) float64 -122.5 -117.5 -112.5 -107.5 ... -77.5 -72.5 -67.5
      * time     (time) int64 0 1 2 3 4 5 6 7 8 ... 388 389 390 391 392 393 394 395



```python
sstpc = xr.open_dataarray('xMCA/examples/data/SSTPac.nc').transpose(*['time', 'lat', 'lon'])
sstpc.name = 'SSTPC'
print(sstpc)
```

    <xarray.DataArray 'SSTPC' (time: 396, lat: 30, lon: 84)>
    [997920 values with dtype=float64]
    Coordinates:
      * lat      (lat) int16 -29 -27 -25 -23 -21 -19 -17 ... 17 19 21 23 25 27 29
      * lon      (lon) uint16 124 126 128 130 132 134 ... 280 282 284 286 288 290
      * time     (time) int64 0 1 2 3 4 5 6 7 8 ... 388 389 390 391 392 393 394 395


### Decompsition and retrieve the first and second loadings and expansion coefficeints 


```python
'''
decomposition, time should be in the first axis
lp is for SSTPC
rp is for USTA
'''

sst_ts = xMCA(sstpc, usta)
sst_ts.solver()
lp, rp = sst_ts.patterns(n=2)
le, re = sst_ts.expansionCoefs(n=2)
frac = sst_ts.covFracs(n=2)
print(frac)
```

    <xarray.DataArray 'frac' (n: 2)>
    array([0.407522, 0.391429])
    Coordinates:
      * n        (n) int64 0 1
    Attributes:
        long_name:  Fractions explained of the covariance matrix between SSTPC an...



```python
fig, (ax1, ax2) = plt.subplots(2, 2, figsize=(12, 5))
lp[0].plot(ax=ax1[0])
le[0].plot(ax=ax1[1])

rp[0].plot(ax=ax2[0])
re[0].plot(ax=ax2[1])
```



![png](https://github.com/Yefee/xMCA/blob/master/xMCA/examples/example_files/example_6_1.png)


### Homogeneous and heterogeneous regression


```python
lh, rh = sst_ts.homogeneousPatterns(n=1)
le, re = sst_ts.heterogeneousPatterns(n=1)
```


```python
fig, (ax1, ax2) = plt.subplots(2, 2, figsize=(12, 5))
lh[0].plot(ax=ax1[0])
rh[0].plot(ax=ax1[1])

le[0].plot(ax=ax2[0])
re[0].plot(ax=ax2[1])
```


![png](https://github.com/Yefee/xMCA/blob/master/xMCA/examples/example_files/example_9_1.png)


### Two-tailed Student-t test
```python
le, re, lphet, rphet = sst_ts.heterogeneousPatterns(n=1, statistical_test=True)
```

```python
fig, (ax1, ax2) = plt.subplots(2, 2, figsize=(12, 5))
le[0].plot(ax=ax1[0])
re[0].plot(ax=ax1[1])

# Only plot where p<0.01
lphet[0].where(lphet[0]<0.01).plot(ax=ax2[0])
rphet[0].where(rphet[0]<0.01).plot(ax=ax2[1])
```
![png](https://github.com/Yefee/xMCA/blob/master/xMCA/examples/example_files/example_statistical_test_1.png)




# Example 2
EOF analysis for US surface air temperature and SST over the Pacific
This example is taken from https://atmos.washington.edu/~breth/classes/AS552/matlab/lect/html/MCA_PSSTA_USTA.html


```python
from xMCA import xMCA
import xarray as xr
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
sstpc = xr.open_dataarray('data/SSTPac.nc').transpose(*['time', 'lat', 'lon'])
sstpc.name = 'SSTPC'
print(sstpc)
```

    <xarray.DataArray 'SSTPC' (time: 396, lat: 30, lon: 84)>
    [997920 values with dtype=float64]
    Coordinates:
      * lat      (lat) int16 -29 -27 -25 -23 -21 -19 -17 ... 17 19 21 23 25 27 29
      * lon      (lon) uint16 124 126 128 130 132 134 ... 280 282 284 286 288 290
      * time     (time) int64 0 1 2 3 4 5 6 7 8 ... 388 389 390 391 392 393 394 395


Decompsition and retrieve the first and second loadings and expansion coefficeints 


```python
'''
decomposition, time should be in the first axis
lp is for SSTPC
rp is for USTA
'''

sst_ts = xMCA(sstpc, sstpc.rename('SSTPC_copy'))
sst_ts.solver()
lp, _ = sst_ts.patterns(n=2)
le, _ = sst_ts.expansionCoefs(n=2)
frac = sst_ts.covFracs(n=2)
print(frac)
```

    <xarray.DataArray 'frac' (n: 2)>
    array([0.873075, 0.04946 ])
    Coordinates:
      * n        (n) int64 0 1
    Attributes:
        long_name:  Fractions explained of the covariance matrix between SSTPC an...



```python
fig, ax1 = plt.subplots(1, 2, figsize=(12, 5))
lp[0].plot(ax=ax1[0])
le[0].plot(ax=ax1[1])

```



![png](https://github.com/Yefee/xMCA/blob/master/xMCA/examples/example_files/example_eof_5_1.png)


Regress PC1 to the original SST field


```python
lh, _ = sst_ts.homogeneousPatterns(n=1)
```



```python
fig, ax1= plt.subplots()
lh[0].plot(ax=ax1)

```



![png](https://github.com/Yefee/xMCA/blob/master/xMCA/examples/example_files/example_eof_8_1.png)

[![Documentation Status](https://readthedocs.org/projects/xmca/badge/?version=latest)](https://xmca.readthedocs.io/en/latest/?badge=latest)

# About xMCA
xMCA, inspired from EOFs https://github.com/ajdawson/eofs,  is Maximum Covariance Analysis (sometimes also called SVD) in xarray. 

API Documentation: http://xmca.readthedocs.io/

# How to install
### via git
```
git clone https://github.com/Yefee/xMCA.git
cd xMCA
python setup.py install
```
# Future Plan
In next version, A Monte Carlo method will be added for statistical test.

# Example 1
MCA analysis for US surface air temperature and SST over the Pacific.
This example is taken from https://atmos.washington.edu/~breth/classes/AS552/matlab/lect/html/MCA_PSSTA_USTA.html


```python
from xMCA import xMCA
import xarray as xr
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
usta = xr.open_dataarray('xMCA/examples/data/USTA.nc').transpose(*['time', 'lat', 'lon'])
usta.name = 'USTA'
print(usta)
```

    <xarray.DataArray 'USTA' (time: 396, lat: 5, lon: 12)>
    array([[[-0.450303, -0.734848, ..., -4.270303, -2.69697 ],
            [ 1.066061,  2.691515, ..., -4.947273, -3.330303],
            ...,
            [      nan, -0.342424, ...,       nan,       nan],
            [      nan,       nan, ...,       nan,       nan]],
    
           [[ 1.524545,  1.370606, ..., -1.430303,  0.048485],
            [ 1.366364,  2.497273, ..., -0.593939, -0.079697],
            ...,
            [      nan,  0.695455, ...,       nan,       nan],
            [      nan,       nan, ...,       nan,       nan]],
    
           ...,
    
           [[ 1.077879,  0.630303, ..., -1.262727, -1.496364],
            [ 1.020606,  0.114848, ..., -0.786667, -0.573939],
            ...,
            [      nan,  1.65    , ...,       nan,       nan],
            [      nan,       nan, ...,       nan,       nan]],
    
           [[ 1.768182,  2.807879, ...,  0.885758,  0.618182],
            [ 1.555152,  3.435152, ..., -0.416667,  0.185152],
            ...,
            [      nan,  0.012121, ...,       nan,       nan],
            [      nan,       nan, ...,       nan,       nan]]])
    Coordinates:
      * lat      (lat) float64 47.5 42.5 37.5 32.5 27.5
      * lon      (lon) float64 -122.5 -117.5 -112.5 -107.5 ... -77.5 -72.5 -67.5
      * time     (time) int64 0 1 2 3 4 5 6 7 8 ... 388 389 390 391 392 393 394 395



```python
sstpc = xr.open_dataarray('xMCA/examples/data/SSTPac.nc').transpose(*['time', 'lat', 'lon'])
sstpc.name = 'SSTPC'
print(sstpc)
```

    <xarray.DataArray 'SSTPC' (time: 396, lat: 30, lon: 84)>
    [997920 values with dtype=float64]
    Coordinates:
      * lat      (lat) int16 -29 -27 -25 -23 -21 -19 -17 ... 17 19 21 23 25 27 29
      * lon      (lon) uint16 124 126 128 130 132 134 ... 280 282 284 286 288 290
      * time     (time) int64 0 1 2 3 4 5 6 7 8 ... 388 389 390 391 392 393 394 395


### Decompsition and retrieve the first and second loadings and expansion coefficeints 


```python
'''
decomposition, time should be in the first axis
lp is for SSTPC
rp is for USTA
'''

sst_ts = xMCA(sstpc, usta)
sst_ts.solver()
lp, rp = sst_ts.patterns(n=2)
le, re = sst_ts.expansionCoefs(n=2)
frac = sst_ts.covFracs(n=2)
print(frac)
```

    <xarray.DataArray 'frac' (n: 2)>
    array([0.407522, 0.391429])
    Coordinates:
      * n        (n) int64 0 1
    Attributes:
        long_name:  Fractions explained of the covariance matrix between SSTPC an...



```python
fig, (ax1, ax2) = plt.subplots(2, 2, figsize=(12, 5))
lp[0].plot(ax=ax1[0])
le[0].plot(ax=ax1[1])

rp[0].plot(ax=ax2[0])
re[0].plot(ax=ax2[1])
```



![png](https://github.com/Yefee/xMCA/blob/master/xMCA/examples/example_files/example_6_1.png)


### Homogeneous and heterogeneous regression


```python
lh, rh = sst_ts.homogeneousPatterns(n=1)
le, re = sst_ts.heterogeneousPatterns(n=1)
```


```python
fig, (ax1, ax2) = plt.subplots(2, 2, figsize=(12, 5))
lh[0].plot(ax=ax1[0])
rh[0].plot(ax=ax1[1])

le[0].plot(ax=ax2[0])
re[0].plot(ax=ax2[1])
```


![png](https://github.com/Yefee/xMCA/blob/master/xMCA/examples/example_files/example_9_1.png)


### Two-tailed Student-t test
```python
le, re, lphet, rphet = sst_ts.heterogeneousPatterns(n=1, statistical_test=True)
```

```python
fig, (ax1, ax2) = plt.subplots(2, 2, figsize=(12, 5))
le[0].plot(ax=ax1[0])
re[0].plot(ax=ax1[1])

# Only plot where p<0.01
lphet[0].where(lphet[0]<0.01).plot(ax=ax2[0])
rphet[0].where(rphet[0]<0.01).plot(ax=ax2[1])
```
![png](https://github.com/Yefee/xMCA/blob/master/xMCA/examples/example_files/example_statistical_test_1.png)




# Example 2
EOF analysis for US surface air temperature and SST over the Pacific
This example is taken from https://atmos.washington.edu/~breth/classes/AS552/matlab/lect/html/MCA_PSSTA_USTA.html


```python
from xMCA import xMCA
import xarray as xr
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
sstpc = xr.open_dataarray('data/SSTPac.nc').transpose(*['time', 'lat', 'lon'])
sstpc.name = 'SSTPC'
print(sstpc)
```

    <xarray.DataArray 'SSTPC' (time: 396, lat: 30, lon: 84)>
    [997920 values with dtype=float64]
    Coordinates:
      * lat      (lat) int16 -29 -27 -25 -23 -21 -19 -17 ... 17 19 21 23 25 27 29
      * lon      (lon) uint16 124 126 128 130 132 134 ... 280 282 284 286 288 290
      * time     (time) int64 0 1 2 3 4 5 6 7 8 ... 388 389 390 391 392 393 394 395


Decompsition and retrieve the first and second loadings and expansion coefficeints 


```python
'''
decomposition, time should be in the first axis
lp is for SSTPC
rp is for USTA
'''

sst_ts = xMCA(sstpc, sstpc.rename('SSTPC_copy'))
sst_ts.solver()
lp, _ = sst_ts.patterns(n=2)
le, _ = sst_ts.expansionCoefs(n=2)
frac = sst_ts.covFracs(n=2)
print(frac)
```

    <xarray.DataArray 'frac' (n: 2)>
    array([0.873075, 0.04946 ])
    Coordinates:
      * n        (n) int64 0 1
    Attributes:
        long_name:  Fractions explained of the covariance matrix between SSTPC an...



```python
fig, ax1 = plt.subplots(1, 2, figsize=(12, 5))
lp[0].plot(ax=ax1[0])
le[0].plot(ax=ax1[1])

```



![png](https://github.com/Yefee/xMCA/blob/master/xMCA/examples/example_files/example_eof_5_1.png)


Regress PC1 to the original SST field


```python
lh, _ = sst_ts.homogeneousPatterns(n=1)
```



```python
fig, ax1= plt.subplots()
lh[0].plot(ax=ax1)

```



![png](https://github.com/Yefee/xMCA/blob/master/xMCA/examples/example_files/example_eof_8_1.png)


# Publications using xMCA

He, C., Liu, Z., Otto-Bliesner, B. L., Brady, E. C., Zhu, C., Tomas, R., Clark, P. U., Zhu, J., Jahn, A., &Gu, S. (2021), Hydroclimate footprint of pan-Asian monsoon water isotope during the last deglaciation, Science Advances, 7(4), eabe2611.

He, C., Liu, Z., Otto-Bliesner, B. L., Brady, E. C., Zhu, C., Tomas, R., Buizert, C., &Severinghaus, J. P. (2021), Abrupt Heinrich Stadial 1 cooling missing in Greenland oxygen isotopes, Science Advances, 7(25), eabh1007.

He, C., Liu, Z., Otto-Bliesner, B. L., Brady, E. C., Zhu, C., Tomas, R., Gu, S., Han, J., &Jin, Y. (2021), Deglacial variability of South China hydroclimate heavily contributed by autumn rainfall, Nature communications, 12(1), 1-9.
