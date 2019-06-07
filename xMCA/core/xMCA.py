# -*- coding: utf-8 -*-
"""
xMCA
=============

"""

import numpy as np
import xarray as xr

class xMCA:
    """xMCA class
    This class is employed to do the MCA and return assocaited results.

    Define the class amd do decomposition as:

    >>> mca = xMCA(leftField, rightField)
    >>> mca.solver()
    """

    def __init__(self, left, right):

        """Create an MCA object.
        The MCA solution is computed after triggering solver method. Method
        calls are used to retrieve computed quantities.

        **Arguments:**
        *dataset*
            Two `xarray.DataArray` (left and right filed) with two or more
            dimensions containing the data to be analysed. Missing values
            must be np.nan, and they are automatically removed during computation.

        **Returns:**
        *solver*
            An `MCA` instance.

        **Examples:**
        MCA analysis::
        >>> from xMCA import xMCA
        >>> mca = xMCA(left, right)
        """

        if not isinstance(left, xr.DataArray) or not isinstance(right, xr.DataArray):
            raise  TypeError('Left and Right field must be xarray DataArray.')

        self._leftField = self._dropSoloCoords(self._checkFieldName(left, 'left'))
        self._rightField = self._dropSoloCoords(self._checkFieldName(right, 'right'))
        self._timeCoords = self._findTimeCoords()

        if self._leftField[self._timeCoords.name]\
                .identical(self._rightField[self._timeCoords.name]):
            pass
        else:
            raise ValueError('Left and Right field must have same time coords.')

    def _dropSoloCoords(self, dsarray):
        """
        drop sole coords in array

        **Return:**
        *array*
            An `xarray.DataArray` with no solo coordinates.

        """
        coordNames = [name for name, coord in dsarray.coords.items()]
        soloCoords = [n for n in coordNames if dsarray[n].size <= 1]
        return dsarray.drop(soloCoords)

    def _findTimeCoords(self):

        """
        find time coordinates of left and right field.
        Time coordinates satisfy one or more of:
        * Values have a dtype of `numpy.datetime64`.
        * Name of the coordinate is 'time'.
        * The coordinate has an attribute 'axis' with value 'T'.

        **Return:**
        *time_coords*
            time coordinates in xr.DataArray.
        """

        for name in self._leftField.dims:
            coord = self._leftField.coords[name]
            is_time = (np.issubdtype(coord.dtype, np.datetime64) or
                       coord.name == 'time' or
                       coord.attrs.get('axis') == 'T')
            if is_time:
                return coord
            else:
                raise ValueError('cannot find a time coordinate, please make sure '''
                                 'the time coord has name "time" or np.datetime64 dtype')

    def _checkFieldName(self, dsarray, name=''):

        """
        check the name of dsarray, if no name, assign accordingly.

        **Arguments:**
        *array*
        An `xarray.DataArray`.

        **Return:**
        *array*
            An `xarray.DataArray` with name assigned if it has no name.

        """

        if dsarray.name is None:
            dsarray.name = name

        return dsarray


    def _findCoordnames(self, dsarray):

        """
        find additional coordinates of left and right field except time.

        **Return:**
        *other_coords*
            a list of additional coords name except time.
        """

        coordNames = [name for name, coord in dsarray.coords.items()]
        coordNamesWithNoTime = [n for n in coordNames if n != self._timeCoords.name]
        return coordNamesWithNoTime

    def _concatenateExtraDims(self, dsarray):
        """
        concatenate additional coordinates of dsarray such that dsarray only has
        two coordinates. One is time and the other one is concatenated dim.

        **Arguments:**
        *array*
        An `xarray.DataArray`.


        **Return:**
        *array*
            An `xarray.DataArray` with time and concatenated dim as coordinates.
        """

        CoordsWithNoTime = self._findCoordnames(dsarray)
        ConField = dsarray.stack(conDim = tuple(CoordsWithNoTime))

        return ConField

    # def dropna(self, dsarray):
    #     return dsarray.dropna('conDim')

    def _remove_time_mean(self, dsarray):
        """
        remove the time average of dsarray.

        **Arguments:**
        *array*
        An `xarray.DataArray`.

        **Return:**
        *array*
            An `xarray.DataArray` with time mean removed.

        """

        return dsarray - dsarray.mean(dim=self._timeCoords.name)

    def _svd(self, left, right):

        """
        Singular Value Decomposition of left and right field.

        **Arguments:**
        *array*
        left field `xarray.DataArray`;
        right field `xarray.DataArray`;

        **Return:**
        *array*
           U: collumns are sigular vectors of the left field;
           V: collumns are sigular vectors of the right field;
           s: sigular values of covariance matirx of left and right field.

        """

        Cxy = np.dot(left.T, right) / (len(self._timeCoords) - 1.)
        U, s, V = np.linalg.svd(Cxy, full_matrices=False)
        # U column is eigenvectors for left, V column is for right
        V = V.T

        self._U = U
        self._V = V
        self._s = s


    def solver(self):
        """
        solver of the MCA analysis. After calling the solver, the
        Singular Value Decomposition is done at background, then users can
        call *patterns*, *expansionCoefs*, *CovFrac* to retrive the
        sigular vectors, expansition coefficients and explained covariance
        fractions.

        **Examples:**

        Initiate the instance::

        >>> test = xMCA(left, right)

        Call the solver::

        >>> test.solver()
        """

        leftConField = self._concatenateExtraDims(self._leftField)
        rightConField = self._concatenateExtraDims(self._rightField)

        left = self._remove_time_mean(leftConField.dropna('conDim'))
        right = self._remove_time_mean(rightConField.dropna('conDim'))

        left = left.transpose(*[self._timeCoords.name, 'conDim'])
        right = right.transpose(*[self._timeCoords.name, 'conDim'])

        # used for construct patterns
        self._leftConField = left
        self._rightConField = right

        self._svd(left, right)


    def patterns(self, n=1, scaling=True):
        """Sigular vectors of left and right fields (SVs).

        **Optional arguments:**
        *scaling*
            Scaling of the SVs by multiplying the standard deviation of
            corresponding expansion coefficients. Default is True.
        *n*
            Number of SVs to retrive. Defaults to the first SVs.

        **Returns:**
        *SVs*
           Two `xarray.DataArray` containing the SVs. The SVs will be reshaped
           to the same as left and right spatial domains.

        **Examples:**

        Initiate the instance::
            
        >>> mca = xMCA(left, right)

        Call the solver::

        >>> mca.solver()

        Retrive the first two SVs of left and right fields with scaling::
        
        >>> lp, rp = mca.patterns(n=2)

        """

        leftPatterns = self._leftConField[0:n].copy()
        leftPatterns.values = self._U[:,0:n].T
        leftPatterns = leftPatterns.unstack().rename({self._timeCoords.name:'n'})
        leftPatterns['n'] = range(n)
        leftPatterns = leftPatterns.rename('leftPattern')

        rightPatterns = self._rightConField[0:n].copy()
        rightPatterns.values = self._V[:,0:n].T
        rightPatterns = rightPatterns.unstack().rename({self._timeCoords.name:'n'})
        rightPatterns['n'] = range(n)
        rightPatterns = rightPatterns.rename('rightPattern')

        if scaling:
            le, re = self.expansionCoefs(n=n, scale=False)
            leftPatterns = leftPatterns * le.std(self._timeCoords.name)
            rightPatterns = rightPatterns * le.std(self._timeCoords.name)

        # name the patterns and attrs
        leftPatterns.name = 'leftPattern'
        leftPatterns.attrs['long_name'] = 'Signular Vectors for ' + self._leftField.name + '.'
        rightPatterns.name = 'rightPattern'
        rightPatterns.attrs['long_name'] = 'Signular Vectors for ' + self._rightField.name + '.'

        return leftPatterns, rightPatterns

    def expansionCoefs(self, n=1, scale=True):

        """Expansion coefficients of left and right fields (PCs).

        **Optional arguments:**
        *scaling*
            Scaling of the PCs to unit variance by deviding the standard deviation of
            corresponding expansion coefficients. Default is True.
        *n*
            Number of PCs to retrive. Defaults to the first PCs.

        **Returns:**
        *PCs*
           Two `xarray.DataArray` containing the PCs.

        **Examples:**
        Initiate the instance::
        
        >>> mca = xMCA(left, right)

        Call the solver::
        
        >>> mca.solver()

        Retrive the first two PCs of left and right fields with scaling::
        
        >>> le, re = mca.expansionCoefs(n=2)

        """


        le = np.dot(self._U[:,0:n].T,self._leftConField.T)
        le = xr.DataArray(le, coords=[range(n), self._timeCoords], dims=['n', self._timeCoords.name])

        re = np.dot(self._V[:,0:n].T,self._rightConField.T)
        re = xr.DataArray(re, coords=[range(n), self._timeCoords], dims=['n', self._timeCoords.name])

        if scale:
            le, re = le / le.std(self._timeCoords.name), re / re.std(self._timeCoords.name)

            le.name = 'leftExpansionCoefs'
            le.attrs['long_name'] = 'Time expansion coefficients for ' + self._leftField.name + '.'
            re.name = 'rightExpansionCoefs'
            re.attrs['long_name'] = 'Time expansion coefficients for ' + self._rightField.name + '.'


        return le, re

    def covFracs(self, n=1):
        """Fractions of covariance explained by first n modes (FCs).

        **Optional arguments:**
        *n*
            Number of FCs to retrive. Defaults to the first FCs.

        **Returns:**
        *FCs*
           Two `xarray.DataArray` containing the FCs.

        **Examples:**
        Initiate the instance::
        
        >>> mca = xMCA(left, right)

        Call the solver::
        
        >>> mca.solver()

        Retrive the first two FCs with scaling::
            
        >>> le, re = mca.covFracs(n=2)

        """

        frac = self._s[0:n] ** 2 /np.sum(self._s ** 2)
        return xr.DataArray(frac,
                            coords=[range(n)],
                            dims=['n'],
                            name='frac',
                            attrs={'long_name': 'Fractions explained of the covariance matrix between ' + \
                                                self._leftField.name + ' and ' + self._rightField.name + '.'}
                            )

    def sigValues(self, n=10):
        """Singular values of the covariance matrix between left and right field (SGs).

        **Optional arguments:**
        *n*
            Number of SGs to retrive. Defaults to the first ten SGs.

        **Returns:**
        *SGs*
           Two `xarray.DataArray` containing the SGs.

        **Examples:**
        Initiate the instance::
        
        >>> mca = xMCA(left, right)

        Call the solver::
            
        >>> mca.solver()

        Retrive the first two SGs with scaling::
        
        >>> sgs = mca.sigValues(n=2)

        """
        return xr.DataArray(self._s[0:n],
                            coords=[range(n)],
                            dims=['n'],
                            name='sigs',
                            attrs={'long_name': 'Singular values of the covariance matrix between '+\
                                   self._leftField.name + ' and ' + self._rightField.name +'.'})


    def _correlationCalc(self, x, y, correlating=True):

        """
        Calculate the correlations between two fields or time series.

        **Optional arguments:**
        *x*
            time series of x in `xarray.DataArray`.

        *y*
            time series of y in `xarray.DataArray`.

        *correlating*
            Calculating correlations or regression coefficients.
            Default is True (correlations). Otherwise, regression coefficients.

        **Returns:**
        *r*
           An `xarray.DataArray` correlations between x and y.

        """
        x = x - x.mean(dim=self._timeCoords.name)
        y = y - y.mean(dim=self._timeCoords.name)
        xy = x * y  # cov(x,y)
        xy = xy.mean(dim=self._timeCoords.name)
        xx = x.std(dim=self._timeCoords.name)
        yy = y.std(dim=self._timeCoords.name)

        if correlating:
            r = xy / xx / yy
        else:
            r = xy / xx ** 2

        return r


    def homogeneousPatterns(self, n=1, correlating=True):
        """
        Sigular modes expressed as the
        correlation between the expansion coefficient time series (PCs)
        and the corresponding time series of the left or right field at each grid
        point.

        **Optional arguments:**
        *n*
            Number of homogeneous Patternss(HPs) to retrive. Defaults to the first HPs.

        *correlating*
            Expressed as correlation maps. Default is True. Otherwise, expressed as
            regression coefficient maps.

        **Returns:**
        *HPs*
           Two `xarray.DataArray` containing the HPs.

        **Examples:**
        Initiate the instance::
        
        >>> mca = xMCA(left, right)

        Call the solver::
        
        >>> mca.solver()

        Retrive the first two PCs expressed as correlations::
        
        >>> sgs = mca.homogeneousPatterns(n=2)

        """
        le, re = self.expansionCoefs(n=n)

        for f in ['l', 'r']:
            r = []
            for i in range(n):
                if f == 'l':
                    tmp = self._correlationCalc(self._leftField, le.sel(n=i))
                else:
                    tmp = self._correlationCalc(self._rightField, re.sel(n=i))

                tmp['n'] = i
                r.append(tmp)

            if f == 'l':
                lHP = xr.concat(r, dim='n')
                # lHP.name = 'lHP'
            else:
                rHP = xr.concat(r, dim='n')
                # rHP.name = 'rHP'

        lHP.name = 'leftHomoPatterns'
        lHP.attrs['long_name'] = 'Homogeneous Patterns for ' + self._leftField.name + '.'
        rHP.name = 'rightHomoPatterns'
        rHP.attrs['long_name'] = 'Homogeneous Patterns for ' + self._rightField.name + '.'

        return lHP, rHP



    def heterogeneousPatterns(self, n=1, correlating=True):
        """
        Sigular modes expressed as the
        correlation between the expansion coefficient time series (PCs)
        and the corresponding time series of the left or right field at each grid
        point.

        **Optional arguments:**
        *n*
            Number of homogeneous Patternss(HPs) to retrive. Defaults to the first HPs.

        *correlating*
            Expressed as correlation maps. Default is True. Otherwise, expressed as
            regression coefficient maps.

        **Returns:**
        *HPs*
           Two `xarray.DataArray` containing the HPs.

        **Examples:**
        Initiate the instance::
        
        >>> mca = xMCA(left, right)

        Call the solver::
        
        >>> mca.solver()

        Retrive the first two PCs expressed as correlations::
        
        >>> sgs = mca.homogeneousPatterns(n=2)

        """
        le, re = self.expansionCoefs(n=n)


        for f in ['l', 'r']:
            r = []
            for i in range(n):
                if f == 'l':
                    tmp = self._correlationCalc(self._leftField, re.sel(n=i))
                else:
                    tmp = self._correlationCalc(self._rightField, le.sel(n=i))

                tmp['n'] = i
                r.append(tmp)

            if f == 'l':
                lHP = xr.concat(r, dim='n')
                # lHP.name = 'lHP'
            else:
                rHP = xr.concat(r, dim='n')
                # rHP.name = 'rHP'

        lHP.name = 'leftHeteroPatterns'
        lHP.attrs['long_name'] = 'Heterogeneous Patterns for ' + self._leftField.name + '.'
        rHP.name = 'rightHeteroPatterns'
        rHP.attrs['long_name'] = 'Heterogeneous Patterns for ' + self._rightField.name + '.'

        return lHP, rHP

