# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2020 Daniel Thompson
# Copyright (C) 2022 John C Peterson

"""
Photoplethysmogram (PPG) Signal Processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Algorithms and signal processing primatives that can be used to convert
the raw PPG samples into heart rate estimates, and other useful things.
"""

import array
import micropython
import watch
import gc

# HRS in 14-bit mode, and ALS in 14-bit
#_POLLING_PERIOD_MS = 24.0
# See make_filters.py for filter coefficients
#_B_HP = 
#_A_HP = 
#_B_HP = 
#_A_HP = 

# HRS in 15-bit mode, and ALS in 15-bit
_POLLING_PERIOD_MS = 48.0
# Type II Tchebyshev HP filter, order 5, > 40 dB supression in stop band
_B_HP = [ 0.84342327297782, -4.20986531397658,  8.41249205009943,
         -8.41249205009943,  4.20986531397658, -0.84342327297782 ]
_A_HP = [ 1.00000000000000, -4.65147387589849,  8.67327978498012,
         -8.10293359342740,  3.79251120240145, -0.71136281740019 ]
# Type II Tchebyshev LP filter, order 3, > 20 dB supression in stop band
_B_LP = [ 0.31811996128135,  0.70768843054069,
          0.70768843054069,  0.31811996128135 ]
_A_LP = [ 1.00000000000000,  0.47477514136428,
          0.48905358959275,  0.08778805268705 ]


_LEN_DEBUG_BUFFER = const(256)	# length of the debug, raw data buffer (samples)
_LEN_HRM_BUFFER_SEC = 10.0	# length of the y buffer for HR estimates (sec)
_INIT_HARDWARE_SEC = 1.5	# time to wait for the hardware to init (sec)
_INIT_FILTER_SEC = 1.0		# time to wait for the filters to converge (sec)

_MIN_HR_BPM = 30.0		# minimum expected heart rate (bpm)
_MAX_HR_BPM = 200.0		# maximum expected heart rate (bpm)
_STOP_ANAL_PCT = 0.50		# fraction of this peak to stop peak analysis
_MAX_ITER_KMEANS = 2		# number of iterations of the K-means algorithm
_ACF_THRESH = 0.25		# Don't return anything for ACF estimates < this


class IIR_Filter():
    """
    General IIR digital filter class
    """

    def __init__(self, b_coeff, a_coeff):
        self._b_coeff = array.array('d', b_coeff)
        self._a_coeff = array.array('d', a_coeff)
        self._M = len(b_coeff)
        self._N = len(a_coeff)
        # _x and _y are "ring" buffers
        self._x = array.array('d', [0.0]*self._M)
        self._y = array.array('d', [0.0]*self._N)
        self._last_x = -1
        self._last_y = 0

    @micropython.native
    def step(self, x_n):
        """
        Compute the output of the IIR filter for the new sample x_n

        a[0]*y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[M]*x[n-M]
                              - a[1]*y[n-1] - ... - a[N]*y[n-N]
        """

        b_coeff = memoryview(self._b_coeff)
        a_coeff = memoryview(self._a_coeff)
        _M = self._M
        _N = self._N

        x = self._x
        y = self._y
        last_x = self._last_x
        last_y = self._last_y

        # If this is the first call, initialize the x, y history buffers

        if last_x < 0:
            for j in range(0, _M):
                x[j] = x_n
                last_x = 0
            for j in range(0, _N):
                y[j] = 1.0
                last_y = 0

        # Update the x history with the most recent sample x_n (ring buffer)

        last_x = (last_x+1) % _M

        x[last_x] = x_n

        # Compute the new filter output y_n

        y_n = 0.0

        for j in range(0, _M):
            k = (last_x+_M-j) % _M
            y_n += b_coeff[j] * x[k]

        for j in range(0, _N-1):
            k = (last_y+_N-j) % _N
            y_n -= a_coeff[j+1] * y[k]

        if a_coeff[0] != 1.0:
            y_n /= a_coeff[0]

        # Update the y history with the most recent value (ring buffer)

        last_y = (last_y+1) % _N

        y[last_y] = y_n

        self._last_x = last_x
        self._last_y = last_y

        return y_n


class Heart_Rate():
    """
    Heart rate estimator class

    The sample history buffer is partitioned into two blocks: the first
    block (of length _min_lag) contains older samples that are used for
    computing the (partial) sample auto correlation function, followed by
    a block that contains the most recent samples (of length buffer_size_sec).
    """

    def __init__(self, buffer_size_sec):
        polling_period_sec = _POLLING_PERIOD_MS / 1000.0
        self._tau_s = polling_period_sec
        self._min_lag = round((60.0/_MAX_HR_BPM)/polling_period_sec)
        self._max_lag = round((60.0/_MIN_HR_BPM)/polling_period_sec)
        self._init_hard = round(_INIT_HARDWARE_SEC/polling_period_sec)
        self._init_filt = round(_INIT_FILTER_SEC/polling_period_sec)
        self._y_buf_size = round(buffer_size_sec/polling_period_sec)
        self._t_buf_size = self._y_buf_size + self._max_lag
        self._min_lag_h = round(0.5*self._min_lag)
        self._max_lag_h = round(0.5*self._max_lag)
        self._initialization = True
        self._spinner_n = -1

        self._wave_list = [ ]
        self._w_min = array.array('f', [0.0]*self._max_lag_h)
        self._autocorr = array.array('f', [0.0]*self._max_lag)

        self._y_buf = array.array('f', [0.0]*self._t_buf_size)
        self._y_buf_len = 0		# length of the buffer (inc. past data)
        self._y_buf_slen = 0.0		# arc length of the y_n samples
        self._y_buf_min = 0.0		# minimum of the y_n samples
        self._y_buf_max = 0.0		# maximum of the y_n samples
        self._y_buf_sum = 0.0		# summation of the y_n samples
        self._y_buf_pwr = 0.0		# summation of the y_n^2 samples
        self._no_samples = 0

        self._i_mean = 0.0	# mean of the buffer indices of the active data
        self._i_sum = 0		# sum of the buffer indices of the active data
        self._i2_sum = 0	# sum of the buffer indices^2 of the active data
        self._iy_sum = 0.0	# sum of the product of buffer index * y_n

        # Initialize some of the (one time) calculations for linear detrending
        self._i_mean = 0.5*(self._max_lag + self._t_buf_size - 1)
        self._i_sum =      (self._max_lag + self._t_buf_size - 1) \
                          * self._y_buf_size // 2

        for j in range(self._max_lag, self._t_buf_size):
            self._i2_sum += j**2


    def step(self, x_n_raw, y_n):
        """
        Update the internal data structures with a new (filtered) sample y_n
        """

        _max_lag = self._max_lag
        _min_lag_h = self._min_lag_h
        _max_lag_h = self._max_lag_h
        _y_buf_len = self._y_buf_len
        _t_buf_size = self._t_buf_size

        _wave_list = self._wave_list
        _autocorr = self._autocorr
        _y_buf = memoryview(self._y_buf)
        _w_min = memoryview(self._w_min)

        # Keep track of how many samples have been processed (for init phase)
        self._no_samples += 1

        # Update the DC offset at the end of hardware initialization phase
        if self._initialization and self._no_samples == self._init_hard:
            self._dc_offset = int(x_n_raw)

        # Save the new filtered sample in the y "ring" buffer
        _y_buf[_y_buf_len] = y_n

        # Update the sample statistics
        if _y_buf_len > _max_lag:
            self._y_buf_slen += abs(y_n - _y_buf[_y_buf_len-1])
            if y_n < self._y_buf_min:
                self._y_buf_min = y_n
            if y_n > self._y_buf_max:
                self._y_buf_max = y_n
            self._y_buf_sum += y_n
            self._y_buf_pwr += y_n * y_n
            self._iy_sum += _y_buf_len*y_n
        elif _y_buf_len == _max_lag:
            # initialize the sample statistics
            self._y_buf_slen = 0.0
            self._y_buf_min = y_n
            self._y_buf_max = y_n
            self._y_buf_sum = y_n
            self._y_buf_pwr = y_n * y_n
            self._iy_sum = _y_buf_len*y_n

        # Update the buffer length
        _y_buf_len += 1

        # Update the sample auto-correlation summations for all lags
        if _y_buf_len > _max_lag:
            for j in range(_max_lag):
                _autocorr[j] += y_n * _y_buf[_y_buf_len-j-2]

        # Check for a Systolic peak, a more thorough happens analysis later
        if (_y_buf_len > _max_lag) and (_y_buf_len < _t_buf_size):
            # We check for possible Systolic peaks by looking for a change
            # in the sign of the first derivative that marks the peak. We
            # then characterize the upslope, rising wave into that point.
            # The following quantities are identified;
            #   w_beg - buffer index of the first point of the up wave
            #   w_end - buffer index of the last point of the up wave
            #   w_min - y value at the first point of the wave (min, valley)
            #   w_max - y value at the last point of the wave (max, peak)
            #
            # When the buffer is full, the heart rate estimator method
            # will do a lot of computationally demanding work and we
            # don't want to do so much work that the next sample gets
            # dropped. So we skip these calculations when the latest
            # sample completely fills the buffer. We can't validate a
            # candidate Systolic peak near the end of the buffer anyway,
            # as there isn't enough data to verify the values drop off
            # after the candidate peak as they should. So nothing lost.
            #
            w_end = _y_buf_len-3
            w_max = _y_buf[w_end]
            if (_y_buf[w_end-2] <  w_max) and \
               (_y_buf[w_end-1] <  w_max) and \
               (_y_buf[w_end+1] <= w_max) and \
               (_y_buf[w_end+2] <= w_max):
                   # The 1-st derivative changes sign at the 3-rd to last
                   # sample in the buffer. Compute the stats characterizing
                   # the rising wave preceding that point.
                   w_beg = 0	# wave start relative to end
                   no_pnts = _max_lag_h	# number of points analyzed
                   _w_min[0] = w_max
                   # working backward, characterize the rising wave
                   for j in range(1, _max_lag_h):
                       # update the minimum of y in wave from here to peak
                       if _y_buf[w_end-j] < _w_min[j-1]:
                          # a new low point in the wave
                          _w_min[j] = _y_buf[w_end-j]
                          w_beg = j
                       else:
                          # no change, use last min
                          _w_min[j] = _w_min[j-1]
                       # decide if we have looked back far enough in time
                       if j > _min_lag_h:
                           thresh = (1.0 - _STOP_ANAL_PCT)*_w_min[j] \
                                         + _STOP_ANAL_PCT*w_max;
                           # check for a significant rise off the lowest point
                           if _y_buf[w_end-j] > thresh:
                               no_pnts = j + 1
                               break
                   # extract info from the point at the start of the wave
                   w_min = _w_min[w_beg]
                   # compute buffer index of the point at the start of wave
                   w_beg = w_end - w_beg

                   # check if any previous up waves are contained in this one
                   for j in range(len(_wave_list)):
                       [w_beg_l, w_end_l, w_min_l, w_max_l] = _wave_list[-1]
                       if w_beg_l >= w_beg:
                           # yes, remove the last wave 
                           _wave_list.pop()
                       else:
                           # that's all
                           break
                   # check for a reasonable number of samples in the up wave
                   if (w_end - w_beg) >= _min_lag_h:
                       # add this up wave to the candidate Systolic peak list
                       self._wave_list.append([w_beg, w_end, w_min, w_max])

        self._y_buf_len = _y_buf_len 


    def k_means(self, period_list):
        """
        Detect beat periods (from missing peaks) that are multiples of the
        fundamental using the K-means algorithm with k=2 and appropriate
        initial conditions. If the converged centroids are far enough
        apart, there are probably multiple beat periods in the set.
        """

        # Check for enough periods to continue. When there are a very
        #   small number of periods, there can be multiple missed beats
        #   making them dominant over the fundamental beats.
        if len(period_list) < 4:
            return period_list

        # Use the median for the 1-st centroid, max for the 2-nd centroid

        period_list.sort(reverse=False)
        centroid_1 = period_list[int(0.5*len(period_list))]
        centroid_2 = max(period_list)

        # Run a few iterations of the k=2 model

        for j_iter in range(_MAX_ITER_KMEANS):
            period_list_1 = [ ]
            period_list_2 = [ ]
            # assign each period to the closest centroid
            for j_period in range(len(period_list)):
                period = period_list[j_period]
                if abs(period - centroid_2) > abs(period - centroid_1):
                    # centroid 1 is closer
                    period_list_1.append(period)
                else:
                    # centroid 2 is closer
                    period_list_2.append(period)
            # check for empty lists (this can happen when they are all equal)
            if len(period_list_1) == 0 or len(period_list_2) == 0:
                # one cluster won, nothing to do
                break
            # recompute the centroids (use the median value for both)
            centroid_1 = period_list_1[int(0.5*len(period_list_1))]
            centroid_2 = period_list_2[int(0.5*len(period_list_2))]

        # Look for significant separation of the centroids
        if len(period_list_1) > len(period_list_2) and len(period_list_2) > 0:
            if min(period_list_2) > 1.5*max(period_list_1):
                # there is a good chance that there are multiple beats
                return period_list_1

        return period_list


    def get_heart_rate(self):
        """
        Compute a new heart rate estimate when enough data is available
        """

        # Wait here for the initialization phases to complete
        if self._initialization:
            if self._y_buf_len == self._t_buf_size:
                self.reset_buffer(False)
            if self._no_samples < self._init_hard:
                return None
            elif self._no_samples == self._init_hard:
                # end of hardware init phase, totally clear the buffer
                self.reset_buffer(True)
                return None
            elif self._no_samples > (self._init_hard + self._init_filt):
                # digital filters should have converged by now
                self.reset_buffer(False)
                self._initialization = False
                return None

        # If the buffer is not full, then nothing to do, nothing to report
        if self._y_buf_len < self._t_buf_size:
            return None

        # Default return values in case of any failed computations

        hr_min = None
        hr_max = None
        hr_mean = None
        hr_median = None
        hr_acf = None
        acf_val = None
        d_qual = None

        ## Compute the data quality metrics for this buffer

        # Compute the sample mean and variance of y

        no_samples = self._y_buf_size
        y_mean = self._y_buf_sum / no_samples
        biased_var = (self._y_buf_pwr / no_samples) - y_mean**2
        y_var = biased_var * no_samples / (no_samples-1)

        # Compute the quality index that quantifies the "filling" of the y range

        y_fill = self._y_buf_slen / (self._y_buf_max - self._y_buf_min)

        d_qual = [ y_mean, y_var, y_fill ]

        # Compute the slope, intercept of the lsq fit of a straight line to y_n
        #   We don't detrend all the data as that wastes resources, but just
        #   compute stuff "on the fly" when it's needed to save CPU and memory.

        b_0 = (no_samples*self._iy_sum - self._i_sum*self._y_buf_sum) \
            / (no_samples*self._i2_sum - self._i_sum*self._i_sum)

        a_0 = y_mean - b_0*self._i_mean

        ## Estimate the heart rate using the Systolic peak picker approach

        _tau_s = self._tau_s
        _min_lag = self._min_lag
        _max_lag = self._max_lag
        _y_buf_len = self._y_buf_len

        _wave_list = self._wave_list
        _autocorr = memoryview(self._autocorr)
        _y_buf = memoryview(self._y_buf)	# many max() calls with slices!

        use_wave = bytearray(len(_wave_list))

        for j in range(len(_wave_list)):
            # default is to accept the up wave and peak
            use_wave[j] = True

        no_peaks = [len(_wave_list)]

        # Discard rising waves that have a small range (relative to buffer var)

        for j in range(len(_wave_list)):
            if use_wave[j]:
                [w_beg, w_end, w_min, w_max] = _wave_list[j]
                if (w_max - w_min)**2 < 1.500*y_var:
                    use_wave[j] = False

        no_peaks.append(sum(use_wave))

        # Discard local peaks that are not above the linear trend line

        for j in range(len(_wave_list)):
            [w_beg, w_end, w_min, w_max] = _wave_list[j]
            # the maximum value should be above the linear trend line
            if w_max < (a_0 + b_0*w_end):
                use_wave[j] = False

        no_peaks.append(sum(use_wave))

        # Discard rising waves that don't fall off *after* the peak

        for j in range(len(_wave_list)):
            if use_wave[j]:
                # samples following the (local) peak should be smaller
                [w_beg, w_end, w_min, w_max] = _wave_list[j]
                no_pnts = w_end - w_beg
                if (w_end+no_pnts) < _y_buf_len:
                    if max(_y_buf[w_end+1:w_end+no_pnts+1]) > w_max:
                        use_wave[j] = False
                else:
                    # not enough points to really be sure, discard it to be safe
                    use_wave[j] = False

        no_peaks.append(sum(use_wave))

        # Append the peak counts to the data quality metrics list:
        #   initial number of local peaks from the step() method (w/full buffer)
        #   after discarding rising waves that have a small range (rel to var)
        #   after discarding local peaks less than the lsq linear trend line
        #   after discarding rising waves that don't drop off *after* the peak

        d_qual.append(no_peaks)

        # Extract the accepted peaks

        peak_list = [ ]
        for j in range(len(_wave_list)):
            if use_wave[j]:
                [w_beg, w_end, w_min, w_max] = _wave_list[j]
                peak_list.append(w_end)

        ## Process the peaks into beat period estimates

        # First check for a meaningful number of local maxima (avoid GIGO)

        if len(peak_list) > 4:

            # compute the time intervals, periods between the accepted peaks

            period_list = [ ]
            for j in range(1,len(peak_list)):
                diff = peak_list[j] - peak_list[j-1]
                if diff > _min_lag and diff < _max_lag:
                    period_list.append(diff)

            # discard multiple beats (from missed peaks) with K-means algorithm

            period_list = self.k_means(period_list)

            # final processing of the accepted periods

            if len(period_list) > 1:
                # find the min, max, median of the periods
                hr_min = 60.0 / (max(period_list)*_tau_s)
                hr_max = 60.0 / (min(period_list)*_tau_s)
                # compute the median heart rate
                period_list.sort(reverse=False)
                s_median = period_list[int(0.5*len(period_list))]
                hr_median = 60.0 / (s_median*_tau_s)
                # compute the mean heart rate
                s_mean = sum(period_list) / len(period_list)
                hr_mean = 60.0 / (s_mean*_tau_s)

        ## Estimate the heart rate using the sample auto-correlation function

        # The sample auto-correlation function will have a local maximum
        #   at a lag associated with any non-trivial low frequency noise,
        #   (such as from motion). So it can produce confusing results
        #   in the absence of other information. But, if we have a good
        #   beat estimate from the peak picker, we can look for a local
        #   maximum in the ACF near that lag and report that. It can
        #   provide better accuracy for faster beat rates, by looking for
        #   local maximum associated with integer multiples (2,3...) of
        #   the fundamental lag (when the fundamental lag << max_lag).

        # Find the first local minimum (after the self correlation peak)

        first_min = _max_lag-1
        for j in range(2, _max_lag-1):
            if _autocorr[j-1] < _autocorr[j-2] and \
               _autocorr[j-1] < _autocorr[j]:
               first_min = j - 1
               break

        # Find all the local maxima in the sample auto-correlation function

        acf_maxima = [ ]

        for j in range(first_min, _max_lag-1):
            # look for all local maximum of the sample auto-correlation function
            if _autocorr[j] > _autocorr[j-1] and \
               _autocorr[j] > _autocorr[j+1]:
                # append this lag to the list (the acf array index + 1)
                acf_maxima.append(j+1)

        if len(acf_maxima) > 0:
            # use information from the peak picker if it's available
            if hr_median is not None:
                # we favor later local peaks as those can give more precision
                for j in range(len(acf_maxima)):
                    # verify the peak lag is a multiple of the fundamental rate
                    ratio = acf_maxima[-j-1] / s_median
                    multiplier = round(ratio)
                    if multiplier > 0 and abs(ratio - multiplier) < 0.25:
                        j_peak = acf_maxima[-j-1]
                        hr_acf = 60.0 * multiplier / (j_peak * _tau_s)

            if hr_acf is not None:
                # This calculation of the scaled sample auto-correlation
                #   function is not 100% precise. We don't have the true
                #   sample mean of the lagged series, but the sample mean
                #   of the unlagged series is a decent approximation.
                acf_val = \
                    ((_autocorr[j_peak]/no_samples) - y_mean**2) / biased_var

                # Don't return a value when the correlation is below _ACF_THRESH
                if acf_val < _ACF_THRESH:
                    # Don't bother returning junk
                    hr_acf = None
                    acf_val = None

        # Update the display activity "spinner" index (modulo 4)

        self._spinner_n = (self._spinner_n + 1) % 4

        # Reset the buffer

        self.reset_buffer(False)

        # Return

        return [ hr_min, hr_max, hr_mean, hr_median, \
                         hr_acf, acf_val, d_qual, self._spinner_n]


    def reset_buffer(self, hard_reset):
        """
        # Reset the filtered sample history buffer
        """

        _max_lag = self._max_lag

        if hard_reset:

            # Perform a "hard" reset if requested
            self._y_buf_len = 0
            for j in range(_max_lag):
                self._autocorr[j] = 0.0
            self._wave_list.clear()

        else:

            # Retain old data so we can seamlessly compute auto-correlations
            if self._y_buf_len > _max_lag:
                # copy the _max_lag most recent samples
                self._y_buf[0:_max_lag] = \
                    self._y_buf[self._y_buf_len-_max_lag:self._y_buf_len]
                self._y_buf_len = _max_lag
            for j in range(_max_lag):
                self._autocorr[j] = 0.0
            self._wave_list.clear()


class PPG():
    """
    PPG processor class
    """

    def __init__(self, x_0_raw):
        self._dc_offset = int(x_0_raw)
        self.debug = None

        # Type II Tchebyshev HP filter, order 5, >40 dB suppression
        self._hpf = IIR_Filter(_B_HP, _A_HP)

        # Type II Tchebyshev LP filter, order 3, >20 dB suppression
        self._lpf = IIR_Filter(_B_LP, _A_LP)

        # Heart rate estimator with specified sample history buffer size
        self._hrm = Heart_Rate(_LEN_HRM_BUFFER_SEC)


    def preprocess(self, x_n_raw):
        """
        Preprocess a PPG sample.
        """

        if self.debug is not None:
            self.debug.append(x_n_raw)

        # Demean by the saved DC offset
        y_n = x_n_raw - self._dc_offset

        # Apply the digital filter(s)
        y_n = self._hpf.step(y_n)
        y_n = self._lpf.step(y_n)

        # Pass the filtered sample to the heart rate estimator
        self._hrm.step(x_n_raw, y_n)


    def dump_debug(self):
        """
        Dump the accumulated raw sample data buffer if debug mode is enabled
        """

        if self.debug is not None:
            if len(self.debug) >= _LEN_DEBUG_BUFFER:
                with open('hrs.data', 'ab') as f:
                    # write a re-sync marker
                    f.write(b'\xff\xff')
                    now = watch.rtc.get_localtime()
                    f.write(array.array('H', now[:6]))
                    f.write(self.debug)
                self.debug = array.array('H')


    def enable_debug(self):

        if self.debug is not None:
            self.debug = array.array('H')

