# Generate IIR filter coefficients for PPG filtering
# Copyright (C) 2022 John C Peterson

from scipy.signal import butter, bessel, cheby2, lfilter

import numpy as np

# PPG sample rate
#tau_s = 0.024
tau_s = 0.048
fs = 1.0 / tau_s
fn = fs / 2.0		# Nyquist frequency

# Heart rate pass band frequencies (BPM)
min_bpm = 30
max_bpm = 300		# include some of the higher harmonics of ~200 BPM

f_1 = (min_bpm/60.0) / fn
f_2 = (max_bpm/60.0) / fn

# Type II Tchebyshev HP filter, order 5
Rs = 40.0
#hp = cheby2(5, Rs, 0.43*f_1, 'highpass')	# tau_s = 24 ms
hp = cheby2(5, Rs, 0.55*f_1, 'highpass')	# tau_s = 48 ms

np.set_printoptions(precision=14)

print('hp = ', hp)

# Type II Tchebyshev LP filter, order 3
Rs = 20.0
#lp = cheby2(3, Rs, 1.80*f_2, 'lowpass');		# tau_s = 24 ms
lp = cheby2(3, Rs, 1.40*f_2, 'lowpass');		# tau_s = 48 ms

print('lp = ', lp)

