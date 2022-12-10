% Generate IIR filter coefficients for PPG filtering
% Copyright (C) 2022 John C Peterson

% PPG sample rate
%tau_s = 0.024;
tau_s = 0.048;
fs = 1.0 / tau_s;
fn = fs / 2.0;		% Nyquist frequency

% Heart rate pass band frequencies (BPM)
min_bpm = 30;
max_bpm = 300;		% include some of the higher harmonics of ~220 BPM, 5Hz

% Minimum attenuation in the stop bands (dB > 0)
%   More attenuation increases the width of the transition bands...

f_1 = (min_bpm/60.0) / fn;
f_2 = (max_bpm/60.0) / fn;
f_HF = 0.40 / fn;	% normal breathing freq

% Type II Tchebyshev HP filter, order 5, high pass filter
Rs = 40.0;
%[b_hp, a_hp] = cheby2(5, Rs, 0.43*f_1, 'high');	% tau_s = 24 ms
[b_hp, a_hp] = cheby2(5, Rs, 0.55*f_1, 'high');		% tau_s = 48 ms

clf;
fvtool(b_hp, a_hp, 'MagScale','Logarithmic'); hold on;
plot((40.0/60.0)/fn, 0, 'ro');
plot(f_1, 0, 'ko');
plot(f_HF, 0, 'go');
set(gca, 'FontSize', 14)
set(gca, 'FontWeight', 'bold')
pause()

printf('b_hp = %.14g %.14g %.14g %.14g %.14g %.14g\n', b_hp);
printf('a_hp = %.14g %.14g %.14g %.14g %.14g %.14g\n', a_hp);

% Type II Tchebyshev LP filter, order 3, low pass filter
Rs = 20.0;
%[b_lp, a_lp] = cheby2(3, Rs, 1.80*f_2);	% tau_s = 24 ms
[b_lp, a_lp] = cheby2(3, Rs, 1.40*f_2);		% tau_s = 48 ms

clf;
fvtool(b_lp, a_lp, 'MagScale','Logarithmic'); hold on;
plot((200.0/60.0)/fn, 0, 'ro');
plot(f_2, 0, 'ko');
set(gca, 'FontSize', 14)
set(gca, 'FontWeight', 'bold')
pause()

printf('b_lp = %.14g %.14g %.14g %.14g %.14g\n', b_lp);
printf('a_lp = %.14g %.14g %.14g %.14g %.14g\n', a_lp);

