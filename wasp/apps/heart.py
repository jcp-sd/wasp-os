# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2020 Daniel Thompson
# Copyright (C) 2022 John C Peterson

"""Heart rate monitor
~~~~~~~~~~~~~~~~~~~~~

A graphing heart rate monitor using a PPG sensor.

.. figure:: res/HeartApp.png
    :width: 179

This program also implements some (entirely optional) debug features to
store the raw PPG samples to the filesystem so that the samples can be used
to further refine the heart rate detection algorithm.

To enable the logging feature you MUST first select the heart rate
application and enable it using the watch UI! Then run the following
command via wasptool:

.. code-block:: sh

    ./tools/wasptool --eval 'wasp.system.app.debug = True'

Once debug mode has been enabled, the watch will automatically log
raw samples from the PPG sensor, (not to be confused with heart rate
estimates). This is only when the heart rate application is running.
Setting the debug flag to False will disable the logging when the
heart rate monitor next exits.

Finally to download the logs for analysis try:

.. code-block:: sh

    ./tools/wasptool --pull hrs.data
"""

import wasp
import machine
import ppg
import gc

# HRS in 14-bit mode, and ALS in 14-bit mode
#_POLLING_PERIOD_MS = 24.0
#_SAMPLES_PER_TICK = const(8)

# HRS in 15-bit mode, and ALS in 15-bit mode
_POLLING_PERIOD_MS = 48.0
_SAMPLES_PER_TICK = const(4)

COLOR565_R = const(0xf800)
COLOR565_Y = const(0xffe0)
COLOR565_G = const(0x07e0)

class HeartApp():
    """Heart rate monitor application."""
    NAME = 'Heart'

    def __init__(self):
        self._debug = False
        self._ppg = None
        self._timer = None
        self._t_timelast = None
        self._t_overflow = None
        self._t_nextsample = None
        self._no_dropped = None
        self._last_dropped = None
        self._last_hr = None
        self._draw_state = None
        self._polling_period_us = round(1000*_POLLING_PERIOD_MS)

    def foreground(self):
        """Activate the application."""

        # Enable the sensor
        wasp.watch.hrs.enable()

        # There is no delay after the hardware enable because the
        #  draw below should take long enough that no delay is needed
        draw = wasp.watch.drawable
        draw.fill()
        draw.set_color(wasp.system.theme('bright'))
        draw.string('Heart Rate Monitor', 0,  5, width=240)
        draw.string('Peak Picker',        0, 40, width=240)
        draw.string('Median:', 0,  70, width=100, right=True)
        draw.string('Mean:',   0,  95, width=100, right=True)
        draw.string('Max:',    0, 120, width=100, right=True)
        draw.string('Min:',    0, 145, width=100, right=True)
        draw.string('---.-', 100,  70, width=80, right=True)
        draw.string('bpm',   190,  70, width=50, right=True)
        draw.string('---.-', 100,  95, width=80, right=True)
        draw.string('bpm',   190,  95, width=50, right=True)
        draw.string('---.-', 100, 120, width=80, right=True)
        draw.string('bpm',   190, 120, width=50, right=True)
        draw.string('---.-', 100, 145, width=80, right=True)
        draw.string('bpm',   190, 145, width=50, right=True)

        draw.string('Sample ACF',   0, 180, width=240)
        draw.string('---.-',       35, 210, width=80, right=True)
        draw.string('bpm',        125, 210, width=50)

        # Start a hardware timer to accurately poll for the PPG samples
        self._timer = machine.Timer(id=1, period=0x7fffff)
        self._t_nextsample = self._polling_period_us
        self._t_overflow = 0
        self._t_timelast = 0
        self._no_dropped = 0
        self._last_dropped = None
        self._last_hr = None
        self._draw_state = 0
        self._timer.start()

        self._ppg = ppg.PPG(wasp.watch.hrs.read_hrs())

        # Get multiple samples per call, wake-up a little early to not be late
        wasp.system.request_tick(
            int((_SAMPLES_PER_TICK-0.25)*_POLLING_PERIOD_MS) )

        if self._debug:
            self._ppg.enable_debug()


    def background(self):
        wasp.watch.hrs.disable()
        self._timer.stop()
        self._timer = None
        self._ppg = None


    def _time32(self):
        """A timer that handles the wrap around of the hardware.Timer method"""

        # Get the current time
        time_now = self._timer.time()

        # Check for timer wrap around
        if time_now < self._t_timelast:
            self._t_overflow += 0x800000

        # Update the last returned time
        self._t_timelast = time_now

        return (self._t_overflow + time_now)


    def _subtick(self, ticks):
        """Notify the application that its periodic tick is due."""

        draw = wasp.watch.drawable

        # Read a new sample, feed it to the filter chain, heart rate estimator
        self._ppg.preprocess(wasp.watch.hrs.read_hrs())

        # Heart rate estimates (updates when the y buffer is full)
        hr = self._ppg._hrm.get_heart_rate()

        # Next action depends on the display drawing state. Calls to
        #   draw.string() are really slooo...ooooooow, so the task of
        #   updating the display is spread out over several subticks
        #   to (hopefully) avoid dropping PPG samples.
        if self._draw_state == 0:
            # watch for a non-empty heart rate estimate, save it if one is found
            if hr is not None:
                self._last_hr = hr
                # advance to the next display state
                self._draw_state += 1
        elif self._draw_state == 1:
            # update the median heart rate estimate (if meaningful)
            [hr_min, hr_max, hr_mean, hr_median, \
                             hr_acf, acf_val, d_qual, spinner_n] = self._last_hr
            txt = '{0:5.1f}'
            if hr_median is not None:
                draw.string(txt.format(hr_median),100,  70, width=80,right=True)
            # advance to the next display state
            self._draw_state += 1
        elif self._draw_state == 2:
            # update the mean heart rate estimate (if meaningful)
            [hr_min, hr_max, hr_mean, hr_median, \
                             hr_acf, acf_val, d_qual, spinner_n] = self._last_hr
            txt = '{0:5.1f}'
            if hr_mean is not None:
                draw.string(txt.format(hr_mean),  100,  95, width=80,right=True)
            # advance to the next display state
            self._draw_state += 1
        elif self._draw_state == 3:
            # update the maximum heart rate estimate (if meaningful)
            [hr_min, hr_max, hr_mean, hr_median, \
                             hr_acf, acf_val, d_qual, spinner_n] = self._last_hr
            txt = '{0:5.1f}'
            if hr_max is not None:
                draw.string(txt.format(hr_max),   100, 120, width=80,right=True)
            # advance to the next display state
            self._draw_state += 1
        elif self._draw_state == 4:
            # update the minimum heart rate estimate (if meaningful)
            [hr_min, hr_max, hr_mean, hr_median, \
                             hr_acf, acf_val, d_qual, spinner_n] = self._last_hr
            txt = '{0:5.1f}'
            if hr_min is not None:
                draw.string(txt.format(hr_min),   100, 145, width=80,right=True)
            # advance to the next display state
            self._draw_state += 1
        elif self._draw_state == 5:
            # update the sample ACF heart rate estimates (if meaningful)
            [hr_min, hr_max, hr_mean, hr_median, \
                             hr_acf, acf_val, d_qual, spinner_n] = self._last_hr
            if hr_acf is not None:
                txt = '{0:5.1f}'
                draw.string(txt.format(hr_acf),    35, 210, width=80,right=True)
            # advance to the next display state
            self._draw_state += 1
        elif self._draw_state == 6:
            [hr_min, hr_max, hr_mean, hr_median, \
                             hr_acf, acf_val, d_qual, spinner_n] = self._last_hr
            # update the spinner to show progress
            draw.string('|/-\\'[spinner_n], 0, 40)
            # update the heart rate estimate quality, status characters
            if hr_median is None:
                draw.set_color(COLOR565_R)
                draw.string('X', 220, 40, width=20, right=True)
            elif hr_max > 2.0*hr_min:
                draw.set_color(COLOR565_Y)
                draw.string('o', 220, 40, width=20, right=True)
            else:
                draw.set_color(COLOR565_G)
                draw.string('*', 220, 40, width=20, right=True)
            if hr_acf is None:
                draw.set_color(COLOR565_R)
                draw.string('X', 220, 180, width=20, right=True)
            elif acf_val < 0.75:
                draw.set_color(COLOR565_Y)
                draw.string('o', 220, 180, width=20, right=True)
            else:
                draw.set_color(COLOR565_G)
                draw.string('*', 220, 180, width=20, right=True)
            # set the color back to the default
            draw.set_color(wasp.system.theme('bright'))
            # advance to the next display state
            self._draw_state += 1
        elif self._draw_state == 7:
            # DEBUG show info on the watch display (TRASHES THE HEADINGS!!!)
            #[hr_min, hr_max, hr_mean, hr_median, \
            #                 hr_acf, acf_val, d_qual, spinner_n] = self._last_hr
            #no_peaks = d_qual[-1]
            #txt = '{}'.format(no_peaks)
            #txt = txt.replace('[','')
            #txt = txt.replace(']','')
            #draw.string(txt.replace(' ',''), 0, 5, width=240)
            #
            #txt = '{}'.format([round(d_qual[0],1),
            #                   round(d_qual[2],1),
            #                   self._no_dropped,
            #                   self._last_dropped])
            #txt = txt.replace('[','')
            #draw.string(txt.replace(' ',''), 20, 40)
            #
            # show just the number of dropped samples
            draw.string('{} '.format(self._no_dropped), 20, 40)
            # return to the initial state (watching for a heart rate estimate)
            self._draw_state = 0

        # Dump the debug raw data buffer when it gets full
        self._ppg.dump_debug()


    def tick(self, ticks):
        """
        This is an outrageous hack. At present, the RTC can not wake
        us up reliably at higher polling, sample rates. We use a hardware
        timer to ensure the polling is done as accurately as possible.
        """

        wasp.system.keep_awake()

        # Get multiple samples on this call
        for j in range(_SAMPLES_PER_TICK):
            while self._time32() < self._t_nextsample:
                pass
            # save the current time to check for dropped samples later
            time_now = self._time32()
            # sample the PPG ADC, preprocess
            self._subtick(1)
            # compute time to read the next sample
            self._t_nextsample += self._polling_period_us
            # check for dropped, missed samples
            if time_now > self._t_nextsample:
                # fell asleep at the wheel, samples were dropped!
                self._no_dropped += 1
                self._last_dropped = self._ppg._hrm._y_buf_len
                # adjust time to read the next sample to avoid repeated samples
                self._t_nextsample = time_now + self._polling_period_us


    @property
    def debug(self):
        return self._debug

    @debug.setter
    def debug(self, value):
        self._debug = value
        if value and self._ppg:
            self._ppg.enable_debug()

