# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2020 Daniel Thompson
# Copyright (C) 2022 John C Peterson

"""HRS3300 driver
~~~~~~~~~~~~~~~~~

"""

from micropython import const

_I2CADDR = const(0x44)

_ID = const(0x00)
_ENABLE = const(0x01)
_ENABLE_HEN = const(0x80)
_ENABLE_HWT = const(0x70)
_ENABLE_PDRIVE1 = const(0x08)
_C1DATAM = const(0x08)
_C0DATAM = const(0x09)
_C0DATAH = const(0x0a)
_PDRIVER = const(0x0c)
_PDRIVER_PDRIVE0 = const(0x40)
_C1DATAH = const(0x0d)
_C1DATAL = const(0x0e)
_C0DATAL = const(0x0f)
_RES = const(0x16)
_HGAIN = const(0x17)

# Wait time between ADC conversions = 0.0 ms, LED drive durrent = 20mA (01b)
#_ENABLE_INIT = const(0x70)	# HRS disabled, PDRIVE[1] = 0
#_PDRIVER_INIT = const(0x6e)	# PDRIVE[0] = 1, low nibble = 0xe is used
#				# by at least two other HRS3300 drivers

# Wait time between ADC conversions = 0.0 ms, LED drive durrent = 30mA (10b)
#_ENABLE_INIT = const(0x78)	# HRS disabled, PDRIVE[1] = 1
#_PDRIVER_INIT = const(0x2e)	# PDRIVE[0] = 0, low nibble = 0xe is used
#				# by at least two other HRS3300 drivers

# Wait time between ADC conversions = 0.0 ms, LED drive durrent = 40mA (11b)
_ENABLE_INIT = const(0x78)	# HRS disabled, PDRIVE[1] = 1
_PDRIVER_INIT = const(0x6e)	# PDRIVE[0] = 1, low nibble = 0xe is used
				# by at least two other HRS3300 drivers

#_RES_INIT = const(0x66)	# HRS in 14-bit mode, and ALS in 14-bit mode
_RES_INIT = const(0x77)		# HRS in 15-bit mode, and ALS in 15-bit mode

_HGAIN_INIT = const(0x10)	# 64x ADC gain (the datasheet recommended value)
#_HGAIN_INIT = const(0x0d)	#  8x ADC gain (taken from reference C code)
#_HGAIN_INIT = const(0x08)	#  4x ADC gain
#_HGAIN_INIT = const(0x04)	#  2x ADC gain
#_HGAIN_INIT = const(0x00)	#  1x ADC gain


class HRS3300:
    def __init__(self, i2c):
        self._i2c = i2c

    def init(self):
        w = self.write_reg

        # ENABLE register (0x01)
        w(_ENABLE, _ENABLE_INIT)

        # HRS LED driver set register (0x0c)
        w(_PDRIVER, _PDRIVER_INIT)

        # RESOLUTION register (0x16)
        w(_RES, _RES_INIT)

        # HGAIN register (0x17)
        w(_HGAIN, _HGAIN_INIT)

    def read_reg(self, addr):
        return self._i2c.readfrom_mem(_I2CADDR, addr, 1)[0]

    def write_reg(self, addr, val):
        self._i2c.writeto_mem(_I2CADDR, addr, bytes((val,)))

    def enable(self):
        self.init()

        enable = self.read_reg(_ENABLE)
        enable |= _ENABLE_HEN
        self.write_reg(_ENABLE, enable)

    def disable(self):
        enable = self.read_reg(_ENABLE)
        enable &= ~_ENABLE_HEN
        self.write_reg(_ENABLE, enable)

    def read_hrs(self):
        # TODO: Try fusing the read of H & L
        m = self.read_reg(_C0DATAM)
        h = self.read_reg(_C0DATAH)
        l = self.read_reg(_C0DATAL)

        return (m << 8) | ((h & 0x0f) << 4) | (l & 0x0f) | ((l & 0x30) << 12)

    def read_als(self):
        # TODO: Try fusing the read of H & L
        m = self.read_reg(_C1DATAM)
        h = self.read_reg(_C1DATAH)
        l = self.read_reg(_C1DATAL)

        return (m << 3) | ((h & 0x3f) << 11) | (l & 0x07)

    def set_gain(self, gain):
        if gain > 64:
            gain = 64
        hgain = 0
        while (1 << hgain) < gain:
            hgain += 1
        self.write_reg(_HGAIN, hgain << 2)

    def set_drive(self, drive):
        """
        Set LED drive current

        Parameters:
            drive (int) LED drive current
                0 = 12.5 mA
                1 = 20   mA
                2 = 30   mA
                3 = 40   mA
        """
        en = self.read_reg(_ENABLE)
        pd = self.read_reg(_PDRIVER)
       
        en = (en & ~_ENABLE_PDRIVE1 ) | ((drive & 2) << 2)
        pd = (pd & ~_PDRIVER_PDRIVE0) | ((drive & 1) << 6)

        self.write_reg(_ENABLE, en)
        self.write_reg(_PDRIVER, pd)

    def set_hwt(self, t):
        """
        Set wait time between each conversion cycle

        Parameters:
            t (int) Wait time between each conversion cycle
                0 = 800   ms
                1 = 400   ms
                2 = 200   ms
                3 = 100   ms
                4 =  75   ms
                5 =  50   ms
                6 =  12.5 ms
                7 =   0   ms
        """
        en = self.read_reg(_ENABLE)
        en = (en & ~_ENABLE_HWT) | (t << 4)
        self.write_reg(_ENABLE, en)

