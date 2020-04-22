# Copyright (C) 2018 Alex Nitz
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""This module provides model classes that assume the noise is Gaussian.
"""

import numpy
import scipy.special

from pycbc import filter as pyfilter
from pycbc.waveform import get_fd_waveform, get_fd_waveform_from_td
from pycbc.detector import Detector

from .gaussian_noise import BaseGaussianNoise

# In this model we only calculate terms up to a constant.
# We are primarily interested in the posterior result


class SingleTemplate(BaseGaussianNoise):
    r"""Model that assumes we know all the intrinsic parameters.

    This model assumes we know all the intrinsic parameters, and are only
    maximizing over the extrinsic ones. We also assume a dominant mode waveform
    approximant only and non-precessing.


    Parameters
    ----------
    variable_params : (tuple of) string(s)
        A tuple of parameter names that will be varied.
    data : dict
        A dictionary of data, in which the keys are the detector names and the
        values are the data (assumed to be unwhitened). All data must have the
        same frequency resolution.
    low_frequency_cutoff : dict
        A dictionary of starting frequencies, in which the keys are the
        detector names and the values are the starting frequencies for the
        respective detectors to be used for computing inner products.
    sample_rate : int, optional
        The sample rate to use. Default is 32768.
    \**kwargs :
        All other keyword arguments are passed to
        :py:class:`BaseGaussianNoise`; see that class for details.
    """
    name = 'single_template'

    def __init__(self, variable_params, data, low_frequency_cutoff,
                 sample_rate=1, **kwargs):
        super(SingleTemplate, self).__init__(
            variable_params, data, low_frequency_cutoff, **kwargs)

        # Generate template waveforms
        df = data[self.detectors[0]].delta_f
        p = self.static_params.copy()
        if 'distance' in p:
            _ = p.pop('distance')
        if 'inclination' in p:
            _ = p.pop('inclination')
        try:
            hp, _ = get_fd_waveform(delta_f=df, distance=1, inclination=0, **p)
        except:
            hp, _ = get_fd_waveform_from_td(delta_f=df, distance=1, inclination=0, **p)

        # Extend template to high sample rate
        flen = int(int(sample_rate) / df) / 2 + 1
        hp.resize(flen)
        hp.save('hp_waveform.hdf')

        # Calculate high sample rate SNR time series
        self.sh = {}
        self.hh = {}
        self.det = {}
        for ifo in self.data:
            flow = self.kmin[ifo] * df
            fhigh = self.kmax[ifo] * df
            # Extend data to high sample rate
            self.data[ifo].resize(flen)
            self.det[ifo] = Detector(ifo)
            self.data[ifo].save('data_{}.hdf'.format(ifo))
            snr, _, norm = pyfilter.matched_filter_core(
                hp, self.data[ifo],
                psd=self.psds[ifo],
                low_frequency_cutoff=flow,
                high_frequency_cutoff=fhigh)

            self.sh[ifo] = 4 * df * snr
            self.hh[ifo] = -0.5 * pyfilter.sigmasq(
                hp, psd=self.psds[ifo],
                low_frequency_cutoff=flow,
                high_frequency_cutoff=fhigh)
                        
            self.sh[ifo].save('snr_{}.hdf'.format(ifo))
        self.time = None
        self.dref = pycbc.detector.Detector('Z1')

        #self.logging = open('logging.txt', 'w')

    def _loglr(self):
        r"""Computes the log likelihood ratio

        Returns
        -------
        float
            The value of the log likelihood ratio.
        """
        # calculate <d-h|d-h> = <h|h> - 2<h|d> + <d|d> up to a constant
        p = self.current_params.copy()
        p.update(self.static_params)

        if self.time is None:
            self.time = p['tc']

        shloglr = hhloglr = 0
        for ifo in self.sh:
            fp, fc = self.det[ifo].antenna_pattern(p['ra'], p['dec'],
                                                   p['polarization'],
                                                   self.time)
            dt = self.det[ifo].time_delay_from_detector(self.dref, p['ra'],
                                                        p['dec'],
                                                        self.time)
            ip = numpy.cos(p['inclination'])
            ic = 0.5 * (1.0 + ip * ip)
            htf = (fp * ip + 1.0j * fc * ic) / p['distance']

            sh = self.sh[ifo].at_time(p['tc'] + dt) * htf
            #print("DEBUGGING SNRs", ifo, int((p['tc'] + dt-self.sh[ifo].start_time)*self.sh[ifo].sample_rate), self.sh[ifo].at_time(p['tc'] + dt), htf, file=self.logging)
            #print("DEBUGGING SNRs", p['tc'] + dt, abs(self.sh[ifo].data).argmax(), self.sh[ifo].data[abs(self.sh[ifo].data).argmax()], self.time, dt, file=self.logging)
            shloglr += sh.real
            hhloglr += self.hh[ifo] * abs(htf)**2

        vloglr = shloglr + hhloglr
        #print ("LL EVAL", p['ra'], p['dec'], p['polarization'], p['tc'], p['inclination'], p['distance'], shloglr, hhloglr, file=self.logging)

        return float(vloglr)
