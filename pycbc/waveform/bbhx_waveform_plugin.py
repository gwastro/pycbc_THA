import numpy as np
from pycbc.types import FrequencySeries, Array
from pycbc import pnutils, conversions
from .waveform import props
from .parameters import location_params

from bbhx.waveformbuild import BBHWaveformFD

def BBHXWaveformFDInterface(run_phenomd=True, nyquist_freq=0.1,
                            sample_points=None, **params):
    # Some of this could go into waveform.py eventually.
    tmp_params = props(None, **params)
    params = location_params.default_dict().copy()
    params.update(tmp_params)


    # FIXME: I need to take an "epoch" and/or tc argument somehow??

    # Is it slow to do this every time?? Does it need caching??
    wave_gen = BBHWaveformFD(amp_phase_kwargs=dict(run_phenomd=run_phenomd)) 

    m1 = params['mass1']
    m2 = params['mass2']
    a1 = params['spin1z']
    a2 = params['spin2z']
    dist = pnutils.megaparsecs_to_meters(params['distance'])
    phi_ref = params['coa_phase']
    f_ref = 0 # This is now NOT standard LAL convention!
    inc = params['inclination'] # Convention here may not match. PLEASE CHECK!
    lam = params['EclipticLongitude'] # Convention here almost certainly does not match.
    beta = params['EclipticLatitude'] # Convention here almost certainly does not match.
    psi = params['polarization'] # Convention here may not match.
    t_ref = 0 # FIXME: This does need to be set somehow!!
    if sample_points is None:
        freqs = np.arange(0, nyquist_freq, params['delta_f'])
    else:
        freqs = sample_points
        # FIXME: Must not hardcode this here!!
        params['delta_f'] = 1E-8
    modes = [(2,2)] # More modes if not phenomd
    direct = False # See the BBHX documentation
    fill = True # See the BBHX documentation
    squeeze = True # See the BBHX documentation
    length = 1024 # An internal generation parameter, not an output parameter

    # FIXME: It would be good to generate the waveform from f_lower, but this
    #        seems difficult to then line up merger times, without using an
    #        expensive "roll" operation. Would like to resolve this, but for
    #        now we just generate from t_0 with merger at the end. This will
    #        result in some undesireable wraparound effects.
    shift_t_limits = False # Times are relative to merger
    #t_obs_start = pnutils.get_imr_duration(m1, m2, a1, a2, params['f_lower'],
    #                                       approximant='IMRPhenomD')
    #t_obs_start = conversions.sec_to_year(t_obs_start)
    t_obs_start = conversions.sec_to_year(1. / params['delta_f'])
    t_obs_end = 0.0 # Generates ringdown as well!

    wave = wave_gen(m1, m2, a1, a2,
                    dist, phi_ref, f_ref, inc, lam,
                    beta, psi, t_ref, freqs=freqs,
                    modes=modes, direct=direct, fill=fill, squeeze=squeeze,
                    length=length, t_obs_start=t_obs_start,
                    t_obs_end=t_obs_end,
                    shift_t_limits=shift_t_limits)[0]

    # Convert outputs to PyCBC arrays
    if sample_points is None:
        output = [FrequencySeries(wave[i], delta_f=params['delta_f'],
                                  epoch=params['tc']- 1/params['delta_f'])
                  for i in range(3)]
    else:
        output = [Array(wave[i]) for i in range(3)]
    return output

