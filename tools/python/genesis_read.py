
import os
import time
import logging
from copy import *
import re
from globals import m_e_eV, epsilon_0, speed_of_light, q_e, h_eV_s, lambda_C_r, I_Alfven, ro_e
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt

_logger = logging.getLogger(__name__) 


def calc_ph_sp_dens(spec, freq_ev, n_photons, spec_squared=1):
    """
    calculates number of photons per electronvolt
    """
    # _logger.debug('spec.shape = {}'.format(spec.shape))
    if spec.ndim == 1:
        axis = 0
    else:
        if spec.shape[0] == freq_ev.shape[0]:
            spec = spec.T
        axis = 1
        #     axis=0
        # elif spec.shape[1] == freq_ev.shape[0]:
        #     axis=1
        # else:
        #     raise ValueError('operands could not be broadcast together with shapes ', spec.shape, ' and ', freq_ev.shape)
    # _logger.debug('spec.shape = {}'.format(spec.shape))

    if spec_squared:
        spec_sum = np.trapz(spec, x=freq_ev, axis=axis)
    else:
        spec_sum = np.trapz(abs(spec) ** 2, x=freq_ev, axis=axis)

    if np.size(spec_sum) == 1:
        if spec_sum == 0:  # avoid division by zero
            spec_sum = np.inf
    else:
        spec_sum[spec_sum == 0] = np.inf  # avoid division by zero

    if spec_squared:
        norm_factor = n_photons / spec_sum
    else:
        norm_factor = np.sqrt(n_photons / spec_sum)

    if spec.ndim == 2:
        norm_factor = norm_factor[:, np.newaxis]
    # _logger.debug('spec.shape = {}'.format(spec.shape))
    # _logger.debug('norm_factor.shape = {}'.format(norm_factor.shape))
    spec = spec * norm_factor
    if axis == 1:
        spec = spec.T
    # _logger.debug('spec.shape = {}'.format(spec.shape))
    return spec
    
    
class GenesisOutput:
    '''
    Genesis output *.out files storage object
    '''

    def __init__(self):
        self.z = []
        self.s = []
        self.I = []
        self.n = []
        self.zSlice = []
        self.E = []
        self.aw = []
        self.qfld = []
        self.leng = []
        self.power_z_max = []
        self.power_z_mean = []
        self.power = []
        self.sliceKeys = []
        self.sliceValues = {}
        self.sliceKeys_used =[]
        self.parameters = {}
        self.filePath = ''

    def fileName(self):
        return filename_from_path(self.filePath)

    def __call__(self, name):

        if name not in self.parameters.keys():
            return None
        else:
            p, = self.parameters[name]
            return float(p.replace('D', 'E'))
            
    def calc_spec(self, mode='mid', npad=0):
        '''
        calculates the on-axis spectrum at every position along the undulator and writes it into "spec" attirube
        
        if mode = "mid" then on-axis power with on-axis phases is used for calculation
        if mode = "int" then transversely integrated power with on-axis phases is used (strictly speaking inaccurate, but informative)
        npad (integer) if > 0 pads the power with zeros in order to increase resolution of spectrum.
        '''
        if self.nSlices == 1:
            raise AssertionError('Cannot calculate spectrum from steady-state simulation')
        
        if (npad%1 != 0) or npad < 0:
            raise ValueError('npad should be positive integer')
        
        if mode == 'ff':
            try:
                power = self.far_field
            except AttributeError:
                mode = 'mid'
                
        if mode == 'mid':
            power = self.p_mid
        elif mode == 'int':
            power = self.p_int
        elif mode == 'ff':
            pass
        else:
            raise ValueError('mode should be either "mid" or "int"')
            
        _logger.debug('calculating spectrum')
        power = power / (2 * self.leng / self('ncar'))**2
        phi_mid = self.phi_mid
        
        zeros = np.zeros((self.nSlices * npad, self.nZ))
        power = np.vstack((power, zeros))
        phi_mid = np.vstack((phi_mid, zeros))
        
        spec = abs(np.fft.fft(np.sqrt(np.array(power)) * np.exp(1.j * np.array(phi_mid)), axis=0))**2 * self.dt**2 * 1e10
        e_0 = h_eV_s * speed_of_light / self('xlamds')
        freq_ev = h_eV_s * np.fft.fftfreq(len(spec), d=self('zsep') * self('xlamds') * self('ishsty') / speed_of_light) + e_0
        
        spec = np.fft.fftshift(spec, axes=0)
        freq_ev = np.fft.fftshift(freq_ev, axes=0)
        freq_lamd = h_eV_s * speed_of_light * 1e9 / freq_ev
        
        self.spec = spec
        self.freq_ev = freq_ev
        self.freq_lamd = freq_lamd
        self.spec_mode = mode
        self.sliceKeys_used.append('spec')
        
        sum_spec = np.sum(self.spec, axis=0)
        sum_spec[sum_spec == 0] = np.inf
        
        self.freq_ev_mean = np.sum(self.freq_ev[:,np.newaxis] * self.spec, axis=0) / sum_spec
        self.freq_ev_mean[self.freq_ev_mean == 0] = np.inf
        
        self.n_photons = self.pulse_energy / q_e / self.freq_ev_mean
        self.spec_phot_density = calc_ph_sp_dens(self.spec, self.freq_ev, self.n_photons)
        # self.spec_phot_density = self.spec #tmp
        self.sliceKeys_used.append('spec_phot_density')
        # print ('        done')
        
    def phase_fix(self, wav=None, s=None, **kwargs):
        '''
        the way to display the phase, without constant slope caused by different radiation wavelength from xlamds. phase is set to 0 at maximum power slice
        '''
        _logger.debug('rewrapping phase')
        
        if 'spec' not in self.sliceKeys_used:
            raise AssertionError('first spectrum should be calculated')
        
        self.phi_mid_disp = deepcopy(self.phi_mid)
        
        if 'phen' in kwargs:
            wav = (h_eV_s * speed_of_light) / kwargs['phen'] * 1e9
        if wav == None:
            spec_idx = np.argmax(self.spec[:, -1])
        else:
            spec_idx = find_nearest_idx(self.freq_lamd, wav)
        
        _logger.debug(ind_str + 'from {}m to {}m carrier'.format(self('xlamds'), wav))
        
        if s == None:
            pow_idx = np.argmax(self.p_mid[:, -1])
        else:
            pow_idx = find_nearest_idx(self.s,s)
            
        for zi in range(np.shape(self.phi_mid_disp)[1]):
            # if debug > 1:
                # print ('      fixing phase display: ' + str(zi) + ' of ' + str(range(shape(self.phi_mid_disp)[1])))

            maxspectrum_wavelength = self.freq_lamd[spec_idx] * 1e-9
            phase = np.unwrap(self.phi_mid[:, zi])
            phase_cor = np.arange(self.nSlices) * (maxspectrum_wavelength - self('xlamds')) / self('xlamds') * self('zsep') * 2 * pi
            phase_fixed = phase + phase_cor
            phase_fixed -= phase_fixed[pow_idx]
            n = 1
            phase_fixed = (phase_fixed + n * pi) % (2 * n * pi) - n * pi
            self.phi_mid_disp[:, zi] = phase_fixed
        self.sliceKeys_used.append('phi_mid_disp')
        # print ('        done')
        
    def calc_radsize(self, weigh_transv=1):
        '''
        weigh_transv = True to average the transverse radiation size over slices with radiation power as a weight
        '''
        if weigh_transv and self.nSlices != 1:
            _logger.debug('calculating the weighted transverse radiation size')
            if np.amax(self.power) > 0:
                weight = self.power + np.amin(self.power[self.power != 0]) / 1e6
            else:
                weight = np.ones_like(self.power)
            self.rad_t_size_weighted = np.average(self.r_size * 1e6, weights=weight, axis=0)
            self.sliceKeys_used.append('rad_t_size_weighted')
        # print ('        done')
        
    def wig(self,z=np.inf,*args,**kwargs):
        return wigner_out(self, z=z, method='mp', *args, **kwargs)
        
    def re_read(self, read_level=2):
        return read_out_file(self.filePath, read_level=read_level)


def read_out_file(filePath, read_level=2, precision=float,debug=1):
    speed_of_light = 299792458
    out = GenesisOutput()
    out.filePath = filePath

    ind_str=''
    precision=float
  #  debug=0
    read_level=2
    chunk = ''
    output_unsorted = []
    nSlice = 0

    wait_attempt = 3
    wait_time = 0.5
    while os.path.isfile(out.filePath) != True:
        _logger.warning(ind_str + 'waiting for "' + out.fileName() + '" ' + str(wait_time) + 's [' + str(wait_attempt) + ']')
        time.sleep(wait_time)  # wait for the .out file to be assembled
        wait_attempt -= 1
        if wait_attempt == 0:
            _logger.error(ind_str + 'file "' + out.filePath + '" not found')
            raise IOError('File ' + out.filePath + ' not found')

    if os.path.getsize(out.filePath) == 0:
        _logger.error(ind_str + 'file "' + out.filePath + '" has zero size')
        raise IOError('File ' + out.filePath + ' has zero size')
    
    start_time = time.time()
    f = open(out.filePath, 'r')

    null = f.readline()
    for line in f:
        tokens = line.strip().split()

        if len(tokens) < 1:
            continue

        if tokens[0] == '**********':
            if read_level == 0:
                break
            chunk = 'slices'
            nSlice = int(tokens[3])
            _logger.log(5, ind_str + 'reading slice # ' + str(nSlice))

        if tokens[0] == 'power':
            chunk = 'slice'
            if len(out.sliceKeys) == 0:  # to record the first instance
                out.sliceKeys = list(copy(tokens))
                _logger.debug(ind_str + 'reading slice values ')
            continue

        if tokens[0] == '$newrun':
            chunk = 'input1'
            _logger.debug(ind_str + 'reading input parameters')
            continue

        if tokens[0] == '$end':
            chunk = 'input2'
            continue

        if tokens == ['z[m]', 'aw', 'qfld']:
            chunk = 'magnetic optics'
            _logger.debug(ind_str + 'reading magnetic optics ')
            continue

        if chunk == 'magnetic optics':
            z, aw, qfld = list(map(precision, tokens))
            out.z.append(z)
            out.aw.append(aw)
            out.qfld.append(qfld)

        if chunk == 'input1':
            tokens = line.replace('=', '').strip().split()
            out.parameters[tokens[0]] = tokens[1:]
            #out.parameters[tokens[0]] = tokens[0:]
            # print 'input:', tokens
        if chunk == 'input2':
            tokens = line.replace('=', '').strip().split()
            out.parameters['_'.join(tokens[1:])] = [tokens[0]]
            #out.parameters[tokens[0]] = tokens[0:]
            # print 'input:', tokens
#
        if chunk == 'slice' and read_level >= 2:

            # tokens_fixed=re.sub(r'([0-9])\-([0-9])',r'\g<1>E-\g<2>',' '.join(tokens))
            # tokens_fixed=re.sub(r'([0-9])\+([0-9])',r'\g<1>E+\g<2>',tokens_fixed)
            # tokens=tokens_fixed.split()
            try:
                vals = list(map(precision, tokens))
            except ValueError:
                _logger.log(5, ind_str + 'wrong E value, fixing')
                _logger.log(5, ind_str + str(tokens))
                tokens_fixed = re.sub(r'([0-9])\-([0-9])', r'\g<1>E-\g<2>', ' '.join(tokens))
                tokens_fixed = re.sub(r'([0-9])\+([0-9])', r'\g<1>E+\g<2>', tokens_fixed)
                tokens_fixed = tokens_fixed.split()
                _logger.log(5, ind_str + str(tokens_fixed))
                vals = list(map(precision, tokens_fixed))

            output_unsorted.append(vals)

        if chunk == 'slices':
            if len(tokens) == 2 and tokens[1] == 'current':
                # print tokens[1]
                out.I.append(float(tokens[0]))
                out.n.append(nSlice)
            elif len(tokens) == 3 and tokens[1] == 'scan':
                out.I.append(float(tokens[0]))
                out.n.append(nSlice)

    #check for consistency
    if chunk == '':
        _logger.error(ind_str + 'File "' + out.filePath + '" has no genesis output information or is corrupted')
        raise ValueError('File "' + out.filePath + '" has no genesis output information or is corrupted')
    
    for parm in ['z', 'aw', 'qfld', 'I', 'n']:
        exec('out.' + parm + ' = np.array(out.' + parm + ')')

    if out('dgrid') == 0:
        rbeam = np.sqrt(out('rxbeam')**2 + out('rybeam')**2)
        ray = np.sqrt(out('zrayl') * out('xlamds') / np.pi * (1 + (out('zwaist') / out('zrayl'))**2))
        out.leng = out('rmax0') * (rbeam + ray)
    else:
        out.leng = 2 * out('dgrid')
    out.ncar = int(out('ncar'))  # number of mesh points
    
    #universal solution?
    out.leng=out('meshsize')*(out.ncar-1)
    
    
    if read_level == 0:
        _logger.debug(ind_str + 'read_level=0, returning *.out header')
        _logger.debug(ind_str + 'done in %.2f seconds' % (time.time() - start_time))
        return out

    out.nSlices = len(out.n)
    if out('entries_per_record') is None:
        _logger.error(ind_str + 'In file "' + out.filePath + '" file header is missing')
        raise ValueError('In file "' + out.filePath + '" file header is missing')
    
    out.nZ = int(out('entries_per_record'))  # number of records along the undulator
    _logger.debug(ind_str + 'nSlices ' + str(out.nSlices))
    _logger.debug(ind_str + 'nZ ' + str(out.nZ))

    if nSlice == 0:
        _logger.error(ind_str + 'In file "' + out.filePath + '" number of recorded slices is zero')
        raise ValueError('In file "' + out.filePath + '" number of recorded slices is zero')

    n_missing = (out.n[-1] - out.n[0]) - (len(out.n) - 1) * out('ishsty')
    if n_missing != 0:
        _logger.error(ind_str + 'File "' + out.filePath + '" is missing at least ' + str(n_missing) + ' slices')
        raise ValueError('File "' + out.filePath + '" is missing at least ' + str(n_missing) + ' slices')
    
    if read_level >= 2:
        output_unsorted = np.array(output_unsorted)  # .astype(precision)
        _logger.debug('output_unsorted.shape = ' + str(output_unsorted.shape))
        _logger.debug(ind_str + 'out.sliceKeys' + str(out.sliceKeys))
        for i in range(len(out.sliceKeys)):
            key = out.sliceKeys[int(i)]
            _logger.debug(ind_str + 'key = ' + str(key))
            if key[0].isdigit():
                hn = key[0]
                if 'bunch' in key:
                    key = 'bunching'
                elif 'phase' in key:
                    key = 'phi_mid'
                elif 'p-mid' in key:
                    key = 'p_mid'
                elif 'power' in key:
                    key = 'power'
                else:
                    pass
                key='h{:}_'.format(hn)+key
                _logger.debug(2*ind_str + 'key_new = ' + str(key))
                
            _logger.log(5, ind_str + 'assembling') 
            # _logger.debug('output_unsorted[:,' + str(i) + '].shape = ' + str(output_unsorted[:,i].shape))
            command = 'out.' + key.replace('-', '_').replace('<', '').replace('>', '') + ' = output_unsorted[:,' + str(i) + '].reshape((' + str(int(out.nSlices)) + ',' + str(int(out.nZ)) + '))'
            _logger.debug(ind_str + command)
            exec(command)
        if hasattr(out, 'energy'):
            out.energy += out('gamma0')
        out.power_z_max = np.max(out.power, 0)
        out.power_z_mean = np.mean(out.power, 0)
        out.sliceKeys_used = out.sliceKeys

    if out('itdp') == True:
        out.s = out('zsep') * out('xlamds') * (out.n - out.n[0])  # np.arange(0,out.nSlices)
        out.t = out.s / speed_of_light * 1.e+15
        out.dt = (out.t[1] - out.t[0]) * 1.e-15
        # out.dt=out('zsep') * out('xlamds') / speed_of_light
        out.beam_charge = np.sum(out.I * out.dt)
        out.sn_Imax = np.argmax(out.I)  # slice number with maximum current
        
        if read_level >= 2:
            out.pulse_energy = np.sum(out.power * out.dt, axis=0)
        
        if read_level >= 3:
            out.calc_spec(npad=2)
            # out.phase_fix()
            # out.calc_radsize(weigh_transv=1)
            
    else:
        out.s = [0]
    
    if out('iscan') != 0:
        out.scv = out.I  # scan value
        out.I = np.linspace(1, 1, len(out.scv))  # because used as a weight
    
    # tmp for back_compatibility
    if read_level >= 2:
        #out.power_int = out.power[:, -1]  # only output the exit power profile
        out.power_int = out.power          #output all positions power profile
        #out.max_power = np.amax(out.power_int)  # remove?
        for parm in [['power', 'p_int'],
                     ['energy', 'el_energy'],
                     ['e_spread', 'el_e_spread'],
                     ]:
            if hasattr(out, parm[0]):
                setattr(out, parm[1], getattr(out, parm[0]))
            for index, parm_key in enumerate(out.sliceKeys_used):
                if parm_key == parm[0]:
                    out.sliceKeys_used[index] = parm[1]
    #             delattr(out,parm[0])
        out.power = out.p_mid[:, -1]
        out.phi = out.phi_mid[:, -1]
    
    _logger.debug(ind_str + 'done in %.2f seconds' % (time.time() - start_time))
    return out


















