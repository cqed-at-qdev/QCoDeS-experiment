# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 13:31:14 2021
​
@author: Michaela
"""
​
#Start with initialization 
​
# SAMPLE_NAME = "M07-19-19.1-6Q2-D-v2-qubit"
SAMPLE_NAME = "M10-30-20.1-6Q5-SE"
EXPERIMENT_NAME = "2nd load"
FIRST_INITIALISATION = 0
​
​
import qcodes as qc
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import time
import yaml
import os
from functools import partial
from qcodes import Monitor
from qcodes import ManualParameter
from qdev_wrappers.file_setup import CURRENT_EXPERIMENT, my_init, close_station
from qdev_wrappers.station_configurator import StationConfigurator
from qcodes.dataset.plotting import plot_by_id
from qcodes.dataset.data_set import load_by_id
from qcodes.dataset.data_export import get_data_by_id
from qcodes import new_experiment, Parameter
from qcodes.dataset.measurements import Measurement
from qdev_wrappers.sweep_functions import do0d, do1d, do2d
from qcodes import ParamSpec, new_data_set
from qdev_wrappers.export import export_by_id
from local_instruments.qdev_fitter_QT1 import qdev_fitter
​
if FIRST_INITIALISATION:
    exp = new_experiment(EXPERIMENT_NAME, sample_name=SAMPLE_NAME)
​
scriptfolder = qc.config.user.scriptfolder
sys.path.append(scriptfolder)
​
mpl.rcParams['figure.subplot.bottom'] = 0.15
mpl.rcParams['font.size'] = 10
mpl.rcParams['image.cmap'] = 'hot'
​
if __name__ == '__main__':
​
    if qc.Station.default:
        close_station(qc.Station.default)
​
    STATION = qc.Station()
​
   
    my_init(SAMPLE_NAME, STATION,
        pdf_folder=True, png_folder=True, analysis_folder=True,
        waveforms_folder=True, calib_config=True,
        annotate_image=False, mainfolder=None, display_pdf=False,
        display_individual_pdf=False, qubit_count=1,
        plot_x_position=0.2)
    
    
# load instruments using station configuration
    scfg = StationConfigurator(station=STATION)
    
     
    vna = scfg.load_instrument('rs_vna_spec', timeout  = 600)#vna for  SPEC
     
    awg5208 = scfg.load_instrument('awg5208') #the AWG that generates our waveforms
     
    switch = scfg.load_instrument('minicircuits_switch')   #allows switching between VNA setup and time domain setup; below we only initiliaze the time domain setup
    
    #Rohde and Schwarz boxes that generate high frequency signals that can be modulated by the AWG    
    cavity = scfg.load_instrument('rs_cavity')  
    localos = scfg.load_instrument('rs_localos')  #second output of the cavity that is not modulated by the AWG
    qubit = scfg.load_instrument('rs_qubit')

​
    Monitor(*scfg.monitor_parameters.values())
​
#%% Then we load the Alazar and the pulse_builder
​
alazar = scfg.load_instrument('alazar')
alazar_ctrl = scfg.load_instrument('alazar_ctrl')
# alazar_ctrl.int_time(2e-7)
# alazar_ctrl.int_delay(2e-7)
​
sample_mag = scfg.load_instrument('sample_mag',parent=alazar_ctrl)
sample_phase = scfg.load_instrument('sample_phase',parent=alazar_ctrl)
sample_rec_mag = scfg.load_instrument('sample_rec_mag',parent=alazar_ctrl)
sample_rec_phase = scfg.load_instrument('sample_rec_phase',parent=alazar_ctrl)
​
avg_f1_mag = scfg.load_instrument('avg_f1_mag',parent=alazar_ctrl)
avg_f1_phase = scfg.load_instrument('avg_f1_phase',parent=alazar_ctrl)
​
avg_f2_mag = scfg.load_instrument('avg_f2_mag',parent=alazar_ctrl)
avg_f2_phase = scfg.load_instrument('avg_f2_phase',parent=alazar_ctrl)
​
avg_f3_mag = scfg.load_instrument('avg_f3_mag',parent=alazar_ctrl)
avg_f3_phase = scfg.load_instrument('avg_f3_phase',parent=alazar_ctrl)
​
avg_f4_mag = scfg.load_instrument('avg_f4_mag',parent=alazar_ctrl)
avg_f4_phase = scfg.load_instrument('avg_f4_phase',parent=alazar_ctrl)
​
rec_f1_mag = scfg.load_instrument('rec_f1_mag',parent=alazar_ctrl)  # we usually measure channel 12 and 13 which corresponds to I (rec_f1_mag) and (Qrec_f1_phase) ...we can also measure channels 2 and 3 but it's difficult to check what the channels are reading out specifically and if averaging/integration is on or off
rec_f1_phase = scfg.load_instrument('rec_f1_phase',parent=alazar_ctrl)
​
rec_f2_mag = scfg.load_instrument('rec_f2_mag',parent=alazar_ctrl)
rec_f2_phase = scfg.load_instrument('rec_f2_phase',parent=alazar_ctrl)
​
rec_f3_mag = scfg.load_instrument('rec_f3_mag',parent=alazar_ctrl)
rec_f3_phase = scfg.load_instrument('rec_f3_phase',parent=alazar_ctrl)
​
rec_f4_mag = scfg.load_instrument('rec_f4_mag',parent=alazar_ctrl)
rec_f4_phase = scfg.load_instrument('rec_f4_phase',parent=alazar_ctrl)
​
rec_f1_complex = scfg.load_instrument('rec_f1_complex',parent=alazar_ctrl)
alazar_ctrl.channels.extend([sample_mag, sample_phase, 
                             sample_rec_mag, sample_rec_phase, 
                             avg_f1_mag, avg_f1_phase, avg_f2_mag, avg_f2_phase,
                             avg_f3_mag, avg_f3_phase, avg_f4_mag, avg_f4_phase,
                             rec_f1_mag, rec_f1_phase, rec_f2_mag, rec_f2_phase,
                             rec_f3_mag, rec_f3_phase, rec_f4_mag, rec_f4_phase, rec_f1_complex])
alazar.sync_settings_to_card()
pulse_builder = scfg.load_instrument('pulse_builder',awg=awg5208,alazar=alazar,alazar_ctrl=alazar_ctrl,qubit=qubit,cavity=cavity)
​
​
#%% Extras
# define extra functions
Qdevfit = qdev_fitter()
​
#%%
vna.add_channel(channel_name='S21_{}'.format(0+1),vna_parameter='S21')
vna.channels.S21_1.format('Linear Magnitude')
cal_traces = [vna.S21_1.trace]
​
vna.channels.S21_1.format('Linear Magnitude')
vna.rf_on()
​
​
​
​
#%% Load Pulsebuilder setup   (feels like a blackbox to me and I am not sure what is used specifically during the measurments..I assume only the Alazar_SingleCavset function ?)
from functools import partial
import numpy as np
from scipy.optimize import curve_fit
​
def CavPrep(switch, qubit, cavity):
  qubit.status('off')
  cavity.status('off')
  switch.all(1)
​
# Skewed Lorentzian, see p. 161 of
# http://web.physics.ucsb.edu/~bmazin/Papers/2008/Gao/Caltech%20Thesis%202008%20Gao.pdf
def SkewLorentzian(f,f0,A1,A2,A3,A4,Q):
    return (A3+A4*(f-f0))/(1+4*Q**2*((f-f0)/f0)**2) + A1 + A2*(f-f0)
​
# Note "must be between 100000.0 and 500000000.0 inclusive; Parameter: alazar_ctrl_sample_rec_mag.demod_freq'"
​
def Alazar_Cavset(trace_list,pulsebuilder,switch,qubit,cavity,detuning=-0.1e6):
    for i, trace in enumerate(trace_list):
        mag_array = trace.get_latest()
        f = trace.setpoints[0]
        try:
            p = [f[np.argmin(mag_array)],mag_array.max(),0,-mag_array.max(),0,2000]
            # b = ((cavity.frequency() - detuning - 500000000, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf), (cavity.frequency() - detuning - 100000, np.inf, np.inf, np.inf, np.inf, np.inf))
            popt, _ = curve_fit(SkewLorentzian, f, mag_array, p0=p)
            f0 = popt[0]
        except:
            print('Fit not converged. Find Resonance by minimum instead.')
            f0 = f[np.argmin(mag_array)]
        getattr(pulsebuilder,'readout_freq_{}'.format(i+1))(f0 + detuning)
    pulsebuilder.update_readout_freqs()
    switch.all(2)
    qubit.status('on')
    cavity.status('on')
   
​
    
def Alazar_SingleCavset(trace,readout_freq_param,switch,qubit,cavity,detuning=+0.3e6):
    mag_array = trace.get_latest()
    f = trace.setpoints[0]
    try:
        popt, _ = curve_fit(SkewLorentzian, f, mag_array, p0=[f[np.argmin(mag_array)],mag_array.max(),0,-mag_array.max(),0,2000])
        f0 = popt[0]
    except:
        print('Fit not converged. Find Resonance by minimum instead.')
        f0 = f[np.argmin(mag_array)]
    demod_freq = readout_freq_param() - cavity.frequency()
    readout_freq_param(f0 + detuning)
    cavity.frequency(f0 + detuning - demod_freq)
    switch.a(1)
    switch.d(2)
    switch.c(2)
    switch.b(2)
    qubit.status('on')
    cavity.status('on')
​
​
# rec_f2_mag.data.setpoints = (tuple(pulse_builder.x_val()-qubit.frequency()+qf_new),)
def QubitShift(pulse_builder,qubit,freq,pwr):
    qubit.frequency(freq)
    qubit.power(pwr)
    #pulse_builder.averages(pulse_builder.averages())
​
def QubitPower(pulse_builder,qubit,pwr):
    qubit.power(pwr)    
    
def set_high_bias(qdac):
    qdac.Q3Bias(-2.538e-3)
    
def set_zero_bias(qdac):
    qdac.Q3Bias(-0.538e-3)
    
def cd_compensate(disp_gate, increment):
    disp_gate(disp_gate() + increment)
​
def pb(diff):
    pulse_builder.averages(int(pulse_builder.averages() + diff))
    print(pulse_builder.averages())
    
def QubitPower(pulse_builder,qubit,pwr):
    qubit.power(pwr)    
​
def compensate_cutter(cutter_gate, plunger_gate, slope, offset):
    current_plunger_voltage = plunger_gate()
    compensated_cutter_voltage = -slope * current_plunger_voltage + offset
    cutter_gate(compensated_cutter_voltage)
   
# def set_localos(localos, cavity, USE_DEMOD_FREQ):
#     localos.frequency(cavity.frequency()+USE_DEMOD_FREQ)
​
Alazar_Cavset_Task = qc.Task(partial(Alazar_Cavset,cal_traces,pulse_builder,switch,qubit,cavity,detuning=0.3e6))
Alazar_SingleCavset_Task = qc.Task(partial(Alazar_SingleCavset,vna.S21_1.trace,pulse_builder.readout_freq_1,switch,qubit,cavity,detuning=0.3e6))
CavPrep_Task = qc.Task(partial(CavPrep,switch,qubit,cavity))
# LocalosPrep_taks = qc.Task(partial(set_localos,localos, cavity, USE_DEMOD_FREQ))
​
​
​
​
​
​
#%%
#Example measurement: Non overlapping spectroscopy with Alazar and AWG; we assume that we did spectroscopy with the VNA beforehand and know the resonator frequency and a range for the qubit frequency
#In this measurement a qubit pulse is sent out first and directly afterwards a readout tone. If the qubit lifetime is sufficiently high we should see a peak in I and Q and knowing that the qubit didn't decay immediately back to the groundstate after the qubit drive
​
#prepare switches and instruments
​
switch.a(2)
switch.b(2)
​
#Turn Rohde and Schwarz boxes on and allow IQ modulation
cavity.status(1)  
cavity.IQ_state(1)
qubit.status(1)
qubit.IQ_state(1)
​
​
# set cavity and qubit power manually
USE_READOUT_FREQ = 6.6104e9
USE_DEMOD_FREQ = 20e6
USE_READOUT_POWER = -15
cavity.power(USE_READOUT_POWER)
qubit.power(-40)
​
cavity.frequency(USE_READOUT_FREQ - USE_DEMOD_FREQ) 
pulse_builder.readout_freq_1(USE_READOUT_FREQ)
​
#prepare channels 
alazar_ctrl.channels[12].demod_freq(USE_DEMOD_FREQ)
alazar_ctrl.channels[13].demod_freq(USE_DEMOD_FREQ)
alazar_ctrl.channels[20].demod_freq(USE_DEMOD_FREQ) # Complex channel
​
​
#set pulse builder settings
pulse_builder.int_time(0.6e-6); #integration time ...Is the limit really 400ns here ? It would be great if we could integrate for even shorter times
pulse_builder.marker_offset(0E-9)  #offset between readout pulse and marker that tells the Alazar to start reading out/integrating
pulse_builder.int_delay(0E-9)  #start integration some time after Alazar started reading out 
pulse_builder.readout_dur(2.5e-6);  #length of the resonator pulse (we refer to this pulse also as readout pulse)
pulse_builder.cycle_time(4e-6); #the time for one full pulse sequence 
pulse_builder.averages(80000); #averaging 
​
​
#update alazar card and choose measurement type..in this case MultiQubit SingleSideBand Spectroscopy Non Overlapping
pulse_builder.MultiQ_SSB_Spec_NoOverlap(4.5e9, 5.0e9, 101, pulse_length=0.6e-6)  #we define the frequency range in which we expect the qubit, here between 4.6GHz and 5GHz, the number of steps (101) and the pulse length of the qubit pulse 
alazar.sync_settings_to_card()
pulse_builder.update_readout_freqs()
alazar.seq_mode('on')
​
do0d(alazar_ctrl.channels[12:14].data)  #measures I and Q, which is channel 12 and 13 of the Alazar 
​
​
​
​
#%% Load and plot data
    
from qdev_wrappers.show_num import show_num, show_meta
from qdev_wrappers.file_setup import CURRENT_EXPERIMENT, my_init
# from qdev_wrappers.transmon.file_helpers import get_latest_counter
dummy_time = qc.ManualParameter('dummy_time')
​
def slow_import(mid):
    
    plot, data = show_num(mid); plt.close()
    
    xvals = data[0].traces[0]['config']['x'][:]
    yvals = data[0].traces[0]['config']['y'][:]
    
    return(np.array([xvals, yvals]))
​
​
#quick plotting
​
plot, data = show_num(3501); #insert number of measurement scan to load data
    
Q= data[0].traces[0]['config']['y'][:]   #sometimes Q and I are mixed up and we always have to double check if I is stored in data[0]... or data[1]...
I = data[1].traces[0]['config']['y'][:]
​
amplitude = np.sqrt(I**2+Q**2)
​
fq = data[0].traces[0]['config']['x'][:]
​
plt.clf()
plt.plot(fq, amplitude, ls = '-', color = 'firebrick')
plt.show()
