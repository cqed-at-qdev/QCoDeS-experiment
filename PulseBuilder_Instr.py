import numpy as np
import math
from qcodes.instrument.base import Instrument
from qcodes import validators as vals


#%% Create funciton for two-qubit readout and spectroscopy

class MultiQ_PulseBuilder(Instrument):
    def __init__(self,name,number_read_freqs,alazar,alazar_ctrl,awg,qubit,cavity,**kwargs):
        super().__init__(name, **kwargs)
        self.awg = awg
        self.qubit = qubit
        self.cavity = cavity
        self.alazar = alazar
        self.alazar_ctrl = alazar_ctrl
        self.filename = 'Waveforms'
        self.x_val = lambda: np.linspace(0,1,10)
        self.SR = 2.5e9
        self.number_read_freqs = number_read_freqs
        
        self.add_parameter('cycle_time',
                      label='Pulse Cycle Time',
                      unit='s',
                      set_cmd= lambda x : x,
                      vals=vals.Numbers(0,10e-3))
        self.add_parameter('int_time',
                      label='Integration time',
                      unit='s',
                      set_cmd= lambda x : x,
                      vals=vals.Numbers(0,0.2e-3))
        self.add_parameter('int_delay',
                      label='Readout Delay',
                      unit='s',
                      set_cmd= lambda x : x,
                      vals=vals.Numbers(0,0.2e-3))
        self.add_parameter('readout_dur',
                      label='Readout Duration',
                      unit='s',
                      set_cmd= lambda x : x,
                      vals=vals.Numbers(0,0.2e-3))
        self.add_parameter('marker_offset',
                      label='Marker Offset',
                      unit='s',
                      set_cmd= lambda x : x,
                      vals=vals.Numbers(-1e-5,1e-5))
        self.add_parameter('averages',
                      label='Averages',
                      unit='',
                      set_cmd=self.num_averages,
                      vals=vals.Numbers(1,1e5))
        for i in range(number_read_freqs):
            self.add_parameter('readout_freq_{}'.format(i+1),
                        label='Readout Frequency {}'.format(i+1),
                        unit='Hz',
                        set_cmd= lambda x : x,
                        vals=vals.Numbers(0,12.5e9))

        
    def MultiQ_SSB_Spectroscopy(self, start, stop, npts, custom_name = None):
        qubitf = np.mean([start,stop])
        self.qubit.frequency(qubitf)
        self.x_val = lambda: np.linspace(start - qubitf,stop - qubitf,npts) + self.qubit.frequency()

        # Here we decide whether we want to use "default" names or a user-specified custom name
        if custom_name == None:
            readout_seq_name = 'Readout_Seq'
            this_seq_name = self.filename
            self.awg.clearSequenceList()
            self.awg.clearWaveformList()
        else:
            this_seq_name = 'Waveforms_' + custom_name
            readout_seq_name = 'Readout_Seq_' + custom_name  
        
        N = int((self.cycle_time()*self.SR+64) - self.cycle_time()*self.SR%64)
        N_offset = int(self.marker_offset()*self.SR)
        time = np.linspace(N/self.SR-self.readout_dur(), N/self.SR, int(self.readout_dur()*self.SR), endpoint=False)
        
        # Create triggers
        ZerosMarker = np.zeros(int(N))
        TriggerMarker = np.zeros(int(N))
        TriggerMarker[-int(self.readout_dur()*self.SR-N_offset):-int(self.readout_dur()*self.SR-N_offset-500e-9*self.SR)] = 1
        
        # Create Drive tones
        wfms = [[], []]
        for i , f in enumerate(self.x_val()-qubitf):
            # SSB drive tone
            #sine_signal = np.concatenate((0.5*np.sin(f*2*np.pi*time),np.zeros(N-len(time))))
            #cosine_signal = np.concatenate((0.5*np.cos(f*2*np.pi*time),np.zeros(N-len(time))))
            sine_signal = np.concatenate((np.zeros(N-len(time)), 0.5*np.sin(f*2*np.pi*time)))
            cosine_signal = np.concatenate((np.zeros(N-len(time)), 0.5*np.cos(f*2*np.pi*time)))
            if i == 0:
                wfm_ch1 = np.array([cosine_signal,TriggerMarker,TriggerMarker])
                wfm_ch2 = np.array([sine_signal,TriggerMarker,TriggerMarker])
            else:
                wfm_ch1 = np.array([cosine_signal,ZerosMarker,ZerosMarker])
                wfm_ch2 = np.array([sine_signal,ZerosMarker,ZerosMarker])
            wfms[0].append(wfm_ch1)
            wfms[1].append(wfm_ch2)

        trig_waits = [0 for _ in range(npts)] 
        nreps = [1 for _ in range(npts)] 
        event_jumps = [0 for _ in range(npts)]
        event_jump_to = [0 for _ in range(npts)]
        go_to = [0 for _ in range(npts)]
        go_to[-1] = 1 # Make the sequence loop back to first step
        
        seqx = self.awg.makeSEQXFile(trig_waits,
                                nreps,
                                event_jumps,
                                event_jump_to,
                                go_to,
                                wfms,
                                [1, 1],
                                this_seq_name, custom_name = custom_name)

        self.awg.sendSEQXFile(seqx, this_seq_name + '.seqx')
        self.awg.loadSEQXFile(this_seq_name + '.seqx')
        
        # Create Readout tones
        self.update_readout_freqs()

        # Create sequence for readout tones SLISt:SEQuence:NEW <sequence_name>,<number_of_steps> [,<number_of_tracks>]
        #self.awg.write('SLISt:SEQuence:NEW \"Readout_Seq\", {}, 2'.format(npts))
        #self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"Readout_Seq\", 1'.format(npts))
        self.awg.write('SLISt:SEQuence:NEW \"' + readout_seq_name + '\", {}, 2'.format(npts))
        self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"'.format(npts) + readout_seq_name + '\", 1')

        # Fill Readout waveforms into sequence
        for i in range(npts):
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet1:WAVeform \"'.format(str(i+1)) + readout_seq_name + '\", \"Readout_I\"')
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet2:WAVeform \"'.format(str(i+1)) + readout_seq_name + '\", \"Readout_Q\"')

        # Assign waveforms
        #self.awg.ch1.setSequenceTrack(self.filename, 1)
        #self.awg.ch2.setSequenceTrack(self.filename, 2)
        #self.awg.ch3.setSequenceTrack('Readout_Seq', 1)
        #self.awg.ch4.setSequenceTrack('Readout_Seq', 2)

        self.awg.ch1.setSequenceTrack(this_seq_name, 1)
        self.awg.ch2.setSequenceTrack(this_seq_name, 2)
        self.awg.ch3.setSequenceTrack(readout_seq_name, 1)
        self.awg.ch4.setSequenceTrack(readout_seq_name, 2)

        self.awg.ch2.state(1)
        self.awg.play()
        
        # Alazar labels
        for ala_chan in self.alazar_ctrl.channels[2:4]:
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('SSB Drive frequency (Overlap)', ala_chan.data.setpoint_labels[1])
            ala_chan.data.setpoint_units = ('Hz',ala_chan.data.setpoint_units[1])

        for n, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]): #AK: 
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('SSB Drive frequency (Overlap)',)
            ala_chan.data.setpoint_units = ('Hz',)

        # prepare channels
        self.num_averages(self._averages)

    def MultiQ_SSB_Spec_NoOverlap(self, start, stop, npts, pulse_length = 2e-6, custom_name = None):

        qubitf = np.mean([start,stop])
        self.qubit.frequency(qubitf)
        self.x_val = lambda: np.linspace(start - qubitf,stop - qubitf,npts) + self.qubit.frequency()

        # Here we decide whether we want to use "default" names or a user-specified custom name
        if custom_name == None:
            readout_seq_name = 'Readout_Seq'
            this_seq_name = self.filename
            self.awg.clearSequenceList()
            self.awg.clearWaveformList()
        else:
            this_seq_name = 'Waveforms_' + custom_name
            readout_seq_name = 'Readout_Seq_' + custom_name  
        
        N = int((self.cycle_time()*self.SR+64) - self.cycle_time()*self.SR%64)
        N_offset = int(self.marker_offset()*self.SR)
        time = np.linspace(0, pulse_length, int(pulse_length*self.SR), endpoint=False)
        
        # Create triggers
        ZerosMarker = np.zeros(int(N))
        TriggerMarker = np.zeros(int(N))
        TriggerMarker[-int(self.readout_dur()*self.SR-N_offset):-int(self.readout_dur()*self.SR-N_offset-500e-9*self.SR)] = 1
        
        # Create Drive tones
        wfms = [[], []]
        for i , f in enumerate(self.x_val()-qubitf):
            # SSB drive tone
            sine_signal = np.concatenate((np.zeros(N-len(time)-int(self.readout_dur()*self.SR)), 0.5*np.sin(f*2*np.pi*time),np.zeros(int(self.readout_dur()*self.SR))))
            cosine_signal = np.concatenate((np.zeros(N-len(time)-int(self.readout_dur()*self.SR)), 0.5*np.cos(f*2*np.pi*time),np.zeros(int(self.readout_dur()*self.SR))))
            if i == 0:
                wfm_ch1 = np.array([cosine_signal,TriggerMarker,TriggerMarker])
                wfm_ch2 = np.array([sine_signal,TriggerMarker,TriggerMarker])
            else:
                wfm_ch1 = np.array([cosine_signal,ZerosMarker,ZerosMarker])
                wfm_ch2 = np.array([sine_signal,ZerosMarker,ZerosMarker])
            wfms[0].append(wfm_ch1)
            wfms[1].append(wfm_ch2)

        trig_waits = [0 for _ in range(npts)] 
        nreps = [1 for _ in range(npts)] 
        event_jumps = [0 for _ in range(npts)]
        event_jump_to = [0 for _ in range(npts)]
        go_to = [0 for _ in range(npts)]
        go_to[-1] = 1 # Make the sequence loop back to first step
        
        seqx = self.awg.makeSEQXFile(trig_waits,
                                nreps,
                                event_jumps,
                                event_jump_to,
                                go_to,
                                wfms,
                                [1, 1],
                                this_seq_name, custom_name = custom_name)

        self.awg.sendSEQXFile(seqx, this_seq_name + '.seqx')
        self.awg.loadSEQXFile(this_seq_name + '.seqx')
        
        # Create Readout tones
        self.update_readout_freqs()

        # Create sequence for readout tones SLISt:SEQuence:NEW <sequence_name>,<number_of_steps> [,<number_of_tracks>]
        #self.awg.write('SLISt:SEQuence:NEW \"Readout_Seq\", {}, 2'.format(npts))
        #self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"Readout_Seq\", 1'.format(npts))
        self.awg.write('SLISt:SEQuence:NEW \"' + readout_seq_name + '\", {}, 2'.format(npts))
        self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"'.format(npts) + readout_seq_name + '\", 1')

        # Fill Readout waveforms into sequence
        for i in range(npts):
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet1:WAVeform \"'.format(str(i+1)) + readout_seq_name + '\", \"Readout_I\"')
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet2:WAVeform \"'.format(str(i+1)) + readout_seq_name + '\", \"Readout_Q\"')

        # Assign waveforms
        #self.awg.ch1.setSequenceTrack(self.filename, 1)
        #self.awg.ch2.setSequenceTrack(self.filename, 2)
        #self.awg.ch3.setSequenceTrack('Readout_Seq', 1)
        #self.awg.ch4.setSequenceTrack('Readout_Seq', 2)

        self.awg.ch1.setSequenceTrack(this_seq_name, 1)
        self.awg.ch2.setSequenceTrack(this_seq_name, 2)
        self.awg.ch3.setSequenceTrack(readout_seq_name, 1)
        self.awg.ch4.setSequenceTrack(readout_seq_name, 2)

        self.awg.ch2.state(1)
        self.awg.play()
        
        # Alazar labels
        for ala_chan in self.alazar_ctrl.channels[2:4]:
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('SSB Drive frequency (Non-overlap)',ala_chan.data.setpoint_labels[1])
            ala_chan.data.setpoint_units = ('Hz',ala_chan.data.setpoint_units[1])

        for n, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]):
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('SSB Drive frequency (Non-overlap)',)
            ala_chan.data.setpoint_units = ('Hz',)

        # prepare channels
        self.num_averages(self._averages)




    def MultiQ_SSB_Two_Seq_Pi(self, npts, carrier, freq1, dur1, freq2, dur2):

        # What is this definition above, taken from NoOverlap spectroscopy.
        # Why does it have to be a lambda function, and why is `qubitf` subtracted and then added back?
        # self.x_val = lambda: np.linspace(start - qubitf,stop - qubitf,npts) + self.qubit.frequency()

        # Set carrier frequency of `qubit` RF source, create `x_val` array, for looping through        
        self.qubit.frequency(carrier)
        self.x_val = lambda: np.arange(npts)

        # Clear AWG
        self.awg.clearSequenceList()
        self.awg.clearWaveformList()
        
        # Determine number of points in a waveform
        wf_npts = int((self.cycle_time()*self.SR+64) - self.cycle_time()*self.SR%64)
        N_offset = int(self.marker_offset() * self.SR) # Not sure what this one is. Not looking into it for now.
        
        # Create triggers
        ZerosMarker = np.zeros(int(wf_npts))
        TriggerMarker = np.zeros(int(wf_npts))
        TriggerMarker[-int(self.readout_dur()*self.SR - N_offset):-int(self.readout_dur()*self.SR-N_offset-500e-9*self.SR)] = 1
        
        # Create Drive tones
        wfms = [[], []]
        freq1_mod = freq1 - carrier
        freq2_mod = freq2 - carrier
        freq1_time_array = np.linspace(0, dur1, int(dur1 * self.SR), endpoint=False)
        freq2_time_array = np.linspace(0, dur2, int(dur2 * self.SR), endpoint=False)

        # for i , f in enumerate(self.x_val() - qubitf):
        # Below is from the poisoning code. We don't even need this. No need for loop at all.
        # for i in np.arange(2):

        wait_init = np.zeros(wf_npts - len(freq1_time_array) - len(freq2_time_array) - int(self.readout_dur() * self.SR))
        pipulse1_sin = 0.5 * np.sin(freq1_mod * 2 * np.pi * freq1_time_array)
        pipulse2_sin = 0.5 * np.sin(freq2_mod * 2 * np.pi * freq2_time_array)
        pipulse1_cos = 0.5 * np.cos(freq1_mod * 2 * np.pi * freq1_time_array)
        pipulse2_cos = 0.5 * np.cos(freq2_mod * 2 * np.pi * freq2_time_array)
        wait_readout = np.zeros(int(self.readout_dur() * self.SR))

        sine_signal = np.concatenate((wait_init, pipulse1_sin, pipulse2_sin, wait_readout))
        cosine_signal = np.concatenate((wait_init, pipulse1_cos, pipulse2_cos, wait_readout))

        wfm_ch1 = np.array([cosine_signal,TriggerMarker,TriggerMarker])
        wfm_ch2 = np.array([sine_signal,TriggerMarker,TriggerMarker])

        # Previously, we had double entries below, i.e.
        # wfms = [[wfm_ch1, wfm_ch1], [wfm_ch2, wfm_ch2]]
        # so that one of them can loop (npts-1) times. Now, however, that we are using `seq_mode = False` for these "historgram" measurements,
        # I think it should be OK to simply loop one element npts times.
        # Or actually, do we even need to loop?
        # We can just let the sequence repeat ... With `seq_mode = False`, the measurement stops because, IIUC, the `records_per_buffer` has been met ...
        
        wfms = [[wfm_ch1], [wfm_ch2]]

        trig_waits = [0] 
        nreps = [1] 
        event_jumps = [0]
        event_jump_to = [0]
        go_to = [1]

        # Setting the names for the AWG sequences
        # This follows a modification I made previously, to enable storing waveforms on the AWG.
        # But here, we stick with the standard sequence name, which is the value of `self.filename`
        this_seq_name = self.filename
        readout_seq_name = 'Readout_Seq'

        # The [1, 1] below is the `amplitudes` argument, of `makeSEQXFile`.
        # My impression is that we need [1, 1], not [1] because we are using both Ch1 and Ch2,
        # i.e., our Sequence has two tracks (and we play Track 1 on Ch1 and Track 2 on Ch2).
        seqx = self.awg.makeSEQXFile(trig_waits,
                                nreps,
                                event_jumps,
                                event_jump_to,
                                go_to,
                                wfms,
                                [1, 1],
                                this_seq_name)

        self.awg.sendSEQXFile(seqx, this_seq_name + '.seqx')
        self.awg.loadSEQXFile(this_seq_name + '.seqx')
        
        # Create Readout tones
        self.update_readout_freqs()

        # Create sequence for readout tones SLISt:SEQuence:NEW <sequence_name>,<number_of_steps> [,<number_of_tracks>]
        self.awg.write('SLISt:SEQuence:NEW \"' + readout_seq_name + '\", {}, 2'.format(1))

        # The below line is saying, that for our 1 waveforms long sequence, after Waveform 1, the instrument should go (back) to Waveform 1
        self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"'.format(1) + readout_seq_name + '\", 1')
        
        # Fill Readout waveforms into sequence
        # Here, we only put one (not 2, or npts many) Readout waveforms into the sequence
        # They are all the same, either way.
        self.awg.write('SLISt:SEQuence:STEP{}:TASSet1:WAVeform \"'.format(str(1)) + readout_seq_name + '\", \"Readout_I\"')
        self.awg.write('SLISt:SEQuence:STEP{}:TASSet2:WAVeform \"'.format(str(1)) + readout_seq_name + '\", \"Readout_Q\"')

        # Assign waveforms
        self.awg.ch1.setSequenceTrack(this_seq_name, 1)
        self.awg.ch2.setSequenceTrack(this_seq_name, 2)
        self.awg.ch3.setSequenceTrack(readout_seq_name, 1)
        self.awg.ch4.setSequenceTrack(readout_seq_name, 2)
        self.awg.ch2.state(1)
        self.awg.play()
        
        # Alazar labels
        for ala_chan in self.alazar_ctrl.channels[2:4]:
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Records',ala_chan.data.setpoint_labels[1])
            ala_chan.data.setpoint_units = ('',ala_chan.data.setpoint_units[1])

        for n, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]):
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Records',)
            ala_chan.data.setpoint_units = ('',)

        # Prepare channels
        self.num_averages(self._averages)


    def MultiQ_Rabi(self, start, stop, npts, custom_name = None):

        self.x_val = lambda: np.linspace(start ,stop ,npts)

        self.awg.ch2.state(0)

        # Here we decide whether we want to use "default" names or a user-specified custom name
        if custom_name == None:
            readout_seq_name = 'Readout_Seq'
            this_seq_name = self.filename
            self.awg.clearSequenceList()
            self.awg.clearWaveformList()
        else:
            this_seq_name = 'Waveforms_' + custom_name
            readout_seq_name = 'Readout_Seq_' + custom_name  
        
        N = int((self.cycle_time()*self.SR+64) - self.cycle_time()*self.SR%64)
        N_offset = int(self.marker_offset()*self.SR)
        
        # Create triggers
        ZerosMarker = np.zeros(int(N))
        TriggerMarker = np.zeros(int(N))
        TriggerMarker[-int(self.readout_dur()*self.SR-N_offset):-int(self.readout_dur()*self.SR-N_offset-500e-9*self.SR)] = 1
        
        # Create Drive tones
        wfms = [[]]
        for i , t in enumerate(self.x_val()):
            # SSB drive tone
            drive = np.zeros(int(N))
            drive[-int((self.readout_dur() + t)*self.SR):-int(self.readout_dur()*self.SR)] = 0.5
            if i == 0:
                wfm_ch1 = np.array([drive,TriggerMarker,TriggerMarker])
            else:
                wfm_ch1 = np.array([drive,ZerosMarker,ZerosMarker])
            wfms[0].append(wfm_ch1)
            
        trig_waits = [0 for _ in range(npts)] 
        nreps = [1 for _ in range(npts)] 
        event_jumps = [0 for _ in range(npts)]
        event_jump_to = [0 for _ in range(npts)]
        go_to = [0 for _ in range(npts)]
        go_to[-1] = 1 # Make the sequence loop back to first step
        
        seqx = self.awg.makeSEQXFile(trig_waits,
                                nreps,
                                event_jumps,
                                event_jump_to,
                                go_to,
                                wfms,
                                [1, 1],
                                this_seq_name, custom_name = custom_name)

        self.awg.sendSEQXFile(seqx, this_seq_name + '.seqx')
        self.awg.loadSEQXFile(this_seq_name + '.seqx')
        
        # Create Readout tones
        self.update_readout_freqs()

        # Create sequence for readout tones SLISt:SEQuence:NEW <sequence_name>,<number_of_steps> [,<number_of_tracks>]
        #self.awg.write('SLISt:SEQuence:NEW \"Readout_Seq\", {}, 2'.format(npts))
        #self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"Readout_Seq\", 1'.format(npts))
        self.awg.write('SLISt:SEQuence:NEW \"' + readout_seq_name + '\", {}, 2'.format(npts))
        self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"'.format(npts) + readout_seq_name + '\", 1')

        # Fill Readout waveforms into sequence
        for i in range(npts):
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet1:WAVeform \"'.format(str(i+1)) + readout_seq_name + '\", \"Readout_I\"')
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet2:WAVeform \"'.format(str(i+1)) + readout_seq_name + '\", \"Readout_Q\"')

        # Assign waveforms
        self.awg.ch1.setSequenceTrack(this_seq_name, 1)
        self.awg.ch3.setSequenceTrack(readout_seq_name, 1)
        self.awg.ch4.setSequenceTrack(readout_seq_name, 2)
        self.awg.play()
        
        # Alazar labels
        for ala_chan in self.alazar_ctrl.channels[2:4]:
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Drive time',ala_chan.data.setpoint_labels[1])
            ala_chan.data.setpoint_units = ('s',ala_chan.data.setpoint_units[1])

        for n, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]):
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Drive time',)
            ala_chan.data.setpoint_units = ('s',)

        # prepare channels
        self.num_averages(self._averages)

   


    def MultiQ_Lifetime(self, start, stop, npts, pipulse = 10e-9):

        self.x_val = lambda: np.linspace(start ,stop ,npts)

        # Clear AWG
        self.awg.ch2.state(0)
        self.awg.clearSequenceList()
        self.awg.clearWaveformList()
        
        N = int((self.cycle_time()*self.SR+64) - self.cycle_time()*self.SR%64)
        N_offset = int(self.marker_offset()*self.SR)
        
        # Create triggers
        ZerosMarker = np.zeros(int(N))
        TriggerMarker = np.zeros(int(N))
        TriggerMarker[-int(self.readout_dur()*self.SR-N_offset):-int(self.readout_dur()*self.SR-N_offset-500e-9*self.SR)] = 1
        
        # Create Drive tones
        wfms = [[]]
        for i , t in enumerate(self.x_val()):
            # SSB drive tone
            drive = np.zeros(int(N))
            drive[-int((self.readout_dur() + t + pipulse)*self.SR):-int((self.readout_dur() + t)*self.SR)] = 0.5
            if i == 0:
                wfm_ch1 = np.array([drive,TriggerMarker,TriggerMarker])
            else:
                wfm_ch1 = np.array([drive,ZerosMarker,ZerosMarker])
            wfms[0].append(wfm_ch1)
            
        trig_waits = [0 for _ in range(npts)] 
        nreps = [1 for _ in range(npts)] 
        event_jumps = [0 for _ in range(npts)]
        event_jump_to = [0 for _ in range(npts)]
        go_to = [0 for _ in range(npts)]
        go_to[-1] = 1 # Make the sequence loop back to first step

        seqx = self.awg.makeSEQXFile(trig_waits,
                                nreps,
                                event_jumps,
                                event_jump_to,
                                go_to,
                                wfms,
                                [1],
                                self.filename)
        self.awg.sendSEQXFile(seqx, self.filename + '.seqx')
        self.awg.loadSEQXFile(self.filename + '.seqx')
        
        # Create Readout tones
        self.update_readout_freqs()

        # Create sequence for readout tones SLISt:SEQuence:NEW <sequence_name>,<number_of_steps> [,<number_of_tracks>]
        self.awg.write('SLISt:SEQuence:NEW \"Readout_Seq\", {}, 2'.format(npts))
        self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"Readout_Seq\", 1'.format(npts))
        # Fill Readout waveforms into sequence
        for i in range(npts):
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet1:WAVeform \"Readout_Seq\", \"Readout_I\"'.format(str(i+1)))
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet2:WAVeform \"Readout_Seq\", \"Readout_Q\"'.format(str(i+1)))

        # Assign waveforms
        self.awg.ch1.setSequenceTrack(self.filename, 1)
        self.awg.ch3.setSequenceTrack('Readout_Seq', 1)
        self.awg.ch4.setSequenceTrack('Readout_Seq', 2)
        self.awg.play()
        
        # ALazar labels
        for ala_chan in self.alazar_ctrl.channels[2:4]:
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Wait time',ala_chan.data.setpoint_labels[1])
            ala_chan.data.setpoint_units = ('s',ala_chan.data.setpoint_units[1])

        for n, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]):
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Wait time',)
            ala_chan.data.setpoint_units = ('s',)

        # prepare channels
        self.num_averages(self._averages)

    def MultiQ_Ramsey(self, start, stop, npts, pi_half_pulse = 5e-9):

        self.x_val = lambda: np.linspace(start ,stop ,npts)

        # Clear AWG
        self.awg.ch2.state(0)
        self.awg.clearSequenceList()
        self.awg.clearWaveformList()
        
        N = int((self.cycle_time()*self.SR+64) - self.cycle_time()*self.SR%64)
        N_offset = int(self.marker_offset()*self.SR)
        
        # Create triggers
        ZerosMarker = np.zeros(int(N))
        TriggerMarker = np.zeros(int(N))
        TriggerMarker[-int(self.readout_dur()*self.SR-N_offset):-int(self.readout_dur()*self.SR-N_offset-500e-9*self.SR)] = 1
        
        # Create Drive tones
        wfms = [[]]
        for i , t in enumerate(self.x_val()):
            # SSB drive tone
            drive = np.zeros(int(N))
            drive[-int((self.readout_dur() + pi_half_pulse)*self.SR):-int(self.readout_dur()*self.SR)] = 0.5
            drive[-int((self.readout_dur() + t + 2*pi_half_pulse)*self.SR):-int((self.readout_dur() + t + pi_half_pulse)*self.SR)] = 0.5
            if i == 0:
                wfm_ch1 = np.array([drive,TriggerMarker,TriggerMarker])
            else:
                wfm_ch1 = np.array([drive,ZerosMarker,ZerosMarker])
            wfms[0].append(wfm_ch1)
            
        trig_waits = [0 for _ in range(npts)] 
        nreps = [1 for _ in range(npts)] 
        event_jumps = [0 for _ in range(npts)]
        event_jump_to = [0 for _ in range(npts)]
        go_to = [0 for _ in range(npts)]
        go_to[-1] = 1 # Make the sequence loop back to first step

        seqx = self.awg.makeSEQXFile(trig_waits,
                                nreps,
                                event_jumps,
                                event_jump_to,
                                go_to,
                                wfms,
                                [1],
                                self.filename)
        self.awg.sendSEQXFile(seqx, self.filename + '.seqx')
        self.awg.loadSEQXFile(self.filename + '.seqx')
        
        # Create Readout tones
        self.update_readout_freqs()

        # Create sequence for readout tones SLISt:SEQuence:NEW <sequence_name>,<number_of_steps> [,<number_of_tracks>]
        self.awg.write('SLISt:SEQuence:NEW \"Readout_Seq\", {}, 2'.format(npts))
        self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"Readout_Seq\", 1'.format(npts))
        # Fill Readout waveforms into sequence
        for i in range(npts):
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet1:WAVeform \"Readout_Seq\", \"Readout_I\"'.format(str(i+1)))
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet2:WAVeform \"Readout_Seq\", \"Readout_Q\"'.format(str(i+1)))

        # Assign waveforms
        self.awg.ch1.setSequenceTrack(self.filename, 1)
        self.awg.ch3.setSequenceTrack('Readout_Seq', 1)
        self.awg.ch4.setSequenceTrack('Readout_Seq', 2)
        self.awg.play()
        
        # ALazar labels
        for ala_chan in self.alazar_ctrl.channels[2:4]:
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Wait time',ala_chan.data.setpoint_labels[1])
            ala_chan.data.setpoint_units = ('s',ala_chan.data.setpoint_units[1])

        for n, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]):
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Wait time',)
            ala_chan.data.setpoint_units = ('s',)

        # prepare channels
        self.num_averages(self._averages)

    def MultiQ_Lambda(self, start, stop, npts, qubitf, delta, Delta, cont_drive_enabled = False, pidrive = 500e-9, wiggles = True, keep_old_seqs = False, histogram = False):

        if histogram == False:
            self.x_val = lambda: np.linspace(start, stop, npts)
        else:
            self.x_val = lambda: np.arange(npts)
            self.awg.ch2.state(0) # Do not need Q modulation for `histogram` mode

        self.qubit.frequency(qubitf)

        TAU_DRIVE_AMPLITUDE = 0.1
        PIPULSE_AMPLITUDE = 0.1

        if wiggles == True:
            pass
        else:
            TAU_DRIVE_AMPLITUDE = 0
            PIPULSE_AMPLITUDE = 0.5

        if keep_old_seqs == False:
            # Clear AWG
            self.awg.clearSequenceList()
            self.awg.clearWaveformList()
        
        N = int((self.cycle_time()*self.SR+64) - self.cycle_time()*self.SR%64)
        N_offset = int(self.marker_offset()*self.SR)
        
        # Create triggers
        ZerosMarker = np.zeros(int(N))
        TriggerMarker = np.zeros(int(N))
        TriggerMarker[-int(self.readout_dur()*self.SR-N_offset):-int(self.readout_dur()*self.SR-N_offset-500e-9*self.SR)] = 1
        
        # Create Drive tones
        wfms = [[], []]

        if cont_drive_enabled == False:
            if histogram == False:
                for i, t in enumerate(self.x_val()):
                    # SSB drive tone
                    sine_drive = np.zeros(int(N))
                    cosine_drive = np.zeros(int(N))
                    temp = np.linspace(0,t,int(self.SR*t))
                    l = len(temp)

                    cosine_drive[-l-int(self.readout_dur()*self.SR):-int(self.readout_dur()*self.SR)] = TAU_DRIVE_AMPLITUDE * np.cos(delta*2*np.pi*temp)+ TAU_DRIVE_AMPLITUDE * np.cos((Delta+delta)*2*np.pi*temp)
                    cosine_drive[-(int((self.readout_dur() + t+pidrive)*self.SR)):-int((self.readout_dur() + t)*self.SR)] = PIPULSE_AMPLITUDE

                    sine_drive[-l-int(self.readout_dur()*self.SR):-int(self.readout_dur()*self.SR)] = TAU_DRIVE_AMPLITUDE * np.sin(delta*2*np.pi*temp)+ TAU_DRIVE_AMPLITUDE * np.sin((Delta+delta)*2*np.pi*temp)
                    sine_drive[-(int((self.readout_dur() + t+pidrive)*self.SR)):-int((self.readout_dur() + t)*self.SR)] = 0
                    
                    if i == 0:
                        wfm_ch1 = np.array([cosine_drive, TriggerMarker, TriggerMarker])
                        wfm_ch2 = np.array([sine_drive, TriggerMarker, TriggerMarker])
                    else:
                        wfm_ch1 = np.array([cosine_drive, ZerosMarker, ZerosMarker])
                        wfm_ch2 = np.array([sine_drive, ZerosMarker, ZerosMarker])
                    wfms[0].append(wfm_ch1)
                    wfms[1].append(wfm_ch2)
            else:
                # This is the `histogram = True` case; The idea is simply to just readout after a \pi pulse.
                # So I am setting l = t = 0, and ignoring the second channel (we do not need Q modulation here).
                wfms = [[]]
                for i in np.arange(2):
                    # SSB drive tone
                    drive = np.zeros(int(N))
                    drive[-int(self.readout_dur()*self.SR):-int(self.readout_dur()*self.SR)] = 0
                    drive[-(int((self.readout_dur() + pidrive)*self.SR)):-int((self.readout_dur())*self.SR)] = PIPULSE_AMPLITUDE
                    if i == 0:
                        wfm_ch1 = np.array([drive,TriggerMarker,TriggerMarker])
                    else:
                        wfm_ch1 = np.array([drive,TriggerMarker,TriggerMarker])
                    wfms[0].append(wfm_ch1)                    
                
        else:
            for i, t in enumerate(self.x_val()):
                # SSB drive tone
                sine_drive = np.zeros(int(N))
                cosine_drive = np.zeros(int(N))
                temp = np.linspace(0,t,int(self.SR*t))
                l = len(temp)

                sine_drive[-l-int(self.readout_dur()*self.SR):-int(self.readout_dur()*self.SR)] = PIPULSE_AMPLITUDE + TAU_DRIVE_AMPLITUDE * np.sin(delta*2*np.pi*temp)+ TAU_DRIVE_AMPLITUDE * np.sin((Delta+delta)*2*np.pi*temp)
                cosine_drive[-l-int(self.readout_dur()*self.SR):-int(self.readout_dur()*self.SR)] = PIPULSE_AMPLITUDE + TAU_DRIVE_AMPLITUDE * np.cos(delta*2*np.pi*temp)+ TAU_DRIVE_AMPLITUDE * np.cos((Delta+delta)*2*np.pi*temp)
                
                if i == 0:
                    wfm_ch1 = np.array([cosine_drive, TriggerMarker, TriggerMarker])
                    wfm_ch2 = np.array([sine_drive, TriggerMarker, TriggerMarker])
                else:
                    wfm_ch1 = np.array([cosine_drive, ZerosMarker, ZerosMarker])
                    wfm_ch2 = np.array([sine_drive, ZerosMarker, ZerosMarker])
                wfms[0].append(wfm_ch1)
                wfms[1].append(wfm_ch2)            
            
        if histogram == False:
            trig_waits = [0 for _ in range(npts)] 
            nreps = [1 for _ in range(npts)] 
            event_jumps = [0 for _ in range(npts)]
            event_jump_to = [0 for _ in range(npts)]
            go_to = [0 for _ in range(npts)]
            go_to[-1] = 1 # Make the sequence loop back to first step
        else:
            trig_waits = [0, 0] 
            nreps = [1, npts - 1] 
            event_jumps = [0, 0]
            event_jump_to = [0, 0]
            go_to = [0, 1]

        this_seq_name = self.filename
        
        if histogram == False:
            amplitudes = [1, 1]
        else:
            amplitudes = [1]
        
        seqx = self.awg.makeSEQXFile(trig_waits,
                                nreps,
                                event_jumps,
                                event_jump_to,
                                go_to,
                                wfms,
                                amplitudes,
                                this_seq_name)
        self.awg.sendSEQXFile(seqx, this_seq_name + '.seqx')
        self.awg.loadSEQXFile(this_seq_name + '.seqx')
        
       # Create Readout tones
        self.update_readout_freqs()

        # Create sequence for readout tones SLISt:SEQuence:NEW <sequence_name>,<number_of_steps> [,<number_of_tracks>]
        # OE: Changing npts -> 2 below, not sure if correct
        self.awg.write('SLISt:SEQuence:NEW \"Readout_Seq\", {}, 2'.format(2))
        self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"Readout_Seq\", 1'.format(2))
        # Fill Readout waveforms into sequence
        for i in range(2):
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet1:WAVeform \"Readout_Seq\", \"Readout_I\"'.format(str(i+1)))
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet2:WAVeform \"Readout_Seq\", \"Readout_Q\"'.format(str(i+1)))

        # Assign waveforms
        self.awg.ch1.setSequenceTrack(this_seq_name, 1)
        if histogram == False:
            self.awg.ch2.setSequenceTrack(this_seq_name, 2)
            self.awg.ch2.state(1)
        self.awg.ch3.setSequenceTrack('Readout_Seq', 1)
        self.awg.ch4.setSequenceTrack('Readout_Seq', 2)
        self.awg.play()
        
        # Alazar labels
        if histogram == False:
            for ala_chan in self.alazar_ctrl.channels[2:4]:
                ala_chan.records_per_buffer(npts)
                ala_chan.data.setpoint_labels = ('Lambda Drive time', ala_chan.data.setpoint_labels[1])
                ala_chan.data.setpoint_units = ('s',ala_chan.data.setpoint_units[1])

            for n, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]):
                ala_chan.records_per_buffer(npts)
                ala_chan.data.setpoint_labels = ('Lambda Drive time',)
                ala_chan.data.setpoint_units = ('s',)
        else:
            for ala_chan in self.alazar_ctrl.channels[2:4]:
                ala_chan.records_per_buffer(npts)
                ala_chan.data.setpoint_labels = ('Records', ala_chan.data.setpoint_labels[1])
                ala_chan.data.setpoint_units = ('',ala_chan.data.setpoint_units[1])

            for n, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]):
                ala_chan.records_per_buffer(npts)
                ala_chan.data.setpoint_labels = ('Records',)
                ala_chan.data.setpoint_units = ('',)            

        # prepare channels
        self.num_averages(self._averages)    

    def MultiQ_Poisoning_Loop(self, npts, keep_old_seqs = False):

        self.x_val = lambda: np.arange(npts)

        self.awg.ch2.state(0)

        # Clear AWG
        if keep_old_seqs == False:
            self.awg.clearSequenceList()
            self.awg.clearWaveformList()
        
        N = int((self.cycle_time()*self.SR+64) - self.cycle_time()*self.SR%64)
        N_offset = int(self.marker_offset()*self.SR)
        
        # Create triggers
        ZerosMarker = np.zeros(int(N))
        TriggerMarker = np.zeros(int(N))
        TriggerMarker[-int(self.readout_dur()*self.SR-N_offset):-int(self.readout_dur()*self.SR-N_offset-500e-9*self.SR)] = 1
        
        # Create Drive tones
        wfms = [[]]
        for i in np.arange(2):
            # SSB drive tone
            drive = np.zeros(int(N))
            if i == 0:
                wfm_ch1 = np.array([drive,TriggerMarker,TriggerMarker])
            else:
                wfm_ch1 = np.array([drive,TriggerMarker,TriggerMarker])
            wfms[0].append(wfm_ch1)
            
        trig_waits = [0, 0] 
        nreps = [1, npts - 1] 
        event_jumps = [0, 0]
        event_jump_to = [0, 0]
        go_to = [0, 1]

        seqx = self.awg.makeSEQXFile(trig_waits,
                                nreps,
                                event_jumps,
                                event_jump_to,
                                go_to,
                                wfms,
                                [1],
                                self.filename)
        self.awg.sendSEQXFile(seqx, self.filename + '.seqx')
        self.awg.loadSEQXFile(self.filename + '.seqx')
        
        # Create Readout tones
        self.update_readout_freqs()

        # Create sequence for readout tones SLISt:SEQuence:NEW <sequence_name>,<number_of_steps> [,<number_of_tracks>]
        # OE: Changing npts -> 2 below, not sure if correct
        # --> This means that he readout waveforms loop within a pair, we do not include `npts` number of waveforms in the Readout Seq.
        #     Since they all look the same, this should be fine. In fact, they might not need to be uploaded in the first place.
        self.awg.write('SLISt:SEQuence:NEW \"Readout_Seq\", {}, 2'.format(2))
        self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"Readout_Seq\", 1'.format(2))
        # Fill Readout waveforms into sequence
        for i in range(2):
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet1:WAVeform \"Readout_Seq\", \"Readout_I\"'.format(str(i+1)))
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet2:WAVeform \"Readout_Seq\", \"Readout_Q\"'.format(str(i+1)))

        # Assign waveforms
        self.awg.ch1.setSequenceTrack(self.filename, 1)
        self.awg.ch3.setSequenceTrack('Readout_Seq', 1)
        self.awg.ch4.setSequenceTrack('Readout_Seq', 2)
        self.awg.play()
        
        # ALazar labels
        for ala_chan in self.alazar_ctrl.channels[2:4]:
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Records',ala_chan.data.setpoint_labels[1])
            ala_chan.data.setpoint_units = ('',ala_chan.data.setpoint_units[1])

        for n, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]):
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Records',)
            ala_chan.data.setpoint_units = ('',)

        # prepare channels
        self.num_averages(self._averages)

    def set_readout_freqs(self,readout_frequencies):
        if len(readout_frequencies) != self.number_read_freqs:
            raise ValueError('Number of given readout frequencies has to \
                                be {}.'.format(self.number_read_freqs))
        for i in range(self.number_read_freqs):
            getattr(self,'readout_freq_{}'.format(i+1))(readout_frequencies[i])
        self.update_readout_freqs()

    def get_readout_freqs(self):
        ret = []
        for i in range(self.number_read_freqs):
            ret.append(getattr(self,'readout_freq_{}'.format(i+1))())
        return np.array(ret)

    def update_readout_freqs(self):
        readout_freqs = self.get_readout_freqs()
        N = int((self.cycle_time()*self.SR+64) - self.cycle_time()*self.SR%64)

        N_offset = int(self.marker_offset()*self.SR)
        # Create trigger
        TriggerMarker = np.zeros(N)
        TriggerMarker[-int(self.readout_dur()*self.SR-N_offset):-int(self.readout_dur()*self.SR-N_offset-500e-9*self.SR)] = 1
    
        # Readout tones
        time = np.linspace(0, self.readout_dur(), int(self.readout_dur()*self.SR), endpoint=False)
        cosine_readout = np.zeros(N)
        sine_readout = np.zeros(N)
        for fr in readout_freqs:
            fi = fr - self.cavity.frequency()
            cosine_readout += (0.5/len(readout_freqs))*np.concatenate((np.zeros(N-len(time)),np.cos(fi*2*np.pi*time)))
            sine_readout += (0.5/len(readout_freqs))*np.concatenate((np.zeros(N-len(time)),np.sin(fi*2*np.pi*time)))
        wfm_ch3 = np.array([cosine_readout,TriggerMarker,TriggerMarker])
        wfm_ch4 = np.array([sine_readout,TriggerMarker,TriggerMarker])
        
        state = self.awg.run_state()
        self.awg.stop()
        wfm_ch3_file = self.awg.makeWFMXFile(wfm_ch3, 1)
        wfm_ch4_file = self.awg.makeWFMXFile(wfm_ch4, 1)
        self.awg.sendWFMXFile(wfm_ch3_file, 'Readout_I.wfmx')
        self.awg.sendWFMXFile(wfm_ch4_file, 'Readout_Q.wfmx')
        self.awg.loadWFMXFile('Readout_I.wfmx')
        self.awg.loadWFMXFile('Readout_Q.wfmx')
        # Only start play if running to begin with
        if state == 'Running':
            self.awg.play()
        
        # Set demod frequencies
        for ala_chan in self.alazar_ctrl.channels[2:4]:
            ala_chan.demod_freq(abs(readout_freqs[0]-self.cavity.frequency()))
            
        for n, ala_chan in enumerate(self.alazar_ctrl.channels[4:12]):
            try:
                ala_chan.demod_freq(abs(readout_freqs[n//2]-self.cavity.frequency()))
            except:
                pass
        for n, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]):#changed to 12:24 from 12:20 to add rec5 and rec6, AK 8/10
            try:
                ala_chan.demod_freq(abs(readout_freqs[n//2]-self.cavity.frequency()))
            except:
                pass

    def num_averages(self,value):
        self._averages = value
        self.alazar_ctrl.int_time(self.int_time())
        self.alazar_ctrl.int_delay(self.int_delay())
        for ala_chan in self.alazar_ctrl.channels[2:4]:
            ala_chan.num_averages(value)
            ala_chan.prepare_channel()
            ala_chan.data.setpoints = (tuple(self.x_val()),ala_chan.data.setpoints[1])
            
        for ala_chan in self.alazar_ctrl.channels[4:12]:
            ala_chan.num_averages(value)
            ala_chan.prepare_channel()
            
        for i, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]): #Changed to 24, AK 8/10
            if i != 8:
                ala_chan.num_averages(value)
            ala_chan.prepare_channel()
            ala_chan.data.setpoints = (tuple(self.x_val()),)


    #define new functions for to have a duty cycle measurement and partly overlapping spec and rabis. Albert Hertel fall 2020 to January 2021
    def MultiQ_Duty_cycle(self, tmin, tmax, npts):

        # Define array with different drive times
        self.x_val = lambda: np.geomspace(tmin ,tmax , num = npts)
        
        # I do not understand why, but we have to add this to the qubit drive to sync drive and readout
        offset = self.cycle_time()-self.readout_dur()

        # Calculate how many duty cycles on can fit into readout_dur
        M = np.zeros(len(self.x_val()), dtype = int)
        for k, l in enumerate(self.x_val()): 
            M[k] = math.floor(self.readout_dur()/l)

        # Clear AWG
        self.awg.ch2.state(0)
        self.awg.clearSequenceList()
        self.awg.clearWaveformList()
        
        N = int((self.cycle_time()*self.SR+64) - self.cycle_time()*self.SR%64)
        N_offset = int(self.marker_offset()*self.SR)
        
        # Create triggers
        ZerosMarker = np.zeros(int(N))
        TriggerMarker = np.zeros(int(N))
        TriggerMarker[-int(self.readout_dur()*self.SR-N_offset):-int(self.readout_dur()*self.SR-N_offset-500e-9*self.SR)] = 1
        
        # Create Drive tones
        wfms = [[]]
        for i , t in enumerate(self.x_val()):
            # SSB drive tone
            drive = np.zeros(int(N))
            for j in range(M[i]):
                    drive[int((j*t+offset)*self.SR):int((j*t+t/2+offset)*self.SR)] = 0.5    
            # print(np.where(drive!=0))
            if i == 0:
                wfm_ch1 = np.array([drive,TriggerMarker,TriggerMarker])
            else:
                wfm_ch1 = np.array([drive,ZerosMarker,ZerosMarker])
            wfms[0].append(wfm_ch1)
            
        trig_waits = [0 for _ in range(npts)] 
        nreps = [1 for _ in range(npts)] 
        event_jumps = [0 for _ in range(npts)]
        event_jump_to = [0 for _ in range(npts)]
        go_to = [0 for _ in range(npts)]
        go_to[-1] = 1 # Make the sequence loop back to first step

        seqx = self.awg.makeSEQXFile(trig_waits,
                                nreps,
                                event_jumps,
                                event_jump_to,
                                go_to,
                                wfms,
                                [1],
                                self.filename)
        self.awg.sendSEQXFile(seqx, self.filename + '.seqx')
        self.awg.loadSEQXFile(self.filename + '.seqx')
        
        # Create Readout tones
        self.update_readout_freqs()

        # Create sequence for readout tones SLISt:SEQuence:NEW <sequence_name>,<number_of_steps> [,<number_of_tracks>]
        self.awg.write('SLISt:SEQuence:NEW \"Readout_Seq\", {}, 2'.format(npts))
        self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"Readout_Seq\", 1'.format(npts))
        # Fill Readout waveforms into sequence
        for i in range(npts):
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet1:WAVeform \"Readout_Seq\", \"Readout_I\"'.format(str(i+1)))
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet2:WAVeform \"Readout_Seq\", \"Readout_Q\"'.format(str(i+1)))

        # Assign waveforms
        self.awg.ch1.setSequenceTrack(self.filename, 1)
        self.awg.ch3.setSequenceTrack('Readout_Seq', 1)
        self.awg.ch4.setSequenceTrack('Readout_Seq', 2)
        self.awg.play()
        
        # ALazar labels
        for ala_chan in self.alazar_ctrl.channels[2:4]:
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Cycle time tau',ala_chan.data.setpoint_labels[1])
            ala_chan.data.setpoint_units = ('s',ala_chan.data.setpoint_units[1])

        for n, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]):
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Cycle time tau',)
            ala_chan.data.setpoint_units = ('s',)

        # prepare channels
        self.num_averages(self._averages)            
        
        
        
    def MultiQ_SSB_Spec_SomeOverlap(self, start, stop, npts, pulse_length = 2e-6, overlap = 1e-6,custom_name = None):

        qubitf = np.mean([start,stop])
        self.qubit.frequency(qubitf)
        self.x_val = lambda: np.linspace(start - qubitf,stop - qubitf,npts) + self.qubit.frequency()

        # Here we decide whether we want to use "default" names or a user-specified custom name
        if custom_name == None:
            readout_seq_name = 'Readout_Seq'
            this_seq_name = self.filename
            self.awg.clearSequenceList()
            self.awg.clearWaveformList()
        else:
            this_seq_name = 'Waveforms_' + custom_name
            readout_seq_name = 'Readout_Seq_' + custom_name  
        
        N = int((self.cycle_time()*self.SR+64) - self.cycle_time()*self.SR%64)
        N_offset = int(self.marker_offset()*self.SR)
        time = np.linspace(0, pulse_length, int(pulse_length*self.SR), endpoint=False)
        #overlap to be added to drive tone. Keep the trigger after the drive tone
        overlap_offset = int(overlap*self.SR)
        overlap_array = np.linspace(0, overlap, int(overlap*self.SR), endpoint=False)

        # Create triggers
        ZerosMarker = np.zeros(int(N))
        TriggerMarker = np.zeros(int(N))
        TriggerMarker[-int(self.readout_dur()*self.SR-N_offset):-int(self.readout_dur()*self.SR-N_offset-500e-9*self.SR)] = 1

        
        # Create Drive tones
        wfms = [[], []]
        for i , f in enumerate(self.x_val()-qubitf):
            # SSB drive tone
            sine_signal = np.concatenate((np.zeros(N-len(time)+len(overlap_array)-int(self.readout_dur()*self.SR)), 0.5*np.sin(f*2*np.pi*time),np.zeros(int(self.readout_dur()*self.SR-len(overlap_array)))))
            cosine_signal = np.concatenate((np.zeros(N-len(time)+len(overlap_array)-int(self.readout_dur()*self.SR)), 0.5*np.cos(f*2*np.pi*time),np.zeros(int(self.readout_dur()*self.SR-len(overlap_array)))))
            if i == 0:
                wfm_ch1 = np.array([cosine_signal,TriggerMarker,TriggerMarker])
                wfm_ch2 = np.array([sine_signal,TriggerMarker,TriggerMarker])
            else:
                wfm_ch1 = np.array([cosine_signal,ZerosMarker,ZerosMarker])
                wfm_ch2 = np.array([sine_signal,ZerosMarker,ZerosMarker])
            wfms[0].append(wfm_ch1)
            wfms[1].append(wfm_ch2)

        trig_waits = [0 for _ in range(npts)] 
        nreps = [1 for _ in range(npts)] 
        event_jumps = [0 for _ in range(npts)]
        event_jump_to = [0 for _ in range(npts)]
        go_to = [0 for _ in range(npts)]
        go_to[-1] = 1 # Make the sequence loop back to first step
        
        seqx = self.awg.makeSEQXFile(trig_waits,
                                nreps,
                                event_jumps,
                                event_jump_to,
                                go_to,
                                wfms,
                                [1, 1],
                                this_seq_name, custom_name = custom_name)

        self.awg.sendSEQXFile(seqx, this_seq_name + '.seqx')
        self.awg.loadSEQXFile(this_seq_name + '.seqx')
        
        # Create Readout tones
        self.update_readout_freqs()

        # Create sequence for readout tones SLISt:SEQuence:NEW <sequence_name>,<number_of_steps> [,<number_of_tracks>]
        #self.awg.write('SLISt:SEQuence:NEW \"Readout_Seq\", {}, 2'.format(npts))
        #self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"Readout_Seq\", 1'.format(npts))
        self.awg.write('SLISt:SEQuence:NEW \"' + readout_seq_name + '\", {}, 2'.format(npts))
        self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"'.format(npts) + readout_seq_name + '\", 1')

        # Fill Readout waveforms into sequence
        for i in range(npts):
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet1:WAVeform \"'.format(str(i+1)) + readout_seq_name + '\", \"Readout_I\"')
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet2:WAVeform \"'.format(str(i+1)) + readout_seq_name + '\", \"Readout_Q\"')

        # Assign waveforms
        #self.awg.ch1.setSequenceTrack(self.filename, 1)
        #self.awg.ch2.setSequenceTrack(self.filename, 2)
        #self.awg.ch3.setSequenceTrack('Readout_Seq', 1)
        #self.awg.ch4.setSequenceTrack('Readout_Seq', 2)

        self.awg.ch1.setSequenceTrack(this_seq_name, 1)
        self.awg.ch2.setSequenceTrack(this_seq_name, 2)
        self.awg.ch3.setSequenceTrack(readout_seq_name, 1)
        self.awg.ch4.setSequenceTrack(readout_seq_name, 2)

        self.awg.ch2.state(1)
        self.awg.play()
        
        # Alazar labels
        for ala_chan in self.alazar_ctrl.channels[2:4]:
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('SSB Drive frequency (Some-overlap)',ala_chan.data.setpoint_labels[1])
            ala_chan.data.setpoint_units = ('Hz',ala_chan.data.setpoint_units[1])

        for n, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]):
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('SSB Drive frequency (Some-overlap)',)
            ala_chan.data.setpoint_units = ('Hz',)

        # prepare channels
        self.num_averages(self._averages)


    def MultiQ_Rabi_overlap(self, start, stop, npts, overlap = 0 ,custom_name = None):

        self.x_val = lambda: np.linspace(start ,stop ,npts)

        self.awg.ch2.state(0)

        # Here we decide whether we want to use "default" names or a user-specified custom name
        if custom_name == None:
            readout_seq_name = 'Readout_Seq'
            this_seq_name = self.filename
            self.awg.clearSequenceList()
            self.awg.clearWaveformList()
        else:
            this_seq_name = 'Waveforms_' + custom_name
            readout_seq_name = 'Readout_Seq_' + custom_name  
        
        N = int((self.cycle_time()*self.SR+64) - self.cycle_time()*self.SR%64)
        N_offset = int(self.marker_offset()*self.SR)
        overlap_array = np.linspace(0, overlap, int(overlap*self.SR), endpoint=False)
        
        # Create triggers
        ZerosMarker = np.zeros(int(N))
        TriggerMarker = np.zeros(int(N))
        TriggerMarker[-int(self.readout_dur()*self.SR-N_offset):-int(self.readout_dur()*self.SR-N_offset-500e-9*self.SR)] = 1
        
        # Create Drive tones
        wfms = [[]]
        for i , t in enumerate(self.x_val()):
            # SSB drive tone
            drive = np.zeros(int(N))
            drive[-int((self.readout_dur() + t-overlap)*self.SR):-int((self.readout_dur()-overlap)*self.SR)] = 0.5
            if i == 0:
                wfm_ch1 = np.array([drive,TriggerMarker,TriggerMarker])
            else:
                wfm_ch1 = np.array([drive,ZerosMarker,ZerosMarker])
            wfms[0].append(wfm_ch1)
            
        trig_waits = [0 for _ in range(npts)] 
        nreps = [1 for _ in range(npts)] 
        event_jumps = [0 for _ in range(npts)]
        event_jump_to = [0 for _ in range(npts)]
        go_to = [0 for _ in range(npts)]
        go_to[-1] = 1 # Make the sequence loop back to first step
        
        seqx = self.awg.makeSEQXFile(trig_waits,
                                nreps,
                                event_jumps,
                                event_jump_to,
                                go_to,
                                wfms,
                                [1, 1],
                                this_seq_name, custom_name = custom_name)

        self.awg.sendSEQXFile(seqx, this_seq_name + '.seqx')
        self.awg.loadSEQXFile(this_seq_name + '.seqx')
        
        # Create Readout tones
        self.update_readout_freqs()

        # Create sequence for readout tones SLISt:SEQuence:NEW <sequence_name>,<number_of_steps> [,<number_of_tracks>]
        #self.awg.write('SLISt:SEQuence:NEW \"Readout_Seq\", {}, 2'.format(npts))
        #self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"Readout_Seq\", 1'.format(npts))
        self.awg.write('SLISt:SEQuence:NEW \"' + readout_seq_name + '\", {}, 2'.format(npts))
        self.awg.write('SLISt:SEQuence:STEP{}:GOTO \"'.format(npts) + readout_seq_name + '\", 1')

        # Fill Readout waveforms into sequence
        for i in range(npts):
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet1:WAVeform \"'.format(str(i+1)) + readout_seq_name + '\", \"Readout_I\"')
            self.awg.write('SLISt:SEQuence:STEP{}:TASSet2:WAVeform \"'.format(str(i+1)) + readout_seq_name + '\", \"Readout_Q\"')

        # Assign waveforms
        self.awg.ch1.setSequenceTrack(this_seq_name, 1)
        self.awg.ch3.setSequenceTrack(readout_seq_name, 1)
        self.awg.ch4.setSequenceTrack(readout_seq_name, 2)
        self.awg.play()
        
        # Alazar labels
        for ala_chan in self.alazar_ctrl.channels[2:4]:
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Drive time',ala_chan.data.setpoint_labels[1])
            ala_chan.data.setpoint_units = ('s',ala_chan.data.setpoint_units[1])

        for n, ala_chan in enumerate(self.alazar_ctrl.channels[12:24]):
            ala_chan.records_per_buffer(npts)
            ala_chan.data.setpoint_labels = ('Drive time',)
            ala_chan.data.setpoint_units = ('s',)

        # prepare channels
        self.num_averages(self._averages)
        