import numpy as np
import broadbean as bb
import time
from qcodes.instrument.base import Instrument
from qcodes import validators as vals
from back_of_beans import BagOfBeans



ramp = bb.PulseAtoms.ramp  # args: start, stop
sine = bb.PulseAtoms.sine  # args: freq, ampl, off, phase

class SequenceBuilder(BagOfBeans):
    """
    Class for generating Sequences with predefined patterns

        Attributes
        ----------
        cycle_time : float
            total time of each cycle 
        pulse_time : float
            duration of the pulse
        readout_time : float
            duration of the readout
        marker_offset : float
            releative offset with respect to the readout_time
        SR: float
            sampling rate 

        Methods
        -------
        MultiQ_SSB_Spec_NoOverlap
            sequence of two channels with orthogonal sine/cosine pulses
        MultiQ_Lifetime_overlap
            One channels containing a pi-pulse
            varying the time between the end of the pi-pulse and the readout
    """

    def __init__(self,name:str,number_read_freqs:int = 1, alazar = None, alazar_ctrl = None , awg = None,qubit  = None,cavity = None,**kwargs):
        super().__init__(name, **kwargs)
        self.awg = awg
        self.qubit = qubit
        self.cavity = cavity
        self.alazar = alazar
        self.alazar_ctrl = alazar_ctrl

        self.add_parameter('cycle_time',
                      label='Pulse Cycle Time',
                      unit='s',
                      set_cmd= lambda x : x,
                      vals=vals.Numbers(0,10e-3))
        self.add_parameter('pulse_time',
                      label='Pulse Time',
                      unit='s',
                      set_cmd= lambda x : x,
                      vals=vals.Numbers(0,10e-3))                      
        self.add_parameter('readout_time',
                      label='Readout Duration',
                      unit='s',
                      set_cmd= lambda x : x,
                      vals=vals.Numbers(0,0.2e-3))
        self.add_parameter('marker_offset',
                      label='Marker Offset',
                      unit='s',
                      set_cmd= lambda x : x,
                      vals=vals.Numbers(-1e-5,1e-5))
        for i in range(number_read_freqs):
            self.add_parameter('readout_freq_{}'.format(i+1),
                        label='Readout Frequency {}'.format(i+1),
                        unit='Hz',
                        set_cmd= lambda x : x,
                        vals=vals.Numbers(0,12.5e9))

    def MultiQ_SSB_Spec_NoOverlap(self, start:float, stop:float, npts:int) -> None:
        """ 
        Updates the broadbean sequence so it contains two channels with orthogonal sine/cosine pulses for an array of  frequencies
        and two channels for the readout for IQ mixing 

            args:
            start (float): Starting point of the frequency interval
            stop (float): Endpoint point of the frequency interval
            npts (int): Number of point in the frequency interval
        """
        self.seq.empty_sequence()
        freq_interval = np.linspace(start,stop,npts)
        readout_freq = self.readout_freq_1.get() #- self.cavity.frequency()

        for i,f in enumerate(freq_interval):
            self.elem = bb.Element()
            if i == 0:
                seg_sin = self.seg_sine(frequency = f,marker=True)
            else:
                seg_sin = self.seg_sine(frequency = f, marker=False)
            seg_cos = self.seg_sine(frequency = f, phase=np.pi/2)
            self.elem.addBluePrint(1, seg_sin)
            self.elem.addBluePrint(2, seg_cos)
            self.elem_add_readout_pulse(readout_freq)
            self.seq.seq.addElement(i+1, self.elem)
            self.seq_settings_infinity_loop(i+1,npts)
        self.seq.seq.setSR(self.SR.get())

        self.seq.set_all_channel_amplitude_offset(amplitude=1, offset=0)
   

    def MultiQ_Lifetime_overlap(self, start:float, stop:float, npts:int) -> None:
        """ 
        Updates the broadbean sequence so it contains one channels containing a pi-pulse
        varying the time between the end of the pi-pulse and the readout
        and two channels for the readout for IQ mixing 
        
            args:
            start (float): Starting point of the delta time
            stop (float): Endpoint point of the delta time
            npts (int): Number of point in the time interval
        """
        
        self.seq.empty_sequence()
        readout_freq = self.readout_freq_1.get()
        pulse_to_readout_time = np.linspace(start,stop,npts)
        readout_freq = self.readout_freq_1.get() #- self.cavity.frequency()
        for i,delta_time in enumerate(pulse_to_readout_time):
            self.elem = bb.Element()
            seg_pi = self.seg_pi(delta_time)
            self.elem.addBluePrint(1, seg_pi)
            self.elem_add_readout_pulse(readout_freq)
            self.seq.seq.addElement(i+1,self.elem)

            self.seq_settings_infinity_loop(i+1,npts)
        self.seq.seq.setSR(self.SR.get())
      
        self.seq.set_all_channel_amplitude_offset(amplitude=1, offset=0)
        
    def test_station(self, start:float, stop:float, npts:int,channel: int) -> None:
        """ 
        Updates the broadbean sequence so it containsone channel with a sine pulse for an array of  frequencies

            args:
            start (float): Starting point of the frequency interval
            stop (float): Endpoint point of the frequency interval
            npts (int): Number of point in the frequency interval
            channel (int): The Channel of the seq/AWG
        """
        self.seq.empty_sequence()
        freq_interval = np.linspace(start,stop,npts)

        for i,f in enumerate(freq_interval):
            elem = bb.Element()
            seg_sin = self.seg_sine(frequency = f)
            elem.addBluePrint(channel, seg_sin)
            self.seq.seq.addElement(i+1, elem)
            self.seq_settings_infinity_loop(i+1,npts)
        self.seq.seq.setSR(self.SR.get())

        self.seq.set_all_channel_amplitude_offset(amplitude=1, offset=0)


    def seg_sine(self,
                frequency:float,
                phase:float = 0,
                marker:bool = False) -> bb.BluePrint:
        """
        Returns a broadbean BluePrint contaning a flat segment, sine segment and a flat segment for readout

        args:
        frequency (float): frequency of the sine 
        phase (float): phase of the sine 
        """
        
        first_time = self.cycle_time-self.pulse_time-self.readout_time 
        
        seg_sin = bb.BluePrint()
        seg_sin.insertSegment(0, ramp, (0, 0), name='first', dur=first_time)
        seg_sin.insertSegment(1, sine, (frequency, 1e-3, 0, phase), name='pulse', dur=self.pulse_time)
        seg_sin.insertSegment(2, ramp, (0, 0), name='read', dur=self.readout_time)
        if marker:
            seg_sin.marker1 = [(first_time+self.pulse_time+self.marker_offset, self.cycle_time)]
        seg_sin.setSR(self.SR.get())
        
        return seg_sin

    def seg_pi(self,
                pulse_to_readout_time:float = 0) -> bb.BluePrint:
        """
        Returns a broadbean BluePrint of a PI pulse 

        args:
        pulse_to_readout_time (float): time between the end of the PI pulse and the readout  
        """
        
        first_time = self.cycle_time-self.pulse_time-self.readout_time-pulse_to_readout_time 
        end_time = self.readout_time+pulse_to_readout_time
        
        seg_sin = bb.BluePrint()
        seg_sin.insertSegment(0, ramp, (0, 0), name='first', dur=first_time)
        seg_sin.insertSegment(1, ramp, (0.05, 0.05), name='pulse', dur=self.pulse_time)
        seg_sin.insertSegment(2, ramp, (0, 0), name='read', dur=end_time)
        seg_sin.marker1 = [(first_time+self.pulse_time+self.marker_offset+pulse_to_readout_time, self.cycle_time)]
        seg_sin.setSR(self.SR.get())
        
        return seg_sin        

    def uploadToAWG(self, awg_amp: list = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]) -> None:
        if '5014' in str(self.awg.__class__):
            #for i,  chan in enumerate(self.seq.get().channels):
            #    self.awg.channels[chan].AMP(float(chbox[chan-1].text()))
            self.awg.ch1_amp(awg_amp[0])
            self.awg.ch2_amp(awg_amp[1])
            self.awg.ch3_amp(awg_amp[2])
            self.awg.ch4_amp(awg_amp[3])
            package = self.seq.get().outputForAWGFile()
            start_time = time.time()
            self.awg.make_send_and_load_awg_file(*package[:])
            print("Sequence uploaded in %s seconds" %(time.time()-start_time));
        elif '5208' in str(self.awg.__class__):
            self.seq.get().name = 'sequence_from_gui'
            self.awg.mode('AWG')
            for chan in self.seq.get().channels:
                self.awg.channels[chan-1].resolution(12)
                self.awg.channels[chan-1].awg_amplitude(awg_amp[chan-1])
                self.seq.get().setChannelAmplitude(chan, self.awg.channels[chan-1].awg_amplitude())
            self.awg.clearSequenceList()
            self.awg.clearWaveformList()
            self.awg.sample_rate(self.seq.get().SR)
            self.awg.sample_rate(self.seq.get().SR)
            
            seqx_input = self.seq.get().outputForSEQXFile()
            start_time=time.time()
            seqx_output = self.awg.makeSEQXFile(*seqx_input)
            # transfer it to the awg harddrive
            self.awg.sendSEQXFile(seqx_output, 'sequence_from_gui.seqx')
            self.awg.loadSEQXFile('sequence_from_gui.seqx')
            #time.sleep(1.300)
            for i,  chan in enumerate(self.seq.get().channels):       
                self.awg.channels[chan-1].setSequenceTrack('sequence_from_gui', i+1)
                self.awg.channels[chan-1].state(1)
            print("Sequence uploaded in %s seconds" %(time.time()-start_time))
 
        else:
            print('Choose an AWG model')

    def runAWG(self):
        if '5014' in str(self.awg.__class__):
            self.awg.run()
        else:
            seq_chan = self.seq.get().channels
            for i, chan in enumerate(self.awg.channels):
                if i+1 in seq_chan:
                    chan.state(1)
                else:
                    chan.state(0)
            self.awg.play()


    def seq_settings_infinity_loop(self, elem_nr:int, last_elem_nr:int) -> None:
        """
        Play element 1 time and go to the next,
        except if you are the last element, then play 1 time and go to the first Element.

        args:
        elem_nr (int): the number of the element
        last_elem_nr (int): the number of the last element in the sequence  
        """
        self.seq.seq.setSequencingTriggerWait(elem_nr, 0)
        self.seq.seq.setSequencingNumberOfRepetitions(elem_nr, 1)
        self.seq.seq.setSequencingEventJumpTarget(elem_nr, 0)
        if elem_nr == last_elem_nr:
            self.seq.seq.setSequencingGoto(elem_nr, 1)
        else:
            self.seq.seq.setSequencingGoto(elem_nr, 0)


    def elem_add_readout_pulse(self, frequency:float, amplitude:float = 0.05):
        seg_sin_readout = self.seg_sine_readout(frequency=frequency, amplitude=amplitude, marker=False)
        seg_cos_readout = self.seg_sine_readout(frequency=frequency, amplitude=amplitude, phase=np.pi/2 ,marker=True)
        self.elem.addBluePrint(3,seg_sin_readout)
        self.elem.addBluePrint(4,seg_cos_readout)


    def seg_sine_readout(self,
                frequency:float,
                phase:float = 0,
                amplitude:float = 1e-3,
                marker:bool = True ) -> bb.BluePrint:
        """
        Returns a broadbean BluePrint contaning a flat segment, sine segment and a flat segment for readout

        args:
        frequency (float): frequency of the sine 
        phase (float): phase of the sine 
        """
        
        first_time = self.cycle_time-self.readout_time 
        
        seg_sin = bb.BluePrint()
        seg_sin.insertSegment(0, ramp, (0, 0), name='first', dur=first_time)
        seg_sin.insertSegment(1, sine, (frequency, amplitude, 0, phase), name='pulse', dur=self.readout_time)
        if marker:
            seg_sin.marker1 = [(first_time+self.marker_offset, self.cycle_time)]
        seg_sin.setSR(self.SR.get())
        
        return seg_sin

