import numpy as np
import broadbean as bb
import time
from qcodes.instrument.base import Instrument
from qcodes import validators as vals
from back_of_beans import BagOfBeans


# The pulsebuilding module comes with a (small) collection of functions appropriate for being segments.
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
    """

    def __init__(self,name:str, awg:Instrument = None,**kwargs):
        super().__init__(name, **kwargs)
        self.awg = awg

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

    def MultiQ_SSB_Spec_NoOverlap(self, start:float, stop:float, npts:int) -> bb.Sequence():
        """ 
        Updates the broadbean sequence so it contains two channels with orthogonal sine/cosine pulses for an array of  frequencies

            args:
            start (float): Starting point of the frequency interval
            stop (float): Endpoint point of the frequency interval
            npts (int): Number of point in the frequency interval
        """
        
        seqtemp = bb.Sequence()
        
        freq_interval = np.linspace(start,stop,npts)

        for i,f in enumerate(freq_interval):
            elem = bb.Element()
            seg_sin = self.seg_sine(frequency = f)
            seg_cos = self.seg_sine(frequency = f, phase=np.pi/2)
            elem.addBluePrint(1, seg_sin)
            elem.addBluePrint(2, seg_cos)
            seqtemp.addElement(i+1, elem)
            seqtemp.setSequencingTriggerWait(i+1, 0)
            seqtemp.setSequencingNumberOfRepetitions(i+1, 0)
            seqtemp.setSequencingEventJumpTarget(i+1, 0)
            if i == npts-1:
                seqtemp.setSequencingGoto(i+1, -1)
            else:
                seqtemp.setSequencingGoto(i+1, 0)
        seqtemp.setSR(self.SR.get())


    def MultiQ_Lifetime_overlap(self, start:float, stop:float, npts:int) -> bb.Sequence():
        """ 
        Updates the broadbean sequence so it contains one channels containing a pi-pulse
        varying the time between the end of the pi-pulse and the readout
        
            args:
            start (float): Starting point of the delta time
            stop (float): Endpoint point of the delta time
            npts (int): Number of point in the time interval
        """
        
        seqtemp = bb.Sequence()
        
        pulse_to_readout_time = np.linspace(start,stop,npts)

        for i,delta_time in enumerate(pulse_to_readout_time):
            elem = bb.Element()
            seg_pi = self.seg_pi(delta_time)
            elem.addBluePrint(1, seg_pi)
            seqtemp.addElement(i+1, elem)
            seqtemp.setSequencingTriggerWait(i+1, 0)
            seqtemp.setSequencingNumberOfRepetitions(i+1, 0)
            seqtemp.setSequencingEventJumpTarget(i+1, 0)
            if i == npts-1:
                seqtemp.setSequencingGoto(i+1, -1)
            else:
                seqtemp.setSequencingGoto(i+1, 0)
        seqtemp.setSR(self.SR.get())
      
        for chan in seqtemp.channels:
            seqtemp.setChannelAmplitude(chan,1)
            seqtemp.setChannelOffset(chan,1)
        self.seq.set(seqtemp)  

    def seg_sine(self,
                frequency:float,
                phase:float = 0) -> bb.BluePrint:
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
        seg_sin.insertSegment(1, ramp, (0.5, 0.5), name='pulse', dur=self.pulse_time)
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
                    chan[i].state(1)
                else:
                    chan[i].state(0)
            self.awg.play()


    def repeat_settings_infinity_loop(seqtemp: bb.Sequence,npts:int) -> bb.Sequence: 
        seqtemp.setSequencingTriggerWait(i+1, 0)
        seqtemp.setSequencingNumberOfRepetitions(i+1, 0)
        seqtemp.setSequencingEventJumpTarget(i+1, 0)
        if i == npts-1:
            seqtemp.setSequencingGoto(i+1, -1)
        else:
            seqtemp.setSequencingGoto(i+1, 0)
        return seqtemp
