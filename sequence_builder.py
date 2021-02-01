import numpy as np
import broadbean as bb
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

    def __init__(self,name,**kwargs):
        super().__init__(name, **kwargs)

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
            seqtemp.setSequencingGoto(i+1, 0)
        seqtemp.setSR(self.SR.get())
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
        seg_sin.insertSegment(2, ramp, (0, 0), name='read', dur=first_time)
        seg_sin.marker1 = [(first_time+self.pulse_time+self.marker_offset, self.cycle_time)]
        seg_sin.setSR(self.SR.get())
        
        return seg_sin