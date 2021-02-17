import qcodes as qc
from qcodes.instrument.base import Instrument, Parameter
import broadbean as bb
from broadbean.plotting import plotter
import numpy as np
from typing import Dict, Callable, List, Optional, Sequence

class BagOfBeans(Instrument):
    """
    Class that turns a broadbean Sequence into a QCoDeS Intrument 

        Parameters
        ----------
        seq (Parseq): Broadbean Sequence tured into a QCoDeS Parameter
        seq_path (str): Path to sequence file
        SR (float): sampling rate          
    """
    def __init__(self, name, **kwargs):
        super().__init__(name, **kwargs)
        self.add_parameter(name='seq',
                           label='Sequence',
                           parameter_class=ParSeq)        
        self.add_parameter(name='seq_path',
                           label='Path to sequence file',
                           get_cmd=None,
                           set_cmd=None) 
        self.add_parameter(name='SR',
                            label='sample rate',
                            unit='Hz',
                            parameter_class=sample_rate,
                            snapshot_exclude = True)

class ParSeq(Parameter):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.seq = bb.Sequence() 
    
    def snapshot_base(
        self, update: bool = False,
        params_to_skip_update: Optional[Sequence[str]] = None) -> Dict:
        return self.seq.description 
    
    def get_raw(self):
        return self.seq

    def set_raw(self,newseq):
        self.seq = newseq

    def from_file(self):
        seq_path = self.root_instrument.seq_path
        self.seq = bb.Sequence.init_from_json(seq_path)
    
    def to_file(self):
        seq_path = self.root_instrument.seq_path
        self.seq.write_to_json(seq_path)

    def empty_sequence(self):
        SR = self.seq.SR
        self.seq = bb.Sequence()
        self.seq.setSR(SR)

    def set_all_channel_amplitude_offset(self,amplitude:float = 1,offset:float = 0) -> None:
        for chan in self.seq.channels:
            self.seq.setChannelAmplitude(chan,amplitude)
            self.seq.setChannelOffset(chan,offset)

    def plot(self):
        plotter(self.seq)

    def plot_elem_nr(self,elem_nr):
        plotter(self.seq.element(elem_nr))

class sample_rate(Parameter):
    
    def get_raw(self):
        seq = self.root_instrument.seq()
        return seq.SR
    
    def set_raw(self,SR):
        seq = self.root_instrument.seq()
        seq.setSR(SR)

