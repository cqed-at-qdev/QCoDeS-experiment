# Configure all settings in the Alazar card
def alazarconfig(alazar):
    """
    function for config of alazar
    """
    with alazar.syncing():    
        alazar.clock_source('INTERNAL_CLOCK')
        alazar.sample_rate(1_000_000_000)
        alazar.clock_edge('CLOCK_EDGE_RISING')
        alazar.decimation(1)
        alazar.coupling1('DC')
        alazar.coupling2('DC')
        alazar.channel_range1(.4)
        alazar.channel_range2(.4)
        alazar.impedance1(50)
        alazar.impedance2(50)
        alazar.trigger_operation('TRIG_ENGINE_OP_J')
        alazar.trigger_engine1('TRIG_ENGINE_J')
        alazar.trigger_source1('EXTERNAL')
        alazar.trigger_slope1('TRIG_SLOPE_POSITIVE')
        alazar.trigger_level1(160)
        alazar.trigger_engine2('TRIG_ENGINE_K')
        alazar.trigger_source2('DISABLE')
        alazar.trigger_slope2('TRIG_SLOPE_POSITIVE')
        alazar.trigger_level2(128)
        alazar.external_trigger_coupling('DC')
        alazar.external_trigger_range('ETR_2V5')
        alazar.trigger_delay(0)
        alazar.timeout_ticks(0)
        alazar.aux_io_mode('AUX_IN_AUXILIARY') # AUX_IN_TRIGGER_ENABLE for seq mode on
        alazar.aux_io_param('NONE') # TRIG_SLOPE_POSITIVE for seq mode on


def alazarconfigttwo(alazar):
    """
    function for config of alazar
    """
    with alazar.syncing():    
        alazar.clock_source('EXTERNAL_CLOCK_10MHz_REF')
        alazar.external_sample_rate(500_000_000)
        alazar.clock_edge('CLOCK_EDGE_RISING')
        alazar.decimation(1)
        alazar.coupling1('DC')
        alazar.coupling2('DC')
        #alazar.channel_range1(.4)
        #alazar.channel_range2(.4)
        #alazar.impedance1(50)
        #alazar.impedance2(50)
        alazar.trigger_operation('TRIG_ENGINE_OP_J')
        alazar.trigger_engine1('TRIG_ENGINE_J')
        alazar.trigger_source1('EXTERNAL')
        alazar.trigger_slope1('TRIG_SLOPE_POSITIVE')
        alazar.trigger_level1(140)
        alazar.trigger_engine2('TRIG_ENGINE_K')
        alazar.trigger_source2('DISABLE')
        alazar.trigger_slope2('TRIG_SLOPE_POSITIVE')
        alazar.trigger_level2(140)
        alazar.external_trigger_coupling('DC')
        alazar.external_trigger_range('ETR_2V5')
        alazar.trigger_delay(0)
        alazar.timeout_ticks(0)
        alazar.aux_io_mode('AUX_IN_TRIGGER_ENABLE') # AUX_IN_TRIGGER_ENABLE for seq mode on
        alazar.aux_io_param('TRIG_SLOPE_POSITIVE') # TRIG_SLOPE_POSITIVE for seq mode on
        #alazar.aux_io_mode('AUX_IN_AUXILIARY') # AUX_IN_TRIGGER_ENABLE for seq mode on
        #alazar.aux_io_param('NONE') # TRIG_SLOPE_POSITIVE for seq mode on


         #parameters:
      #clock_source: {initial_value: 'EXTERNAL_CLOCK_10MHz_REF'}
      #external_sample_rate: {initial_value: 500_000_000} # We normally use 500_000_000
      #coupling1: {initial_value: 'DC'}
      #coupling2: {initial_value: 'DC'}
      #trigger_operation: {initial_value: 'TRIG_ENGINE_OP_J'}
      #trigger_engine1: {initial_value: 'TRIG_ENGINE_J'}
      #trigger_source1: {initial_value: 'EXTERNAL'}
      #trigger_slope1: {initial_value: 'TRIG_SLOPE_POSITIVE'}
      #trigger_level1: {initial_value: 140}
      #trigger_source2: {initial_value: 'DISABLE'}
      #trigger_slope2: {initial_value: 'TRIG_SLOPE_POSITIVE'}
      #trigger_level2: {initial_value: 140}
      #external_trigger_coupling: {initial_value: 'DC'}
      #external_trigger_range: {initial_value: 'ETR_2V5'}
      #trigger_delay: {initial_value: 0}
      #timeout_ticks: {initial_value: 0}
      #aux_io_mode: {initial_value: AUX_IN_AUXILIARY}
      #aux_io_param: {initial_value: NONE}
      ##seq_mode: {monitor: false}       
        
        
        
        

        