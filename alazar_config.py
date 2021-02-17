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
        
        
        
        

        