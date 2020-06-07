function sm_pulleys_xytable_cross_configmotion(mdlname,config)
% Copyright 2018 The MathWorks, Inc.

if(strcmpi(config,'XY Position'))
    set_param([mdlname '/Platform/Prescribe Motion'],...
        'OverrideUsingVariant','On');
    set_param([mdlname '/Pulleys/Revolute P2'],...
        'OverrideUsingVariant','Measure');
    set_param([mdlname '/Pulleys/Revolute P6'],...
        'OverrideUsingVariant','Measure');
    set_param([mdlname '/Pulleys'],...
        'MaskDisplay','image(''sm_pulleys_xytable_cross_pulleys_IMAGE.png'');');
    set_param([mdlname '/Platform'],...
        'MaskDisplay','image(''sm_pulleys_xytable_cross_platformXYMotion_IMAGE.png'');');
    set_param([mdlname '/Motion'],...
        'MaskDisplay','image(''sm_pulleys_xytable_cross_motionPulleyAngles_IMAGE.png'');');
    set_param([mdlname '/Motion'],...
        'MaskDisplay','image(''sm_pulleys_xytable_cross_motionXYMotion_IMAGE.png'');');
    
elseif(strcmpi(config,'Pulley Angles'))
    set_param([mdlname '/Platform/Prescribe Motion'],...
        'OverrideUsingVariant','Off');
    set_param([mdlname '/Pulleys/Revolute P2'],...
        'OverrideUsingVariant','Actuate');
    set_param([mdlname '/Pulleys/Revolute P6'],...
        'OverrideUsingVariant','Actuate');
    set_param([mdlname '/Pulleys'],...
        'MaskDisplay','image(''sm_pulleys_xytable_cross_pulleysPulleyAngles_IMAGE.png'');');
    set_param([mdlname '/Platform'],...
        'MaskDisplay','image(''sm_pulleys_xytable_cross_platform_IMAGE.png'');');
    set_param([mdlname '/Motion'],...
        'MaskDisplay','image(''sm_pulleys_xytable_cross_motionPulleyAngles_IMAGE.png'');');
end
    
    
    