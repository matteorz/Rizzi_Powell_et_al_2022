function fn = tCSF_flash_v2_0_3_annulus_spatialfrequency_gabor(metaParams, graphicParams, stimParams, psychometricParams, eyetrackParams, audioParams)    
% Script for measuring temporal contrast sensitivity function (tCSF), using
% an eyetracker to ensure fixed eccentricity / steady fixation
%
% Requires:         ivis v1.3
%                   palamedes v1.7.0
%                   PTB-3
%                   Tobii SDK 3.0
%                   PsychTestRig
%                   CRS Visage Toolbox
%
% Matlab:           v2012 onwards
%
% Author(s):    	Pete R Jones <petejonze@gmail.com>
% 

% Version History:  1.0.0	PJ  29/03/2016 	Initial build.
%                   1.1.0	PJ  22/08/2016	Added option to specify
%                                           multiple cone types (L, M, LM, S, LS, etc.) and also to
%                                           do a YesNo version of the task.
%                   1.2.0	PJ  01/10/2016	Add screen masking options.
%                                           Fixed bugs
%                                               
%
% Copyright 2016 : P R Jones
% *********************************************************************
% 

% @todo
% add calibration check?
% add full calibration? (probably not necessary)

USE_EYETRACKER = false;
                
                
	%%%%%%%%%
    %% 1.1 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Basic init %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %-------------Check OS/Matlab version------------------------------
        if ~strcmpi(computer(),'PCWIN')
            error('This code has only been tested on Windows 7 running 32-bit Matlab\n  Detected architecture: %s', computer());
        end
        
       	%-------------Ready workspace-------------------------------------- 
        tmp = javaclasspath('-dynamic');
        clearJavaMem();
        close all;
        if length(tmp) ~= length(javaclasspath('-dynamic'))
            % MATLAB calls the clear java command whenever you change
            % the dynamic path. This command clears the definitions of
            % all Java classes defined by files on the dynamic class
            % path, removes all variables from the base workspace, and
            % removes all compiled scripts, functions, and
            % MEX-functions from memory.
            error('clearJavaMem:MemoryCleared','clearJavaMem has modified the java classpath (any items in memory will have been cleared)\nWill abort, since this is highly likely to lead to errors later.\nTry running again, or see ''help PsychJavaTrouble'' for a more permenant solution\n\ntl;dr: Try running again.');
        end

        %-------------Check for requisite toolkits-------------------------
        AssertOpenGL();             % PTB-3 correctly installed? Abort otherwise.
    	PAL_version();              % check we have palamedes installed
        assertPTRversion(0.7);      % check we're using the right PsychTestRig version
        IvMain.assertVersion(1.5);  % check we're using the right ivis version

        %-------------Check classpath--------------------------------------
        ivis.main.IvMain.checkClassPath();

        %-------------Hardcoded User params--------------------------------  
%         IN_DEBUG_MODE = false; % true;
%         % media files
%         RESOURCES_DIR   = fullfile('..', 'resources');
%         SND_DIR         = fullfile(RESOURCES_DIR, 'audio', 'wav');
%         IMG_DIR         = fullfile(RESOURCES_DIR, 'images');
        % eye-tracking log files
        LOG_RAW_DIR     = fullfile('..', '..', 'data', '__EYETRACKING', 'raw');
        LOG_DAT_DIR     = fullfile('..', '..', 'data', '__EYETRACKING', 'data');
        
        %-------------Add any requisite paths------------------------------ 
        import ivis.main.* ivis.classifier.* ivis.video.*  ivis.broadcaster.* ivis.math.* ivis.graphic.* ivis.audio.* ivis.log.* ivis.calibration.*;

        
	%%%%%%%%%
    %% 1.2 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Validate user inputs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fprintf('\nVaidating inputs...\n');

        %-------------graphicParams----------------------------------------
        p = inputParser; p.StructExpand = true;
        p.addParameter('screenNum',                   	@isNonNegativeInt);
        p.addParameter('Fr',                         	@isPositiveInt);
        p.addParameter('screenWidth_px',              	@isPositiveInt);
        p.addParameter('screenHeight_px',            	@isPositiveInt);
        p.addParameter('screenWidth_cm',              	@isPositiveNum);
        p.addParameter('screenHeight_cm',              	@isPositiveNum);
        p.addParameter('assumedViewingDistance_cm',    	@isPositiveNum);
        p.addParameter('COMMENT',                     	@ischar);
        p.parse(graphicParams);
                
        %-------------stimParams-------------------------------------------
        p = inputParser; p.StructExpand = true;
        p.addParameter('eccentricity_x_deg',          	@(x)all(isPositiveNum(x) | isnan(x)));
        p.addParameter('eccentricity_y_deg',          	@(x)all(isPositiveNum(x) | isnan(x)));
        p.addParameter('temporalfreq_hz',           	@(x)all(isPositiveNum(x)));
        p.addParameter('sensors_fn',                    @ischar);
        p.addParameter('deviceSPD_fn',                  @ischar);
        p.addParameter('coneType',                    	@(x)ismember(x,{'L','M','S'}) );
        p.addParameter('diameter_deg',                  @isnumeric);  
        p.addParameter('useOverlay',                    @(x)ismember(x,[0 1 2]));
        p.addParameter('apperture_x_deg',               @(x)all(isPositiveNum(x) | isnan(x)));
        p.addParameter('apperture_y_deg',               @(x)all(isPositiveNum(x) | isnan(x)));
        p.addParameter('apperture_deg',                 @isnumeric);
        p.addParameter('cosRamp_sec',                   @(x)isPositiveNum(x));
        p.addParameter('COMMENT',                     	@ischar);
        p.parse(stimParams);

        %-------------psychometricParams-----------------------------------
        p = inputParser; p.StructExpand = true;
        p.addParameter('paradigm',                      @(x)ismember(x,{'QUEST','staircase','MCS'}) );
     	p.addParameter('useLogScale',                   @islogical);
     	p.addParameter('qparams',                       @(x)1==1);
        p.addParameter('scparams',                      @(x)1==1);
        p.addParameter('mcsparams',                     @(x)1==1);
        p.addParameter('instructions',                 	@ischar );
        p.addParameter('breakAfterNBlocks',            	@(x)all(isPositiveNum(x)) );
        p.addParameter('conditionSequence',           	@(x)ismember(x,{'inOrder','inRandOrder','interleaved'}) );
        p.addParameter('preInt_sec',                  	@isPositiveNum);
        p.addParameter('stimInt_sec',                  	@isPositiveNum);
        p.addParameter('fbackInt_sec',                 	@isPositiveNum);
        p.addParameter('allowRespBeforeEndOfStim',  	@islogical);
     	p.addParameter('useStaticWNoisePostMask',       @islogical);
        p.addParameter('COMMENT',                     	@ischar);
        p.parse(psychometricParams);

        %-------------eyetrackParams---------------------------------------
        p = inputParser; p.StructExpand = true;
        p.addParameter('ivisVersion',                  	@isPositiveNum);
        p.addParameter('type',                        	@(x)ismember(x,{'tobii','mouse'}));
        p.addParameter('id',                            @ischar);
        p.addParameter('Fs',                            @(x)ismember(x,[60 120]));
        p.addParameter('eye',                        	@(x)ismember(x,[0 1 2]));
        p.addParameter('useGUI',                      	@islogical);
        p.addParameter('gaze_nOutliersPermitted',     	@isNonNegativeInt);
        p.addParameter('gaze_outlierCriterion_px',     	@isNonNegativeInt);
        p.addParameter('gaze_nMissingSamplesPermitted',	@isNonNegativeInt);
        p.addParameter('COMMENT',                       @ischar);
        p.parse(eyetrackParams);
        
        %-------------audioParams------------------------------------------
        p = inputParser; p.StructExpand = true;
        p.addParameter('devID',                       	@(x)isempty(x) || isNonNegativeInt(x));
        p.addParameter('COMMENT',                    	@ischar);
        p.parse(audioParams);
    
        
	%%%%%%%%%
    %% 1.3 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Further validation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Check psychometric params
        switch lower(psychometricParams.paradigm)
            case 'quest'
                tmp = QUESTwrapper(psychometricParams.qparams); %#ok
            case 'staircase'
                tmp = AdaptiveTrack(psychometricParams.scparams); %#ok
            case 'mcs'
                tmp = ConstantStimuli(psychometricParams.mcsparams.levels, psychometricParams.mcsparams.nTrialsPerLevel); %#ok
                if any(psychometricParams.mcsparams.levels>100 | psychometricParams.mcsparams.levels<=0)
                    error('all contrast levels must be between 0 and 100 (0 < x <= 100)');
                end
            otherwise % defensive
                error('Unknown paradigm: %s', psychometricParams.paradigm)
        end
        clear tmp;
        
        % set nConditions
        psychometricParams.nConditions = numel(stimParams.temporalfreq_hz);
                
        % check sequence
        switch lower(psychometricParams.conditionSequence)
            case 'inorder'
                % do nothing
            case 'inrandorder'
                stimParams.temporalfreq_hz = Shuffle(stimParams.temporalfreq_hz);
            case 'interleaved'
                error('functionality not yet written');
            otherwise % defensive
                error('Unknown sequence option: %s', psychometricParams.conditionSequence)
        end
        
        % check all temporal refresh rates involve whole-numbers of frames
        if any(abs(graphicParams.Fr./stimParams.temporalfreq_hz - round(graphicParams.Fr./stimParams.temporalfreq_hz)) > 0.0001)
            error('all temporal refresh rates must involve whole-numbers of frames');
        elseif any(graphicParams.Fr./stimParams.temporalfreq_hz ~= round(graphicParams.Fr./stimParams.temporalfreq_hz))
            % correct for any minor rounding errors
            stimParams.temporalfreq_hz = graphicParams.Fr./round(graphicParams.Fr./stimParams.temporalfreq_hz);
        end
        
        
	%%%%%%%%%
    %% 2.2 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Hardware: Set up the Tobii X120 eyetracker %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if USE_EYETRACKER
            % initialise, and launch the ivis toolbox
            IvMain.initialise(IvParams.getSimpleConfig('graphics.testScreenNum',0, 'GUI.useGUI',false, 'eyetracker.type',eyetrackParams.type, 'eyetracker.id',eyetrackParams.id, 'keyboard.handlerClass',NaN, 'eyetracker.fixationMarker','cursor', 'log.raw.dir',LOG_RAW_DIR, 'log.data.dir',LOG_DAT_DIR));
            [eyetracker, logs] = IvMain.launch();
        end
                
	%%%%%%%%%
    %% 2.3 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Hardware: Set up the CRS VSG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %-------------Set up the CRS VSG----------------------------------- 
        global CRS;
        errCode = vsgInit;
        if errCode<0
            error('VSG not working???');
        end
        crsSetSpatialUnits(CRS.DEGREEUNIT);
        
        % Colour Definitions
        % ==================
        RedRGB      = [1.0,0.0,0.0];
        GreenRGB	= [0.0,1.0,0.0];
        BlueRGB     = [0.0,0.0,1.0];
        BlackRGB    = [0.0,0.0,0.0];
        WhiteRGB    = [1.0,1.0,1.0];
        BackgroundRGB = [0.5,0.5,0.5];

        midpoint = 0.5; %around what should we vary the contrasts
        
        MaskRGB    = [0.5,0.5,0.5]; %Set background of screen
        AnnulusRGB=[0,0,0]; %Set annulus contrast
        controlRGB = [midpoint,midpoint,midpoint]; %set control oval contrast (want to keep this grey)
        
        Contrast_prcnt = 10; %maintain stimulus contrast at 100%
        
        
        %Set section of pallet to use for static components
        AnnulusLevel = 1;
        BackgroundLevel= 2;
        controlLutLevel = 3;
        %SpatialFrequency = 0.2/1; %cycles per pixel
        
        padval = 30; %how many px to pad everything by
          
        % Silent Substitution: Parameters
        %=========================================================
        % LMS
        Sensors = stimParams.sensors_fn; % e.g., 'ConeSensitivities_SS_2degELin3908301.mat';
        % Device SPD
        deviceSPD = stimParams.deviceSPD_fn; % e.g., 'RizziInterp.mat';
        
        % Silent Substitution: Set cone modulation
        %=========================================================
        % Define which cone will be modulated.
        % The order is: cone_modulation = [L M S];
        % The value different from 0 indicates which cone-class will be isolated.
        % The value equal to 0 indicates which cone-classes will generate a constant response.
        % In this demo, L-cones will be isolated.
        % The modulation value represent the contrast from the background.
        % Make sure that the modulation value you want to use will generate RGB
        % values in the monitor gamut.
        %[~,coneIdx] = ismember(stimParams.coneType, {'L','M','S'});
        coneIdx = nan(1,3);
        coneIdx(1) = ~isempty(regexp(stimParams.coneType,'L','once'));
        coneIdx(2) = ~isempty(regexp(stimParams.coneType,'M','once'));
        coneIdx(3) = ~isempty(regexp(stimParams.coneType,'S','once'));
        coneIdx = coneIdx==1;
        
        %=========================================================
        % DEFINE PIXEL LEVELS.
        %=========================================================
%         stimulusLutLevel = 1;
%         crsPaletteSetPixelLevel(stimulusLutLevel, [0.5,0.5,0.5]);
%         controlLutLevel=2;
%         crsPaletteSetPixelLevel(controlLutLevel, [0.5,0.5,0.5]);
        minsine = 0;
        maxsine = 1;
        sinepal = minsine:(maxsine-minsine)./length(4:(253-1)):maxsine;

        %=========================================================
        % Compute background LMS color
        %=========================================================
        fromCS               = 'CS_RGB';
        toCS                 = 'CS_LMS';
        fromColour           = BackgroundRGB;
        [toColour ErrorCode] = ctGetColourTrival(fromCS,toCS,fromColour,deviceSPD,Sensors);
        BackgroundLMS        = toColour;

        

        %=========================================================
        % Check that max contrast is within range
        %=========================================================

        % determine max contrast that yeilds a valid RGB vector
        cone_modulation = [0 0 0]';
        cone_modulation(coneIdx) = 1;
        N = 100; % increase for greater fidelity
        stimulusRGB = nan(N, 3);
        contrasts = linspace(0,1,N);
        for i = 1:N
            stimulusLMS  = BackgroundLMS.*(1 + cone_modulation*contrasts(i)); % !!!!!!! Set Contrast
            stimulusRGB(i,:)  = ctGetColourTrival('CS_LMS','CS_RGB',stimulusLMS',deviceSPD,Sensors);
        end
        isValid = all(stimulusRGB>=0 & stimulusRGB<=1, 2);
        idx = find(isValid,1,'last');
        %maxContrast_prcnt = 100*contrasts(idx);
        if idx < length(isValid)
            firstInvalidContrast_prcnt = 100*contrasts(idx+1);
            % throw error if user input is invalid
            %if psychometricParams.qparams.deltaMax >= firstInvalidContrast_prcnt
             %   error('Specified maximum contrast (%1.2f) is greater than the expected maximum possible contrast (%1.2f) for cone receptor: %s', psychometricParams.qparams.deltaMax, firstInvalidContrast_prcnt, stimParams.coneType)
            %end
        end
        
        
        
        
        % ColourMode
        % ==========
        % The ColourMode parameter selects how the contrast value is interpreted
        % and how it will affect the current object's colours.
        %
        % The ColourMode should be CRS.UNIPOLAR for Gaussian patches, and should
        % be CRS.BIPOLAR for Gratings and Gabors.
        %
        % The reason for this is that both Sine and Gabor functions use the notion
        % of a mid-level background, half way between the function extrema.
        %
        % Gaussian patches, on the other hand, sit on a background that is coincident
        % with the function minimum. These functions will therefore behave
        % differently when the contrast is changed.
        ColourMode = CRS.UNIPOLAR;
        
        % Object pixel levels
        % ===================
        % The Object is animated by modifying a section of the palette. Which
        % section to modify is given by the following parameters:
        MinValue  = 4;
        MaxValue  = 253;
        NumValues = (MaxValue - MinValue);
        
        % System colours
        % ==============
        % Some palette entries are reserved for system colours. Here, we are using
        % a background colour and a fixation colour.
        %
        % For a CRS.BIPOLAR stimulus, the background colour should be half way
        % between the object minimum and maximum values. For a CRS.UNIPOLAR stimulus,
        % the background colour should be the same as the object minimum value.
        FixRGB   = [0 0 0];
        
        % Video Pages
        % ===========
        BlankPage    = 1;
        StimulusPage = 2;
        PrestimPage    = 3;
        PostStimPage = 4;
        BlahPage = 5;
        
        
      
        crsPaletteSetPixelLevel(BackgroundLevel, MaskRGB);
        crsPaletteSetPixelLevel(controlLutLevel, controlRGB);
        crsPaletteSetPixelLevel(AnnulusLevel, AnnulusRGB);
        crsPaletteSetPixelLevel(CRS.BACKGROUND, BackgroundRGB);
        crsPaletteSetPixelLevel(CRS.FIXATION, FixRGB);
        for k=4:253
            crsPaletteSetPixelLevel(k, repmat(sinepal(k-3),[1,3]));
        end
        
%         
        Buff_new = zeros(256,3);
        Buff_new(AnnulusLevel,1) = AnnulusRGB(1);
        Buff_new(AnnulusLevel,2) = AnnulusRGB(2);
        Buff_new(AnnulusLevel,3) = AnnulusRGB(3);
        Buff_new(controlLutLevel,1) = controlRGB(1);
        Buff_new(controlLutLevel,2) = controlRGB(2);
        Buff_new(controlLutLevel,3) = controlRGB(3);
        Buff_new(4:253,1) = sinepal;
        Buff_new(4:253,2) = sinepal;
        Buff_new(4:253,3) = sinepal;
        Buff_new(BackgroundLevel,1) = MaskRGB(1);
        Buff_new(BackgroundLevel,2) = MaskRGB(2);
        Buff_new(BackgroundLevel,3) = MaskRGB(3);
        %Buff_new(1,1) = MaskRGB(1); % ensure background remain visible/invariant
        %Buff_new(1,2) = MaskRGB(2);
        %Buff_new(1,3) = MaskRGB(3);
        %                 Buff_new(2,1) = WhiteRGB(1); % ensure background remain visible/invariant
        %                 Buff_new(2,2) = WhiteRGB(2);
        %                 Buff_new(2,3) = WhiteRGB(3);
        Buff_new(CRS.BACKGROUND,1) = BackgroundRGB(1); % ensure background remain visible/invariant
        Buff_new(CRS.BACKGROUND,2) = BackgroundRGB(2);
        Buff_new(CRS.BACKGROUND,3) = BackgroundRGB(3);
        %                 Buff_new(3,1) = AnnulusContrastRGB(1); % ensure annulus remains visible/invariant
        %                 Buff_new(3,2) = AnnulusContrastRGB(2);
        %                 Buff_new(3,3) = AnnulusContrastRGB(3);
        Buff_new(CRS.FIXATION,1) = FixRGB(1); % ensure fixation marker remains visible/invariant
        Buff_new(CRS.FIXATION,2) = FixRGB(2);
        Buff_new(CRS.FIXATION,3) = FixRGB(3);
        
        crsLUTBUFFERWrite(1,Buff_new);
        

        % Stimulus Parameters (some of these will get set later)
        % ===================
        ScreenSize     = crsGetScreenSizeDegrees;
        MinDimension   = min(ScreenSize(:));
        
        % Fixation Point Parameters (some of these will get overwritten later)
        % =========================
        fix_xy_deg       = [0,0]; % Degrees of visual angle (x,y)
        fix_diameter_deg = 11.5;   % Degrees of visual angle (diameter)
%         fix_diameter_deg = 1.3;
        % Check hardware running as expected
        % =========================
        % check expected viewing distance matches what the VSG thinks
        if crsGetViewDistMM ~= 10*graphicParams.assumedViewingDistance_cm
            error('VSG view distance (%1.2f mm) does not match the value expected in the config file (%1.2f cm)', crsGetViewDistMM, graphicParams.assumedViewingDistance_cm)
        end
        % check frame rate matches what the VSG thinks
        if crsGetFrameRate ~= graphicParams.Fr
            error('VSG refresh rate (%1.2f hz) does not match the value expected in the config file (%1.2f hz)', crsGetFrameRate, graphicParams.Fr)
        end
        
	%%%%%%%%%
    %% 2.3 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Hardware: Set up the keyboard input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        InH = MyInputHandler();
        
    
	%%%%%%%%%
    %% 3.1 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Final checks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        % check expected viewing distance matches the empirical data from
        % the Tobii
warning('add distance checks');        
% <do me>        
  
	try %#ok
        
	%%%%%%%%%
    %% 3.1 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Instructions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        crsSetSpatialUnits(CRS.PIXELUNIT);
        BackLevel = MinValue; % CRS.BACKGROUND;
        StimLevel = MaxValue;
        
        % Display a blank page while we work.
        crsSetDrawPage(BlankPage);
        crsClearPage(BlankPage,BackLevel);
        crsSetDisplayPage(BlankPage);
        
        % set text page
        crsSetDrawPage(StimulusPage);
        crsClearPage(StimulusPage,BackLevel);

        % Set the drawing origin to the top left of the screen
        crsSetDrawOrigin([100,100]);

        % Fill the palette with the colours that we have chosen
        %MinPixelLevel  = BackLevel;
        %NumPixelLevels = (StimLevel - MinPixelLevel) + 1;
        %crsPaletteSetSection(MinPixelLevel,linspace(BlackRGB(1),WhiteRGB(1),NumPixelLevels),...
        %    linspace(BlackRGB(2),WhiteRGB(2),NumPixelLevels),...
        %    linspace(BlackRGB(3),WhiteRGB(3),NumPixelLevels));

        % Set antialiasing on. (Needed to make strings look good).
        oldMode = crsGetDrawMode;
        Mode = oldMode;
        Mode.AntiAliasing = true;
        crsSetDrawMode(Mode);

        % Set our string parameters
        LineHeight  = 28;
        LineSpacing = 8;
        StringSize = [10,LineHeight-LineSpacing];
        Horizontal = CRS.ALIGNLEFTTEXT;
        Vertical   = CRS.ALIGNBOTTOMTEXT;
        Angle      = 0;
        Flags      = CRS.FONTNORMAL;
        crsSetStringMode(StringSize,Horizontal,Vertical,Angle,Flags);

        % As well as the font that we will use
        crsSetTrueTypeFont('Arial');

        % Select the section of the palette that we have chosen
        crsSetPen1(StimLevel);
        crsSetPen2(BackLevel);

        % *************************************

        % Draw the text (ASCII)
        crsDrawString([0,1 * LineHeight],psychometricParams.instructions);

        % *************************************

        % Finally, Display the stimulus that we have drawn
        crsSetDisplayPage(StimulusPage);

        % wait for appropriate input before continuing
        InH.waitForInput([InH.INPT_MIDDLE InH.INPT_LEFTARROW])
    
        % restore settings
        crsSetDrawMode(oldMode);
        crsSetSpatialUnits(CRS.DEGREEUNIT);
        % crsSetDrawOrigin(CRS.SCREENWIDTH/2, CRS.SCREENHEIGHT/2)  % these units are not right for some reason?
        crsSetDrawOrigin(ScreenSize(1)/2, ScreenSize(2)/2);
        
	%%%%%%%%%
    %% 3.2 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calibrating the eye-tracker %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        if USE_EYETRACKER
            eyetrackerIsCalibrated = false;
            xyFixation_px = [NaN NaN];

            % Display a blank page with fixation cross
            crsSetDrawPage(BlankPage);
            crsClearPage(BlankPage,CRS.BACKGROUND);
            % Draw a fixation point
            crsSetPen1(CRS.FIXATION);
            crsSetPen2(CRS.FIXATION);
            crsSetSpatialUnits(CRS.DEGREEUNIT);
            crsDrawOval(fix_xy_deg,[fix_diameter_deg,fix_diameter_deg]);
            crsSetSpatialUnits(CRS.PIXELUNIT);
            % ###
            crsSetFixationColour(  WhiteRGB);
            crsSetBackgroundColour(BackgroundRGB);

            % display page
            crsSetDisplayPage(BlankPage);
            crsPresent();

            % run
            fprintf('Running eyetracker calibration. Ask participant to direct gaze towards fixation cross. Press SPACE when ready\n');
            InH.keyboard.getInput()
            while ~eyetrackerIsCalibrated
                if any(InH.keyboard.getInput() == InH.keyboard.INPT_SPACE.code)
                    fprintf('Evaluating if eyes currently tracked...\n')
                    eyetracker.refresh(true); % false to suppress data logging
                    if logs.data.getN()<10
                        fprintf('FAILED: No eyetracking data found\n')
                    else
                        [xy, t] = logs.data.getLastKnownXY(10, true, false); % [useRaw, allowNan]
                        timeNow = GetSecs();
                        if (timeNow-t(1)) > 0.4 % must have 10 recent non-NaN data
                            fprintf('FAILED: No recent eyetracking data found (Current time: %1.2f; Last known time: %1.2f\n', timeNow, t(1));
                        else
                            xyFixation_px = mean(xy);
                            eyetrackerIsCalibrated = true;
                        end
                    end

                end
                WaitSecs(0.01);
            end
        else
            % display page
            crsSetDisplayPage(BlankPage);
            crsPresent();
        end

        
	%%%%%%%%%
    %% 4.1 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Run the experiment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

        % will primarily export data to block-by-block .csv files. But will
        % also save a matlab .mat file for completeness
        matExportData = struct('metaParams',metaParams, 'blockID',[], 'trialID',[], 'targRespCode',[], 'respCode',[], 'anscorrect',[], 'eccentricity_x_deg',[], 'eccentricity_y_deg',[], 'temporalFrequency_hz',[], 'Spatfreq_prcnt',[], 'respLatency_sec',[], 'currentGuess',[]);
            
        nBlocks = psychometricParams.nConditions;
        for blockID = 1:nBlocks
            
            % establish eyetracking log filename
            eyetrackLogFn = sprintf('%s-%i-%i-%i-%s', metaParams.expID, metaParams.partID, metaParams.sessID, getBlockNum(), datestr(now(),30));
                
            % set parameters for this condition. Contrast and orientation
            % will also be varied, on a trial-by-trial basis
            temporalFrequency_hz = stimParams.temporalfreq_hz(blockID);
            fprintf('====================================================\nCommencing new test block:\n > Temporal freq: %1.2f hz\n====================================================\n\n', temporalFrequency_hz);
            
            % Create a control object, P, for the psychophysical paradigm
            switch lower(psychometricParams.paradigm)
                case 'quest'
                    P = QUESTwrapper(psychometricParams.qparams);
                case 'staircase'
                    P = AdaptiveTrack(psychometricParams.scparams);
                case 'mcs'
                    P = ConstantStimuli(psychometricParams.mcsparams.levels, psychometricParams.mcsparams.nTrialsPerLevel);
                otherwise % defensive
                    error('Unknown paradigm: %s', psychometricParams.paradigm)
            end
            
            % Set stimulus modulation
            %=========================================================
            %nFrames = graphicParams.Fr/temporalFrequency_hz;
            %itime     = linspace(-1, 1,nFrames);
            %mplitude = (sin(pi*(itime)));
            amplitude = getPureTone(temporalFrequency_hz, graphicParams.Fr, psychometricParams.stimInt_sec , stimParams.cosRamp_sec);
            nFrames = length(amplitude);

        
            trialID = 0;
            while ~P.isFinished() % run psychometric block (e.g., 1 staircase)

                % init
                currentGuess = NaN;
                anscorrect = NaN;
                respLatency_sec = NaN;            
                
                % ---------------------------------------------------------
                % 1 Show pre-trial page
                % ---------------------------------------------------------
 
                
                % Display a blank page with fixation cross
                crsSetDrawPage(PrestimPage);
%                 if stimParams.useOverlay==2
%                     crsClearPage(PrestimPage,1);
%                 else
                    %crsClearPage(PrestimPage,CRS.BACKGROUND);
%                 end
                
                %=========================================================
                % DEFINE PIXEL LEVELS.
                %=========================================================
                      
%                 stimulusLutLevel = 1;
%                 crsPaletteSetPixelLevel(stimulusLutLevel, [0.5,0.5,0.5]);
%                 controlLutLevel=2;
%                 crsPaletteSetPixelLevel(controlLutLevel, [0.5,0.5,0.5]);

        Buff_new = zeros(256,3);
        Buff_new(AnnulusLevel,1) = AnnulusRGB(1);
        Buff_new(AnnulusLevel,2) = AnnulusRGB(2);
        Buff_new(AnnulusLevel,3) = AnnulusRGB(3);
        Buff_new(controlLutLevel,1) = controlRGB(1);
        Buff_new(controlLutLevel,2) = controlRGB(2);
        Buff_new(controlLutLevel,3) = controlRGB(3);
        Buff_new(4:253,1) = sinepal;
        Buff_new(4:253,2) = sinepal;
        Buff_new(4:253,3) = sinepal;
        Buff_new(BackgroundLevel,1) = MaskRGB(1);
        Buff_new(BackgroundLevel,2) = MaskRGB(2);
        Buff_new(BackgroundLevel,3) = MaskRGB(3);
        %Buff_new(1,1) = MaskRGB(1); % ensure background remain visible/invariant
        %Buff_new(1,2) = MaskRGB(2);
        %Buff_new(1,3) = MaskRGB(3);
        %                 Buff_new(2,1) = WhiteRGB(1); % ensure background remain visible/invariant
        %                 Buff_new(2,2) = WhiteRGB(2);
        %                 Buff_new(2,3) = WhiteRGB(3);
        Buff_new(CRS.BACKGROUND,1) = BackgroundRGB(1); % ensure background remain visible/invariant
        Buff_new(CRS.BACKGROUND,2) = BackgroundRGB(2);
        Buff_new(CRS.BACKGROUND,3) = BackgroundRGB(3);
        %                 Buff_new(3,1) = AnnulusContrastRGB(1); % ensure annulus remains visible/invariant
        %                 Buff_new(3,2) = AnnulusContrastRGB(2);
        %                 Buff_new(3,3) = AnnulusContrastRGB(3);
        Buff_new(CRS.FIXATION,1) = FixRGB(1); % ensure fixation marker remains visible/invariant
        Buff_new(CRS.FIXATION,2) = FixRGB(2);
        Buff_new(CRS.FIXATION,3) = FixRGB(3);
                
                crsLUTBUFFERWrite(1,Buff_new);
                crsSetSpatialUnits(CRS.PIXELUNIT);
                % Create the LUT Object and set parameters
                screenDims_px = crsGetScreenSizePixels;
                
      
                
                %Extract stimulus parameters
                eccentricity_x_deg = round(crsDegreesToPixels(stimParams.eccentricity_x_deg));
                eccentricity_y_deg = round(crsDegreesToPixels(stimParams.eccentricity_y_deg));
                stim_deg = round(crsDegreesToPixels((stimParams.diameter_deg./2)));
                appeture_x_deg = round(crsDegreesToPixels(stimParams.apperture_x_deg));
                appeture_y_deg = round(crsDegreesToPixels(stimParams.apperture_y_deg));
                appeture_deg = round(crsDegreesToPixels((stimParams.apperture_deg./2)));
                
                
                %Select which side to display stimulus on
                
                palletgabor =ones(screenDims_px(2),screenDims_px(1))*BackgroundLevel;
                %make sine wave
 
                %draw annulus
                palletgabor = padarray(palletgabor,[1 1],BackgroundLevel,'pre');
                c_mask = zeros(size(palletgabor))*BackgroundLevel; 
                c_mask2 = zeros(size(palletgabor))*BackgroundLevel; 
                for side = 1:2
                    middle_px = screenDims_px./2;
                    %if side == 1
                    % middle_px = middle_px+1;
                    %end
                    cx = appeture_x_deg(side)+middle_px(1);
                    cy= appeture_y_deg(side)+middle_px(2);
                    %     ix=screenDims_px(1)+50;
                    %     iy=screenDims_px(2)+50;
                    ix=screenDims_px(1);
                    iy=screenDims_px(2);
                    %f = 0.
                    r=appeture_deg;
                    %Plot the annulus ring
                    [x,y] = meshgrid(-cx:(ix-cx),-cy:(iy-cy));
                    c_mask = c_mask + ((x.^2 + y.^2) <= r^2);
                    
                    
                    %Plot the flashstim ring
                    %                      cxf = round(crsDegreesToPixels(stimParams.eccentricity_x_deg(1)))+middlepx(1);
                    %                      cyf = round(crsDegreesToPixels(stimParams.eccentricity_y_deg(1)))+middlepx(2);
                    %                      flashr = round(crsDegreesToPixels((stimParams.diameter_deg./2)));
                    xc=eccentricity_x_deg(side)+middle_px(1);
                    yc=eccentricity_y_deg(side) +middle_px(2);
                    [x,y] = meshgrid(-xc:(ix-xc),-yc:(iy-yc));
                    c_mask2 = c_mask2 + ((x.^2 + y.^2) <= stim_deg^2);
                end
                palletgabor(~c_mask) = BackgroundLevel;
                palletgabor(logical(c_mask - c_mask2)) = AnnulusLevel;
                
                middle_px = screenDims_px./2;
                % if targRespCode == 1
                %     middle_px = middle_px+1;
                % end
                for controlside=1:2
                    xc=eccentricity_x_deg(controlside)+middle_px(1);
                    yc=eccentricity_y_deg(controlside) +middle_px(2);
                    [x,y] = meshgrid(-xc:(ix-xc),-yc:(iy-yc));
                    controlmask = ((x.^2 + y.^2) <= stim_deg^2);
                    palletgabor(controlmask==1) = controlLutLevel;
                end
                
                %palletgabor = palletgabor(2:end,2:end);
                palletgabor = padarray(palletgabor,[padval,padval],BackgroundLevel,'both');
                
                
                
                
                crsSetDrawPage(PrestimPage);
                crsDrawMatrixPalettised(0,0,palletgabor);   
                
                % Draw a fixation point
                crsSetPen1(CRS.FIXATION);
                crsSetPen2(CRS.FIXATION);
                crsSetSpatialUnits(CRS.DEGREEUNIT);
                crsDrawOval(fix_xy_deg,[fix_diameter_deg,fix_diameter_deg]);
                crsSetSpatialUnits(CRS.PIXELUNIT);

                % display page
                crsSetDisplayPage(PrestimPage);
                crsPresent();
                tmp_sec = GetSecs();


                % ---------------------------------------------------------
                % 2 Set Parameters
                % --------------------------------------------------------- 
                trialID = trialID + 1;
                targRespCode = Randi(2);
                %eccentricity_x_deg = stimParams.eccentricity_x_deg(targRespCode);
                %eccentricity_y_deg = stimParams.eccentricity_y_deg(targRespCode);
                
                % Set spatial frequency
                % ========
                % in cycles per pixel
                %Spatfreq_prcnt = 99; % Percent of maximum range
                Spatfreq_raw = P.getDelta();

                % map to log-spaced value if so requested
                if psychometricParams.useLogScale
                    %rng_linscale = psychometricParams.qparams.deltaMin:psychometricParams.qparams.beta:psychometricParams.qparams.deltaMax
                    %rng_logscale = logspace(log10(psychometricParams.qparams.deltaMin), log10(psychometricParams.qparams.deltaMax), length(rng_linscale))
                    lin_min = psychometricParams.qparams.deltaMin;
                    lin_max = psychometricParams.qparams.deltaMax;
                    log_min = log10(lin_min);
                    log_max = log10(lin_max);
                    %
                    targ_lin = Spatfreq_raw;
                    targ_lin_prcnt = (targ_lin-lin_min)/(lin_max - lin_min);
                    targ_log = 10^((log_max - log_min)*targ_lin_prcnt);
                    %
                    Spatfreq_prcnt = targ_log * psychometricParams.qparams.deltaMin;
                else
                   Spatfreq_prcnt = Spatfreq_raw;
                end

                % ---------------------------------------------------------
                % 3 Create stimuli
                % ---------------------------------------------------------
                
                FixRGB=[0 0 0];
                
                % set contrast
                cone_modulation = [0 0 0];
                cone_modulation(coneIdx) = Contrast_prcnt/100;
                
                %Select which side to display stimulus on
                side = (targRespCode);
                
                % report
                fprintf('Block %i, Trial %i, Eccentricity %1.2f, SpatFreq=%1.10f, cone_modulation=[%1.2f %1.2f %1.2f]  [targRespCode=%i]...', blockID, trialID, eccentricity_x_deg(side), Spatfreq_prcnt, cone_modulation, targRespCode);
                %fprintf(['\nContrast is: ', num2str(Contrast_prcnt), ' percent'])
                %fprintf(['\nEccentricity is: ', num2str(eccentricity_x_deg(side))])
                %=========================================================
                % Write LUT.
                %=========================================================
         Buff_new = zeros(256,3);
        Buff_new(AnnulusLevel,1) = AnnulusRGB(1);
        Buff_new(AnnulusLevel,2) = AnnulusRGB(2);
        Buff_new(AnnulusLevel,3) = AnnulusRGB(3);
        Buff_new(controlLutLevel,1) = controlRGB(1);
        Buff_new(controlLutLevel,2) = controlRGB(2);
        Buff_new(controlLutLevel,3) = controlRGB(3);
        Buff_new(4:253,1) = sinepal;
        Buff_new(4:253,2) = sinepal;
        Buff_new(4:253,3) = sinepal;
        Buff_new(BackgroundLevel,1) = MaskRGB(1);
        Buff_new(BackgroundLevel,2) = MaskRGB(2);
        Buff_new(BackgroundLevel,3) = MaskRGB(3);
        %Buff_new(1,1) = MaskRGB(1); % ensure background remain visible/invariant
        %Buff_new(1,2) = MaskRGB(2);
        %Buff_new(1,3) = MaskRGB(3);
        %                 Buff_new(2,1) = WhiteRGB(1); % ensure background remain visible/invariant
        %                 Buff_new(2,2) = WhiteRGB(2);
        %                 Buff_new(2,3) = WhiteRGB(3);
        Buff_new(CRS.BACKGROUND,1) = BackgroundRGB(1); % ensure background remain visible/invariant
        Buff_new(CRS.BACKGROUND,2) = BackgroundRGB(2);
        Buff_new(CRS.BACKGROUND,3) = BackgroundRGB(3);
        %                 Buff_new(3,1) = AnnulusContrastRGB(1); % ensure annulus remains visible/invariant
        %                 Buff_new(3,2) = AnnulusContrastRGB(2);
        %                 Buff_new(3,3) = AnnulusContrastRGB(3);
        Buff_new(CRS.FIXATION,1) = FixRGB(1); % ensure fixation marker remains visible/invariant
        Buff_new(CRS.FIXATION,2) = FixRGB(2);
        Buff_new(CRS.FIXATION,3) = FixRGB(3); % ensure annulus remains visible/invariant
                %Buff_new(6,2) = AnnulusContrast_One_RGB(2);
                %Buff_new(6,3) = AnnulusContrast_One_RGB(3);
                %Buff_new(7,1) = AnnulusContrast_Two_RGB(1); % ensure annulus remains visible/invariant
                %Buff_new(7,2) = AnnulusContrast_Two_RGB(2);
                %Buff_new(7,3) = AnnulusContrast_Two_RGB(3);
                crsLUTBUFFERWrite(1,Buff_new);
                screenDims_px = crsGetScreenSizePixels;
           
                
                crsSetSpatialUnits(CRS.PIXELUNIT);
                
                %Extract stimulus parameters
                eccentricity_x_deg = round(crsDegreesToPixels(stimParams.eccentricity_x_deg));
                eccentricity_y_deg = round(crsDegreesToPixels(stimParams.eccentricity_y_deg));
                stim_deg = round(crsDegreesToPixels((stimParams.diameter_deg./2)));
                appeture_x_deg = round(crsDegreesToPixels(stimParams.apperture_x_deg));
                appeture_y_deg = round(crsDegreesToPixels(stimParams.apperture_y_deg));
                appeture_deg = round(crsDegreesToPixels((stimParams.apperture_deg./2)));
                
                %Set section of pallet to use for gabor
                MinValue = 4;
                MaxValue = 253;
   
                f = Spatfreq_prcnt;
                
                t = -stim_deg:stim_deg;
                wv = (sin(2*pi*f*t)) .* (Contrast_prcnt/100) ;
                sinewave = repmat(wv,[length(t),1]);
                
                [X,Y] = meshgrid(-stim_deg:stim_deg,-stim_deg:stim_deg);
                sigma = stim_deg./0.1; %was 2
                exponent = (((X-(median(X(1,:)))).^2)+((Y-median(Y(:,1))).^2))/(2*sigma.^2);
                amp = 1;
                gauss = amp *exp(-exponent);
                
                gabor = sinewave.*gauss;
                gabor = (gabor +1)*midpoint;
              
                
                %Put the gabor into the masks
                cx=0;
                cy=0;
                ix=screenDims_px(1);
                iy=screenDims_px(2);
                %f = 0.
                r=stim_deg;
                
                %Plot the annulus ring
                [x,y] = meshgrid(-cx:(ix-cx),-cy:(iy-cy));
                stim_mask = ((X.^2 + Y.^2) <= r^2);
                
                
                if round(mean(gabor(stim_mask)),10) ~= 0.5
                    warning(['The stimulus does not average to grey (0.5)! Instead, it averages to ',...
                        num2str(mean(palletgabor(stim_mask))),...
                        '. This might not be a fair comparison between stimulus and control!']);
                    
                end
                
                %Make palletised gabor wave
                palletgabor = (gabor.*(MaxValue-MinValue))+MinValue;
                
                palletgabor(~stim_mask) = BackgroundLevel;
                
                stimscreen = ones(screenDims_px(2)+1,screenDims_px(1)+1)*BackgroundLevel;
                c_mask = ones(size(stimscreen))*BackgroundLevel;
                c_mask2 = ones(size(stimscreen))*BackgroundLevel;
                for side = 1:2
                    middle_px = screenDims_px./2;
                    %if side == 1
                    % middle_px = middle_px+1;
                    %end
                    cx = appeture_x_deg(side)+middle_px(1);
                    cy= appeture_y_deg(side)+middle_px(2);
                    %     ix=screenDims_px(1)+50;
                    %     iy=screenDims_px(2)+50;
                    ix=screenDims_px(1);
                    iy=screenDims_px(2);
                    %f = 0.
                    r=appeture_deg;
                    %Plot the annulus ring
                    [x,y] = meshgrid(-cx:(ix-cx),-cy:(iy-cy));
                    c_mask = c_mask + ((x.^2 + y.^2) <= r^2);
                    
                    
                    %Plot the flashstim ring
                    %                      cxf = round(crsDegreesToPixels(stimParams.eccentricity_x_deg(1)))+middlepx(1);
                    %                      cyf = round(crsDegreesToPixels(stimParams.eccentricity_y_deg(1)))+middlepx(2);
                    %                      flashr = round(crsDegreesToPixels((stimParams.diameter_deg./2)));
                    xc=eccentricity_x_deg(side)+middle_px(1);
                    yc=eccentricity_y_deg(side) +middle_px(2);
                    [x,y] = meshgrid(-xc:(ix-xc),-yc:(iy-yc));
                    c_mask2 = c_mask2 + ((x.^2 + y.^2) <= stim_deg^2);
                end
                stimscreen(~c_mask) = BackgroundLevel;
                stimscreen(logical(c_mask - c_mask2)) = AnnulusLevel;
                
                if targRespCode == 1
                    controlside = 2;
                else
                    controlside =1;
                end
                % Fill in the control side
                middle_px = screenDims_px./2;
                xc=eccentricity_x_deg(controlside)+middle_px(1);
                yc=eccentricity_y_deg(controlside) +middle_px(2);
                [x,y] = meshgrid(-xc:(ix-xc),-yc:(iy-yc));
                controlmask = ((x.^2 + y.^2) <= stim_deg^2);
                stimscreen(controlmask==1) = controlLutLevel;
                
                % Fill in the stimulus size
                middle_px = screenDims_px./2;
                xc=eccentricity_x_deg(targRespCode)+middle_px(1);
                yc=eccentricity_y_deg(targRespCode) +middle_px(2);
                [x,y] = meshgrid(-xc:(ix-xc),-yc:(iy-yc));
                
                targmask = ((x.^2 + y.^2) <= stim_deg^2);
                stimscreen(targmask==1) = palletgabor(stim_mask);
%                 stimscreen(targmask==1) = controlLutLevel;
%                  [x,y] = meshgrid(-cx:(ix-cx),-cy:(iy-cy));
%                 stim_mask = ((X.^2 + Y.^2) <= r^2);

                palletgabor = stimscreen;
 
                palletgabor = padarray(palletgabor,[padval,padval],BackgroundLevel,'both');
                
                crsSetDrawPage(StimulusPage);
                disp(['size stim: ',num2str(size(palletgabor))])
                crsDrawMatrixPalettised(0,0,palletgabor);
                
                % Draw a fixation point
                crsSetPen1(CRS.FIXATION);
                crsSetPen2(CRS.FIXATION);
                
                crsSetSpatialUnits(CRS.DEGREEUNIT);
                crsDrawOval(fix_xy_deg,[fix_diameter_deg,fix_diameter_deg]);
                crsSetSpatialUnits(CRS.PIXELUNIT);
                
                % Set the system colours
                crsSetFixationColour(  FixRGB);
                %                 crsSetBackgroundColour(BlackRGB);
                
                
                % ---------------------------------------------------------
                % 4 Pause for pre-stimulus interval - wait for eyes if necessary
                % ---------------------------------------------------------
                WaitSecs(psychometricParams.preInt_sec - (GetSecs()-tmp_sec)); % don't include time already spent doing the above
                
                % ---------------------------------------------------------
                % 5 Present stimuli - abort trial if gaze deviates from
                % fixation spot
                % ---------------------------------------------------------
                % The system colours and Object parameters will
                % change in synchrony when crsPresent is called.
                % The drawing and display pages are also swapped
                % over at this point.
                crsSetDisplayPage(StimulusPage)
                crsPresent();
                
% crsResponseBoxBuzzer(CRS.respSEC05,CRS.respTONE1);                
%                 beep(); % give auditory cue that the signal/noise has appeared


                nOutliers = 0;
                nMissingSamples = 0;
                isTrialAborted = false;
                
                respCode = NaN;
             	InH.getInput(); % flush
                if USE_EYETRACKER
                    eyetracker.refresh(false); % flush
                end
                t = GetSecs();
                while (GetSecs()-t)<psychometricParams.stimInt_sec
                    % check for eye-gaze deviation
                    if USE_EYETRACKER
                        n = eyetracker.refresh(true); % false to supress logging
                        if n > 0 % Update classifier
                            % check for deviant samples (e.g., not fixating fixation cue)
                            gaze_xy_px = logs.data.getSinceT(t, 1:2);
                            deviation_px = sqrt(sum(bsxfun(@minus, gaze_xy_px, xyFixation_px).^2,2));
                            nOutliers = sum(deviation_px>eyetrackParams.gaze_outlierCriterion_px);
                            if nOutliers > eyetrackParams.gaze_nOutliersPermitted
                                crsResponseBoxBuzzer(CRS.respSEC05,CRS.respTONE1);
                                fprintf('Gaze not on target\n');
                                isTrialAborted = true;
                                break
                            end
                            % check for missing samples (e.g., closed eyes or turned away from the screen)
                            nMissingSamples = sum(any(isnan(deviation_px),2));
                            if nMissingSamples > eyetrackParams.gaze_nMissingSamplesPermitted
                                crsResponseBoxBuzzer(CRS.respSEC05,CRS.respTONE2);
                                fprintf('Eyes not detected\n');
                                isTrialAborted = true;
                                break
                            end
                        end
                    end
                    % check for user input (if allowed to abort
                    % stimulus prematurely)
                    if psychometricParams.allowRespBeforeEndOfStim
                        usrInput = InH.getInput();
                        usrInput = usrInput(1);
                        if ismember(usrInput, [InH.INPT_LEFTARROW InH.INPT_RIGHTARROW])
                            respCode = usrInput;
                            break;
                        end
                    end
                    
                    WaitSecs(0.001);
                end
               	endOfStimTime_sec = GetSecs();
                %crsSetCommand(CRS.CYCLELUTDISABLE); % stop LUT cycling
                
                
                % #############
                crsSetDrawPage(BlahPage);
%                 if stimParams.useOverlay>=1
%                     crsClearPage(  BlahPage,1);
%                 else
                    %crsClearPage(  BlahPage);
%                 end
             crsSetSpatialUnits(CRS.PIXELUNIT);
                % Create the LUT Object and set parameters
                screenDims_px = crsGetScreenSizePixels;
             
                
                %Extract stimulus parameters
                eccentricity_x_deg = round(crsDegreesToPixels(stimParams.eccentricity_x_deg));
                eccentricity_y_deg = round(crsDegreesToPixels(stimParams.eccentricity_y_deg));
                stim_deg = round(crsDegreesToPixels((stimParams.diameter_deg./2)));
                appeture_x_deg = round(crsDegreesToPixels(stimParams.apperture_x_deg));
                appeture_y_deg = round(crsDegreesToPixels(stimParams.apperture_y_deg));
                appeture_deg = round(crsDegreesToPixels((stimParams.apperture_deg./2)));
                
                
                %Select which side to display stimulus on
                
                palletgabor = ones(screenDims_px(2),screenDims_px(1))*BackgroundLevel;
                %make sine wave
 
                %draw annulus
                palletgabor = padarray(palletgabor,[1 1],BackgroundLevel,'pre');
                c_mask = zeros(size(palletgabor))*BackgroundLevel;
                c_mask2 = zeros(size(palletgabor))*BackgroundLevel; 
                for side = 1:2
                    middle_px = screenDims_px./2;
                    %if side == 1
                    % middle_px = middle_px+1;
                    %end
                    cx = appeture_x_deg(side)+middle_px(1);
                    cy= appeture_y_deg(side)+middle_px(2);
                    %     ix=screenDims_px(1)+50;
                    %     iy=screenDims_px(2)+50;
                    ix=screenDims_px(1);
                    iy=screenDims_px(2);
                    %f = 0.
                    r=appeture_deg;
                    %Plot the annulus ring
                    [x,y] = meshgrid(-cx:(ix-cx),-cy:(iy-cy));
                    c_mask = c_mask + ((x.^2 + y.^2) <= r^2);
                    
                    
                    %Plot the flashstim ring
                    %                      cxf = round(crsDegreesToPixels(stimParams.eccentricity_x_deg(1)))+middlepx(1);
                    %                      cyf = round(crsDegreesToPixels(stimParams.eccentricity_y_deg(1)))+middlepx(2);
                    %                      flashr = round(crsDegreesToPixels((stimParams.diameter_deg./2)));
                    xc=eccentricity_x_deg(side)+middle_px(1);
                    yc=eccentricity_y_deg(side) +middle_px(2);
                    [x,y] = meshgrid(-xc:(ix-xc),-yc:(iy-yc));
                    c_mask2 = c_mask2 + ((x.^2 + y.^2) <= stim_deg.^2);
                end
                palletgabor(~c_mask) = BackgroundLevel;
                palletgabor(logical(c_mask - c_mask2)) = AnnulusLevel;
                
                middle_px = screenDims_px./2;
                % if targRespCode == 1
                %     middle_px = middle_px+1;
                % end
                for controlside=1:2
                    xc=eccentricity_x_deg(controlside)+middle_px(1);
                    yc=eccentricity_y_deg(controlside) +middle_px(2);
                    [x,y] = meshgrid(-xc:(ix-xc),-yc:(iy-yc));
                    controlmask = ((x.^2 + y.^2) <= stim_deg.^2);
                    palletgabor(controlmask==1) = controlLutLevel;
                end
                
                palletgabor = padarray(palletgabor,[padval,padval],BackgroundLevel,'both');
                %
                
                
                
                
                
                crsSetDrawPage(BlahPage);
                crsDrawMatrixPalettised(0,0,palletgabor);   
                
                % Draw a fixation point
                crsSetPen1(CRS.FIXATION);
                crsSetPen2(CRS.FIXATION);
                crsSetSpatialUnits(CRS.DEGREEUNIT);
                crsDrawOval(fix_xy_deg,[fix_diameter_deg,fix_diameter_deg]);
                crsSetSpatialUnits(CRS.PIXELUNIT);
                crsSetDisplayPage(BlahPage);
                

                % Draw a fixation point
%                 crsSetPen1(CRS.FIXATION);
%                 crsSetPen2(CRS.FIXATION);
%                 crsDrawOval(fix_xy_deg,[fix_diameter_deg,fix_diameter_deg]);
%                 crsSetDisplayPage(BlankPage);
                
                % clean up memory
                crsObjDestroyAll();
                
                if ~isTrialAborted
                    % ---------------------------------------------------------
                    % 6 Get & evaluate response
                    % ---------------------------------------------------------
                    if isnan(respCode)
                        respCode = InH.waitForInput([InH.INPT_LEFTARROW InH.INPT_RIGHTARROW]); %1 for left-arrow, 2 for right-arrow
                    end
                    respLatency_sec = GetSecs() - endOfStimTime_sec;
                    anscorrect = targRespCode==respCode;
                    
                    % ---------------------------------------------------------
                    % 7 Update psychophysical params
                    % ---------------------------------------------------------
                    P.update(anscorrect);
                end
                
                % ---------------------------------------------------------
                % 8 Record trial data
                % ---------------------------------------------------------
                if strcmpi(psychometricParams.paradigm, 'quest')
                    currentGuess = P.estThresh();
                    % map to log-spaced value if so requested
                    if psychometricParams.useLogScale
                        targ_lin = currentGuess;
                        targ_lin_prcnt = (targ_lin-lin_min)/(lin_max - lin_min);
                        targ_log = 10^((log_max - log_min)*targ_lin_prcnt);
                        currentGuess = targ_log;
                    end
                end
                
                eccentricity_x_deg_side = eccentricity_x_deg(targRespCode);
                eccentricity_y_deg_side = eccentricity_y_deg(targRespCode);
                
                writeData(blockID, trialID, targRespCode, respCode, anscorrect, eccentricity_x_deg_side, eccentricity_y_deg_side, temporalFrequency_hz, Spatfreq_prcnt, respLatency_sec, currentGuess, nOutliers, nMissingSamples, isTrialAborted, eyetrackLogFn);

                matExportData.blockID(end+1)                    = blockID;
                matExportData.trialID(end+1)                    = trialID;
                matExportData.targRespCode(end+1)               = targRespCode;
                matExportData.respCode(end+1)                   = respCode;
                matExportData.anscorrect(end+1)                 = anscorrect;
                matExportData.eccentricity_x_deg(end+1)        	= eccentricity_x_deg_side;
                matExportData.eccentricity_y_deg(end+1)      	= eccentricity_y_deg_side;
                matExportData.temporalFrequency_hz(end+1)       = temporalFrequency_hz;
                matExportData.Spatfreq_prcnt(end+1)          	= Spatfreq_prcnt;
                matExportData.respLatency_sec(end+1)            = respLatency_sec;
                matExportData.currentGuess(end+1)               = currentGuess;
                matExportData.eyetrackLogFn                     = eyetrackLogFn;
                
                % ---------------------------------------------------------
                % 9 Provide feedback
                % ---------------------------------------------------------
                fprintf('  anscorrect=%i\n', anscorrect);
% crsResponseBoxBuzzer(CRS.respSEC05,CRS.respTONE3);                
% warning('Provide feedback')                 
    
                if isTrialAborted
                    crsSetFixationColour(WhiteRGB);   
                    FixRGB = WhiteRGB;
                elseif anscorrect
                    crsSetFixationColour(  BlackRGB);
                    FixRGB = BlackRGB;
                    
                else
                    crsSetFixationColour(  BlackRGB);
                    FixRGB = BlackRGB;
                    
                end
               	
                % init page
                crsSetDrawPage(PostStimPage);
%                 if stimParams.useOverlay>=1
                   % crsClearPage(  PostStimPage,1);
%                 else
                   % crsClearPage(  PostStimPage);
%                 end
                
                
                      crsSetSpatialUnits(CRS.PIXELUNIT);
                % Create the LUT Object and set parameters
                screenDims_px = crsGetScreenSizePixels;
                
          
                %Extract stimulus parameters
                eccentricity_x_deg = round(crsDegreesToPixels(stimParams.eccentricity_x_deg));
                eccentricity_y_deg = round(crsDegreesToPixels(stimParams.eccentricity_y_deg));
                stim_deg = round(crsDegreesToPixels((stimParams.diameter_deg./2)));
                appeture_x_deg = round(crsDegreesToPixels(stimParams.apperture_x_deg));
                appeture_y_deg = round(crsDegreesToPixels(stimParams.apperture_y_deg));
                appeture_deg = round(crsDegreesToPixels((stimParams.apperture_deg./2)));
                
                
                %Select which side to display stimulus on
                
                palletgabor = ones(screenDims_px(2),screenDims_px(1))*BackgroundLevel;
                %make sine wave
 
                %draw annulus
                palletgabor = padarray(palletgabor,[1 1],BackgroundLevel,'pre');
                c_mask = zeros(size(palletgabor))*BackgroundLevel;
                c_mask2 = zeros(size(palletgabor))*BackgroundLevel; 
                for side = 1:2
                    middle_px = screenDims_px./2;
                    %if side == 1
                    % middle_px = middle_px+1;
                    %end
                    cx = appeture_x_deg(side)+middle_px(1);
                    cy= appeture_y_deg(side)+middle_px(2);
                    %     ix=screenDims_px(1)+50;
                    %     iy=screenDims_px(2)+50;
                    ix=screenDims_px(1);
                    iy=screenDims_px(2);
                    %f = 0.
                    r=appeture_deg;
                    %Plot the annulus ring
                    [x,y] = meshgrid(-cx:(ix-cx),-cy:(iy-cy));
                    c_mask = c_mask + ((x.^2 + y.^2) <= r^2);
                    
                    
                    %Plot the flashstim ring
                    %                      cxf = round(crsDegreesToPixels(stimParams.eccentricity_x_deg(1)))+middlepx(1);
                    %                      cyf = round(crsDegreesToPixels(stimParams.eccentricity_y_deg(1)))+middlepx(2);
                    %                      flashr = round(crsDegreesToPixels((stimParams.diameter_deg./2)));
                    xc=eccentricity_x_deg(side)+middle_px(1);
                    yc=eccentricity_y_deg(side) +middle_px(2);
                    [x,y] = meshgrid(-xc:(ix-xc),-yc:(iy-yc));
                    c_mask2 = c_mask2 + ((x.^2 + y.^2) <= stim_deg^2);
                end
                palletgabor(~c_mask) = BackgroundLevel;
                palletgabor(logical(c_mask - c_mask2)) = AnnulusLevel;
                
                middle_px = screenDims_px./2;
                % if targRespCode == 1
                %     middle_px = middle_px+1;
                % end
                for controlside=1:2
                    xc=eccentricity_x_deg(controlside)+middle_px(1);
                    yc=eccentricity_y_deg(controlside) +middle_px(2);
                    [x,y] = meshgrid(-xc:(ix-xc),-yc:(iy-yc));
                    controlmask = ((x.^2 + y.^2) <= stim_deg^2);
                    palletgabor(controlmask==1) = controlLutLevel;
                end
                
                palletgabor = padarray(palletgabor,[padval,padval],BackgroundLevel,'both');
                %
                
                
                
                
                
                crsSetDrawPage(PostStimPage);
                %disp(['size 2: ',num2str(size(palletgabor))])
                crsDrawMatrixPalettised(0,0,palletgabor);  
                
                
   
                
                % Draw a fixation point
                crsSetPen1(CRS.FIXATION);
                crsSetPen2(CRS.FIXATION);                
                crsSetSpatialUnits(CRS.DEGREEUNIT);
                crsDrawOval(fix_xy_deg,[fix_diameter_deg,fix_diameter_deg]);
                crsSetSpatialUnits(CRS.PIXELUNIT);
                % display page
                crsSetDisplayPage(PostStimPage);
                crsPresent();
                
                if ~isTrialAborted
                    WaitSecs(psychometricParams.fbackInt_sec);
                else
                   warning('Wait until fixating again') 
                    WaitSecs(1);
                end
                % reset
                crsSetFixationColour(FixRGB);
                crsObjDestroyAll();
            end % end-of-trial

            % end of block - save eye-tracking data and flush log
            if USE_EYETRACKER
                logs.data.getInstance().save(eyetrackLogFn, false);
            end
            
            if blockID < nBlocks
                % manually tell PTR to start a new block
                newBlock();
                
                % compute score
                switch lower(psychometricParams.paradigm)
                    case 'quest'
                        scoretxt = sprintf('%1.1f percent', 100 - P.estThresh());
                    case 'staircase'
                        error('functionality not yet written');
                    case 'mcs'
                        [pc, trialval, anscorrect] = P.getPC();
                        scoretxt = [sprintf('%1.0f, ', P.getPC()*100) 'percent'];
                    otherwise % defensive
                        error('Unknown paradigm: %s', psychometricParams.paradigm)
                end
                
                % show block feedback
                txt = sprintf('Well done. Block %i of %i complete. Your score was: %s', blockID, nBlocks, scoretxt);
                crsSetSpatialUnits(CRS.PIXELUNIT);
                BackLevel = MinValue; % CRS.BACKGROUND;
                StimLevel = MaxValue;
                crsSetDrawPage(BlankPage);
                crsClearPage(BlankPage,BackLevel);
                crsSetDisplayPage(BlankPage);
                crsSetDrawPage(StimulusPage);
                crsClearPage(StimulusPage,BackLevel);
                crsSetDrawOrigin([100,100]);
                MinPixelLevel  = BackLevel;
                NumPixelLevels = (StimLevel - MinPixelLevel) + 1;
                crsPaletteSetSection(MinPixelLevel,linspace(BlackRGB(1),WhiteRGB(1),NumPixelLevels),...
                    linspace(BlackRGB(2),WhiteRGB(2),NumPixelLevels),...
                    linspace(BlackRGB(3),WhiteRGB(3),NumPixelLevels));
                oldMode = crsGetDrawMode;
                Mode = oldMode;
                Mode.AntiAliasing = true;
                crsSetDrawMode(Mode);
                LineHeight  = 28;
                LineSpacing = 8;
                StringSize = [10,LineHeight-LineSpacing];
                Horizontal = CRS.ALIGNLEFTTEXT;
                Vertical   = CRS.ALIGNBOTTOMTEXT;
                Angle      = 0;
                Flags      = CRS.FONTNORMAL;
                crsSetStringMode(StringSize,Horizontal,Vertical,Angle,Flags);
                crsSetTrueTypeFont('Arial');
                crsSetPen1(StimLevel);
                crsSetPen2(BackLevel);
                crsDrawString([0,1 * LineHeight],txt);
                crsSetDisplayPage(StimulusPage);
                InH.waitForInput([InH.INPT_MIDDLE InH.INPT_LEFTARROW]) % wait for appropriate input before continuing
                crsSetDrawMode(oldMode);
                crsSetSpatialUnits(CRS.DEGREEUNIT);
                crsSetDrawOrigin(ScreenSize(1)/2, ScreenSize(2)/2);
        
                % take a break if one is due
                if ismember(blockID, psychometricParams.breakAfterNBlocks)
                    fprintf('\n\n<<<<<<<<<<<< BREAK >>>>>>>>>>>>\n\n');
                    error('Take-a-break functionality not yet written');
                end
                
            end

        end % end-of-block

        % save matlab variable
        fn = fullfile('..', '..', 'data', '__matlabdata', sprintf('%s_p%i_s%i_%s.mat', metaParams.expID, metaParams.partID, metaParams.sessID, datestr(now(),30)));
        save(fn, 'matExportData');
        
        % test complete
        fprintf('Testing complete!')
        
        
	%%%%%%%%%
    %% 5.1 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Debrief %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
warning('add debrief screen');    
    
    
	%%%%%%%%%
    %% 6.1 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Finish up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SingletonManager.clearAll();
        IvMain.finishUp();
        
    catch ME
        % ensure graceful exit
        clear CRS
        SingletonManager.clearAll();
        IvMain.finishUp();
        rethrow(ME);
    end