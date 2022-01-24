%  function demoSilentSubstitutionDisplay
% DEMOSILENTSUBSTITUTIONDISPLAY demonstrates how to modulate a flickering light
% using silent substitution. 
%------------------------------------------------------------------------------
% demoSilentSubstitutionDisplay
% ======================
% This demo illustrates how to modulate a flickering light
% using silent substitution. "Silent Substitution" technique 
% was originally introduced by M. Ishihara at the begining of the century.
% Though the method is known mostly through the work of W.A.H. Rushton et al. (1973).
% [See also  Estevez and Spekreijse, 1974, 1982].
% The technique is a method of isolating cone responses. For example, 
% consider the case in which we wish to isolate L-cones responses. 
% In order to do so we need to provide constant excitations for the S- and
% M-cones. 
% In this demo, a central disk appears in the centre of the screen. 
% The disk is temporally modulated at 1 Hz. The stimulus modulation 
% generates constant M- and S-cone responsess but different L-cone
% responses, i.e. L-cone responses have been isolated.
% The demo uses the function <ctGetColourTrival> to transform LMS into RGB.
% The function <ctGetColourTrival> requires in input the 
% Spectral Power Distributions of the phosphors of the display device and a set of cone sensitivities. 
% We provide a standard set of SPDs measured for a Sony Trinitron and a set
% of cone sensitivities [Stockman and Sharpe (2000)].
% Please note that this function requires a ViSaGe or VSG system to display
% the stimulus.

% 29-Jan -2007   Caterina Ripamonti (c.ripamonti@ucl.ac.uk) Wrote it for Cambridge Research Systems Ltd.
% 17-Sept-2010   CR                 Use Phosphors_Barco_SpectraCal3807801 instead of Sony.
% 27-Nov-2012    CR                 Edit description.

fprintf('\ninitialising VSG ...\n')
global CRS;
vsgInit;


%=========================================================
% Set pages.
%=========================================================
blankPage    = 1;
stimulusPage = 2;
BackgroundRGB = [0.5 0.5 0.5];
crsPaletteSetPixelLevel(CRS.BACKGROUND,BackgroundRGB);

% Display a blank page while we work.
crsSetDrawPage   (blankPage);
crsClearPage     (blankPage,CRS.BACKGROUND);
crsSetDisplayPage(blankPage);
  
% Draw 
crsSetDrawPage(stimulusPage);
crsClearPage(  stimulusPage,CRS.BACKGROUND);

%=========================================================
% Set screen resolution.
%=========================================================
theScreen.height = crsGetScreenHeightPixels;
theScreen.width  = crsGetScreenWidthPixels;

%=========================================================
% LMS
Sensors = 'ConeSensitivities_SS_2degELin3908301.mat';

%=========================================================
% Device SPD
deviceSPD='Phosphors_Barco_SpectraCal3807801.mat';

%=========================================================
% Set temporal frequency and define number of frames per cycle.
%=========================================================
frameRate       = crsGetFrameRate;
frequency       = 0.3; 
nFramesPerCycle = frameRate/frequency;

%=========================================================
% Set cone modulation.
%=========================================================
% Define which cone will be modulated. 
% The order is: cone_modulation = [L M S];
% The value different from 0 indicates which cone-class will be isolated. 
% The value equal to 0 indicates which cone-classes will generate a constant response. 
% In this demo, L-cones will be isolated.
% The modulation value represent the contrast from the background.
% Make sure that the modulation value you want to use will generate RGB
% values in the monitor gamut.

cone_modulation = [0.15 0.0 0.0];

%=========================================================
% DEFINE PIXEL LEVELS.
%=========================================================
stimulusLutLevel = 127;

%=========================================================
% DRAW STIMULUS.
%=========================================================
theStimulus.size = 200;

% Draw the matching stimulus.
crsSetDrawMode(CRS.SOLIDFILL + CRS.CENTREXY);
crsSetPen1(stimulusLutLevel);
[error] = crsDrawOval(-10,0,theStimulus.size,theStimulus.size);

%=========================================================
% Set stimulus colour.
%=========================================================
fromCS               = 'CS_RGB';
toCS                 = 'CS_LMS';
fromColour           = BackgroundRGB;
[toColour ErrorCode] = ctGetColourTrival(fromCS,toCS,fromColour,deviceSPD,Sensors);
BackgroundLMS        = toColour;

%=========================================================
% Set stimulus modulation.
%=========================================================
amplIndex = 1;
itime     = linspace(-1, 1,nFramesPerCycle);
amplitude = (sin(pi*(itime)))/amplIndex;

%=========================================================

fprintf('\nAbout to show stimuli ...\n');

%=========================================================
% Write LUT.
%=========================================================
Buff_new(1:256,1) = BackgroundRGB(1);
Buff_new(1:256,2) = BackgroundRGB(2);
Buff_new(1:256,3) = BackgroundRGB(3);

% clc
for nLut=1:nFramesPerCycle
    % Modulate LMS values.
%     (1 + amplitude(nLut).*cone_modulation (:))
    stimulusLMS  = BackgroundLMS.*(1 + amplitude(nLut).*cone_modulation (:));
  
    % Convert LMS into RGB.
    fromCS      = 'CS_LMS';
    toCS        = 'CS_RGB';
    fromColour            = stimulusLMS';
    [toColour ErrorCode]  = ctGetColourTrival(fromCS,toCS,fromColour,deviceSPD,Sensors);
    stimulusRGB   = toColour

    % Set LUTBuffer.
    Buff_new(stimulusLutLevel,:) = stimulusRGB;
    crsLUTBUFFERWrite(nLut,Buff_new);

    theNumber=1000*nFramesPerCycle;
    theFrameDelay=1;
    theFirstLUT=1;
    theLastLUT=nFramesPerCycle;
    theLUTStoSkip=1;
    theStartLUT=1;
    theTriggerLUT=-1;
    crsLUTBUFFERCyclingSetup(theNumber,theFrameDelay,theFirstLUT,theLastLUT,theLUTStoSkip,theStartLUT,theTriggerLUT);
    
end
crsSetCommand(CRS.CYCLELUTENABLE);
crsSetDisplayPage(stimulusPage);