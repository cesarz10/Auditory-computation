% example_cfg - example measurement configuration file -
%
% This matlab skript is called by afc_main when starting
% the experiment 'example'.
% example_cfg constructs a structure 'def' containing the complete
% configuration for the experiment.
% To design an own experiment, e.g., 'myexperiment'
% make changes in this file and save it as 'myexperiment_cfg.m'.
% The default values of all parameters are defined in 'default_cfg.m'
%
% See also help example_set, example_user, afc_main
%

% Copyright (c) 1999-2013 Stephan Ewert, Universitaet Oldenburg.
% $Revision: 1.30 beta$  $Date: 25.01.2013 10:47 $

sentences_nr = 20; % 10 for simple, 20 for double, 30 for triple lists
def.sentences_nr = sentences_nr;

def.variableCheck = 0; % turn off warning about 'not-recognised sentences_nr field...'

% general measurement procedure
def.measurementProcedure = 'vlmatrix'; % measurement procedure
def.intervalnum = 1;                   % number of intervals
def.repeatnum = 1;                     % number of repeatitions of the experiment

% experimental variable (result of procedure yields dependent variable)
def.startvar = 0; % starting value of the tracking variable
def.expvarunit = 'dB'; % unit of the tracking variable
def.expvardescription = 'Speech reception threshold';	% description of the tracking variable

% limits for tracking variable
def.minvar = -100; % minimum value of the tracking variable
def.maxvar = 20;   % maximum value of the tracking variable
def.terminate = 1; % terminate execution on min/maxvar hit: 0 = warning, 1 = terminate
def.endstop = 3;   % Allows x nominal levels higher/lower than the limits before terminating (if def.terminate = 1)

%def.maxiter = 15; % set below

% interleaved (typical OLSA interleaving converging against 0.2 and 0.8 target intelligibility with split testlist across tracks

def.interleaved = 0;

% One list can be split to measure different SRTs (at different intelligibilities):
def.OLSA_SplitTestlistAcrossTracks = 0; %one single testlist is split across tracks of an interleaved run. If exppar1 must define only one testlist or the same for all tracks!
def.OLSA_TargetIntelligibility = 0.5; % 0.5 for SRT (50% correct) % [0.2; 0.8]; % assign different targetIntelligibility to the two tracks, must be column vector

if def.interleaved == 1
    def.interleavenum = 2;
    def.maxiter = sentences_nr/def.interleavenum; % must be half testlist length if interleaved and interleavenum == 2 FIXME calculate automatically later
else
    def.maxiter = sentences_nr; % must be testlist length if non interleaved FIXME calculate automatically later
end

% experimental variable (independent variable)
def.exppar1 = input('Enter a list number (1, 2, or 4 for this demo): '); % From 1 to 13. vector containing experimental parameters for which the exp is performed
def.exppar1unit = 'Testlist';				% unit of experimental parameter
def.exppar1description = 'Name of testlist';% description of the experimental parameter

% experimental variable 2...N (independent variable)
% add here if required

% interface, feedback and messages
def.afcwin = 'vlmatrix_win'; 		% overload response window, either 'olsa_win' for closed test, or 'olsa_open_win' for operator window (open test)
%def.afcwin = 'vlmatrix_open_win'; % NOT IMPLEMENTED YET, NOT TESTED % overload response window
def.windetail = 0;					% displays additional information in the response window
def.mouse = 0;							% enables mouse/touch screen control (1), or disables (0)
def.feedback = 0;						% no feedback
def.messages = 'vlmatrix'; % message configuration file, if 'autoSelect' AFC automatically selects depending on expname and language setting, fallback is 'default'. If 'default' or any arbitrary string, the respectively named _msg file is used.
def.language = 'EN'; % EN = english, DE = german, FR = french, DA = danish
def.soundmexMark = 0; % disable all button marking
def.markinterval = 0; % disable all button marking
def.OLSA_WordAlternativeMapping = repmat([10:-1:1]',1,5); % Ordering of word buttons in matrix on sceen 

% save paths and save function
def.result_path = vlmatrix_getpaths('result_path'); % where to save results
def.control_path = vlmatrix_getpaths('control_path'); % where to save control files
def.savefcn = 'default';			% function which writes results to disk

% samplerate, number of channels and sound output
def.samplerate = 44100;	% sampling rate in Hz
%def.bits = 32;			% output bit depth: 8 or 16 see def.externSoundCommand for 32 bits
%def.externSoundCommand = 'soundmexprofree';	%this experiment requires soundmexpro
% def.externSoundCommand = 'soundmexpro';
def.outputChannels = 2;
%def.trackMap =[0 0 0 1];
%def.deviceID = 0;

% continuos background signal and mixing (requires soundmex(pro) enabled: def.externSoundCommand = 'soundmexprofree', see above )
def.continuous = 0;			% set to 1 to use continous background noise
if def.continuous == 1
    error('AO, 4-Oct-2019: Not validated yet...')
    % % def.bgloopwav = '.\addons\sounds\pnoise_44k.wav';
    % % def.bgloopwav = '.\addons\sounds\pnoise_44k.wav';
    % % def.bgloopwav = '.\experiments\binolsa\maskers\masker.wav';
    % def.bgloopwav = './addons/sounds/VlMatrixnoise_ltass.wav';
    % def.bgloopwavFullScaleLevel = 80; % for volume setting if using soundmexPro, work.currentCalLevel has to be defined either explicitely _set, or
    %                               % by using def.calScript
end

% computing
def.allowpredict = 0; % not enabled otherwise for speech materials so far