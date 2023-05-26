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
% AM detection thresholds is the minimal modulation depth required to distinguish a modulated tone from a pure tone

% Copyright (c) 1999-2013 Stephan Ewert.

% The results will get saved in the this directory.
dir_where = [fileparts(which('afc.m')) filesep];
result_path = [dir_where 'experiments' filesep 'AM_detection_ptp' filesep 'results' filesep];


% general measurement procedure
def.measurementProcedure = 'transformedUpDown';	% measurement procedure
def.intervalnum = 3;				% number of intervals
def.rule = [1 2];				% [up down]-rule: [1 2] = 1-up 2-down
def.varstep = [10 5 3 1];				% [starting stepsize ... minimum stepsize] of the tracking variable
def.steprule = -1;				% stepsize is changed after each upper (-1) or lower (1) reversal
def.reversalnum = 6;				% number of reversals in measurement phase
def.repeatnum = 1;				% number of repeatitions of the experiment

% experimental variable (result of procedure yields dependent variable)
def.startvar = -6;				% starting value of the tracking variable
def.expvarunit = 'dB';				% unit of the tracking variable
def.expvardescription = 'Modulation index';	% description of the tracking variable

% limits for experimental variable
def.minvar = -100;				% minimum value of the tracking variable
def.maxvar = 0;					% maximum value of the tracking variable
def.terminate = 1;				% terminate execution on min/maxvar hit: 0 = warning, 1 = terminate
def.endstop = 3;				% Allows x nominal levels higher/lower than the limits before terminating (if def.terminate = 1) 

% experimental parameter (independent variable)
def.exppar1 = [120];				% vector containing experimental parameters for which the exp is performed
def.exppar1unit = 'Hz';				% unit of experimental parameter
def.exppar1description = 'Modulation frequency';% description of the experimental parameter

% experimental parameter 2...N (independent variable)
% 1333 is the unprocessed condition
def.exppar2 = [1333];				% vector containing experimental parameters for which the exp is performed
%1001 is the new processing for 
def.exppar1unit = 'Type';				% unit of experimental parameter
def.exppar1description = 'Profile to compensate';% description of the experimental parameter
def.parrand = 1;

% interface, feedback and messages 
def.mouse = 1;					% enables mouse/touch screen control (1), or disables (0) 
def.markinterval = 1;				% toggles visual interval marking on (1), off(0)
def.feedback = 1;				% visual feedback after response: 0 = no feedback, 1 = correct/false/measurement phase
def.messages = 'default';			% message configuration file, if 'autoSelect' AFC automatically selects depending on expname and language setting, fallback is 'default'. If 'default' or any arbitrary string, the respectively named _msg file is used.
def.language = 'EN';				% EN = english, DE = german, FR = french, DA = danish
def.windetail = 1;

% save paths and save function
def.result_path = result_path;				% where to save results
def.control_path = result_path;				% where to save control files
def.savefcn = 'default';			% function which writes results to disk

% samplerate and sound output
def.samplerate = 48000;	% 44100;		% sampling rate in Hz
def.intervallen = 24000; % 22050;			% length of each signal-presentation interval in samples (might be overloaded in 'expname_set')
def.pauselen = 24000; % 22050;  				% length of pauses between signal-presentation intervals in samples (might be overloaded in 'expname_set')
def.presiglen = 10;				% length of signal leading the first presentation interval in samples (might be overloaded in 'expname_set')
def.postsiglen = 10;				% length of signal following the last presentation interval in samples (might be overloaded in 'expname_set')
def.bits = 16;					% output bit depth: 8 or 16 see def.externSoundCommand for 32 bits

% computing
def.allowpredict = 0;				% if 1 generate new stimuli during sound output if def.markinterval disabled

% tweaking
%def.keyboardResponseButtonMapping = {'a','s','d'};
%def.soundmexMark = 1;
%def.markinterval = 1;

% def.internSoundCommand = 'playrec';
%def.internSoundCommand = 'wavplay';
%def.internSoundCommand = 'audioplayer'; % 'sound' or 'audioplayer'
%def.markIntervalDelay = 0.5;        % tweak if audioplayer is used, default is 0

% def.externSoundCommand = 'playrec';

%def.bits = 16;
%def.deviceID = 0;

% eof
