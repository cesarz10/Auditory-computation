% example_set - setup function of experiment 'example' -
%
% This function is called by afc_main when starting
% the experiment 'example'. It defines elements
% of the structure 'setup'. The elements of 'setup' are used 
% by the function 'example_user.m'.
% 
% If an experiments can be started with different experimental 
% conditions, e.g, presentation at different sound preasure levels,
% one might switch between different condition dependent elements 
% of structure 'setup' here.
%
% For most experiments it is also suitable to pregenerate some 
% stimuli here.
% 
% To design an own experiment, e.g., 'myexperiment',
% make changes in this file and save it as 'myexperiment_set.m'.
% Ensure that this function does exist, even if absolutely nothing 
% is done here.
%
% See also help example_cfg, example_user, afc_main

function example_set

global def
global work
global setup

% make condition dependend entries in structure set

% define the calibration level (assume 90 dB SPL for 0 dB FS)
setup.currentCalLevel = 113.5; %113;

% define signals in structure setup
setup.exppar1 = 4000;
setup.modsine = sin([0:def.intervallen-1]'*2*pi*setup.exppar1/def.samplerate);

if iscell(def.exppar1)
    setup.modcosine = cos([0:def.intervallen-1]'*2*pi*def.exppar1{1}/def.samplerate+pi);
else
    setup.modcosine = cos([0:def.intervallen-1]'*2*pi*def.exppar1/def.samplerate+pi);
end

setup.level = 70;
setup.hannlen = round(0.05*def.samplerate);
setup.window = hannfl(def.intervallen,setup.hannlen,setup.hannlen);

% eof