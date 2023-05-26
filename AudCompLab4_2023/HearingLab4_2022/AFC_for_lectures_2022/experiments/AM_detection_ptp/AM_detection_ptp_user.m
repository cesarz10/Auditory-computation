% example_user - stimulus generation function of experiment 'example' -
%
% This function is called by afc_main when starting
% the experiment 'example'. It generates the stimuli which
% are presented during the experiment.
% The stimuli must be elements of the structure 'work' as follows:
%
% work.signal = def.intervallen by 2 times def.intervalnum matrix.
%               The first two columns must contain the test signal
%               (column 1 = left, column 2 = right) ...
% 
% work.presig = def.presiglen by 2 matrix.
%               The pre-signal. 
%               (column 1 = left, column 2 = right).
%
% work.postsig = def.postsiglen by 2 matrix.
%                The post-signal. 
%               ( column 1 = left, column 2 = right).
%
% work.pausesig = def.pausesiglen by 2 matrix.
%                 The pause-signal. 
%                 (column 1 = left, column 2 = right).
% 
% To design an own experiment, e.g., 'myexperiment',
% make changes in this file and save it as 'myexperiment_user.m'.
% Ensure that the referenced elements of structure 'work' are existing.
%
% See also help example_cfg, example_set, afc_main

function AM_detection_ptp_user

global def
global work
global setup

% 1000 ... 4000 Hz bandlimited Gaussian noise
% The experiment is meant to be 3-interval 3-alternative forced choice, so we need 3 signals.
% One signal holds the target amplitude modulation.
% work.expvaract holds the current value of the tracking variable of the experiment.

%% experiment
tref1 = setup.modsine; % * 10^((setup.level - setup.currentCalLevel)/20);
tref2 = setup.modsine; % * 10^((setup.level - setup.currentCalLevel)/20);
tuser = setup.modsine .* (1 + (10^(work.expvaract/20) * setup.modcosine));
% * 10^((setup.level - setup.currentCalLevel)/20) .* ... // holds the target signal

% pure tone will have rms of 70 dB-SPL
scaling_factor = 10^((setup.level - setup.currentCalLevel)/20) / rms(tref1);
%.* (1 + setup.modcosine);    % 100% modulated AM tone - for ptp calibration

tref1 = tref1 * scaling_factor;
tref2 = tref2 * scaling_factor;

%%% for having the same peak
% tuser = tuser * (max(abs(tref1)) / max(abs(tuser)));
%%% for having the same carrier amplitude
tuser = tuser * scaling_factor;

cs_par=[fix(work.exppar2./100),mod(fix(work.exppar2./10),10),mod(work.exppar2,10)];

tuserp = process_speech_nofb_mod(tuser,cs_par,1);
%%% for having the same peak
tuser = tuserp * (max(abs(tuser)) / max(abs(tuserp)));
%%% for having the same rms
% tuser = tuserp * rms(tuser) / rms(tuserp);

%%
% Hanning windows

tref1 = tref1 .* setup.window;
tref2 = tref2 .* setup.window;
tuser = tuser .* setup.window;

% pre-, post- and pausesignals (all zeros)

presig = zeros(def.presiglen,2);
postsig = zeros(def.postsiglen,2);
pausesig = zeros(def.pauselen,2);

% make required fields in work

work.signal = [tuser tuser tref1 tref1 tref2 tref2];	% left = right (diotic) first two columns holds the test signal (left right)
work.presig = presig;											% must contain the presignal
work.postsig = postsig;											% must contain the postsignal
work.pausesig = pausesig;										% must contain the pausesignal

% eof