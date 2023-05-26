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

function exampleOLSA_user

global def
global work
global setup


%if ~isfield( work, 'firstPresentation' ) 
%		work.firstPresentation = 0;
%end
%
%if ( work.stepnum{work.pvind} == 1 )
%	if ( work.firstPresentation == 0 )
%		pause(3)
%		work.firstPresentation = 1;
%	end
%else
%	work.firstPresentation = 0;
%end

% path to official speech material
% has to be received from Hoertech, Daniel Berg
dir_afc = Get_UGent_paths('tb_AFC'); %%% AO
speechMaterialPath = [dir_afc 'procedures' filesep 'olsa' filesep 'wavs' filesep];

if ( def.interleaved & def.OLSA_SplitTestlistAcrossTracks )
	% take every other sentence in interleaved tracks with split testlist
	currentSentence = work.pvind + (work.stepnum{work.pvind}(end)-1)*def.interleavenum;
else
	currentSentence = work.stepnum{work.pvind}(end);
end

%wavName = ['procedures\olsa\OlSa_raw\' setup.testlistStrings{work.pvind}(currentSentence,1:5) '.wav'];
wavName = [speechMaterialPath setup.testlistStrings{work.pvind}(currentSentence,1:5) '.wav'];
[sentence,fs] = wavread(wavName);

% check samplerate
if ( fs ~= def.samplerate )
    error('wav files have wrong samplerate')
end

% DEBUG: print current rms of OLSA sentence 
%20*log10(rms(sentence))
% DO NOT NORMALIZE INDIVIDUALLY, sentences have to be
% assumed to have -24 dB rms re full scale.

% scale the signal and take left channel to make it mono
sentence = sentence(:,1) * 10^(work.expvaract/20);


% we might want to change the track mapping
%def.trackMap =[0 0 0 1];
%afc_sound('bgloopwav_restart');

% pre-, post- and pausesignals (all zeros)

presig = zeros(def.presiglen,def.outputChannels);
postsig = zeros(def.postsiglen,def.outputChannels);
pausesig = zeros(def.pauselen,def.outputChannels);

% make required fields in work

work.signal = [ sentence sentence ];	% left = right (diotic) first two columns holds the test signal (left right)
work.presig = presig;											% must contain the presignal
work.postsig = postsig;											% must contain the postsignal
work.pausesig = pausesig;										% must contain the pausesignal

% eof