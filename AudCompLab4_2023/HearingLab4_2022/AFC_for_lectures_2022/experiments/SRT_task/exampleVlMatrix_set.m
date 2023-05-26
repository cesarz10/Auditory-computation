function exampleVlMatrix_set
% function exampleVlMatrix_set - setup function of experiment 'exampleVlMatrix'
%
% See also help exampleVlMatrix_cfg, exampleVlMatrix_user, afc_main

global def
global work
global setup

% make condition dependend entries in structure set

work.sentences_nr = 20;
work.varnoise = 0;
switch work.condition
case {'00','broadband','original'} % % default
   work.sounds2look = 'vlmatrix00';
case {'01'}
   work.sounds2look = 'vlmatrix01';
case {'02'}
   work.sounds2look = 'vlmatrix02';
case 'test'
    work.sentences_nr = 10;
   work.sounds2look = 'vlmatrix00';
otherwise
   error('Condition not recognised, this should indicate the number of sentences in the lists');
end
work.varspeech = ~work.varnoise;

%%% Loading some locations:
speechMaterialPath = vlmatrix_getpaths(work.sounds2look); 
testListPath       = vlmatrix_getpaths('TestLists'); 

work.speechMaterialPath = speechMaterialPath;
work.noise_filename = 'VlMatrixnoise_ltass.wav';

%%% Loading the noise:
work.noise_path  = speechMaterialPath;
[noise,fs_noise] = audioread([speechMaterialPath work.noise_filename]);
if fs_noise ~= def.samplerate
    error('The noise to be used has a different sampling frequency than the speech samples');
end
work.noise = noise;
%%%

if work.sentences_nr ~= def.sentences_nr
    error('Condition specified in afc_main mismatches the number of lists specified in cfg.m file (changed manually)...')
end

% select 20 or 30 testlists or demo (go to exampleVlMatrix_cfg)
switch work.sentences_nr
    case 10
        testListFolderFilePart = 'vlmatrix10.';
    case 20
        testListFolderFilePart = 'vlmatrix20.';
    otherwise
        error('Only examples for either 10 or 20 sentences')
end

for idx = 1:def.interleavenum
	testListFile = [ testListPath testListFolderFilePart num2str(work.int_exppar1{idx}) ];
	
	tmp = strvcat(textread( testListFile ,'%s','headerlines',4));
	tmp = tmp(1:8:end-1,:);
	setup.testlistStrings{idx} = tmp(:,1:5);
	
	%%%%%%%%%%% sanity check
	if ( def.interleaved & def.OLSA_SplitTestlistAcrossTracks )
		if (size(setup.testlistStrings{idx},1) ~= def.maxiter*def.interleavenum )
			error('Items in testlist do not match def.maxiter. def.OLSA_SplitTestlistAcrossTracks is set to 1, def.maxiter must equal half the number of items in the testlist.');
		end
	else
		if (size(setup.testlistStrings{idx},1) ~= def.maxiter )
			error('Items in testlist do not match def.maxiter');
		end
	end
	%%%%%%%%%%%%%%%%%%%%%%
	
	for i=1:size(setup.testlistStrings{idx},1)
	    for k=1:5
	    		% SE 23.04.2013 17:11 moved testlist{} to structure work
	        work.OLSA_testlist{idx}(i,k) = str2num(setup.testlistStrings{idx}(i,k));
	    end
	end
	
	%setup.testlist
	
	% define signals in structure set
	work.OLSA_ResponseVector{idx} = zeros(1,5);
end

%%% Cosine ramp to be applied to the background noise
dur_ramp     = 50e-3; 
dur_ramp_samples = def.samplerate*dur_ramp;
win_both     = hannfl(2*dur_ramp_samples,dur_ramp_samples,dur_ramp_samples);
setup.win_up = win_both(1:dur_ramp_samples);
setup.win_dn = win_both(dur_ramp_samples:end);
%%%

work.currentCalLevel = 100; 
%soundmexpro('show');
warning('AO: Calibration assumed to be equal to AMT convention...')

% eof
