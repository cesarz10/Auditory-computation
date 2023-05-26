% example_set - setup function of experiment 'example' -
%
% This function is called by afc_main when starting
% the experiment 'example'. It defines elements
% of the structure 'set'. The elements of 'set' are used 
% by the function 'example_user.m'.
% 
% If an experiments can be started with different experimental 
% conditions, e.g, presentation at different sound preasure levels,
% one might switch between different condition dependent elements 
% of structure 'set' here.
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

function exampleOLSA_set

global def
global work
global setup

% make condition dependend entries in structure set

switch work.condition
case 'cd1'
   %set.level=1;
case 'cd2'
   %set.level=0.25;
otherwise
   % error('condition not recognized');
   warning('AO: condition not recognized (check this)');
end

% load test list
%set.testlist = work.exppar1;

%%% AO:
dir_afc = Get_UGent_paths('tb_AFC');

testListPath = [dir_afc 'procedures' filesep 'olsa' filesep 'TestLists' filesep];

% select 20 or 30 testlists or demo
%testListFolderFilePart = '30\olsa30.';
%testListFolderFilePart = '20\olsa20.';
testListFolderFilePart = 'demo/demo30.'; % \ changed by / by AO

for idx = 1:def.interleavenum
	testListFile = [ testListPath testListFolderFilePart num2str(work.int_exppar1{idx}) ];
	
	tmp = strvcat(textread( testListFile ,'%s','headerlines',4));
	tmp = tmp(1:8:end-1,:);
	setup.testlistStrings{idx} = tmp(:,1:5);
	
%if (0)
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
%end
	
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

work.currentCalLevel = 100; %90;
%soundmexpro('show');

warning('AO: Calibration assumed to be equal to AMT convention...')

% eof