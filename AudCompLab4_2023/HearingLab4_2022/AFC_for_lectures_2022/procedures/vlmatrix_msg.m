% olsa_msg - english message definition file -
% Version 1.30.0, last modified 16.04.2013 09:36
%
% ready_msg			displayed when ready for user response
% measure_msg		displayed when entering measurement phase
% correct_msg		displayed after correct response
% false_msg			displayed after false response
% maxvar_msg		displayed when maxvar is reached
% minvar_msg		displayed when minvar is reached
% start_msg			displayed when the experiment starts
% next_msg			displayed when the next parameter is presented
% finished_msg		displayed when the experiment is finished

%------------------------------------------------------------------------------
% AFC for Mathworks MATLAB
%
% Author(s): Stephan Ewert
%
% Copyright (c) 1999-2013, Stephan Ewert. 
% Some rights reserved.
%
% This work is licensed under the 
% Creative Commons Attribution-NonCommercial-NoDerivs 3.0 Unported License (CC BY-NC-ND 3.0). 
% To view a copy of this license, visit
% http://creativecommons.org/licenses/by-nc-nd/3.0/ or send a letter to Creative
% Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
%------------------------------------------------------------------------------

msg=struct(...
    'measure_msg','Beginning measurement',	...
    'correct_msg','--- CORRECT ---',			...
    'false_msg','--- WRONG ---',				...
    'maxvar_msg','Maximum level reached',	...
    'minvar_msg','Minimum level reached' ...
    );

msg.ready_msg    = {'Select the words from', ...
    'the matrix and press OK to continue'};
msg.start_msg    = {'You have started a new measurement.', ...
    'Press any key to continue.'};
msg.next_msg     = {'End of Run.', ...
    'Press "s" for a new run or "e" to end.'};
msg.finished_msg = {'Experiment Done.', ...
    'Press "e" to end.'};

msg.experiment_windetail = 'Experiment: %s';
msg.measurement_windetail = 'Measurement %d of %d';
msg.measurementsleft_windetail = '%d of %d measurements left';

msg.okButtonString = 'OK';
msg.allCorrectButtonString = 'Alles verstanden';

msg.startButtonString = 'Start';
msg.endButtonString = 'End';

% List copied from ../Speech/dutch/Matrix/woordenmatrix.xls
msg.buttonString = {...
'David' ,	'draagt','twee',   'beige', 'bedden'; ...
'Ellen' , 	'heeft', 'drie',   'blauwe','boten'; ...
'Emma'  , 	'kiest', 'vier',   'bruine','doeken'; ...
'Jacob' ,	'koopt', 'vijf',   'gele',  'dozen'; ...
'Jeroen',	'krijgt','zes',	   'grijze','fietsen'; ...
'Johan' ,	'leent', 'acht',   'groene','jassen'; ...
'Lucas' ,	'maakt', 'tien',   'paarse','kousen'; ...
'Sara'  ,	'wint',  'elf',    'rode',	'manden'; ...
'Sofie' ,	'ziet',  'twaalf', 'witte',	'pennen'; ...
'Thomas',	'zoekt', 'veel',   'zwarte','ringen'};

% eof
