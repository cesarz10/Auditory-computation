% vlmatrix_init - init script for VlMatrix procedure
% Version 1.30.0, last modified 11.04.2013 16:36

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

% add procedure specific path if required paths
%addpath([fileparts(which(mfilename)), filesep, 'vlmatrix'])
% including subfolders 22.07.2013 14:27
addpath(genpath([fileparts(which(mfilename)), filesep, 'vlmatrix']))

% % some splash screen
% vlmatrix_splashscreen( 'open' );

% change save function to method specific default if default save function is requested
if ( strcmp( def.savefcn, 'default' ) )
    def.savefcn = 'vlmatrix_srt';
end

% sanity check
if ( def.allowpredict ~= 0 )
	def.allowpredict = 0; % not implemented for OLSA
	warning('def.allowpredict must equal 0 for OLSA or other speech matrials');
end

% prepare for interleaving
if ~iscell(def.OLSA_TargetIntelligibility)
    expvarTmp = def.OLSA_TargetIntelligibility;
    def=rmfield(def,'OLSA_TargetIntelligibility');

    if ( def.interleaved )
        if size( expvarTmp, 1 ) == 1
            expvarTmp = repmat(expvarTmp, def.interleavenum, 1);
        elseif size( expvarTmp, 1 ) ~= def.interleavenum
            error('afc_main: def.expvar dimensions mismatch def.interleavenum');
        end
    end

    for ( i=1:def.interleavenum )
        def.OLSA_TargetIntelligibility{i} = expvarTmp(i);
    end
end

for (i=1:def.interleavenum)

    % initialize starting value for all interleaved tracks
    work.expvarnext{i} = def.startvar{i};				% pre buffer (cell array) for expvaract, copied to expvaract in afc_interleave

    % add procedure specific entries to structure work
    work.OLSA_Direction{i} = [];
    work.OLSA_PercentCorrect{i} = 0;
    work.OLSA_SelectedWords{i} = [];
    work.OLSA_PresentedWords{i} = [];

end

% eof
