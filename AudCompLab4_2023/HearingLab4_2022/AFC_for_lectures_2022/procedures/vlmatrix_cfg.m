% olsa_cfg - configuration file for custom procedure
% Version 1.30.0, last modified 16.04.2013 09:36

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


% Add new configuration variables and their default values for the procedure here.
% Always include isfield check!
if ~isfield(def,'OLSA_TargetIntelligibility')
	def.OLSA_TargetIntelligibility = 0.5;		% Intelligibility of sentence to converge at
end

if ~isfield(def,'OLSA_WordAlternativeMapping')
	def.OLSA_WordAlternativeMapping = repmat([10:-1:1]',1,5); % Ordering of word buttons in matrix on sceen, as default exactly the order as defined in msg
end

if ~isfield(def,'OLSA_SplitTestlistAcrossTracks')
	def.OLSA_SplitTestlistAcrossTracks = 1;		% if 1, a single testlist is split across tracks of an interleaved run. 
																						% exppar1 must define only one testlist or the same for all tracks!
end

%if ~isfield(def,'newConfig')					% new config variable
%	def.newConfig = 0;
%end

% Append all new configuration variables to list
% varname, type ( 0 = not numeric, 1 = numeric, 2 = both possible)
def.configurationVariableList = [def.configurationVariableList; { ...
  'OLSA_TargetIntelligibility',1; ...
  'OLSA_WordAlternativeMapping',1; ...
  'OLSA_SplitTestlistAcrossTracks',1; ...
% 'newConfig', 1; ...
}];
                                                                            
% eof
