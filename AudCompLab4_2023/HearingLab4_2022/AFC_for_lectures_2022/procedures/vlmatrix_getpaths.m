function misc = vlmatrix_getpaths(type) 
% function misc = vlmatrix_getpaths(type) 
%   
%   Get local paths that need to be used for running the Flemish Matrix 
%   material. If 'type' is specified, then 'misc' will be a character 
%   containing the
%   that specific path. See example 2 (below).
% 
% % Example 1: to get all important directories
%       misc = vlmatrix_getpaths;
%
% % Example 2: to get all important directories
%       misc = vlmatrix_getpaths('VlMatrix');
%
% Programmed by Alejandro Osses (ale.a.osses@gmail.com), Hearing Technology,
% WAVES, UGent, Belgium, 2018-2019
% Created on    : 04/10/2019
% Last update on: 04/10/2019
% Last use on   : 27/04/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir_where = [fileparts(which('afc.m')) filesep];

l = length('vlmatrix');
if length(type) >= l
    if strcmp(type(1:l),'vlmatrix')
        switch type
            case 'vlmatrix00'
                misc.VlMatrix  = [dir_where 'sounds' filesep 'vlmatrix' filesep '00-original' filesep];
            case 'vlmatrix01'
                misc.VlMatrix  = [dir_where 'sounds' filesep 'vlmatrix' filesep '01-LP' filesep];
            case 'vlmatrix02'
                misc.VlMatrix  = [dir_where 'sounds' filesep 'vlmatrix' filesep '02-HP' filesep];
            otherwise
                % Just unprocessed
                misc.VlMatrix  = [dir_where 'sounds' filesep 'vlmatrix' filesep '00-original' filesep];
        end
        type = 'VlMatrix'; % for passing final check in this script
    end
end

misc.tb_AFC    = dir_where;
misc.TestLists = [misc.tb_AFC 'procedures' filesep 'vlmatrix' filesep 'TestLists' filesep];

misc.result_path = [dir_where filesep 'experiments' filesep 'SRT_task' filesep 'results' filesep];
misc.control_path = misc.result_path;

if nargin==1
    if (~isfield(misc,type))
        error('Invalid type');
    end
    misc=misc.(type);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
