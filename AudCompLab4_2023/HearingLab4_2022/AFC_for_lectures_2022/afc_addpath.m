function afc_addpath
%------------------------------------------------------------------------------
% AFC for Mathwork's MATLAB
%
% Version 1.40.0
%
% Author(s): Stephan Ewert
%
% Copyright (c) 1999-2014, Stephan Ewert. 
% All rights reserved.
%
% This work is licensed under the 
% Creative Commons Attribution-NonCommercial-NoDerivs 4.0 International License (CC BY-NC-ND 4.0). 
% To view a copy of this license, visit
% http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to Creative
% Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
%------------------------------------------------------------------------------

% last modified 26-04-2014 11:47:24

addpath([fileparts(which(mfilename))])
addpath([fileparts(which(mfilename)), filesep, 'scripts'])
addpath([fileparts(which(mfilename)), filesep, 'base'])
addpath([fileparts(which(mfilename)), filesep, 'experiments'])
addpath(genpath([fileparts(which(mfilename)), filesep, 'experiments', filesep, 'AM_detection_ptp']))
addpath(genpath([fileparts(which(mfilename)), filesep, 'experiments', filesep, 'SRT_task']))
addpath([fileparts(which(mfilename)), filesep, 'models'])
addpath([fileparts(which(mfilename)), filesep, 'calibration'])
addpath([fileparts(which(mfilename)), filesep, 'gui'])
addpath([fileparts(which(mfilename)), filesep, 'soundmexpro', filesep, 'bin'])
addpath([fileparts(which(mfilename)), filesep, 'procedures'])
addpath([fileparts(which(mfilename)), filesep, 'results'])
addpath([fileparts(which(mfilename)), filesep, 'addons'])
addpath([fileparts(which(mfilename)), filesep, 'sessions'])
