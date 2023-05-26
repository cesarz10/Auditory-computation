
%------------------------------------------------------------------------------
% AFC for Mathworks MATLAB
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

% last modified 29-06-2006 15:37:51
% revision 0.93 beta, modified 07/08/00

% default savefcn for psychometric function measurement procedure 
% compare to psyfun_savefcn for afc 

function cs_psyfcn_savefcn

global def
global work

%str=[def.result_path def.expname '_' work.vpname work.userstr '_' work.condition];	% name of save file
str=[def.result_path work.filename];

%%%%%%%%% build headerline %%%%%%%%
headerl = '% run#   track#';

for i = 1:def.expparnum
   eval(['parunit = def.exppar' num2str(i) 'unit;']);
   headerl = [headerl '   exppar' num2str(i) '[' parunit ']'];
end
headerl = [headerl '   expvar[' def.expvarunit ']   n.pres   n.correct'];

if def.modelEnable == 1
    global simwork
    global simdef
    
    if isfield(simwork,'mueact_nonoise')
        bAdd_mue = 1;
    else
        bAdd_mue = 0;
    end
else
    bAdd_mue = 0;
end

if bAdd_mue
    headerl = [headerl '   mue_target   mue_ref   sigma'];
end
%%%%%%%% end headerline %%%%%%%%%

%%%%%%%% get data %%%%%%%%%%%%%%%
tmpinterleavenum = 1; % fixme: remove this in future versions
if def.interleaved > 0
   tmpinterleavenum = def.interleavenum;
end

for i = 1:tmpinterleavenum
    endstepsizeindex = sum(def.practicenum{i});
    restmp = work.expvar{i}(endstepsizeindex + 1 : end); % get expvars in measurement phase
    correcttmp = work.correct{i}(endstepsizeindex + 1 : end);

    res.range{i} = sort(def.expvar{i}); % determine expvar range
    res.n_pres{i} = res.range{i} * 0;
    res.n_correct{i} = res.range{i} * 0;

    for k = 1:length(restmp) % gathering statistics
        tmpindex = find(res.range{i} == restmp(k));
        res.n_pres{i}(tmpindex) = res.n_pres{i}(tmpindex) + 1;
        res.n_correct{i}(tmpindex) = res.n_correct{i}(tmpindex) + correcttmp(k);
    end
end
%%%%%%%% end get data %%%%%%%%%%%%%

if bAdd_mue == 0
    r='   %.8f   %i   %i';
else
    r='   %.8f   %i   %i    %.8f   %.8f   %.12f';
    rextended='   %.8f   %i   %i    %e   %e   %e';
end

ex = exist([str '.dat']);

if ex == 0 
    % If does not exist, then we are writing for the first time
    
    %dat=['% exppar[' def.exppar1unit ']   expvar[' def.expvarunit ']'];
    fid=fopen([str '.dat'],'w');
    fprintf(fid,['%s\n'],headerl);
	for k=1:tmpinterleavenum    
        for l = 1:length(res.range{k})
            fprintf(fid,'%i',work.numrun);			% current run number
            fprintf(fid,'   %i',k);					% track number if interleaved
            for i=1:def.expparnum
                eval(['tmp = work.int_exppar' num2str(i) '{' num2str(k) '};']);
                if def.exppartype(i) == 0
                    fprintf(fid,['   %.8f'],tmp);
                else
                    fprintf(fid,['   %s'],tmp);
                end
            end
            if bAdd_mue == 0 % default in AFC toolbox
                fprintf(fid,[r '\n'],[res.range{k}(l) res.n_pres{k}(l) res.n_correct{k}(l)] );
            else
                extra2print = [work.answer_extra1{k}(l) work.answer_extra2{k}(l) work.answer_extra3{k}(l)];
                % if length( find(extra2print < 1e-8) )>=2
                    fprintf(fid,[rextended '\n'],[res.range{k}(l) res.n_pres{k}(l) res.n_correct{k}(l) extra2print] );
                %else
                %    fprintf(fid,[r         '\n'],[res.range{k}(l) res.n_pres{k}(l) res.n_correct{k}(l) extra2print] );
                %runend
                % fprintf([r '\n'],[res.range{k}(l) res.n_pres{k}(l) res.n_correct{k}(l) work.answer_extra1{k}(l) work.answer_extra2{k}(l) work.answer_extra3{k}(l)] );
            end
        end
    end
	fclose(fid);
else
    % Then this is not the first time we write
    fid=fopen([str '.dat'],'a');
    for k=1:tmpinterleavenum    
        for l = 1:length(res.range{k})
            fprintf(fid,'%i',work.numrun);			% current run number
            fprintf(fid,'   %i',k);					% track number if interleaved
            for i=1:def.expparnum
                eval(['tmp = work.int_exppar' num2str(i) '{' num2str(k) '};']);
                if def.exppartype(i) == 0
                    fprintf(fid,['   %.8f'],tmp);
                else
                	fprintf(fid,['   %s'],tmp);
                end
            end
            if bAdd_mue == 0 % default in AFC toolbox
                fprintf(fid,[r '\n'],[ res.range{k}(l) res.n_pres{k}(l) res.n_correct{k}(l) ] );
            else
                extra2print = [work.answer_extra1{k}(l) work.answer_extra2{k}(l) work.answer_extra3{k}(l)];
                % if length( find(extra2print < 1e-8) )>=2
                fprintf(fid,[rextended '\n'],[res.range{k}(l) res.n_pres{k}(l) res.n_correct{k}(l) extra2print] );
                % fprintf(fid,[r '\n'],[res.range{k}(l) res.n_pres{k}(l) res.n_correct{k}(l) work.answer_extra1{k}(l) work.answer_extra2{k}(l) work.answer_extra3{k}(l)] );
            end
        end
    end
	fclose(fid);
end

if def.debug == 1
    save([str '.mat'], 'work'); 
end

% eof
