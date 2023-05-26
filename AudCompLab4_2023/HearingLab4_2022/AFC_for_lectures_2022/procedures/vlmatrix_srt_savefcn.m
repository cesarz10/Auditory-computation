% olsa_srt_savefcn - custom save function for OLSA procedure
% Version 1.30.0, last modified 22.01.2013 09:36

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

function vlmatrix_srt_savefcn

global def
global work

%global bla

%str=[def.result_path def.expname '_' work.vpname work.userstr '_' work.condition];	% name of save file
str=[def.result_path work.filename];

%%%%%%%%% build headerline %%%%%%%%
if def.interleaved > 0
   headerl = '% track#';
else
   headerl = '%';
end
for i = 1:def.expparnum
   eval(['parunit = def.exppar' num2str(i) 'unit;']);
   headerl = [headerl '   exppar' num2str(i) '[' parunit ']'];
end
headerl = [headerl '   expvar[' def.expvarunit ']'];
%%%%%%%% end headerline %%%%%%%%%

% write experimental variable at all reversals during the measurement phase to disk
%res=[work.exppar work.expvarrev(end - def.reversalnum+1:end)];
%
%r=[];
%for i=1:def.reversalnum
%   r=[r '   %8.8f'];
%end

% write mean and standard deviation to disk  

%res = work.expvarrev(end - def.reversalnum+1:end);
%res=[work.exppar mean(res) std(res,1)];

tmpinterleavenum = 1;
if def.interleaved > 0
   tmpinterleavenum = def.interleavenum;
end

for i = 1:tmpinterleavenum
	 	
    %   res{i} = FitOlLX(work.expvar{i}(1:end), work.answer{work.pvind}(1:end), def.OLSA_TargetIntelligibility{work.pvind});
    %
    %	r='   %8.8f';
    %SE 22.01.2013 18:18 save SNR and target intelligibility
    thresVl = FitVlLX(work.expvar{i}(1:end), work.answer{i}(1:end), def.OLSA_TargetIntelligibility{i});
	res{i}=[thresVl def.OLSA_TargetIntelligibility{i}];
    
	r='   %8.8f   %8.8f';
	
end

% write median, lower and upper quartiles to disk
% not yet implemented

ex = exist([str '.dat']);

if ex == 0
   %dat=['% exppar[' def.exppar1unit ']   expvar[' def.expvarunit ']'];
   fid=fopen([str '.dat'],'w');
   fprintf(fid,['%s\n'],headerl);
   for k=1:tmpinterleavenum
      if def.interleaved > 0
         fprintf(fid,'%i',k);
      end
   	for i=1:def.expparnum
      	eval(['tmp = work.int_exppar' num2str(i) '{' num2str(k) '};']);
      	if def.exppartype(i) == 0
         	fprintf(fid,['   %8.8f'],tmp);
      	else
         	fprintf(fid,['   %s'],tmp);
      	end
      end
      fprintf(fid,[r '\n'],res{k});
   end
   
   
   %fprintf(fid,['%8.8f' r '\n'],res);
	fclose(fid);
else
   fid=fopen([str '.dat'],'a');
   for k=1:tmpinterleavenum
      if def.interleaved > 0
         fprintf(fid,'%i',k);
      end
   	for i=1:def.expparnum
      	eval(['tmp = work.int_exppar' num2str(i) '{' num2str(k) '};']);
      	if def.exppartype(i) == 0
         	fprintf(fid,['   %8.8f'],tmp);
      	else
         	fprintf(fid,['   %s'],tmp);
      	end
      end
      fprintf(fid,[r '\n'],res{k});
   end

   %	fprintf(fid,['%8.8f' r '\n'],res);
	fclose(fid);
end

if def.debug == 1
save([str '.mat'], 'work'); 
end

% eof
