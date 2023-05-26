% olsa_proto - custom protocol file for OLSA procedure
% Version 1.30.0, last modified 16.04.2013 09:36

%------------------------------------------------------------------------------
% AFC for Mathwork’s MATLAB
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

% The protocol file is opened and closed in afc_proto.m,
% here only the writing of the procedure-specific protocol entry is performed.
% The protocol entry must follow AFC's protocol file conventions

% write the header for next entry
fprintf(fid,['%s\n'],'%');
fprintf(fid,['%s\n'],['% new entry (date/experiment/' parstr 'sentence/response/' varstr 'percent correct/procedure)']);

% write date/time and current experiment file name
fprintf(fid,['%s\n'],dateTime);
fprintf(fid,['%s\n'],expFileName);

% prepare and write the data for each track (just one track if not interleaved, def.interleavenum = 1)
for i = 1:def.interleavenum

		% prepare the data for saving
    rev = work.answer{i}(1:end);

    expvarrev{i} = num2str(rev,'%8.8f   '); 									% expvarrev at reversals in measurement phase
    expvar{i} = num2str(work.expvar{i},'%8.8f   ');

    tmpX = work.OLSA_PresentedWords{i}{1};
    for idx=2:length(work.OLSA_PresentedWords{i})
        tmpX = [tmpX '  ' work.OLSA_PresentedWords{i}{idx}];
    end
    pos{i} = tmpX;  % sentence words
    tmpX = work.OLSA_SelectedWords{i}{1};
    for idx=2:length(work.OLSA_SelectedWords{i})
        tmpX = [tmpX '  ' work.OLSA_SelectedWords{i}{idx}];
    end
    res{i} = tmpX;    % pressed words


    %%%%%%%%%%% start writing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % write track number if interleaved
    if def.interleaved > 0
        fprintf(fid,['%s\n'],['% track ' num2str(i)]);
    end

    % current values for exppar(s)
    for k = 1:def.expparnum
        eval(['expparStr = num2str(work.int_exppar' num2str(k) '{i}, ''%8.8f'' );']);
        fprintf(fid,['%s\n'], expparStr );
    end

    % selected data as defined above
    fprintf(fid,['%s\n'],pos{i});
    fprintf(fid,['%s\n'],res{i});
    fprintf(fid,['%s\n'],expvar{i});
    fprintf(fid,['%s\n'],expvarrev{i});

    % append name of measurement procedure
    fprintf(fid,['%s\n'], def.measurementProcedure );
end

% eof
