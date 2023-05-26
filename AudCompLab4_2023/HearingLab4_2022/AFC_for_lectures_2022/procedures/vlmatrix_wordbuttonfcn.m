% vlmatrix_wordbuttonfcn - custom function for VlMatrix procedure
% Adapted from OLSA_wordbuttonfcn, Version 1.30.0, last modified 21.03.2014 14:05

%------------------------------------------------------------------------------
% AFC for Mathworks MATLAB
%
% Author(s): Stephan Ewert
%
% Copyright (c) 1999-2014, Stephan Ewert. 
% Some rights reserved.
%
% This work is licensed under the 
% Creative Commons Attribution-NonCommercial-NoDerivs 3.0 Unported License (CC BY-NC-ND 3.0). 
% To view a copy of this license, visit
% http://creativecommons.org/licenses/by-nc-nd/3.0/ or send a letter to Creative
% Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
%------------------------------------------------------------------------------

function vlmatrix_wordbuttonfcn( wordNum, wordAlternative )

global work
global def
global msg

h = findobj('Tag','afc_win');
y = get(h,'UserData');		% ready flag: 0 = not ready, 1 = accept commands,
% 2 = accept only end and restart command,
% -1 = proceed in afc_work if terminate and end flags are not set
% -2 = accept only end
% 4 = any keyboard and window button (waitfor in afc_main)

% early out if not ready for any action yet SE 21.06.2007 16:45
if ( y ~= 1 )
    return;
end

% current step % FIXME doublicated code
if ( def.interleaved & def.OLSA_SplitTestlistAcrossTracks )
    % take every other sentence in interleaved tracks with split testlist
    currentSentence = work.pvind + (work.stepnum{work.pvind}(end)-1)*def.interleavenum;
else
    currentSentence = work.stepnum{work.pvind}(end);
end

% update response vector
% if ( wordNum ~= -1 & wordAlternative ~= -1 ) % -1 no update, 0 operator window, > 0 selection from matrix
if ( wordNum > 0 & wordAlternative > 0 )
    % write in response vector
    work.OLSA_ResponseVector{work.pvind}(wordNum)=wordAlternative;
elseif ( wordAlternative == 0 )
	% SE 20.03.2014 10:03 open test
	if ( wordNum ~= 0 )
		% did we already give this response? then toggle and set to not
		% selected again
		if ( work.OLSA_ResponseVector{work.pvind}(wordNum) == (work.OLSA_testlist{work.pvind}(currentSentence,wordNum) + 1) )
			work.OLSA_ResponseVector{work.pvind}(wordNum) = 0;
		else
			% write correct one in response vector
    	work.OLSA_ResponseVector{work.pvind}(wordNum) = work.OLSA_testlist{work.pvind}(currentSentence,wordNum) + 1;
    end
	else
		% all correct
		for i=1:5
			work.OLSA_ResponseVector{work.pvind}(i) = work.OLSA_testlist{work.pvind}(currentSentence,i) + 1;
		end 
	end
end

% write selected word in respective text box

% initialize percent correct
work.OLSA_PercentCorrect{work.pvind} = 0;

for i=1:5
    if (work.OLSA_ResponseVector{work.pvind}(i)~= 0 )
        h=findobj('Tag',['afc_wordSelection' num2str(i)]);

        if (~isempty(h))
            set(h,'string', msg.buttonString{work.OLSA_ResponseVector{work.pvind}(i),i});
        end

        if work.OLSA_ResponseVector{work.pvind}(i) == work.OLSA_testlist{work.pvind}(currentSentence,i) + 1
            work.OLSA_PercentCorrect{work.pvind} = work.OLSA_PercentCorrect{work.pvind} + 1/5;
        end
    % SE clear text field again if empty
    elseif (work.OLSA_ResponseVector{work.pvind}(i)== 0 )
        h=findobj('Tag',['afc_wordSelection' num2str(i)]);

        if (~isempty(h))
            set(h,'string', ' ');
        end
    end
end

%currentSentence
%sprintf('%i%i%i%i%i',work.OLSA_ResponseVector{work.pvind}-1)
%work.OLSA_PercentCorrect{work.pvind}

%strTmp

%sprintf('%i%i%i%i%i',afc_set.testlist(currentSentence,:))
%afc_set.testlistStrings(currentSentence,:)


% track the selected and presented words
strTmp = [];
for idx=1:5
    if ( work.OLSA_ResponseVector{work.pvind}(idx) == 0 )
        strTmp = [ strTmp '-' ];
    else
        strTmp = [ strTmp num2str(work.OLSA_ResponseVector{work.pvind}(idx)-1) ];
    end
end

work.OLSA_SelectedWords{work.pvind}{work.stepnum{work.pvind}(end)} = strTmp;
work.OLSA_PresentedWords{work.pvind}{work.stepnum{work.pvind}(end)} = sprintf('%i%i%i%i%i',work.OLSA_testlist{work.pvind}(currentSentence,:));
