
%------------------------------------------------------------------------------
% AFC for Mathwork�s MATLAB
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

% last modified 26-03-2014 12:46:45
% version 0.92, modified 23/02/00

function afc_pressfcn(num,x)

global work
global def

h = findobj('Tag','afc_win');
y = get(h,'UserData');		% ready flag: 0 = not ready, 1 = accept commands,
% 2 = accept only end and restart command,
% -1 = proceed in afc_work if terminate and end flags are not set
% -2 = accept only end
% 4 = any keyboard and window button (waitfor in afc_main)

% early out if not ready for any action yet SE 21.06.2007 16:45
if ( y == 0 )
    return;
end

% 02-06-2005 11:26 give focus back to an invisible dummy button, un**$** matlab 7 persistent button focus.
% This requires a version check, since Matlab 5.x or 6.x would crash on this line
if ( work.matlabVersion > 6 )
    uicontrol(findobj('Tag','afc_focusdummy'));
end

% remember pressed button
pressedbutton = x;

% 08.11.2006 09:56 check for start/end buttom
if ~isnumeric(x)	% called with a string
    %disp('here we go')
elseif x == 0		% called with zero, retrieve currentcharacter

    %% nice str compare 13/03/03, don't need that, checked for numeric x below
    %tmp = get(h,'currentcharacter');
    %for (i=1:num)
    %	if (tmp == num2str(i))
    %	x = i;
    %	end
    %end
    %if x == 0

    % new keyboardResponseButtonMapping 07.03.2007 10:36
    if ( ~isempty( def.keyboardResponseButtonMapping ) & ~strcmp(get(h,'currentcharacter'),'s') & ~strcmp(get(h,'currentcharacter'),'e') )
        % 16.07.2010 13:52
        % is it a numeric table of Ascii Dec codes
        if ( isnumeric( def.keyboardResponseButtonMapping ) )
            % this returns empty x if not mapped
            x = min( find( def.keyboardResponseButtonMapping == double( get(h,'currentcharacter') ) ));
        else
            % this returns empty x if not mapped
            x = min( find( strcmp( def.keyboardResponseButtonMapping, get(h,'currentcharacter') )));
        end
    else
        % old code
        x=str2num(get(h,'currentcharacter'));
        if isempty(x) == 1
            x=get(h,'currentcharacter');
        end
    end

else		% called with number, mark button
    hb=findobj('Tag',['afc_button' num2str(x)]);
    set(hb,'backgroundcolor',[0.5 0.5 0.5])
    set(hb,'backgroundcolor',[0.75 0.75 0.75])
end



if y==1 & isnumeric(x)
    % 13/03/03 x securely numeric now

    % 08.09.2005 14:26 new acceptButton
    if ( ~isempty( def.acceptButton ) )
        acceptButton = def.acceptButton;
    else
        acceptButton = [1:num];
    end

    if ( ismember(x, acceptButton ) ) %sum(x==[1:num])==1
        %work = afc_control(x,work);				% experimental run

        % SE 21.07.2010 15:09 set to "not ready" right now or afc_control might be called multiple times on mad button tabbing resulting in an error
        % set(h,'Userdata',0);
        % sanity check is done afc_control

        % 05.04.2012 15:03 WARNING DIRTY
        if isfield(work,'correctResponse') % Only TRUE for 1I-2AFC Method
            % SE check whether correctResponse is in accepted range
            if ( ismember(work.correctResponse, acceptButton ) )
                work.position{work.pvind}(end) = work.correctResponse;
            end
        end

        afc_control(x);
    end

elseif y == 2 & x == 's'
    work.terminate=1;				% abort and restart
    set(h,'Userdata',-1);
elseif y == 2 & x == 'e'
    work.terminate=1;				% abort and end
    work.abortall=1;
    set(h,'Userdata',-1);
elseif y == -2 & x == 'e'
    feval(def.afcwin, 'close');% end
    if ( def.showEnable > 0 )
        afc_show('close');
    end
    %afc_end							% end
elseif y == 4
    set(h,'Userdata',0);
end

% se 14.09.2006 11:48
% ahhhh, there is no other way to get rid of the R14 button focus when the user moves the mouse fast off
% the window (like I do). Help.
if ( (work.matlabVersion > 6 ) & (pressedbutton > 0) & isnumeric(pressedbutton))

    if (y == 4)
        % waitfor in main, clear focus from all buttons
        clearjob = [1 def.intervalnum];

        % dont ask why, otherwise 1st button cant get cleared
        uicontrol(findobj('Tag','afc_focusdummy'));
    else
        clearjob = [pressedbutton pressedbutton];
    end

    for ( clearb = clearjob(1):clearjob(2) )
        %SE 26.03.2014 11:35 only clear enabled buttons
        if ( ~(def.mouse == 0) & ( isempty(def.acceptButton) | ismember(clearb, def.acceptButton ) ) )
            hb=findobj('Tag',['afc_button' num2str(clearb)]);
            set(hb,'enable','off');
            drawnow;
            set(hb,'enable','on');
            drawnow;
        end
    end
end
