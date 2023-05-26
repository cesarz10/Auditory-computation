% olsa_procedure - custom measurement procedure
% Version 1.30.0, last modified 10.04.2013 16:36
%
% called by afc_control
%
% answer - holds the input passed to afc_control by _pressfcn (in most cases the number of the button pressed)
%
% must update:
% work.expvarnext{work.pvind} - updated value of experimental variable
% work.stepnum{work.pvind} - history of steps, append last step number
% work.trackfinished{work.pvind}
% work.finishedcount

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


% Append the current value of the experimental variable to the history
% FIXME this should go to afc_control, FIXED 23.04.2013 16:39
%work.expvar{work.pvind}=[work.expvar{work.pvind} work.expvarnext{work.pvind}];

% Calculate the next value of the experimental variable dependent on answer
if ( isempty(work.OLSA_Direction{work.pvind} ) )
    % first (initialization) call of level control
    % work.firstPresentation = 1;
    [work.expvarnext{work.pvind}, work.OLSA_Direction{work.pvind}, work.OLSA_SpeedFactor{work.pvind}] = LevelControlVlLX(work.expvarnext{work.pvind},answer,def.OLSA_TargetIntelligibility{work.pvind});
else
    
    [work.expvarnext{work.pvind}, work.OLSA_Direction{work.pvind}, work.OLSA_SpeedFactor{work.pvind}] = LevelControlVlLX(work.expvarnext{work.pvind},answer,def.OLSA_TargetIntelligibility{work.pvind},work.OLSA_Direction{work.pvind}, work.OLSA_SpeedFactor{work.pvind});
end


% Sanity check lower and upper bounds
if work.expvarnext{work.pvind} > def.maxvar
    work.expvarnext{work.pvind} = def.maxvar;
    work.minmaxcount{work.pvind} = work.minmaxcount{work.pvind} + 1;
    if ( work.predict == 0 & def.afcwinEnable > 0 )
        feval(def.afcwin, 'maxvar');
    end
elseif work.expvarnext{work.pvind} < def.minvar
    work.expvarnext{work.pvind} = def.minvar;
    work.minmaxcount{work.pvind} = work.minmaxcount{work.pvind} + 1;
    if ( work.predict == 0 & def.afcwinEnable > 0 )
        feval(def.afcwin, 'minvar');
    end
end


% If we are in prediction state the procedure is not actually run. Time to leave.
if work.predict == 1
    return;
end


%----------------------- experimental run -----------------------------
%%%%%%% track already terminated but was called again (due to interleaving) %%%%%%%%
if work.trackfinished{work.pvind}==1	% undo changes in work
    
    % changed 4/17/01 undo only if tracks are not to be continued
    if ~def.holdtrack
        work = worksave;
    else		% otherwise fake continue
        work.stepnum{work.pvind}=[work.stepnum{work.pvind} max(work.stepnum{work.pvind})+1];
    end
else
    %%%%%%% track termination on minmaxcount %%%%%%%%%
    if work.minmaxcount{work.pvind} * def.terminate >= def.endstop

        if ~def.allterminate
            work.trackfinished{work.pvind} = 1;						% was work.terminate before
            work.finishedcount = work.finishedcount + 1;

            % changed 4/17/01 remove track if hold off
            if ~def.holdtrack
                work.removetrack = work.pvind;
            else		% otherwise fake continue
                work.stepnum{work.pvind}=[work.stepnum{work.pvind} max(work.stepnum{work.pvind})+1];
            end

            %disp('track skipped');
            work.expvarrev{work.pvind} = ones(1,def.reversalnum) * 0.12345678;

            %afc_result;					% moved to main
            %work.writeresult =1;

            % else part from allterminate
        else
            for i = 1:def.interleavenum
                work.trackfinished{i} = 1; % was work.terminate before
                work.finishedcount = work.finishedcount + 1;
                work.removetrack = i;
                %disp('track skipped');
                work.expvarrev{i} = ones(1,def.reversalnum) * 0.12345678;
            end
        end


        %%%%%%% track termination on end of testlist %%%%%%%%%
    elseif work.stepnum{work.pvind}(end) >= def.maxiter
        work.trackfinished{work.pvind}=1;						% end of threshold run
        work.finishedcount = work.finishedcount + 1;

        % changed 4/17/01
        work.reachedendreversalnum{work.pvind} = work.stepnum{work.pvind}(end);

        % changed 4/17/01 remove track if hold off
        if ~def.holdtrack
            work.removetrack = work.pvind;
        else		% otherwise fake continue
            work.stepnum{work.pvind}=[work.stepnum{work.pvind} max(work.stepnum{work.pvind})+1];
        end

        %afc_result;					% moved to main
        %work.writeresult =1;

        %%%%%%% go on %%%%%%%%
    else
        %if work.stepnum{work.pvind}(end) ~= 1 | ((work.firstPresentation == 1) & (answer >=0.6))
            work.stepnum{work.pvind}=[work.stepnum{work.pvind} max(work.stepnum{work.pvind})+1];
            % work.firstPresentation = 0;
            
        %else  % then answer < 0.6 for the first presentation
            % work.stepnum{work.pvind}=[work.stepnum{work.pvind} max(work.stepnum{work.pvind})]; % then we will repeat the same sentence
            % if get(h,'Userdata')==1
            %   set(h,'Userdata',1);
            % end
        %end
        
        %set(h,'Userdata',-1);

    end
end % if already terminated


% clear response vector
work.OLSA_ResponseVector{work.pvind} = zeros(1,5);
work.OLSA_PercentCorrect{work.pvind} = 0;

% eof
