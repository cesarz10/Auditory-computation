function afc_SRT_extraction(subject,dir)
% Script to extract the SRT results from AFC toolbox
% subject: subject name is given as a character vector
% dir: directory to load the results from - if not defined
% it reads from the ./results/subject subdirectory
%
% Fotis Drakopoulos, UGent 10/2019
% Updated 9/2020
% 

% clearvars -except results1;
close all;

if nargin < 2
    dir = ['experiments' filesep 'SRT_task' filesep 'results' filesep ];
end

% subject='FD_Subject08';
condsp={'exampleVlMatrix'};

cond='varnoise20';

for testi=1:length(condsp)
    test=condsp{testi};
    
fileID = fopen([dir test '_protocol.pro'],'r');
S = textscan(fileID,'%s','delimiter','\n') ;
S=S{1};
fclose(fileID);

cntr=0;
condition=[];
for i=1:length(S)
    if strcmp(S{i},[test '_' subject '_00'])
        cntr=cntr+1;
        list{cntr}=str2num(S{i+1});
        condition{cntr}='BB';
        condition_speech{cntr}=cond;
        results{cntr}=str2num(S{i+4});
    end
        if strcmp(S{i},[test '_' subject '_01'])
        cntr=cntr+1;
        list{cntr}=str2num(S{i+1});
        condition{cntr}='LP';
        condition_speech{cntr}=cond;
        results{cntr}=str2num(S{i+4});
        end
    if strcmp(S{i},[test '_' subject '_02'])
        cntr=cntr+1;
        list{cntr}=str2num(S{i+1});
        condition{cntr}='HP';
        condition_speech{cntr}=cond;
        results{cntr}=str2num(S{i+4});
    end
end

figure
for i=1:length(results)
    legend('autoupdate','on')
    plot(results{i},'-s','DisplayName',condition{i})
    hold on
    if i ~= length(results)
        legend('autoupdate','off')
    end
    SRT(i) = mean(results{i}(end-5:end));
    plot(length(results{i})-5:length(results{i}),mean(results{i}(end-5:end))*ones(1,6),'k--','DisplayName','Mean SRT')
    hold on
    text(length(results{i}),SRT(i),num2str(SRT(i),'%1.2f'))
end
xlim([0 22])
% ylim([-12 6])
xlabel('Trials')
ylabel('SNR [dB]')
grid on
legend('boxoff')
title('VlMatrix SRT tracking '+ string(subject))


fprintf('Mean of last 6 trials - %s - %s is \n',test,subject);

for i=1:length(condition)
    fprintf('Condition %s : %1.2f \n',condition{i},SRT(i));
end

end