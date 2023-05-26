function afc_AM_extraction(subject,dir)
% Script to extract the AM modulation threshold results from AFC toolbox
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
    dir = ['experiments' filesep 'AM_detection_ptp' filesep 'results' filesep ];
end

% subject='FD_Subject08';
condsp={'AM_detection_ptp'};

cond='AM_detection_ptp';

for testi=1:length(condsp)
    test=condsp{testi};
    
fileID = fopen([dir test '_protocol.pro'],'r');
S = textscan(fileID,'%s','delimiter','\n') ;
S=S{1};
fclose(fileID);

cntr=0;
condition=[];
for i=1:length(S)
    if strcmp(S{i},[test '_' subject])
        cntr=cntr+1;
        list{cntr}=str2num(S{i+1}); % mod freq
        condition{cntr}='AM_detection';
        condition_speech{cntr}=cond;
        results{cntr}=str2num(S{i+5});
        mod_index{cntr}=str2num(S{i+7});
    end
end

mod_index = mod_index{1}(1);

figure
for i=1:length(results)
    legend('autoupdate','on','interpreter','none','location','south')
    plot(results{i},'-s','DisplayName','AM detection')
    hold on
    if i ~= length(results)
        legend('autoupdate','off','interpreter','none','location','south')
    end
    plot(length(results{i})-25:length(results{i}),mod_index*ones(1,26),'r--','DisplayName','Mean modulation index')
    hold on
    text(length(results{i})+3,mod_index,num2str(mod_index,'%1.2f'),'color','r');
end
xlabel('Trials')
ylabel('Modulation index [dB]')
grid on
legend('boxoff','location','south')
title('AM detection threshold tracking '+ string(subject))

end