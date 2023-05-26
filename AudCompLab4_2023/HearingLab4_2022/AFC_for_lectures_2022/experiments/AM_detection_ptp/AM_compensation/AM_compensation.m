% Cochlear Synaptopathy and Outer Hair Cell Loss 
% Compensation for AM stimuli
% 
% Fotis Drakopoulos, Ghent University
%
% The AM_compensation folder needs to be inside the 
% 'Verhulstetal2018Model-master' folder
% 
% Define the kind of stimuli in the "declare parameters" section
%
% All synaptopathy profiles are included.
% Included OHC-loss profiles: 
% 'Flat00_Slope05','Flat00_Slope10', 'Flat00_Slope35','Flat35','Slope35'
%

clear;
close all;

save_dir=''; % '/media/robspear/fotisdr/hi_compensation/outputs'
model_dir='/home/fdrakopo/Documents/Verhulstetal2018Model'; 
% model_dir='/media/fdrakopo/EA05C216/Verhulstetal2018Model'; 
% model_dir='/media/robspear/fotisdr/Verhulstetal2018Model';

% model_dir='E:\Windows\fdrakopo\Downloads\hi_compensation\Verhulstetal2018Model-master'; 

mydir=pwd;
addpath(genpath(mydir));
cd(model_dir);

%% declare parameters
cs_arr=[7,0,0];        % synaptopathy profile in [HSR,MSR,LSR] array
poles_str='';           % ohc loss profile name (empty string for NH)

L=[30 50 70 90];        % stimuli levels (one value or array of levels)
f=4000;                 % carrier frequency
fm=120;                 % modulation frequency
m=1;                    % modulation depth percentage (1 is 100%, 0 is 0%)

dur=125e-3;              % duration of stimuli
fsp=44100;              % original stimuli frequency (should be 44100 for proper fb analysis)
fs=100e3;               % frequency for the model
model_fc = 'abr';

%% processing

ts=(0:1/fsp:dur);
p0=2e-5;
bands=32;
fl=50;
cf = flipud(ERBSpace(fl, fsp/2, bands));

for n=1:numel(L)
    stim(n,:)=(1+m.*cos(2*pi*fm*ts+pi)).*sin(2*pi*f*ts);
    stim(n,:)=p0*10^(L(n)/20)*stim(n,:)./rms(stim(n,:));
end
for n=1:numel(L)
%     stimt=fbsynthesis1(fbanalysis1(stim(n,:),fsp,bands,fl,1)); % account for fb modification
%     stimt=stimt.*max(abs(stim(n,:)))./max(abs(stimt));
    stimt=stim(n,:);
    if fs ~= fsp
        stimr(n,:)=resample(stimt,fs,fsp);
    else
        stimr(n,:)=stimt;
    end
end
clearvars stimt;

% NH response
sheraP=load('StartingPoles.dat');
output=model2018(stimr,fs,model_fc,1,'evihmlbw',1,sheraP,0.05,'vel',13,3,3,1,[pwd(),'/']);
fs_an=output.fs_an;
t=[0:size(output(1).anfH,1)-1]/fs_an;
for n=1:numel(L)
    ANn(n,:,:)=output(n).an_summed;
    W1n(n,:)=output(n).w1;
%     W5n(n,:)=output(n).w5;
    vn=output(n).v;
end

if ~isempty(poles_str)
    sheraP=load(['Poles/' poles_str '/StartingPoles.dat']);
end
output=model2018(stimr,fs,model_fc,1,'evihmlbw',1,sheraP,0.05,'vel',cs_arr(1),cs_arr(2),cs_arr(3),1,[pwd(),'/']);
for n=1:numel(L)
    ANh(n,:,:)=output(n).an_summed;
    W1h(n,:)=output(n).w1;
%     W5h(n,:)=output(n).w5;
end
clearvars output fs_an;

stimp=process_speech_doublesided_nofb(stim,cs_arr);
for n=1:numel(L)
    if fs ~= fsp
        stimpr(n,:)=resample(stimp(n,:),fs,fsp);
    else
        stimpr(n,:)=stimp(n,:);
    end
end
stimpr=stimpr(:,1:size(stimr,2));

output=model2018(stimpr,fs,model_fc,1,'evihmlbw',1,sheraP,0.05,'vel',cs_arr(1),cs_arr(2),cs_arr(3),1,[pwd(),'/']);
for n=1:numel(L)
    ANp(n,:,:)=output(n).an_summed;
    W1p(n,:)=output(n).w1;
%     W5p(n,:)=output(n).w5;
    vp=output(n).v;
end

stimp1=process_speech_doublesided_nofb(stim,cs_arr); %process_speech_doublesided_stft(stim,cs_arr,poles_str,fsp);
for n=1:numel(L)
    if fs ~= fsp
        stimpr1(n,:)=resample(stimp1(n,:),fs,fsp);
    else
        stimpr1(n,:)=stimp1(n,:);
    end
end
stimpr1=stimpr1(:,1:size(stimr,2));

output=model2018(stimpr1,fs,model_fc,1,'evihmlbw',1,sheraP,0.05,'vel',cs_arr(1),cs_arr(2),cs_arr(3),1,[pwd(),'/']);
for n=1:numel(L)
%     ANp(n,:,:)=output(n).an_summed;
    W1p1(n,:)=output(n).w1;
%     W5p(n,:)=output(n).w5;
end

%% plotting
cd(mydir);

set(0,'DefaultTextInterpreter','none');
if isempty(poles_str)
    poles_str='NH';
end
headtag=[poles_str ' - ' regexprep(num2str(cs_arr),'\s+',',')];
if numel(L)>1
    spi(1)=2;
else
    spi(1)=1;
end
spi(2)=ceil(numel(L)/spi(1));

tst=(0:1/fsp:(length(stim)-1)/fsp);
figure
for n=1:numel(L)
subplot(spi(1),spi(2),numel(L)-n+1)
plot(1000*tst,stim(n,:),'DisplayName','Original stim')
hold on
plot(1000*tst,stimp(n,:),'DisplayName','Processed stim')
title(['L = ',num2str(L(n)),' dB-SPL'])
xlabel('Time [ms]'),ylabel('Amplitude')
legend('Location','Southeast')
legend('boxoff')
xlim([50 66.66])
% xlim([dur/2*1000 dur*1000])
% xlim([0 2000])
grid on
end
% p=mtit(headtag,'xoff',0,'yoff',.025, 'fontsize',14,'color',[0.1 0.1 0.5]);

% h=gcf;
% set(h,'PaperPositionMode','auto');         
% set(h,'PaperOrientation','landscape');
% set(h,'Position',[50 50 1200 800]);
% print(gcf, '-dpdf', [headtag '-1.pdf'])

figure
for n=1:numel(L)
subplot(spi(1),spi(2),numel(L)-n+1)
plot(1000*t,W1n(n,:),'DisplayName','NH - NoS')
hold on
plot(1000*t,W1h(n,:),'DisplayName','HI - LS & Slope35') % & Slope35
hold on
plot(1000*t,W1p(n,:),'DisplayName','HI - LS & Slope35 - Processed') % & Slope35
if exist('W1p1','var')
    hold on
    plot(1000*t,W1p1(n,:),'DisplayName','HI - LS & Slope35 - Processed1') % & Slope35
end
title(['Summed AN response - L = ',num2str(L(n)),' dB-SPL'])
xlabel('Time [ms]'),ylabel('AN population response [Volts]')
legend('Location','Southeast')
legend('boxoff')
l = legend();
set(l, 'Interpreter', 'none')
% xlim([100 116.67])
% xlim([dur/2*1000 dur*1000-0.2])
% xlim([0 2000])
grid on
end
% p=mtit(headtag,'xoff',0,'yoff',.025, 'fontsize',14,'color',[0.1 0.1 0.5]);

h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'Position',[50 50 600 800]);
print(gcf, '-dpdf', [headtag '-2.pdf'])
