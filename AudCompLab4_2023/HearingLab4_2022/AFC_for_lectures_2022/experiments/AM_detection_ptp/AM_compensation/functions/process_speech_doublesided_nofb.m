function [stim,env,g]=process_speech_doublesided_nofb(stim,sp_par,spl70_flag,seg,height,prominence,distance)

if nargin<3
    spl70_flag=0;
end
if nargin<4
    seg=0;
end
if nargin<5
    height=0;
end
if nargin<6
    prominence=1.35;
end
if nargin<7
    distance=1;
end

if iscolumn(stim)
    stim=stim';
    inv_flag=1;
else
    inv_flag=0;
end

% stimo=stim;
if seg
    for n=1:size(stim,1)
        stimp(n,:)=abs(envelope(stim(n,:),seg,'rms'));
        if exist('smoothdata')
            stimp(n,:)=smoothdata(stimp(n,:),'movmean',seg);
        else
            stimp(n,:)=smooth(stimp(n,:),'moving',seg);
        end
    end
    
    peaks=1;
    minima=[];
    cntr=1;
    for i=3:size(stimp,2) % two-sided prominence?
        if stimp(n,i-1)>height & stimp(n,i-1)>stimp(n,i) & ...
                stimp(n,i-1)>stimp(n,i-2) & i-1>=peaks(cntr)+distance
            [~,pt]=min(stimp(n,peaks(cntr):i-1));
            if stimp(n,i-1)>prominence*stimp(n,peaks(cntr)+pt-1)
                minima=[minima peaks(cntr)+pt-1];
                cntr=cntr+1;
                peaks(cntr)=i-1;
            end
        end
    end
else
    for n=1:size(stim,1)
        stimp(n,:)=abs(envelope(stim(n,:)));
    end
end

minima(1)=1;
if minima(end)~=size(stimp,2)
    minima=[minima,size(stimp,2)];
end

if seg
    peaks=peaks(2:end);
    minimat=1;
    for i=1:length(peaks)
        if stimp(n,peaks(i))>=prominence*stimp(n,minima(i+1))
            minimat=[minimat minima(i+1)];
        end
    end
    if minimat(end)~=size(stimp,2)
        minimat=[minimat size(stimp,2)];
    end
    minima=minimat;
end

% stimp=abs(stim);

%% synaptopathy compensation
if any(sp_par~=[13,3,3])
    scale=1;
    fname = sprintf('parameters_cs/cs_%d_%d_%d.mat',sp_par(1),sp_par(2),sp_par(3));
    load(fname,'b70','pg','pc','gf70','xdatam');
    
    for j=1:length(minima)-1
        index=minima(j):minima(j+1);
        stim_seg=stim(:,index);
        
        for n=1:size(stim_seg,1)
            xdata(n)=max(abs(stim_seg(n,:)));
            gf(n)=xdata(n)./gf70;
            if gf(n)==0
                gf(n)=1e-20;
            end
        end
        
        if exist('pg','var') & any(pg ~= 0)
            for n=1:size(stim_seg,1)
                if spl70_flag % process as if it is 70 dB-SPL
                    if xdata(n)/gf(n)<xdatam
                        gl(n)=10.^(pg(1)+exp(-pg(2).*log10(xdata(n)./gf(n))-pg(3)));
                    else
                        gl(n)=1;
                    end
                else
                    if xdata(n)<xdatam
                        gl(n)=10.^(pg(1)+exp(-pg(2).*log10(xdata(n))-pg(3)));
                    else
                        gl(n)=1;
                    end
                end
            end
        else
            error('gl estimate not found')
        end
        if exist('pc','var') && any(pc ~= 0)
            if spl70_flag % process as if it is 70 dB-SPL
                c=pc(1)+pc(2)./(1+exp(-pc(3).*log10(xdata./gf)-pc(4)));
            else
                c=pc(1)+pc(2)./(1+exp(-pc(3).*log10(xdata)-pc(4)));
            end
        else
            error('c estimate not found')
        end
        
        for n=1:size(stim_seg,1)
            env(n,index)=stimp(n,index).*xdata(n)./max(abs(stimp(n,index)));
            if spl70_flag % preserve the same peak
                gl(n)=1./max(1./(1+exp(scale.*(-b70.*env(n,index)./gf(n)+c(n)))));
            end
            g(n,index)=gl(n)./(1+exp(scale.*(-b70.*env(n,index)./gf(n)+c(n))));
        end
    end
    stim=g.*stim;
end

%% ohc loss compensation
% bands=80; % bands for fb analysis: 80-band fb is used so far
%
% clearvars stimp;
% for n=1:size(stim,1)
%     stima=fbanalysis1(stim(n,:),fsp,bands);
%     stima=fbanalysis1(stim(n,:).*max(abs(stim(n,:)))./max(abs(fbsynthesis1(stima))),fsp,bands);
%     g=ones(size(stima));
%     if ~isempty(poles_str)
%         load(['parameters_ohc/' poles_str '.mat'],'pg');
%         for cbi=1:bands
% %             env=abs(stima(cbi,:)); %abs(envelope(stima(cbi,:),seg,'rms')); %
%             env=abs(envelope(stima(cbi,:)));
% %             env=env.*max(abs(stima(cbi,:)))./max(env);
%             f=(1-(log10(1+exp((-pg(cbi,4)./(env))+pg(cbi,5))))./pg(cbi,3));
%             g(cbi,:)=(pg(cbi,1)+pg(cbi,2).*sign(f).*abs(f).^0.8);
%         end
%
%         % smooth g
%         for ti=1:size(stima,2)
%             for cbi=2:size(stima,1)-1
%                 gs(cbi,ti)=min(g(cbi-1:cbi+1,ti));
%             end
% %             gs(:,ti)=smoothdata(g(:,ti),'movmean',10);
%         end
%         gs(1,:)=g(1,:);
%         gs(bands,:)=g(end,:);
%         g=gs;
%
%     end
%     stimap=g.*stima;
%     stimp(n,:)=fbsynthesis1(stimap);
% end
% stim=stimp;

if inv_flag
    stim=stim';
end
