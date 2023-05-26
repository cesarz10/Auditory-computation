function [stim,env,g]=process_speech(stim,sp_par,poles_str,fsp,seg,height,prominence,distance)

if nargin<5
    seg=0;
end
if nargin<6
    height=0;
end
if nargin<7
    prominence=1.35;
end
if nargin<8
    distance=1;
end

if iscolumn(stim)
    stim=stim';
end

stimo=stim;
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

% stimp=abs(stim);

%% synaptopathy compensation
if any(sp_par~=[13,3,3])
    if sp_par(1)<11
        scale=0.5;
    else
        scale=1;
    end
    fname = sprintf('parameters_cs/cs_%d_%d_%d.mat',sp_par(1),sp_par(2),sp_par(3));
    load(fname,'b70','pg','pc','gf70','pgi');
    
    for j=1:length(minima)-1
        index=minima(j):minima(j+1);
        stim_seg=stim(:,index);
        for n=1:size(stim_seg,1)
            xdata(n)=max(abs(stim_seg(n,:)));
            gf(n)=xdata(n)./gf70;
        end
%         xdata=gf70;
%         gf=1;
        
        if exist('pgi','var')
            for n=1:size(stim_seg,1)
                if xdata(n)<pgi
                    gl(n)=pg(2,1)+exp(-pg(2,2).*log10(xdata(n))-pg(2,3));
                else
                    gl(n)=pg(1,1)+exp(-pg(1,2).*log10(xdata(n))-pg(1,3));
                end
            end
        else
            gl=pg(1)+exp(-pg(2).*log10(xdata)-pg(3));
        end
        c=pc(1)+pc(2)./(1+exp(-pc(3).*log10(xdata)-pc(4)));
        
        for n=1:size(stim_seg,1)
            env(n,index)=stimp(n,index).*xdata(n)./max(abs(stimp(n,index)));
            g(n,index)=gl(n)./(1+exp(scale.*(-b70.*env(n,index)./gf(n)+c(n))));
%             g(n,:)=gl(n)./(1+exp(scale.*(-b70.*envelope(stim(n,:))./gf(n)+c(n))));
        end
    end
    stim=g.*stim;
end

%% ohc loss compensation
bands=80; % bands for fb analysis: 80-band fb is used so far

clearvars stimp;
for n=1:size(stim,1)
    stima=fbanalysis1(stim(n,:),fsp,bands);
    stima=fbanalysis1(stim(n,:).*max(abs(stim(n,:)))./max(abs(fbsynthesis1(stima))),fsp,bands);
    g=ones(size(stima));
    if ~isempty(poles_str)
        load(['parameters_ohc/' poles_str '.mat'],'pg');
        for cbi=1:bands
            env=envelope(stima(cbi,:)); % abs?
            env(env<0)=abs(env(env<0));
            f=(1-(log10(1+exp((-pg(cbi,4)./(env))+pg(cbi,5))))./pg(cbi,3));
            g(cbi,:)=(pg(cbi,1)+pg(cbi,2).*sign(f).*abs(f).^0.8);
        end
        
        % smooth g
%         for ti=1:size(stima,2)
%             for cbi=2:size(stima,1)-1
%                 gs(cbi,ti)=min(g(cbi-1:cbi+1,ti));
%             end
% %             gs(:,ti)=smoothdata(g(:,ti),'movmean',10);
%         end
%         gs(1,:)=g(1,:);
%         gs(bands,:)=g(end,:);
%         g=gs;
        
    end
    stimap=g.*stima;
    stimp(n,:)=fbsynthesis1(stimap);
end
stim=stimp;
