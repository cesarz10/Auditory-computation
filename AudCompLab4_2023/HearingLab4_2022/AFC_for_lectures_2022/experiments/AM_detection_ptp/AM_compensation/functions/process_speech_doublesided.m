function [stim,env,g]=process_speech_doublesided(stim,sp_par,poles_str,fsp,seg,height,prominence,distance)

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
    if sp_par(1)<10
        scale=0.5;
    else
        scale=1;
    end
%     fname = sprintf('parameters_cs/cs_%d_%d_%d.mat',sp_par(1),sp_par(2),sp_par(3));
    fname = sprintf('Synaptopathy-new/cs_%d_%d_%d-1.mat',sp_par(1),sp_par(2),sp_par(3));
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
bands=32; % bands for fb analysis: 80-band fb is used so far

clearvars stimp;
for n=1:size(stim,1)
    stima=fbanalysis1(stim(n,:),fsp,bands);
    stima=fbanalysis1(stim(n,:).*max(abs(stim(n,:)))./max(abs(fbsynthesis1(stima))),fsp,bands);
    g=ones(size(stima));
    if ~isempty(poles_str)
%         load(['parameters_ohc/' poles_str '.mat'],'pg','xdata');
%         for cbi=1:bands
% %             env=abs(stima(cbi,:)); %abs(envelope(stima(cbi,:),seg,'rms')); %
%                         env=abs(envelope(stima(cbi,:)));
%             %             env=env.*max(abs(stima(cbi,:)))./max(env);
%             f=(1-(log10(1+exp((-pg(cbi,4)./(env))+pg(cbi,5))))./pg(cbi,3));
%             g(cbi,:)=(pg(cbi,1)+pg(cbi,2).*sign(f).*abs(f).^0.8);
% %             g(cbi,env<xdata(11))=1;
%         end
        
        load(['parameters_ohc/' poles_str '-32.mat'],'g','xdata');
        gth=0.25;
        for cbi=1:bands
            thdata=g(cbi,:);
            xd=(1:length(thdata))';
            gn(cbi,:)=interp1(xd(thdata>gth),thdata(thdata>gth),xd,'linear','extrap');
            gn(cbi,:)=smoothdata(gn(cbi,:),'movmean',3);
        end
        clearvars g xd thdata gth;
        for cbi=1:bands
            env=abs(envelope(stima(cbi,:)));
            for i=1:length(env)
                g(cbi,i) = (interp1(log10(xdata),gn(cbi,:),log10(env(i)),'linear','extrap'));
            end
%             g(cbi,g(cbi,:)<1)=1;
        end
    
    % smooth g
%     smbands=1;
%     gs=g;
%     for ti=1:size(stima,2)
%         for cbi=smbands+1:size(stima,1)-smbands
% %             gs(cbi,ti)=min(g(cbi-smbands:cbi+smbands,ti));
%             if g(cbi-smbands:cbi+smbands,ti) > 1
%                 gs(cbi-smbands,ti)=1;
%                 gs(cbi+smbands,ti)=1;
%             end
%         end
%         %             gs(:,ti)=smoothdata(g(:,ti),'movmean',5);
%     end
% %     gs(1:smbands,:)=g(1:smbands,:);
% %     gs(bands-smbands+1:bands,:)=g(bands-smbands+1:bands,:);
%     g=gs;

    end
    
    stimap=g.*stima;
    stimp(n,:)=fbsynthesis1(stimap);
end
stim=stimp;
