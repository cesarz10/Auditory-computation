function [stim,env,g]=process_am_tone_nofb(stim,sp_par,spl70_flag,seg,height,prominence,distance)

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

xdata_min(n)=min(stimp(n,100:end-100));
xdata_max(n)=max(stimp(n,:)-xdata_min(n));

stimp(n,:)=(stimp(n,:)-xdata_min(n))./xdata_max(n);
stimp(n,stimp(n,:)<0)=0;

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
%     if sp_par(1)<10
%         scale=0.5;
%     else
        scale=1;
%     end
    fname = sprintf('cs_parameters/cs_%d_%d_%d.mat',sp_par(1),sp_par(2),sp_par(3));
    load(fname,'b70','pg','pc','gf70','xdatam');
    
    for j=1:length(minima)-1
        index=minima(j):minima(j+1);
        stim_seg=stim(:,index);
        
        for n=1:size(stim_seg,1)
            xdata(n)=max(abs(stim_seg(n,:)));
            gf(n)=xdata(n)./gf70;
        end
        
        if spl70_flag % process as if it is 70 dB-SPL
            if exist('pgi','var') & pgi ~= 0
                for n=1:size(stim_seg,1)
                    if xdata(n)<pgi
                        gl(n)=pg(2,1)+exp(-pg(2,2).*log10(xdata(n)./gf(n))-pg(2,3));
                    else
                        gl(n)=pg(1,1)+exp(-pg(1,2).*log10(xdata(n)./gf(n))-pg(1,3));
                    end
                end
            elseif exist('pg','var') & any(pg ~= 0)
%                 gl=pg(1)+exp(-pg(2).*log10(xdata./gf(n))-pg(3));
                for n=1:size(stim_seg,1)
                    if xdata(n)/gf(n)<xdatam
                        gl(n)=10.^(pg(1)+exp(-pg(2).*log10(xdata(n)./gf(n))-pg(3)));
                    else
                        gl(n)=1;
                    end
                end
            else
%                 load(fname,'gl');
%                 if length(gl)>71
%                     gl=gl(71);
%                 end
                error('gl estimate not found')
            end
            if exist('pc','var') && any(pc ~= 0)
                c=pc(1)+pc(2)./(1+exp(-pc(3).*log10(xdata./gf)-pc(4)));
            else
%                 load(fname,'c70');
%                 c=c70;
                error('c estimate not found')
            end
        else
            if exist('pgi','var') & pgi ~= 0
                for n=1:size(stim_seg,1)
                    if xdata(n)<pgi
                        gl(n)=pg(2,1)+exp(-pg(2,2).*log10(xdata(n))-pg(2,3));
                    else
                        gl(n)=pg(1,1)+exp(-pg(1,2).*log10(xdata(n))-pg(1,3));
                    end
                end
                load(fname,'gl');
            elseif exist('pg','var') & any(pg ~= 0)
%                 gl=pg(1)+exp(-pg(2).*log10(xdata)-pg(3));
                for n=1:size(stim_seg,1)
                    if xdata(n)<xdatam
                        gl(n)=10.^(pg(1)+exp(-pg(2).*log10(xdata(n))-pg(3)));
                    else
                        gl(n)=1;
                    end
                end
            else
%                 load(fname,'gl');
%                 if length(gl)>71
% %                     gl=gl(71);
%                     gl=gl(21:10:91);
%                 end
                error('gl estimate not found')
            end
            if exist('pc','var') & any(pc ~= 0)
                c=pc(1)+pc(2)./(1+exp(-pc(3).*log10(xdata)-pc(4)));
            else
%                 load(fname,'c');
% %                 c=c(31:10:91);
%                 c=c70;
                error('c estimate not found')
            end
        end        
        
        for n=1:size(stim_seg,1)
            env(n,index)=stimp(n,index).*xdata(n)./max(abs(stimp(n,index)));
            gl(n)=1./max(1./(1+exp(scale.*(-b70.*env(n,index)./gf(n)+c(n)))));
            g(n,index)=gl(n)./(1+exp(scale.*(-b70.*env(n,index)./gf(n)+c(n))));
            %             g(n,:)=gl(n)./(1+exp(scale.*(-b70.*envelope(stim(n,:))./gf(n)+c(n))));
        end
    end
    
    env1=g.*env.*xdata_max./xdata+xdata_min;
%     env1=env1.*max(abs(env))./max(abs(env1));
    
    g=env1./envelope(stimo);
    
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
