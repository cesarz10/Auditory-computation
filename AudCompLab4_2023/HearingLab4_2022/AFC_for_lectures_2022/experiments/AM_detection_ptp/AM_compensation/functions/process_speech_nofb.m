function [stim,env,g]=process_speech_nofb(stim,sp_par,spl70_flag,seg,height,prominence,distance)

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
stimp=zeros(size(stim));
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

%% synaptopathy compensation
if any(sp_par~=[13,3,3])
    scale=1;
    fname = sprintf('parameters_cs/cs_%d_%d_%d.mat',sp_par(1),sp_par(2),sp_par(3));
    load(fname,'b70','pg','pc','gf70','xdatam');
    
    xdata=zeros(size(stim,1));
    gf=zeros(size(stim,1));
    gl=zeros(size(stim,1));
    env=zeros(size(stim));
    g=zeros(size(stim));
    
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

if inv_flag
    stim=stim';
end
