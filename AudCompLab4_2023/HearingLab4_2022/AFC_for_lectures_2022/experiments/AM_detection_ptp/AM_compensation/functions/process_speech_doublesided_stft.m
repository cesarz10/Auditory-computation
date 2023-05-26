function [stim,env,g]=process_speech_doublesided_stft(stim,sp_par,poles_str,fsp,seg,height,prominence,distance)

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
bands=32; % bands for fb analysis: 80-band fb is used so far

clearvars stimp;
for n=1:size(stim,1)
    stima=fbanalysis1(stim(n,:),fsp,bands);
    stima=fbanalysis1(stim(n,:).*max(abs(stim(n,:)))./max(abs(fbsynthesis1(stima))),fsp,bands);
    g=ones(size(stima));
%     if ~isempty(poles_str)
% %         load(['parameters_ohc/' poles_str '.mat'],'pg','xdata');
% %         for cbi=1:bands
% % %             env=abs(stima(cbi,:)); %abs(envelope(stima(cbi,:),seg,'rms')); %
% %                         env=abs(envelope(stima(cbi,:)));
% %             %             env=env.*max(abs(stima(cbi,:)))./max(env);
% %             f=(1-(log10(1+exp((-pg(cbi,4)./(env))+pg(cbi,5))))./pg(cbi,3));
% %             g(cbi,:)=(pg(cbi,1)+pg(cbi,2).*sign(f).*abs(f).^0.8);
% % %             g(cbi,env<xdata(11))=1;
% %         end
%         
%         load(['parameters_ohc/' poles_str '.mat'],'g','xdata');
%         for cbi=1:bands
% %             if g(cbi,1)>g(cbi,2)
% %                 g(cbi,1)=0.5.*g(cbi,2);
% %             end
% %             for i=2:bands-1
% %                 if g(cbi,i)<0.1
% %                     g(cbi,i)=(g(cbi,i-1)+g(cbi,i+1))/2;
% %                 end
% %             end
%             gn(cbi,:)=smoothdata(g(cbi,:),'movmean',3);
%         end
%         clearvars g;
%         for cbi=1:bands
%             env=abs(envelope(stima(cbi,:)));
%             for i=1:length(env)
%                 g(cbi,i) = (interp1(log10(xdata),gn(cbi,:),log10(env(i)),'linear','extrap'));
%             end
%             g(cbi,g(cbi,:)<1)=1;
%         end
%     
%     % smooth g
% %     smbands=1;
% %     gs=g;
% %     for ti=1:size(stima,2)
% %         for cbi=smbands+1:size(stima,1)-smbands
% % %             gs(cbi,ti)=min(g(cbi-smbands:cbi+smbands,ti));
% %             if g(cbi-smbands:cbi+smbands,ti) > 1
% %                 gs(cbi-smbands,ti)=1;
% %                 gs(cbi+smbands,ti)=1;
% %             end
% %         end
% %         %             gs(:,ti)=smoothdata(g(:,ti),'movmean',5);
% %     end
% % %     gs(1:smbands,:)=g(1:smbands,:);
% % %     gs(bands-smbands+1:bands,:)=g(bands-smbands+1:bands,:);
% %     g=gs;
% 
%     end
    
    stimap=g.*stima;
    stimp(n,:)=fbsynthesis1(stimap);
end
stim=stimp;

%% FFT filterbank

% load(['parameters_ohc/' poles_str '.mat'],'g','gn','xdata','cf');
% gnt=gn;
% gn=20*log10(g);
% gn(26:end,:)=20*log10(gnt(26:end,:));
% clearvars g gnt;
% for i=1:size(gn,1)
%     gn(i,:)=smoothdata(gn(i,:),'movmean',3);
% end
% % gt=20*log10(xdata)+94;
% gt=[0:1:80];
load(['parameters_ohc/' poles_str '-stft.mat'],'g','xdata','cf','gt');
gth=0.5;
for i=1:size(g,1)
    thdata=g(i,:);
    xdata=(1:length(thdata))';
    gn(i,:)=interp1(xdata(thdata>gth),thdata(thdata>gth),xdata,'linear','extrap');
    gn(i,:)=smoothdata(gn(i,:),'movmean',3);
end
% gt=xdata;
clearvars g xdata thdata gth;

W=2048; % frame size & fft size
h=0.125; % overlap percentage 
wnd=hamming(W);
npad=W; % zero-padding size 
pads=(1-h)*W;

% ta=0.02; % attack time constant (20ms)
% tr=0.1; % release time constant (100ms)
spl_norm=0; %26; % dB to add to the spectrogram to convert the scale to dB-SPL
db_offset=200; % just an offset to perform logarithmic interpolation over negative values

% if the gain table provided is not described for each frequency bin,
% compute the rest of them via logarithmic interpolation
f=fsp/2*linspace(0,1,(W+npad)/2+1); % the correct one
erb=hz2erb(f);
cntr=0;
fstep=3;
for i=fstep:fstep:fstep*ceil(ceil(max(erb))./fstep)
    cntr=cntr+1;
    Wet(cntr,:) = (erb>=(i-fstep) & erb<i);
    cf_erb(cntr)=mean(f(Wet(cntr,:)));
end
bands=length(cf_erb);
C=[1;2.*ones(W-1,1);1];
We=Wet.*C';

% Gt=smoothdata(gn,'movmean',3);
Gt=gn;
clearvars gn;
for n=1:size(Gt,2)
    Gt(Gt(:,n)<1,n)=1;
    for k=1:bands
        gn(k,n) = (interp1(log10(cf),(Gt(:,n)),log10(cf_erb(k)),'linear','extrap'));
        if gn(k,n)<0
            gn(k,n)=0;
        end
    end
end
clear Gt f;

for l=1:size(stim,1)
    x=[zeros(pads,1);stim(l,:)';zeros(pads,1)];
    
    % segment signal for FFT (STFT)
    sp=fix(W.*h);
    N=fix((length(x)-W)/sp+1); %number of segments
    Index=(repmat(1:W,N,1)+repmat((0:(N-1))'*sp,1,W))';
    hw=repmat(wnd,1,N);
    xs=x(Index).*hw;
    
    % pad with zeros
    pad=zeros((npad/2),size(xs,2));
    x_pad=[pad' xs' pad'];
    xs=x_pad';
    clear x_pad;
    
    % STFT
    X=fft(xs); % whole STFT
    XPhase=angle(X(1:ceil(size(X,1)/2+1),:)); % Phase
    X1=abs(X(1:ceil(size(X,1)/2+1),:)); % Spectrogram
    X2=abs(X(1:ceil(size(X,1)/2+1),:)).^2; % Spectrogram
    
    frames=size(X,2);
    
%     aa=0; %exp(-1./(ta*fs)); % attack coefficient
%     ar=0; %exp(-1./(tr*fs)); % release coefficient
    for n=1:frames
        for k=1:bands
            Lin(k,n) = 10*log10(sum(abs(X1(:,n).*We(k,:)').^2))+spl_norm;
%             if Xlev(k,n) <= Lin(k,n-1)
%                 Lin(k,n) = aa*Lin(k,n-1)+(1-aa)*Xlev(k,n);
%             else
%                 Lin(k,n) = ar*Lin(k,n-1)+(1-ar)*Xlev(k,n);
%             end
            Gt(k,n) = (interp1(log10(gt+db_offset),(gn(k,:)),log10(Lin(k,n)+db_offset),'linear','extrap'));
            G(Wet(k,:),n)=Gt(k,n);
            % gain to be applied computed with logarithmic interpolation
            % extrapolation used for values outside the given SPL range.
        end
%         G(isinf(G(:,n)),n)=-100;
    end
    
    Y2=10*log10(X2)+G-0.6641; % processed STFT % 0.6641 is the correction for the STFT
    
    y=OverlapAdd2(10.^(Y2/20),XPhase,(W+npad),h*W);
    y=y(pads+npad/2+1:end-pads-npad/2); % processing signal in time domain
    stim(l,:)=y(1:size(stim,2))*2*h;
end

end

function ReconstructedSignal=OverlapAdd2(XNEW,yphase,windowLen,ShiftLen,ifftsym)

%Y=OverlapAdd(X,A,W,S);
%Y is the signal reconstructed signal from its spectrogram. X is a matrix
%with each column being the fft of a segment of signal. A is the phase
%angle of the spectrum which should have the same dimension as X. if it is
%not given the phase angle of X is used which in the case of real values is
%zero (assuming that its the magnitude). W is the window length of time
%domain segments if not given the length is assumed to be twice as long as
%fft window length. S is the shift length of the segmentation process ( for
%example in the case of non overlapping signals it is equal to W and in the
%case of %50 overlap is equal to W/2. if not givven W/2 is used. Y is the
%reconstructed time domain signal.
%Sep-04
%Esfandiar Zavarehei

if nargin<2
    yphase=angle(XNEW);
end
if nargin<3
    windowLen=size(XNEW,1)*2;
end
if nargin<4
    ShiftLen=windowLen/2;
end
if fix(ShiftLen)~=ShiftLen
    ShiftLen=fix(ShiftLen);
    disp('The shift length have to be an integer as it is the number of samples.')
    disp(['shift length is fixed to ' num2str(ShiftLen)])
end

if nargin<5
    ifftsym='symmetric';
end

[~, FrameNum]=size(XNEW);

Spec=XNEW.*exp(1i*yphase);

% if mod(windowLen,2) %if FreqResol is odd
%     Spec=[Spec;flipud(conj(Spec(2:end,:)))];
% else
%     Spec=[Spec;flipud(conj(Spec(2:end-1,:)))];
% end
% sig=zeros((FrameNum-1)*ShiftLen+windowLen,1);
sig=zeros((FrameNum)*ShiftLen+windowLen,1);
% weight=sig;
for i=1:FrameNum
    start=(i-1)*ShiftLen+1;
    spec=Spec(:,i);
    ibuffer=real(ifft(spec,windowLen,ifftsym));
    sig(start:start+windowLen-1)=sig(start:start+windowLen-1)+ibuffer;
end
ReconstructedSignal=sig;

% if max(abs(ReconstructedSignal))>1
%     ReconstructedSignal=ReconstructedSignal./max(abs(ReconstructedSignal));
% end

end
