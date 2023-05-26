function [signaln]=fbanalysis1(signal,fs,bands,fl,w)

if nargin<4
    fl=50;
end
if nargin<5
    w=1;
end

fcoefs = MakeERBFilters(fs,bands,fl,w);
fcoefs = flipud(fcoefs);
signaln = ERBFilterBank(signal,fcoefs);
signaln = fliplr(ERBFilterBank(fliplr(signaln),fcoefs));

end