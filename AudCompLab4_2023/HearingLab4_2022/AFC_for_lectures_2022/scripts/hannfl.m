function h = hannfl(len,h1len,h2len);
% function h = hannfl(len,h1len,h2len);
%
%   len, h1len, h2len are all in samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if len > 0
    h = ones(len,1);

        switch h1len
        case 0
            otherwise
                % h(1:h1len)=(1-cos(pi/h1len*[0:h1len-1]))/2;
                hvar = transpose(0:h1len-1);
                if ~isempty(hvar)
                    h(1:h1len)=(1-cos(pi/h1len*hvar))/2;
                end
        end

        switch h2len
        case 0
            otherwise
                % h(end-h2len+1:end)=(1+cos(pi/h2len*[1:h2len]))/2;
                hvar = transpose(1:h2len);
                if ~isempty(hvar)
                    h(end-h2len+1:end,1)=(1+cos(pi/h2len*hvar))/2;
                end
        end

else
end


