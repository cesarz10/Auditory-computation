% datread.m reads floats from ascii file. Tabs and space characters are ignored.
% Every line with a leading '%' or '//' is skipped. If all lines were skipped,
% datread returns an empty char.
%
% Usage: out = datread('filename');
%
% Typically used to read AFC data ".dat" files
%
% See also PARSEDAT, PSYDATM, PSYDATZ, ALLMEAN


function out=datread_custom(file);

fid=fopen(file,'r');
bla=1;%fgetl(fid);
out=[];

bla2=fgets(fid);

bFirst_read = 0;

while bla ~= -1
   bla=fgets(fid);
   bla2=sscanf(bla2,'%f');
   	
   	% added this 27/03/03 for MATLAB 6.5
   	if ~isempty(bla2)
        switch bFirst_read
            case 0
                out=[out bla2];
                bFirst_read = 1;
                N = length(bla2);
            case 1
                out = [out bla2(1:N)];
        end
	end
	bla2=bla;
end
fclose(fid);

out=out';
