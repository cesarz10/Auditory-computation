function av_data = allmedian(data,percL,percU);
% function av_data = allmedian(data,percL,percU);
%
%allmean(data)
% Sorts data and outputs the means and stds.
% The 4 output columns are: 
%   Col Nr.
%   1       : Parameter value 
%   2       : Median result
%   3       : Interpercentile range (distance between lower and higher percentile)
%   4       : Nr. of presentations
%   5       : Lower percentile (default: percentile 25)
%   6       : Higher percentile (default: percentile 75)
%
% Adapted from allmean.m by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    percU = 75;
end

if nargin < 2
    percL = 25;
end

if nargin < 1
   help allmedian
   return
end   

n_datapoints = size(data,1);

data = sortrows(data);

temp_data = data(1,:);
av_data = [];

for i = 2:n_datapoints
   if data(i,1) == data(i-1,1)
      temp_data = [temp_data;data(i,:)];
   else
      n_points = size(temp_data,1);
      dataL = prctile(temp_data(:,2),percL);
      dataU = prctile(temp_data(:,2),percU);
      av_data = [av_data; data(i-1,1) prctile(temp_data(:,2),50) dataU-dataL n_points dataL dataU];
      temp_data = data(i,:);
   end
end
n_points = size(temp_data,1);
dataL = prctile(temp_data(:,2),percL);
dataU = prctile(temp_data(:,2),percU);
av_data = [av_data; data(n_datapoints,1) prctile(temp_data(:,2),50) dataU-dataL n_points dataL dataU];
