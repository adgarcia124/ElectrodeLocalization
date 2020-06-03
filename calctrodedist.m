function [trodedist,datetable]=calctrodedist(trodenum,querydates,xlsfile)
%CALCTRODEDIST - Determine the electrode distance at each queried date.
%
%Input:
%   trodenum - Nx1 vector. Numerical ID of electrode (repeats are fine)
%   querydates - Nx1 array of strings. Dates for each electrode
%   xlsfile - string or Mx1 struct.  Filepath of excel file containing
%       turning log for each electrode, OR datetable struct (see below).
%
%Output:
%   trodedist - Nx1 vector.  Distance of each electrode at queried date
%   datetable - Mx1 struct.  Distance by date table for each electrode.
%
%
%Jon Rueckemann 2018

if isempty(querydates)
    querydates={date}; %populate with current date
elseif ~iscell(querydates)
    querydates={querydates};
end

if numel(querydates)>1
    assert(numel(trodenum)==numel(querydates),...
        'The same number of elements must be in first and second inputs.');
end

if nargin<3
    xlsfile=[];
end

%Create table of dates and distances for each electrode from turning log
if ~isstruct(xlsfile)
    datetable=trodelog(xlsfile);
else
    datetable=xlsfile;
end
tablechannel=cellfun(@str2num,{datetable.Channel});

%Populate electrode and date information if incomplete
if isempty(trodenum)
    trodenum=tablechannel;
end
if numel(trodenum)~=numel(querydates)
    querydates=repmat(querydates,numel(trodenum),1);
end

%Process each electrode entry
trodedist=nan(size(trodenum));
for t=1:numel(trodenum)
    %Load the relevant date and distance information from the table
    idx=find(trodenum(t)==tablechannel);
    
    %Find the most recent distance measurement for each query day
    if numel(datetable(idx).Dist)>=2
        trodedist(t)=interp1(datetable(idx).NumericDates,...
            datetable(idx).Dist,datenum(querydates(t)),...
            'previous','extrap');
    elseif numel(datetable(idx).Dist)==1
        trodedist(t)=datetable(idx).Dist;
    end
end
end