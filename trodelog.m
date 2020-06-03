function [datetable]=trodelog(xlsfile,sumdist,rmentry)
%Determine the depth of electrodes by date 
%(particular to Rueckemann/GrayMatter conventions)
%
%Input:
%   xlsfile-path of excel file with turn counts
%   sumdist-summate distances if not originally entered as cumulatively
%
%Output:
%   datatable-struct with elements organized by channel number, containing
%   dates (as string), dates (as datevec), and the corresponding turn
%   distance
%
%Jon Rueckemann 2018

if nargin<1 || isempty(xlsfile)
    [xlsfile,path]=uigetfile;
    xlsfile=fullfile(path,xlsfile);
end
if nargin<2 || isempty(sumdist)
    sumdist=false;
end
if nargin<3
    rmentry=[];
end

%List electrode sheets in excel workbook (assumes sheetnames are numerical)
[~,sheetnames]=xlsfinfo(xlsfile);
if all(cellfun(@(x) isempty(str2num(x)),sheetnames)) %#ok<ST2NM>
    warning(['Sheet names are not numerical. '...
        'Electrode number will be extracted from numeric characters.']);
     ChID=cellfun(@(x) str2num(x(isstrprop(x,'digit'))),...
         sheetnames,'uni',0); %#ok<ST2NM>
     
     sheetnames(cellfun(@isempty,ChID))=[];
     ChID(cellfun(@isempty,ChID))=[];
else
    sheetnames(cellfun(@isempty,cellfun(@str2num,sheetnames,'uni',0)))=[];
    ChID=sheetnames;
end


%Process data for each electrode
datetable=struct;
for m=1:numel(sheetnames)
    disp(sheetnames{m});
    [~,~,dates]=xlsread(xlsfile,sheetnames{m},'A4:A1000');
    [~,~,dist]=xlsread(xlsfile,sheetnames{m},'D4:D1000');
    dist=cell2mat(dist);
    
    if ~isempty(rmentry)
        [~,~,experimenter]=xlsread(xlsfile,sheetnames{m},'B4:B1000');
        dropidx=strcmpi(experimenter,rmentry);
        dates(dropidx)=[];
        dist(dropidx)=[];
    end
    
    %Find last row index in current spreadsheet
    lastidx=max(find(~cellfun(@isnumeric,dates),1,'last'),...
        find(dist>0,1,'last'));
    dates=dates(1:lastidx);
    dist=dist(1:lastidx);
    
    %Create a date for each distance entry
    dateidx=cumsum(~cellfun(@isnumeric,dates)); %indices for each date
    dates=dates(~cellfun(@isnumeric,dates)); %sparse date list
    dates=dates(dateidx);    
    
    %Reduce the entries; only include location at the end of the day    
    numericdates=cellfun(@datenum,dates); %Convert dates to numeric codes    
    [numericdates,dateendidx]=unique(numericdates,'last');
    dist=dist(dateendidx);
    if sumdist
        dist=cumsum(dist);
    end
    dates=dates(dateendidx);    
    
    datetable(m).Channel=num2str(ChID{m});
    datetable(m).Dates=dates;
    datetable(m).NumericDates=numericdates;
    datetable(m).Dist=dist;
end
end