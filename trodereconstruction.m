function [epos,trodedist,datetable]=trodereconstruction(trodedata,stlout,trodenum,querydates,xlsfile,stlrad)
%TRODERECONSTRUCTION - Create an STL reconstruction of each electrode at a
%query date.
%
%Input:
%   trodedata - Tx1 struct. Contains 3D start positions and trajectory
%       norms for electrodes, indentified by electrode numerical ID
%   stlout - Path of output STL file
%   trodenum - Nx1 vector, empty.  Channel names to render. Default: empty;
%       uses all available channels
%   querydates - Nx1 string array, string, empty. Dates for each electrode.
%       Default: empty; uses most recent date (today).
%   xlsfile - string or Mx1 struct.  Filepath of excel file containing
%       turning log for each electrode, OR datetable struct (see below).
%
%Output:
%   epos - Nx3 matrix. 3D position of electrode.
%   trodedist - Nx1 vector.  Distance of each electrode at queried date
%   datetable - Mx1 struct.  Distance by date table for each electrode.
%
%Jon Rueckemann 2019


%Extract channel identifiers
ChannelID=[trodedata.ChannelID];
ChannelID=ChannelID(:);

%Populate default values
if nargin<2 || isempty(stlout)
    [stlout,path]=uiputfile;
    stlout=fullfile(path,stlout);
end
if nargin<3 || isempty(trodenum)
    trodenum=ChannelID;
end
if nargin<4
    querydates=[];
end
if nargin<5 || isempty(xlsfile)
    [xlsfile,path]=uigetfile;
    xlsfile=fullfile(path,xlsfile);
end
if nargin<6 || isempty(stlrad)
    stlrad=0.25;
end
    
%Retrieve applicable electrode start locations and norms
[~,trodeidx]=ismember(trodenum,ChannelID);
if any(trodeidx==0)
    warning('Not all trodenum entries matched the ChannelID in trodedata');
    trodeidx=trodeidx(trodeidx~=0);
end
trodestart=[trodedata(trodeidx).StartPos];
trodestart=reshape(trodestart,3,[])';
norms=[trodedata(trodeidx).Norms];
norms=reshape(norms,3,[])';

%Calculate electrode distances based on date of recording
[trodedist,datetable]=calctrodedist(trodenum,querydates,xlsfile);

%Determine 3D position of each electrode
epos=trodestart+norms.*repmat(trodedist(:),1,3);

[~,~,test]=fileparts(stlout);
if ~isempty(test)
    %Create 3D STL of each electrode position
    [x,y,z]=sphere;
    [a,b,c]=deal(repmat({nan(size(x))},numel(trodenum),1));
    for m=1:numel(trodenum)
        a{m}=epos(m,1)+x*stlrad;
        b{m}=epos(m,2)+y*stlrad;
        c{m}=epos(m,3)+z*stlrad;
    end
    a=cell2mat(a);
    b=cell2mat(b);
    c=cell2mat(c);
    
    stlwrite(stlout,surf2patch(a,b,c,'triangles'));
end
end