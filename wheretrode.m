function [inregion,regionidx,epos]=wheretrode(stlfiles,trodedata,datetable,querydates,trodenum)
%WHERETRODE - Determines whether electrodes are within STL models
%
%Jon Rueckemann 2020

if nargin<3 || isempty(datetable)
    %Create table of dates and distances for electrodes from turning log
    datetable=trodelog();
end
if nargin<4
    querydates=[];
end
if nargin<5 || isempty(trodenum)
    trodenum=cellfun(@str2num,{datetable.Channel});
    trodenum=trodenum(:);
end

%Extract channel identifiers
ChannelID=[trodedata.ChannelID];
ChannelID=ChannelID(:);

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

%Determine electrode distance for queried date(s)
trodedist=calctrodedist(trodenum,querydates,datetable);

%Determine whether electrodes are within anatomical region
inregion=false(numel(trodedist),numel(stlfiles));
regionidx=nan(numel(trodedist),1);
for m=1:numel(stlfiles)
    [F,V]=stlread(stlfiles{m});
    triobj=triangulation(F,V); %triangulation object
    
    %Determine if region contains electrodes
    [~,~,~,~,Adist]=triintersect(triobj,trodestart,norms);
    n_int=cellfun(@(x) numel(x),Adist);
    if any(mod(n_int,2)>0)
        warning(['Structure ' num2str(m) ' has an odd number of '...
            'intercepts for electrodes: ' ....
            num2str(reshape(find(n_int>2),1,[]))]);
    end
    thisregion=cellfun(@(x,y) any(x(1:2:end)<y&x(2:2:end)>y),...
        Adist,num2cell(trodedist));
    inregion(:,m)=thisregion;
    regionidx(thisregion)=m; %overwrites previous region identity
end

%Determine 3D position of each electrode
if nargout>2
    epos=trodestart+norms.*repmat(trodedist(:),1,3);
end
end