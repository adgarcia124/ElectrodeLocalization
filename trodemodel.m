function [inregion,regionidx,epos,sh,e,ax]=trodemodel(stlfiles,trodestruct,datetable,trodenum,querydates,makefig,ax,params,distoffset)
%TRODEMODEL - Determines whether electrodes are within STL models and
%optionally creates a 3D figure of recording targets and electrode 
%recording positions.
%
%Jon Rueckemann 2019

%Default parameters
regionalpha=0.45;
regioncolor=[];
regionedgealpha=0.35;
regionedges='--'; %dotted line
facelighting='gouraud';
lightsource=[-15 15 -5];

trodealpha=0.95;
trodecolor=[];
trodeedgealpha=0.01;
trodeedges='none'; %dashed line
troderad=0.15;

plotall=true;
reverseyz=true;  %#ok<*UNRCH>
viewangle=[-108 4];

if nargin<3
    datetable=[];
end
if nargin<4 || isempty(trodenum)
    trodenum=cellfun(@str2num,{datetable.Channel});
    trodenum=trodenum(:);
end
if nargin<5
    querydates=[];
end
if nargin<6
    makefig=false;
end
if nargin<9
    distoffset=0;
end

%Load analysis options from params struct
if nargin==8 && isstruct(params)
    fnames=fieldnames(params);
    for m=1:numel(fnames)
        eval([lower(fnames{m}) '=params.(fnames{m});']);
    end
end

%Update params with workspace variables; effectively adds default values
vars=who; %list variables in function workspace
vars=vars(~strcmp(vars,'stlfiles')&~strcmp(vars,'trodedist')...
    &~strcmp(vars,'trodestart')&~strcmp(vars,'norms')...
    &~strcmp(vars,'ax')&~strcmp(vars,'savefile')&~strcmp(vars,'m')...
    &~strcmp(vars,'fnames')&~strcmp(vars,'params')); %keep analysis params
for m=1:numel(vars)
    eval(['params.(vars{m})=' vars{m} ';']);
end

%Determine electrode distance for queried date(s)
trodedist=calctrodedist(trodenum,querydates,datetable);
trodedist=trodedist+distoffset;

%Get electrode starts and norms for each entry
[~,trodeidx]=ismember(trodenum,[trodestruct.ChannelID]);
trodestart=[trodestruct(trodeidx).StartPos];
trodestart=reshape(trodestart,3,[])';
norms=[trodestruct(trodeidx).Norms];
norms=reshape(norms,3,[])';

%Populate values for each region
n_trode=numel(trodedist);

if numel(regionalpha)==1
    regionalpha=regionalpha.*ones(numel(stlfiles),1);
end
if numel(regionedgealpha)==1
    regionedgealpha=regionedgealpha.*ones(numel(stlfiles),1);
end
if isempty(regioncolor)
    regioncolor=lines(numel(stlfiles));
elseif size(regioncolor,1)==1
    regioncolor=repmat(regioncolor,numel(stlfiles),1);
end

if numel(trodealpha)==1
    trodealpha=trodealpha.*ones(n_trode,1);
end
if numel(trodeedgealpha)==1
    trodeedgealpha=trodeedgealpha.*ones(n_trode,1);
end
if isempty(trodecolor)
    trodecolor=repmat([1 0 0],n_trode,1); %red
elseif numel(trodecolor)==1
    trodecolor=repmat(trodecolor,n_trode,1);
end

%Ensure axes for plotting exists
if makefig
    if nargin<7 || isempty(ax)
        fig=figure;
        ax=axes('parent',fig);
    end
    hold(ax,'on');
end

%Plot each STL surface and tag the electrodes within the region
inregion=false(n_trode,numel(stlfiles));
regionidx=nan(n_trode,1);
sh={}; %#ok<*AGROW>
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
    
    if makefig
        %Render surface of region
        cmap=interpcmap(1:64,[regioncolor(m,:); 0 0 0]);
        cmap=cmap(1:64,:);
        lightdist=pdist2(lightsource, V);
        cval=interpcmap(lightdist,cmap);
        %cval=interpcmap([0 lightdist],cmap);
        %cval=cval(2:end,:);
        if reverseyz
            V=V(:,[1 3 2]);
            V(:,1)=-V(:,1);
            sh{m}=patch('Faces',F,'Vertices',V(:,[1 3 2]),'facevertexcdata',...
                cval,'facecolor','interp','parent',ax);
        else
            sh{m}=patch('Faces',F,'Vertices',V,'facevertexcdata',...
                cval,'facecolor','interp','parent',ax);
        end
        sh{m}.FaceAlpha=regionalpha(m);
        sh{m}.EdgeAlpha=regionedgealpha(m);
        sh{m}.LineStyle=regionedges;
        sh{m}.FaceLighting=facelighting;
    end
end

%Determine 3D position of each electrode and plot
epos=trodestart+norms.*repmat(trodedist,1,3);

if makefig
    [x,y,z]=sphere;
    e={};
    for m=1:n_trode
        if ~isnan(regionidx(m))||plotall
            if reverseyz
                e{m}=patch(surf2patch(epos(m,1)+x*troderad,...
                    epos(m,3)+z*troderad,epos(m,2)+y*troderad,'triangles'),...
                    'FaceColor',trodecolor(m,:),'parent',ax);
            else
                e{m}=patch(surf2patch(epos(m,1)+x*troderad,...
                    epos(m,2)+y*troderad,epos(m,3)+z*troderad,'triangles'),...
                    'FaceColor',trodecolor(m,:),'parent',ax);
            end
            e{m}.FaceAlpha=trodealpha(m);
            e{m}.EdgeAlpha=trodeedgealpha(m);
            e{m}.LineStyle=trodeedges;
        end
    end
    
    %Change view
    view(viewangle);
end
end