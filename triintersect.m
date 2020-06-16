function [idx,intdist,intpt,Aidx,Adist,Aintpt]=triintersect(triobj,pts,norms)
%Find the nearest intersection of rays with a triangulation object
%
%Jon Rueckemann 2017

%Find the normal to the plane that defines the face of each triangle
fn=faceNormal(triobj);

P1=triobj.Points(triobj.ConnectivityList(:,1),:);
P2=triobj.Points(triobj.ConnectivityList(:,2),:);
P3=triobj.Points(triobj.ConnectivityList(:,3),:);

%Iterate through each ray
idx=nan(size(pts,1),1);
intdist=nan(size(pts,1),1);
intpt=nan(size(pts,1),3);
[Aidx,Adist,Aintpt]=deal(cell(size(pts,1),1));
for m=1:size(pts,1)
    curpt=repmat(pts(m,:),size(fn,1),1);
    curnorm=repmat(norms(m,:),size(fn,1),1);
    
    %Find distance from each point to the plane of each triangle
    curdist=(dot(P1'-curpt',fn')./dot(curnorm',fn'))'; %may need to be vectorized
    
    %Get intersection of each ray with the plane of each triangle
    P0=curpt+curdist.*curnorm;
    
    %Determine if planar intersection is within the triangle
    curidx=(dot(cross(P2'-P1',P0'-P1'),fn')>=0 & ...
        dot(cross(P3'-P2',P0'-P2'),fn')>=0 & ...
        dot(cross(P1'-P3',P0'-P3'),fn')>=0)';
    
    %Find the nearest triangle face
    curdist(~curidx)=nan;
    curdist=abs(curdist);
    [intdist(m),idx(m)]=min(curdist);
    if ~isnan(intdist(m))
        intpt(m,:)=P0(idx(m),:);
    else
        idx(m)=nan;
    end
    
    if nargout>3
        %Find all intersecting triangle faces
        [Adist{m},srtidx]=sort(curdist(curidx));
        x=find(curidx);
        Aidx{m}=x(srtidx);
        y=P0(curidx,:);
        Aintpt{m}=y(srtidx,:);
    end
end
end