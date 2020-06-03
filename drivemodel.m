function [drivebottom,drivenorm,drivetop]=drivemodel(topvertices,bottomvertices,toptemplate)
%DRIVEMODEL - Fit template to vertices of the drive STL marking electrodes.
%Note: This algorithm assumes that the electrodes are parallel and the
%template accurately describes the positions at the top and bottom of the
%drive, meaning that the bottom of the drive can be projected onto the top.
%
%topvertices - Nx3 matrix of STL vertices defining the drive top
%bottomvertices - Mx3 matrix of STL vertices defining the drive bottom
%toptemplate - Px3 matrix of electrode coordinates defining the template
%               If toptemplate is Px2, the third column is filled with 0s
%
%Output:
%  drivebottom - Px3 matrix of coordinates where electrodes exit the drive
%  drivenorm - 3x1 unit vector describing the trajectory of the electrodes
%  drivetop - Px3 matrix of electrode coordinates at the top of the drive
%
%Jon Rueckemann 2020

%Find optimal plane (and norm) for the template of the top of the drive
DM=[toptemplate(:,1),toptemplate(:,2),ones(size(toptemplate,1),1)]; %design
B_template=DM\toptemplate(:,3); %fit coefficients via least squares
if all(B_template)==0; B_template=[0 0 1]; end

%Rotate the template to be parallel to the X-Y plane
R=calculatePlaneRotation([0 0 1],B_template);
newtemplate=(R*toptemplate')';

%Translate the center of mass of template to the origin
CoM=sum(newtemplate)/size(newtemplate,1);
newtemplate=newtemplate-CoM;

%Find optimal plane (and norm) for the vertices at the top of the STL
DM=[topvertices(:,1),topvertices(:,2),ones(size(topvertices,1),1)]; %design
B_stl=DM\topvertices(:,3); %fit coefficients via least squares
drivenorm=B_stl/norm(B_stl);

%Rotate the drive top vertices and norm to be parallel to the X-Y plane
R_flat=calculatePlaneRotation([0 0 1],B_stl);
newtopvertices=(R_flat*topvertices')';

%Translate the center of mass of vertices to the origin
topCoM=sum(newtopvertices)/size(newtopvertices,1);
newtopvertices=newtopvertices-topCoM;

%Fit the template to the drive top vertices
mindist=inf(size(topvertices,1),1);
iter=0;
while any(mindist>0.5) %nearest electrode to each vertex should be <0.5mm
    %Rotate the template around the plane norm to align the orientation
    func=@(x) sum(mindistance(...
        rotateDataAroundNorm([0 0 1],x,newtemplate),newtopvertices));
    newang=fminbnd(func,0,pi);
    newtemplate=rotateDataAroundNorm([0 0 1],newang,newtemplate);
    %Minimize the distances between newtemplate & topvertices
    
    %Translate template along X axis
    func=@(x) sum(mindistance(newtemplate+x*[1 0 0],newtopvertices));
    xoffset=fminbnd(func,-100,100); %Normal axis
    newtemplate=newtemplate+xoffset*[1 0 0];
    
    %Translate template along Y axis
    func=@(x) sum(mindistance(newtemplate+x*[0 1 0],newtopvertices));
    yoffset=fminbnd(func,-100,100); %Normal axis
    newtemplate=newtemplate+yoffset*[0 1 0];
    
    %Find the distance of each vertex to nearest electrode
    mindist=mindistance(newtemplate,newtopvertices);
    
    iter=iter+1;
    if iter==21
        error(['Algorithm has not converged to <0.5mm after 20 '...
            'iterations. Either the template is ill-defined or the top '...
            'vertices are not well isolated.']);
    end
end

%Apply rotation and translation to the vertices at the bottom of the drive
newbottomvertices=(R_flat*bottomvertices')'-topCoM;

%Calculate distance of projection of bottom vertices to nearest template
bottomproj=newbottomvertices;
bottomproj(:,3)=0;
[D,idx]=mindistance(newtemplate,bottomproj);
assert(all(D<0.5),['Vertices from drive bottom do not align with top. '...
    'Either the model fit is bad or the vertex isolation is bad.']);

%Find the distance from each vertex to the new template plane
vertex_dist=newbottomvertices(:,3);

%Average the distance of each vertex corresponding to each electrode
drivebottom=newtemplate;
for m=1:size(drivebottom,1)
    drivebottom(m,3)=mean(vertex_dist(idx==m));
end

%Reverse the transforms that recentered the drive top at the origin in X-Y
drivetop=(R_flat\(newtemplate+topCoM)')';
drivebottom=(R_flat\(drivebottom+topCoM)')';
end

function R=calculatePlaneRotation(A,B)
%Rotate A to B
A=A./norm(A);
B=B./norm(B);
if abs(dot(A,B))~=1
    v=cross(A,B);
    ssc=[0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];%skew-symmetric X-matrix
    R=eye(3)+ssc+ssc^2*(1+dot(A,B))^-1; %Rodrigues' rotation formula
else
    R=eye(3);
end
end

function rot_data=rotateDataAroundNorm(A,ang,data)
%Rotate data around vector A by ang radians
A=A./norm(A);
ssc=[0 A(3) -A(2); -A(3) 0 A(1); A(2) -A(1) 0];%skew-symmetric cross matrix
R=eye(3)+sin(ang)*ssc+(1-cos(ang))*ssc^2; %Rodrigues' rotation formula
rot_data=(R*data')';
end

function [d,idx]=mindistance(U,V)
[d,idx]=min(pdist2(U,V),[],1);
end