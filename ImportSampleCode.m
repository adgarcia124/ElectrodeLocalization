%% Spock
load('SpockTrodeInfo.mat')
datetable=trodelog; %Use GUI to select turn log Excel file
stlpath='curSpock.stl';



%Find all electrode positions today and create STL file
epos=trodereconstruction(Spock,stlpath,[],[],datetable);

%Or find all electrode positions for specified query date
qdate={'16-Nov-2019'};
epos=trodereconstruction(Spock,stlpath,[],qdate,datetable);

%Or query individual channels and specific dates (each electrode number is
%individually paired with a query date)
    %This example finds the position electrode 50 (E50) on Nov 10th,
    %E110 on Nov 12th, E50 on Nov 16th, and E120 on Nov 16th.
tnum=[50; 110; 50; 120];
qdate={'10-Nov-2019';'12-Nov-2019';'16-Nov-2019';'16-Nov-2019'};
epos=trodereconstruction(Spock,stlpath,tnum,qdate,datetable);

%% Gromit
load('GromitTrodeInfo.mat')
datetable=trodelog; %Use GUI to select turn log Excel file 
stlpath='curGromit.stl';

%Find all electrode positions today and create STL file
epos=trodereconstruction(Gromit,stlpath,[],[],datetable);