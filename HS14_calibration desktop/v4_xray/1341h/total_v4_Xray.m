%% this is the 4th recon version, made for xrays

%get name of thisscript
[~,name,~]=fileparts(mfilename('fullpath')); %get m file name

%init script
s4_init

% flag for displaying the proofs
M.proofs=1;

% M for main structure, d for data
M.d.filepath=('C:/data/20161216 Campaign 5 Robert Z/13h41 leer 110kV 5mA.seq');
M.d.imsize=[1024 640];
M.d.header_size=2048;         % varian specific header size of .seq files

load('M') % 
%M=f4_loadVideo(M); % memory bug workaround
%save('M','M')


%% find all the relevant frames and dispose the rest

% manual guess
M.d.crop=[30,1386]; % kepp only thoes frames (including)
M.d.crop=[3,1457]; % kepp only thoes frames (including)
M.d.axlim=[70,130,420,460]; % required for the small axes

M=f4_CropFrames(M); % CropFrames

%% correct for the beam non-uniformity
% manual input of the dose areas
M.dose=[[5,427,429];...     % x start value
    [162,630,627];...       % x end value
    [100,100,700];...       % y start value
    [1010,229,1014]];       % y end value

% make the dose areas and check them
M=f4_DoseArea(M);

%% calculate & apply the beam non uniformity
M=f4_FitBeam(M);

% check upon it
f4_FitBeamCheck(M);
M=rmfield(M,'im');
M=rmfield(M,'unibeam');


%% go towards reconstruction
% cut away undesired stuff
M.rpx=[145,435,50,1015]; % reconpix x1 x2 y1 y2
M.recsize=291;

%M = f4_centering(M);










