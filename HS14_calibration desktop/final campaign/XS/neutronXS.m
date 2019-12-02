% make a plot for the neutron xs

fnames = {'H1.txt',...
    'H2.txt',...
    'C12.txt',...
    'C13.txt',...
    'O16.txt',...
    'O17.txt',...
    'O18.txt',...
    'Al27.txt',...
    'Cl35.txt',...
    'Cl37.txt',};

fig=figure(123);clf
ax=gca();
hold on;
clear data;

for fil=1:size(fnames,2) % rad in data, find min & max energy & data point number
    %sprintf(num2str(fil))
    data{fil}=dlmread(fnames{fil},' ',2);
    
    % remove duplicates
    [uniqueA i j] = unique(data{fil}(:,1),'first');
    indexToDupes = find(not(ismember(1:numel(data{fil}(:,1)),i)));
    data{fil}(indexToDupes,:)=[];
    
    
    emin(fil)=data{fil}(1,1);
    emax(fil)=data{fil}(end,1);
    ndata(fil)=size(data{fil},1);
    
end
iconraw=dlmread('icon.spec',',');
iconraw(1:2,:)=[]; %cut away the top 2 points because they behave bad
iconraw(:,1)=iconraw(:,1)*1000000; %MeV to eV
[uniqueA i j] = unique(iconraw(:,1),'first');
indexToDupes = find(not(ismember(1:numel(iconraw(:,1)),i)));
iconraw(indexToDupes,:)=[];

% create energy vector
energy=logspace(log10(max(emin)),log10(min(emax)),2000);
%%
% interpolate xs to energy vector
clear isoxs;
for fil=1:size(fnames,2)
    isoxs(fil,:)=interp1(data{fil}(:,1),data{fil}(:,2),energy,'linear',0);
    icon=interp1(iconraw(:,1),iconraw(:,2),energy,'linear',0);
    plot(energy,isoxs(fil,:),'displayname',fnames{fil}(1:end-4));
    
end
plot(iconraw(:,1),iconraw(:,2),'k','displayname','ICON')
legend()
set(ax,'XScale','log','YScale','log')
grid on
title('neutron xs of isotopes')
savefig('neutron xs isotopes')
print('neutron xs isotopes','-dpng')

%% make elements
eles={'H','D','C','O','Al','Cl'};
% elemtx=[0.9998 0.0002  0      0     0      0      0     0   0    0 ; %H
%         0      0       0.989 0.011  0      0      0     0   0    0; %C
%         0      0       0      0     0.9976 0.0004 0.002 0   0    0; %O
%         0      0       0      0     0      0      0     1   0    0 ; % Al
%         0      0       0      0     0      0      0     0   0.76 0.24]; % Cl

elemtx=[ 0.9998 0.0002   0      0     0      0      0     0   0    0 ; %H
    0.002 0.998    0      0     0      0      0     0   0    0; %D
    0      0       0.989 0.011  0      0      0     0   0    0; %C
    0      0       0      0     0.9976 0.0004 0.002 0   0    0; %O
    0      0       0      0     0      0      0     1   0    0 ; % Al
    0      0       0      0     0      0      0     0   0.76 0.24]; % Cl

% check: sums are 1
disp('sum of the element isotopes:')
disp(sum(elemtx,2))

elexs=elemtx*isoxs;

figure(255);clf;
ax=gca();
hold on
for i=1:size(elexs,1)
    plot(energy,elexs(i,:),'displayname',eles{i});
    
    
end

plot(energy,icon,'k','displayname','ICON')
legend()
set(ax,'XScale','log','YScale','log')
grid on
title('neutron xs of elements (nat)')
savefig('neutron xs elements (nat)')
print('neutron xs elements (nat)','-dpng')


%% mke substances micro xs
subs={'H2O','D2O','CHCl3','Al'};

submtx=[2 0 1 0 0 0; %H2O
    0 2 1 0 0 0; %D2O
    1 0 0 1 0 3; %CHCL3
    0 0 0 0 1 0]; % AL

microxs=submtx*elexs;

fig=figure(2155);clf;
pub.BFfigure(fig,T.F)
ax=gca();
pub.BFaxis(ax,T.F)
hold on
for i=1:size(microxs,1)
    plot(energy,microxs(i,:),'displayname',subs{i});
end

plot(energy,icon*10,'k','displayname','ICON')
set(ax,'XScale','log','YScale','log')

%ylh=ylabel('microscopic cross section \sigma [barn] \n relative ICON spectrum [-]',...
%    'Interpreter','latex')

ylh=ylabel({'microscopic cross section, [barn]','ICON relative spectrum [-]'})
pub.BFylab(ylh,T.F)

lh=legend();
pub.BFlegend(lh,T.F)


grid on
% title('mico xs')
savefig('mico xs')
print('mico xs','-dpng')

%% attenuation coefficient
rho=[0.955 1.1 1.411 2.7]; % g/cm3
atomicweight=[18.01528 20.0286 119.38 26.98]; %g/mol
na=6.022e23; % atoms per mol
barn=1e-24; %cm2
att=na*barn*repmat(rho./atomicweight,[size(microxs,2),1])'.*microxs;

fh=figure(1246);
figfac=1;
pub.BFfigure2(fh,T.F,figfac)

clf;
ax=gca;
pub.BFaxis(ax,T.F)
hold on
col=get(groot,'DefaultAxesColorOrder');
coli=[1,4,2,3];
att(3,1236)=1.2

for i=1:size(microxs,1)
    effxs(i)=dot(icon,att(i,:))/sum(icon)
    p(i)=plot(energy,att(i,:),...
        'displayname',sprintf('%s, %.2f 1/cm',subs{i},effxs(i)),...
        'color',col(coli(i),:));
    pub.BFLine(p(i),T.F)
end
p(i+1)=plot(energy,icon,'k','displayname','ICON spectrum');
pub.BFLine(p(i+1),T.F)

set(ax,'XScale','log','YScale','log')
xlim([1e-5,1e7]);
%xlim(1e3*[20,110])
ylim([1e-2,1e2])




grid on;
xlh=xlabel('Energy, [eV]');
pub.BFxlab(xlh,T.F)
ylh=ylabel({'attenuation coefficient, [cm^{-1}]','relative spectrum'});
pub.BFylab(ylh,T.F)
% th=title('Attenuation');
% pub.BFtitle(th,T.F)
pos=ax.Position;
ax.Position=[pos(1:2), 1-pos(1:2)];
yticks([1e-2 1e-1 1e0 1e1]);
lh=legend('location','northeast');
lpos=lh.Position;
lh.Position=[1-lpos(3:4),lpos(3:4)]
pub.BFlegend(lh,T.F)
fname='NeutAttCoeff';
fpath=sprintf('%s%s',T.F.saveto,fname);
savefig(fh,fpath)
print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
print(sprintf('%s%s',T.F.saveto,fname),'-dpng')

%%
% %% now xrays, comppund data already from nist in g/cm2

% read in files
%first the data from spekcalc
if 1
    
    % **** COMMENT ****
    %
    % **** INPUTS ****
    % kVp [kV] hvMIN [keV] Dhv [keV]
    % 110  11  1
    % Angle [deg.]
    % 12
    % t_AIR t_BE t_AL t_CU t_SN t_W t_Wa [mm]
    % 1000  0  3.4  0.2  0  0  0
    % Nf P
    % 0.68  0.33
    % ****CALCULATED OUTPUTS ****
    % Brem[uGy/mAs@1m] Char[uGy/mAs@1m]
    % 37.77275  4.150026
    % HVL1[cm AL] HVL2[cm AL] HVL1[cmCu] HVL2[cmCu] MeanE[keV] EffEAl[keV] EffECu[keV]
    % 0.7530923  0.9011965  0.0377079  0.0587936  61.55855  52.37369  54.62214
    % **** CALCULATED SPECTRUM ****
    % Energy[keV]  N[keV cm^2 mAs]^-1 @ 1 meter
    
    spec= [11 , 4.235940e-15;
        12 , 1.056486e-10 ;
        13 , 2.567868e-7  ;
        14 , 0.0000852    ;
        15 , 0.0074988    ;
        16 , 0.2509247    ;
        17 , 3.885374     ;
        18 , 34.772       ;
        19 , 206.5956     ;
        20 , 886.4131     ;
        21 , 2924.845     ;
        22 , 7997.122     ;
        23 , 18415.28     ;
        24 , 37476.11     ;
        25 , 67421.13     ;
        26 , 112373.8     ;
        27 , 171939.6     ;
        28 , 252642       ;
        29 , 344636.6     ;
        30 , 451831.1     ;
        31 , 570758.7     ;
        32 , 700179       ;
        33 , 835063.7     ;
        34 , 974990.4     ;
        35 , 1.112241e+6  ;
        36 , 1.246623e+6  ;
        37 , 1.375792e+6  ;
        38 , 1.500727e+6  ;
        39 , 1.613081e+6  ;
        40 , 1.720311e+6  ;
        41 , 1.815014e+6  ;
        42 , 1.902275e+6  ;
        43 , 1.977505e+6  ;
        44 , 2.045455e+6  ;
        45 , 2.098338e+6  ;
        46 , 2.146997e+6  ;
        47 , 2.184832e+6  ;
        48 , 2.214508e+6  ;
        49 , 2.237069e+6  ;
        50 , 2.250773e+6  ;
        51 , 2.259073e+6  ;
        52 , 2.256035e+6  ;
        53 , 2.253790e+6  ;
        54 , 2.244392e+6  ;
        55 , 2.233469e+6  ;
        56 , 2.212175e+6  ;
        57 , 2.190791e+6  ;
        58 , 6.109498e+6  ;
        59 , 9.108508e+6  ;
        60 , 2.110834e+6  ;
        61 , 2.078199e+6  ;
        62 , 2.043164e+6  ;
        63 , 2.005841e+6  ;
        64 , 1.971676e+6  ;
        65 , 1.929713e+6  ;
        66 , 1.891193e+6  ;
        67 , 4.386712e+6  ;
        68 , 1.808032e+6  ;
        69 , 2.442177e+6  ;
        70 , 1.522234e+6  ;
        71 , 1.490227e+6  ;
        72 , 1.459920e+6  ;
        73 , 1.429941e+6  ;
        74 , 1.397765e+6  ;
        75 , 1.366144e+6  ;
        76 , 1.331933e+6  ;
        77 , 1.298882e+6  ;
        78 , 1.265607e+6  ;
        79 , 1.231013e+6  ;
        80 , 1.196406e+6  ;
        81 , 1.161623e+6  ;
        82 , 1.127735e+6  ;
        83 , 1.092687e+6  ;
        84 , 1.056985e+6  ;
        85 , 1.021636e+6  ;
        86 , 987102.2     ;
        87 , 951791.1     ;
        88 , 917056.8     ;
        89 , 881159.1     ;
        90 , 846482.2     ;
        91 , 810705.8     ;
        92 , 775893.3     ;
        93 , 740281.6     ;
        94 , 704619.7     ;
        95 , 669828.8     ;
        96 , 633688.3     ;
        97 , 598350.6     ;
        98 , 561787.3     ;
        99 , 526209       ;
        100,  489580.9    ;
        101,  453083.4    ;
        102,  415478.2    ;
        103,  377931.4    ;
        104,  338111.9    ;
        105,  296507.4    ;
        106,  255760.7    ;
        107,  214721.3    ;
        108,  166169.4    ;
        109,  81787.88    ;
        110,  0           ];
end  % the spectrum goes here

spec(:,1)=spec(:,1)*1000; %keV to eV


xfnames = {'H2O.txt',...
    'CHCL3.txt',...%'Anticorodal.txt',...
    'Al_xray.txt'};
xplotfnames = {'H2O',...
    'CHCl3',...%'Anticorodal.txt',...
    'Al'};


for fil=1:size(xfnames,2) % rad in data, find min & max energy & data point number
    %sprintf(num2str(fil))
    temp=dlmread(xfnames{fil},' ',3);
    xdata(:,fil)=temp(:,2);
end

xenergy=temp(:,1)*1000000; % MeV to eV;
xspec=interp1(spec(:,1),spec(:,2),xenergy,'linear',0);


% combine energy range of the xs and spectrum data
enecomb=sort([xenergy;spec(:,1)]);
% construct energy bin width
denecomb=enecomb(2:end)-enecomb(1:end-1); % forward bin sizes
enebin=0.5*([denecomb(1);denecomb]+[denecomb;denecomb(end)]);

speccomb=interp1(spec(:,1),spec(:,2),enecomb,'linear',0);
for fil=1:size(xfnames,2) % rad in data, find min & max energy & data point number
xdatacomb(:,fil)=interp1(xenergy,xdata(:,fil),enecomb,'linear',0);
end


fh=figure(1231);clf
pub.BFfigure2(fh,T.F,figfac)
ax=gca();
hold on;
clear data p2;
for fil=1:size(xfnames,2) % rad in data, find min & max energy & data point number
    xeffxs(i)=dot(speccomb.*enebin,xdatacomb(:,fil))/sum(speccomb.*enebin);
    p2(fil)=plot(enecomb/1000,xdatacomb(:,fil),...
        'displayname',sprintf('%s, %.2f 1/cm',xplotfnames{fil},xeffxs(i)));
    
end
p2(fil+1)=plot(enecomb/1000,speccomb/1500000,'k','displayname','110kV spectrum');

pub.BFLine(p2,T.F)

%set(ax,'YScale','log')
xlim(1e0*[20,110]);
ylim([1e-1,1e4])
ylim([0 8])
grid on;
xlh=xlabel('Energy, [keV]');
pub.BFxlab(xlh,T.F)
ylh=ylabel({'attenuation coefficient, [cm^{-1}]','relative spectrum'});
pub.BFylab(ylh,T.F)
yticks([0 2 4 6])
pos=ax.Position;
ax.Position=[pos(1:2), 1-pos(1:2)];

lh2=legend('location','best');
%lpos2=lh2.Position

pub.BFlegend(lh2,T.F)


fname='XrayAttCoeff';
fpath=sprintf('%s%s',T.F.saveto,fname);
savefig(fh,fpath)
print(sprintf('%s%s',T.F.saveto,fname),'-dpdf')
print(sprintf('%s%s',T.F.saveto,fname),'-dpng')


% create energy vector











