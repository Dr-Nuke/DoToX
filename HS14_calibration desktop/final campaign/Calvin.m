
% type in the base folder
basefol='E:\20171213 final campaign\damsohn\Liq_Films_Exp\'; 

%type in the individual folders, for best functionality in sets of three
folds={'Sp00_air\','Sp00_c4f8\','Sp00_he\',...
    'Sp05_air\','Sp05_c4f8\','Sp05_he\'}%,...
    %'Sp06_air\','Sp06_c4f8\','Sp06_he\'}; 
fre=10000; % i guess frequency
%clear C % C the Calvin struct
nfolds=length(folds); % number of folders

jgs=[]; % collects all available J_g values
betas=[]; % same for betas
nfolds=length(folds)

% now we first run a loop to scan available stuff, and in the second loop
% we actually read in

for fol=1:nfolds % scanner loop
    folder=folds{fol};
    
    if folder(6)=='a' % define the gas
        gas=1;
    elseif folder(6)=='c'
        gas=2;
    elseif folder(6)=='h'
        gas=3;
    else
        warning('gas is not found')
    end
    oldfolder=cd(strcat(basefol,folder)); % swicht directory
    files=dir('J*b*.cal'); % list files
    cd(oldfolder) % swich back to old directory
    for fil=1:length(files) % scan through stuff
        
        %record some properties
        fname=strcat(basefol,folder,files(fil).name); % file name
        jg=str2num(files(fil).name(2:3));   % J_g
        bet=str2num(files(fil).name(6:9));  % 1-beta
        disp(sprintf('%s %d %d ',folder,jg,bet)) % command line debug
        
        jgs=[jgs,jg];       % dat list of velocities
        betas=[betas,bet];  % same for 1-betas
        
    end
    ujg=unique(jgs);        % list of unique velocites
    ubetas=unique(betas);   % list of unique betas
    
    njg=length(ujg);        % number of....
    nbetas=length(ubetas);  % ....
end
% done scanning

for fol=1:nfolds % reader loop
    folder=folds{fol};
    
    if folder(6)=='a' % define the gas
        gas=1;
    elseif folder(6)=='c'
        gas=2;
    elseif folder(6)=='h'
        gas=3;
    else
        warning('gas is not found')
    end
    oldfolder=cd(strcat(basefol,folder));
    files=dir('J*b*.cal');
    cd(oldfolder)
    
    for fil=1:length(files) % check out each file
        
        % open & read in file; close handle afterwards
        fname=strcat(basefol,folder,files(fil).name);
        file=fopen(fname,'r');
        S=double(fread(file,64*16*fre*10,'int16'));
        fclose(file);
        
        jg=str2num(files(fil).name(2:3)); % J_g
        bet=str2num(files(fil).name(6:9));% 1-beta
        spac=str2num(folder(4))+1;          % spacer number, +1 increment... 
        % because it starts at zero, bad for indexing
        
        ind_jg=find(ujg==jg);       % find the how-many-th J_g of unique
        ind_bet=find(ubetas==bet);  % same for 1-beta
        
        % debug
        disp(sprintf('%s %d %d %d %d %d ',folder,spac,gas,jg,bet))
        
        % transform data into block of x*y*frames, same format as dotox
        S =(reshape(S,[64 16 fre*10]));
        S =permute(S,[2 1 3]); % get x and y in right order
        S =flipdim(S,2); %reverse y dim
        S =flipdim(S,1); % correct x dim
        
        % example plot
        % frame =4567;
        % imshow(S(:,:,frame)',[]);set(gca,'YDir','normal');
        
        % write mean value and some properties to C struct
        exp.mean=mean(S,3);
        exp.folder=folder;
        exp.fil=files(fil).name;
        exp.jg=jg;
        exp.beta=bet;
        exp.spac=spac;
        exp.gas=gas;
        
        C.exp(spac,gas,ind_jg,ind_bet)=exp;
        C.map(spac,gas,ind_jg,ind_bet)=1;
        clear exp;
        
    end
end
disp('done read in')

C.angles=((1:16)-1)/15*180;
C.heights=((1:64)-1)/63*120;
C.Dh=18.85;
C.Dhs=C.heights/C.Dh;
C.jg=ujg;
C.beta=ubetas;

%% make case-overview-plots

figure(11223);clf
spacs=[1,6];
liqs={'air','C4F8','He'};
for spa=1:2 % reader loop
    for gas=1:3
        ind=gas+(spa-1)*3;
        
        subplot(2,3,ind);
        
        imagesc(squeeze(C.map(spacs(spa),gas,:,:))','xdata',C.jg);
        set(gca,'YDir','normal');
        title(sprintf('%d %d %s',spa,gas,liqs{gas}))
        ytx=1:2:length(C.beta);
        yticks(ytx);
        yticklabels(C.beta(ytx)/10000);
        xlabel('J_g')
        ylabel('1-beta')
        xtx=1:2:length(C.jg);
        xticks([C.jg(xtx)]);
        
        %     figure(664+fol);
        %     clf;
        %     ax=gca;
        %     imagesc(squeeze(C.map(fol,:,:))','xdata',ujg);
        %     set(gca,'YDir','normal');
        %     title(folder(1:end-1))
        %     xlabel('J_g [m/s]')
        %     ylabel('1-beta')
    end
    
end

figure(7824);clf
im=cat(3,squeeze(C.map(6,1,:,:))',...
    squeeze(C.map(6,2,:,:))',...
    squeeze(C.map(6,3,:,:))');
imagesc(im,'xdata',C.jg);
set(gca,'YDir','normal');
title('r=air, g=C4F8, b=He')
yticks(ytx);
yticklabels(C.beta(ytx)/10000);
xlabel('J_g')
ylabel('1-beta')
       xtx=1:2:length(C.jg);
        xticks([C.jg(xtx)]);


%%
figure(251465);clf
gases=[3 1 2];
spacstr={'no spac','Sp2-A'};
j=6;
b=3;

for spa=1:2 % reader loop
    for gas=1:3
        ind=(gas-1)*2+spa;
        ax=subplot(3,2,ind);
        gass=gases(gas);
        
        imagesc((C.exp(spacs(spa),gass,j,b).mean)');
        set(gca,'YDir','normal');
        colormap(flipud(gray))
        title(sprintf('%s %s %d m/s, %0.4f',spacstr{spa},liqs{gases(gas)},C.exp(1,1,j,b).jg,...
        C.exp(1,1,j,b).beta/10000))
        caxis([0 500])
   
        
        
    end
end
xlabel('compare with manus diss p.86')

%% no spacer plot


    
