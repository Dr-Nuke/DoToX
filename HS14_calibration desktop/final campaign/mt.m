classdef mt % the MT method test class
    
%% function collection file for the liquid film thickness computation
methods(Static) %evil function scope hack
    
function raw=loadData(T,F,MT) % grab all the relevant war data, i.e. the corrected
                        % added and shifted sinograms
    nplanes=MT.raw.nplanes;     %number of planes to use
    startplane=MT.raw.startplanes; % starting plane. 300 is somewhat over the spacer 
    planes=startplane:(startplane+nplanes-1); % planes list
    raw=zeros(F.Ncases,T.imsize(1),nplanes,T.q360.nFrames,'single'); % preallocation
    cass=[1,2,4];
    for i=1:3 %leave out mid case
        cas=cass(i);
        rfname=sprintf('2_add_%02d.mat',cas);
        d=f.loadSingleVariableMATFile(strcat(T.d.DataPath,rfname));
        raw(i,:,:,:)=d(:,planes,:);
    end    
end

function recon=OldRec(T,F,MT) % old reconstruction method
    cstart=170;
    cstop=488;
    recsize=cstop-cstart+1;
    recon=zeros(3,recsize-1,recsize-1,MT.raw.nplanes,'single'); % block of recons
    cas_old=[1 2 4];
    for cas =1:3
        caso=cas_old(cas);

        d2=squeeze(MT.raw.d(cas,:,:,:));
        if cas>1 % if we have to divide
            d2=-log(d2./squeeze(MT.raw.d(1,:,:,:)));
        end
        d2=d2(cstart:cstop,:,:);
    
        % turn off warning from iradon about NaNs in lines
        id='MATLAB:interp1:NaNstrip';
        warning('off',id);
        
        temprec=zeros(recsize-1,recsize-1,MT.raw.nplanes,'single');
        for j=1:MT.raw.nplanes
            f.f_BoCount(j,5,10,3)
            try
                imc=circshift(squeeze(d2(:,j,:)),200,2); % rotate such that...
                    %  one quadrant per heating rod (for rd finding later)
                imc=f.f3_center(imc,T.c.fitshift(caso,j+MT.raw.startplanes-1));
                imc=imc(2:end-1,:);
                rec=iradon(imc,T.angles,'spline','Hann',1,recsize-1);
                %disp(sprintf('%d %d',j,size(imc,1)));
                temprec(:,:,j)=rec; %hardcode after-recon-cropping
            catch
                disp(sprintf('%4d didnt reconstruct',j))
            end
        end
        recon(cas,:,:,:)=temprec;
    end
end

function recon=NewRec1(T,F,MT) % new reconstruction :substract after reconstruction
    cstart=170;
    cstop=488;
    recsize=cstop-cstart+1;
    recon=zeros(3,recsize-1,recsize-1,MT.raw.nplanes,'single'); % block of recons
    cas_old=[1 2 4];
    for cas =1:3
        caso=cas_old(cas);
        d2=squeeze(MT.raw.d(cas,:,:,:));
        d2=d2(cstart:cstop,:,:);
    
        % turn off warning from iradon about NaNs in lines
        id='MATLAB:interp1:NaNstrip';
        warning('off',id);
        
        temprec=zeros(recsize-1,recsize-1,MT.raw.nplanes,'single');
        for j=1:MT.raw.nplanes
            f.f_BoCount(j,5,10,3)
            try
                imc=circshift(squeeze(d2(:,j,:)),200,2); % rotate such that...
                    %  one quadrant per heating rod (for rd finding later)
                imc=f.f3_center(imc,T.c.fitshift(caso,j+MT.raw.startplanes-1));
                imc=imc(2:end-1,:);
                rec=iradon(imc,T.angles,'spline','Hann',1,recsize-1);
                %disp(sprintf('%d %d',j,size(imc,1)));
                temprec(:,:,j)=rec; %hardcode after-recon-cropping
            catch
                disp(sprintf('%4d didnt reconstruct',j))
            end
        end
        
        
        recon(cas,:,:,:)=temprec;
        if cas >1
            recon(cas,:,:,:)=recon(1,:,:,:)-recon(cas,:,:,:);
        end
    end
end

function Slide2Recons(MT,i)

    crange=[-0.001,0.0015];
    block=squeeze(cat(2,MT.rec1(i,:,:,:),MT.rec2(i,:,:,:)));
    slide_viewer(block,crange,3);

end

function [kernel] = KernelGen(d,rk,w)
    %returns a gaussian ring kernel 
    % with edg size d
    %with radius r and width w
if mod(d,2)~=1;
    disp('bad kernel diameter not odd');
end
%pixel diameter of kernel . inner and outer. outer must be odd

r0=((d-1)/2); % middle pixel
[x,y]=meshgrid(-r0:r0,-r0:r0);
r=sqrt(x.^2+y.^2);
kernel=exp(-(r-rk).^2/(2*w^2));
kernel=kernel/sum(kernel(:));

end

function [kernel] = KernelGenBin(d,do,di)
    % creates a binary mask
    % d = mask size (edge)
    % do = outer diameter
    % di = inner diameter
if mod(d,2)~=1;
    disp('bad kernel diameter not odd');
end
%pixel diameter of kernel . inner and outer. outer must be odd
ro=do/2; %outer radius of kernel
ri=di/2; %inner radius of kernel
r=((d+1)/2); % middle pixel
[x,y]=meshgrid(1:d,1:d);
kernel=zeros(d);
kernel(and((x-r).^2+(y-r).^2>=ri^2,(x-r).^2+(y-r).^2<=ro^2))=1;
end
function SOQuestion(F,recon)
    % generates images for a stack exhchange question
    h=figure(405);clf
    imshow(squeeze(recon(2,:,:,400)),[])
    
    ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
    saveas(h,'Stackexchange raw.png')
    
    
end

function [kernel] = ElliKernelGen(dx,dy,sx,sy)
    %returns a elliptical gaussian kernel 
    % x-size dy
    % y-size dy
    % x-sigma sx
    % y-sigma sy
if any(mod([dx,dy],2)~=1)
    disp('bad kernel diameter not odd');
end

%pixel diameter of kernel . inner and outer. outer must be odd

x0=((dx-1)/2); % middle pixel
y0=((dy-1)/2); % middle pixel

[x,y]=ndgrid(-x0:x0,-y0:y0);

kernel=exp(-((x).^2/(2*sx^2)+(y).^2/(2*sy^2)));
kernel=kernel/sum(kernel(:));

end

end %static
end %class
















