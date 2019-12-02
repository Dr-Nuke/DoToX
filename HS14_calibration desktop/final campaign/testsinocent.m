
%% test of the centering
l=301; % rough width of my sinograms

% phantom:
a=ones(l);
px=120; % impulse position
a(px,px)=0;
shift=4.5 % artificial non-centerdness
ang=360;  % angular range
nang=1357; % number of angles

% angles list
angles=linspace(1,ang,nang+1);
angles(end)=[];

% show phantom
fig=figure(8);clf;
fig.Position(1:2)=[0,0];
imagesc(a');set(gca,'YDir','normal')
setfig(fig)
title('phantom')

% apply the artificial non-centerdness
b=fraccircshift(radon(a,angles),shift);

% show de-centered sinogram
fig=figure(9);clf;
fig.Position(1:2)=[500,0];
imagesc(b');set(gca,'YDir','normal')
setfig(fig)
title(sprintf('sinogramm offset by %.2f',shift))

% find the off-centerdness
[c,lags]=xcov(b(:,1),flipud(b(:,ceil(nang/2))),40,'coeff');
[~,ind]=max(c);
% 3 point gauss fit
centershift=(lags(ind)+0.5*(log(c(ind-1))-log(c(ind+1)))/...
    (log(c(ind-1))-2*log(c(ind))+log(c(ind+1))))/2

% correct for the off-centerdness
d=fraccircshift(b,-centershift);

% plot centered sinogram
fig= figure(10);clf;
fig.Position(1:2)=[0,500];
imagesc(d');set(gca,'YDir','normal')
setfig(fig)
title('shifted back')

% do the iradon
e=iradon(d,angles,'spline','Hann',1,l);

% show result
fig=figure(11);clf;
fig.Position(1:2)=[500,500];
imagesc(e');set(gca,'YDir','normal')
title('recon')
setfig(fig)

% zoom in to the relevant area
xlim([px-10,px+10])
ylim([px-10,px+10])

% show difference to original
fig=figure(12);clf;
fig.Position(1:2)=[1000,500];
imagesc((e-a)');set(gca,'YDir','normal')
title('recon deviation')
xlim([px-10,px+10])
ylim([px-10,px+10])

setfig(fig)


function setfig(fid)
    xlabel('b')
    ylabel('v')
    colormap(gray)
    colorbar
end

function B = fraccircshift(A,shiftsize)
    %fraccircshift expands circshift to fractional shifts values, using linear
    %interpolation. In contrast to other approaches to non-integer shifts
    %of matrices on the base of fft2 or interp2, the number of dimensions of A
    %is not limited, and the syntax of circshift applies. For integer elements
    %of shiftsize, fracircshift and circshift give the same results.
    
    int = floor(shiftsize);     %integer portions of shiftsize
    fra = shiftsize - int;      %fractional portions of shiftsize
    dim = numel(shiftsize);
    B = A;
    for n = 1:numel(shiftsize)  %The dimensions are treated one after another.
        intn = int(n);
        fran = fra(n);
        shift1 = zeros(dim,1);
        shift1(n) = intn;
        shift2 = zeros(dim,1);
        shift2(n) = intn+1;
        %Linear intepolation:
        B = (1-fran)*circshift(B,shift1) + fran*circshift(B,shift2);
    end
end
