shiftmap=main(T,d)

function shiftmap=main(T,d)

%set up stuff
cas=8;T.cas=cas;
rep=3;T.rep=rep;

crange=[-0.2,0.2];
ref=T.Raw.Sino(:,:,1,1);
disl=T.Raw.Sino(:,:,cas,rep);

%make pre-mask

env=15; % 
x0=1;
xmax=640;
y0=20;
ymax=1443;
wr=[x0+env, 195, 420,  1300;... %xstart xstop ystart ystop % window range
    455, xmax-env, 20,    615;...
    455, xmax-env, 1090, ymax-env];
premask=zeros(size(ref),'logical');
for i=1:size(wr,1)
    premask(wr(i,1):wr(i,2),wr(i,3):wr(i,4))=1;
end

diff=ref-disl;
figure(1);clf;
imshow(diff.*(premask+0.5),crange)
colorbar
title('default diff + mask')

% make pin mask
maskthresh=0.6;
dilbox=5;

refmask=imdilate((ref<maskthresh),ones(dilbox)).*premask;
dislmask=imdilate((disl<maskthresh),ones(dilbox)).*premask;

figure(2);clf;
imshow((refmask+dislmask+1).*ref,[])

uni=sum(or(dislmask(:),refmask(:)));
sec=sum(and(dislmask(:),refmask(:)));

overlay=sec/uni;
title(sprintf('overlay factor=%.2f',overlay))



shiftmap=zeros([size(ref),2]);
pointlist=find(dislmask==1);
s=size(ref);

tic
for point=1:10:length(pointlist)
      %generate section
     f.f_BoCount(point,1000,10,6)
     [x,y]=ind2sub(s,pointlist(point));
     refsec=ref(x-env:x+env,y-env:y+env);
     dislsec=disl(x-env:x+env,y-env:y+env);
     
     % find shift
     tshift=dftregistration(fft2(refsec),fft2(dislsec),1000);
     shiftmap(x,y,:)=tshift(3:4);

    
end

toc
figure(3);clf;
imshow(cat(1,shiftmap(:,:,1),shiftmap(:,:,2)),[])


end