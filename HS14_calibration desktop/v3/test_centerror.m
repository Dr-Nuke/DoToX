       

%% step by step cenering to find out that integer errör
k=1030;
sino=squeeze(block(:,k,:));

inte=fix(shift(k));
        frac=shift(k)-fix(shift(k));

        if fix(shift(k))>0 %in case of positive shift value, shifting to the left
            sino1=[sino;repmat(sino(end,:),abs(inte),1)]; %integer shift to center
            bs=size(sino1); %blocksize

            if mod(size(sino1,1),2) % if the number of pixels is odd
                [XI,YI]=ndgrid(1:bs(1),1:bs(2));
                sino2=interpn(XI,YI,sino1,XI,YI+double(frac),'cubic'); %subgrid shift
            else % if it is even, make it odd
                [XI,YI]=ndgrid(1:bs(1)+1,1:bs(2));
                sino2=interpn(XI,YI,[sino1;sino1(end,:)],XI,YI+double(frac)-0.5,'cubic'); %subgrid shift
            end %this is to get always odd sized projections, that preserves centeredness in iradon

        else %in case of negative shift value, shifting axis to the right
            sino1=[repmat(sino(1,:),abs(inte),1);sino]; %integer shift to center
            bs=size(sino1); %blocksize

            if mod(size(sino1,1),2) % if the number of pixels is odd
                [XI,YI]=ndgrid(1:bs(1),1:bs(2));
                sino2=interpn(XI,YI,sino1,XI,YI+double(frac),'cubic'); %subgrid shift
            else % if it is even, make it odd
                [XI,YI]=ndgrid(1:bs(1)+1,1:bs(2));
                sino2=interpn(XI,YI,[sino1;sino1(end,:)],XI,YI+double(frac)-0.5,'cubic'); %subgrid shift
            end %this is to get always odd sized projections, that preserves centeredness in iradon
        end

%%
figure(20)
clf
figure(19)
clf

for i = 1000
    
    disp(i)
        sino=squeeze(block(:,i,:)); %sinogram
        sino(30:31,20:40)=1;
        sino(20:40,30:31)=1;
        sino=f3_center(sino,shift(i));
        imshow(sino(1:50,1:50),[0,1],'InitialMagnification','fit');
        title(i)
        drawnow
        figure(20)
        imshow(sino,[0,1]);
end