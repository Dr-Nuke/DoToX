function [ sino ] = f3_center( sino,shift )
% this function performs the sub pixel shift with a given sinogram
        inte=fix(shift);
        frac=shift-fix(shift);

        if fix(shift)>0 %in case of positive shift value, shifting to the left
            sino=[sino;repmat(sino(end,:),abs(inte),1)]; %integer shift to center
            bs=size(sino); %blocksize

            if mod(size(sino,1),2) % if the number of pixels is odd
                [XI,YI]=ndgrid(1:bs(1),1:bs(2));
                sino=interpn(XI,YI,sino,XI+double(frac),YI,'cubic'); %subgrid shift
            else % if it is even, make it odd
                [XI,YI]=ndgrid(1:bs(1)+1,1:bs(2));
                sino=interpn(XI,YI,[sino;sino(end,:)],XI+double(frac)-0.5,YI,'cubic'); %subgrid shift
            end %this is to get always odd sized projections, that preserves centeredness in iradon

        else %in case of negative shift value, shifting axis to the right
            sino=[repmat(sino(1,:),abs(inte),1);sino]; %integer shift to center
            bs=size(sino); %blocksize

            if mod(size(sino,1),2) % if the number of pixels is odd
                [XI,YI]=ndgrid(1:bs(1),1:bs(2));
                sino=interpn(XI,YI,sino,XI+double(frac),YI,'cubic'); %subgrid shift
            else % if it is even, make it odd
                [XI,YI]=ndgrid(1:bs(1)+1,1:bs(2));
                sino=interpn(XI,YI,[sino;sino(end,:)],XI+double(frac)-0.5,YI,'cubic'); %subgrid shift
            end %this is to get always odd sized projections, that preserves centeredness in iradon
        end


end

