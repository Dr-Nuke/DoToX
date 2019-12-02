function [ detcoords ] = f4_DetCoords( d,nd,l)
% creates the coordinates of the detector pixels
% d = hieght of the detector
% nd= number of pixels
% l =distance from rotation axis

if any([d<0, nd <=0])
    error('f4_DetCoords: detector height is <0 or pixelnumber is <=0 ')
end

detcoords(2,:)=linspace(-d/2,d/2,nd);
detcoords(1,:)=l;

end

