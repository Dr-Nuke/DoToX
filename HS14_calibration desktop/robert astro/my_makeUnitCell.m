function [ Iout ] = my_makeUnitCell( pixelSize, attnPin, attnFlow, noPin )

% takes pixel size in mm and water density around pin to create image of
% unit cell consisting of pin and surrounding water subchannel

% makes an image of a unit cell with empty Zr pin

pitch = 12.9; % mm (rounded from 12.87)
pinOuter = 9.8; % mm
pinInner = 8.5; % mm , from 0.65 mm cladding

% pixel size must evenly divide into pitch
rows = pitch/pixelSize;
cols = rows;

% this gives the center effective row/col for distance from center calculations
if mod(rows,2)==0 %if rows are even
    center = rows/2+0.5;
else % else if rows are odd
    center = rows/2;
end

Iout = zeros(rows,cols);
for i=1:rows
    for j=1:cols
        % distance from center of image
        r = sqrt( (i-center)^2 + (j-center)^2 ) * pixelSize;
        if r<pinInner/2 % if inside pin
            % do nothing
        elseif r<pinOuter/2 % if part of cladding
            Iout(i,j) = attnPin;
        else
            Iout(i,j) = attnFlow;
        end
    end
end

if noPin==1 % overwrite everything with water if no pin
    Iout(:,:) = attnFlow;    
end

end

