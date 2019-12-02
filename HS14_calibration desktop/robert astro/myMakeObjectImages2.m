


% make map of densities, from 0-10% over 100 steps
densityMap = [1:100] * 0.2/100; % absolute fraction
densityMap = reshape(densityMap,10,10);
densityMap = densityMap';
% make 5x5 map of absent pins
noPinMap =        [0 0 0 0 0 0 0 0 0 0 
                   0 0 0 0 0 0 0 0 0 0 
                   0 0 0 0 0 0 0 0 0 0 
                   0 0 0 0 0 0 0 0 0 0 
                   0 0 0 0 1 1 0 0 0 0 
                   0 0 0 0 1 1 0 0 0 0 
                   0 0 0 0 0 0 0 0 0 0 
                   0 0 0 0 0 0 0 0 0 0 
                   0 0 0 0 0 0 0 0 0 0 
                   0 0 0 0 0 0 0 0 0 0]; 
% use filled water when no pin present, overwrite density map
densityMap(noPinMap==1)=1;              

for version = 1:3
    overallImage = [];
    [pinRows,pinCols] = size(densityMap);
    for i=1:pinRows
        thisColumn = [];
        for j=1:pinCols
            density = densityMap(i,j);
            noPin = noPinMap(i,j);
            if version==1,density = 0;end % make empty for 1st case
            if version==3,density = 1;end % make full for 3rd case
            if noPin==1,density = 0;end;
            thisAttnFlow = attnFlow * density;
            thisCell = my_makeUnitCell( pixelSize1, attnPin, thisAttnFlow , noPin);
            thisColumn = [thisColumn;thisCell];
        end
        overallImage = [overallImage thisColumn];
    end
    allImages{version} = overallImage;
end

Iempty = allImages{1}*pixelSize2/10;
Ivaried = allImages{2}*pixelSize2/10;
Ifilled = allImages{3}*pixelSize2/10;

Iempty = imresize(Iempty,pixelSize1/pixelSize2);
Ivaried = imresize(Ivaried,pixelSize1/pixelSize2);
Ifilled = imresize(Ifilled,pixelSize1/pixelSize2);


