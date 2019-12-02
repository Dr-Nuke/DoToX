

[rows,cols] = size(I);
% pixel size must evenly divide into pitch
% this gives the center effective row/col for distance from center calculations
if mod(rows,2)==0 %if rows are even
    center = rows/2+0.5;
else % else if rows are odd
    center = rows/2;
end

c = 16 ; % line intercept +c and -c in x and y
t = 1.0 ; % line thickness in mm
Iout = zeros(rows,cols);
for i=1:rows
    for j=1:cols
        % horizontal position
        x = (i-center) * pixelSize2;
        % vertical position
        y = (j-center) * pixelSize2;
        % distance from center of image
        %r = sqrt( x^2 + y^2 ) * pixelSize;
        if      y <   x+c && y >=  x+c-t && abs(x)<c && abs(y)<c ; I(i,j)=attn3;
        elseif  y <  -x+c && y >= -x+c-t && abs(x)<c && abs(y)<c ; I(i,j)=attn3;
        elseif  y >= -x-c && y <  -x-c+t && abs(x)<c && abs(y)<c ; I(i,j)=attn3;
        elseif  y >=  x-c && y <   x-c+t && abs(x)<c && abs(y)<c ; I(i,j)=attn3;
        end
    end
end

crossThickness = 1; % cross thickness in mm
crossPixels = round(crossThickness/pixelSize2); % get it in pixels
[rows,cols] = size(I);
mp = floor(rows/2); % middle pixel index
% starts top left and goes clockwise
% image shows below where A, B, and C are added in the cross
% Q1 [A] Q2
% [   B   ]
% Q4 [C] Q3
Q1 = I(     1:mp ,    1:mp  ) ; % top left[
[Q1x,Q1y]=size(Q1); % save Q1 size
Q2 = I(     1:mp , mp+1:end ) ; % top right
[Q2x,Q2y]=size(Q2); % save Q2 size
Q3 = I( mp+1:end , mp+1:end ) ; % bottom right
[Q3x,Q3y]=size(Q3); % save Q3 size
Q4 = I( mp+1:end ,    1:mp  ) ; % bottom left
[Q4x,Q4y]=size(Q4); % save Q4 size
cp = c/pixelSize2; % intercept of box in pixels
A = ones(mp,crossPixels) * attn3; % vertical line top part
A(end-cp:end,:)=0; % empty inside box
B = ones(crossPixels,rows+crossPixels)*attn3; % horizontal
B(:,mp-cp:mp+cp)=0; % empty inside box
C = ones(rows-mp,crossPixels) * attn3; % vertical line bottom part
C(1:cp,:)=0; % empty inside box
I = [[Q1 A Q2];B;[Q4 C Q3]]; % combine all the parts
% now add the outer parts
% [    ] [mid] [    ]
% [side] [ I ] [side]
% [    ] [mid] [    ]
to = 1.5; % outer walls thickness in mm
tp = round(to/pixelSize2); % get it in pixels
[rows,cols] = size(I);
Imid = ones(tp,cols)*attn3;
Iside = ones(rows+tp+tp,tp)*attn3;
I = [Iside [Imid;I;Imid] Iside];

Iout = I;
