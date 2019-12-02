function benebeach
A = imread('peppers.png');
 A = A(31:330,1:500,:);
 figure(300);clf
 imshow(A)
 title('Original Image','FontSize',14)
 conformal = maketform('custom', 2, 2, [], @conformalInverse, []);
% type conformalInverse.m
 uData = [ -1.25   1.25];  % Bounds for REAL(w)
vData = [  0.75  -0.75];  % Bounds for IMAG(w)
xData = [ -2.4    2.4 ];  % Bounds for REAL(z)
yData = [  2.0   -2.0 ];  % Bounds for IMAG(z)
B = imtransform( A, conformal, 'cubic', ...
                'UData', uData,'VData', vData,...
                'XData', xData,'YData', yData,...
                'Size', [300 360], 'FillValues', 255 );
figure(301);clf
imshow(B)
title('Transformed Image','FontSize',14)
end