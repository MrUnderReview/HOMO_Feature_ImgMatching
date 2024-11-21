%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CVPR 2025 Submission Paper ID #16689
% This code is only for the purpose of reviewing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = ShiTomasi(img, scale)

hx = [-1,0,1;-2,0,2;-1,0,1];  % 一阶梯度 Sobel算子; first-order gradient Sobel operator
hy = hx';
Gx = imfilter(img, hx, 'replicate');
Gy = imfilter(img, hy, 'replicate');

W = floor(scale/2);  % 窗半径; window radius
dx = -W : W;         % 邻域x坐标; neighborhood x-coordinate
dy = -W : W;         % 邻域y坐标; neighborhood y-coordinate
[dx,dy] = meshgrid(dx,dy);
Wcircle = ((dx.^2 + dy.^2) < (W+1)^2)*1.0;  % 圆形窗; circular window
h = fspecial('gaussian',[scale+1,scale+1], scale/6) .* Wcircle;
Gxx = imfilter(Gx.*Gx, h, 'replicate');
Gyy = imfilter(Gy.*Gy, h, 'replicate');
Gxy = imfilter(Gx.*Gy, h, 'replicate');

value = zeros(size(img,1),size(img,2));
warning off
for i=1:size(img,1)
    for j=1:size(img,2)
        try
            M = [Gxx(i,j),Gxy(i,j); Gxy(i,j),Gyy(i,j)];
            [value(i,j),~] = min(eig(M));
        catch
            value(i,j) = 0;
        end
    end
end