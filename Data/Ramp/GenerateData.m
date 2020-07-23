clc 
clear all
close all
addpath ../AuxFunctions/

% Build syntetic data;
figure('name','Z');

%% A) Generate surface 
m = 256; n = 256;
Z = zeros(m,n);
imagesc(Z), colormap jet, colorbar
[X,Y] = meshgrid([1:n],[1:m]);

% GaussianS2;
x0 = 180; y0 = 60;
aux = 2*50*50*fspecial('gaussian',[30 30],5); 
Z(y0+[1:30],x0+[1:30]) = aux;

x0 = 200; y0 = 30;
aux = 2*50*50*fspecial('gaussian',[30 30],5); 
Z(y0+[1:30],x0+[1:30]) = aux;

x0 = 160; y0 = 30;
aux = 2*50*50*fspecial('gaussian',[30 30],5); 
Z(y0+[1:30],x0+[1:30]) = aux;

% Gaussian1;
x0  = 100; y0 = 100;
%aux = 6*150*150*fspecial('gaussian',[150 150],25); 
aux = 13*150*150*fspecial('gaussian',[150 150],25); 
Z(y0+[1:150],x0+[1:150]) = aux;
clf, imagesc(Z), colormap jet, colorbar

% Ramp;
x0  = 30; y0 = 40;
aux = .3*Y(1:200,1:80); 
Z(y0+[1:200],x0+[1:80]) = aux;
clf, imagesc(Z), colormap jet, colorbar

figure,
p.Type = 'jet';
display3D(Z,p)
clear aux x0 y0 

I = uint8(Z);
imwrite(uint8(cat(3,I,I,I)),['GroundTruth.bmp'],'bmp');
    

K = .4;
D = K*Z;
    

%% aux
%dx = .5;
%dy = .5;
%D = dx*X + dy*Y;
%display3D(D,p);

%% B) Simulate Patterns;
% I) BCS 
h1 = figure('name','plane');
h2 = figure('name','test');

for k = 1:8;
    P = mean(double(imread(['BinaryCodedStripes/Stripes_BinCoded_N8_' num2str(k) '_p.bmp'])),3);
    
    [I] = SimulateProjection(P,D);
    
    imwrite(uint8(cat(3,I,I,I)),...
        ['BinaryCodedStripes/Stripes_BinCoded_N8_' num2str(k) '.bmp'],'bmp');
    
end

% II) Simulate SIN;
for k = 1:3;
    P = mean(double(imread(['Sin/Sin' num2str(k) '_T100_p.bmp'])),3);
    I = SimulateProjection(P,D);
    
    imwrite(uint8(cat(3,I,I,I)),...
        ['Sin/Sin' num2str(k) '_T100.bmp'],'bmp');
    
end

% III) Simulate TAK;
P = mean(double(imread(['Takeda/Stripes_Tak_strW_4_p.bmp'])),3);
I = SimulateProjection(P,D);

imwrite(uint8(cat(3,I,I,I)),...
    ['Takeda/Stripes_Tak_strW_4.bmp'],'bmp');

%% C) Simulate D3D;
% I) Sum,
P = mean(double(imread(['D3D/Stripes_Sum_strW_4_p.bmp'])),3);
I = SimulateProjection(P,D);
imwrite(uint8(cat(3,I,I,I)),['D3D/Stripes_Sum_strW_4.bmp'],'bmp');

% II) Prod
P = mean(double(imread(['D3D/Stripes_Pro_strW_4_p.bmp'])),3);
I = SimulateProjection(P,D);
imwrite(uint8(cat(3,I,I,I)),['D3D/Stripes_Pro_strW_4.bmp'],'bmp');

% III) Color
P = double(imread(['D3D/Stripes_Color_strW_4_p.bmp']));
I = SimulateProjection(P,D);
imwrite(uint8(I),['D3D/Stripes_Color_strW_4.bmp'],'bmp');

%a = ( -dx/sqrt(2) ) / ( 1 - dy/sqrt(2))
%verpendiente









%


