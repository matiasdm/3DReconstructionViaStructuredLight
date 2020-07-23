clc 
clear all
close all
addpath ../AuxFunctions/

% Build syntetic data;
figure('name','Z');

close all
m = 751; 
n = 501;

I = uint8(255*ones(m,n));
imwrite(cat(3,I,I,I),'white.bmp','bmp');
imwrite(cat(3,I,I,I),'white_p.bmp','bmp');
    
%% B) Simulate Patterns;
% I) BCS
h1 = figure('name','plane');
h2 = figure('name','test');

for k = 1:8;
    P = mean(double(imread(['BinaryCodedStripes/Stripes_BinCoded_N8_' num2str(k) '_p.bmp'])),3);
    I = imresize(P,[m n],'bicubic');
    
    imwrite(uint8(cat(3,I,I,I)),...
        ['BinaryCodedStripes/Stripes_BinCoded_N8_' num2str(k) '_p.bmp'],'bmp');
    
end

% II) Simulate SIN;
for k = 1:3;
    P = mean(double(imread(['Sin/Sin' num2str(k) '_T100_p.bmp'])),3);
    I = imresize(P,[m n],'bicubic');
    
    imwrite(uint8(cat(3,I,I,I)),...
        ['Sin/Sin' num2str(k) '_T100_p.bmp'],'bmp');
    
end

% III) Simulate TAK;
P = mean(double(imread(['Takeda/Stripes_Tak_strW_4_p.bmp'])),3);
I = imresize(P,[m n],'bicubic');

imwrite(uint8(cat(3,I,I,I)),...
    ['Takeda/Stripes_Tak_strW_4_p.bmp'],'bmp');

%% C) Simulate D3D;
% I) Sum,
P = mean(double(imread(['D3D/Stripes_Sum_strW_4_p.bmp'])),3);
I = imresize(P,[m n],'bicubic');
imwrite(uint8(cat(3,I,I,I)),['D3D/Stripes_Sum_strW_4_p.bmp'],'bmp');

% II) Prod
P = mean(double(imread(['D3D/Stripes_Pro_strW_4_p.bmp'])),3);
I = imresize(P,[m n],'bicubic');
imwrite(uint8(cat(3,I,I,I)),['D3D/Stripes_Pro_strW_4_p.bmp'],'bmp');

% III) Color
P = double(imread(['D3D/Stripes_Color_strW_4_p.bmp']));
I = imresize(P,[m n],'bicubic');
imwrite(uint8(I),['D3D/Stripes_Color_strW_4_p.bmp'],'bmp');

%a = ( -dx/sqrt(2) ) / ( 1 - dy/sqrt(2))
%verpendiente



