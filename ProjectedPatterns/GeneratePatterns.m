%% Generate patterns to be projected.
%
% matias di martino, matiasdm@fing.edu.uy                              2014    

clc 
close all
clear all

% Set variables, 
W =  1024; % set projector width resolution 
H =  768 ; % set projector heihgt resolution 

%% 0) Center, used for calibration
P      = zeros(H,W);
A      = zeros(H,W,3);
center = [round(W/2) round(H/2)];
[X,Y]  = meshgrid([1:W],[1:H]);
r = 3;
P(sqrt((X-center(1)).^2+(Y-center(2)).^2)<r) = 255;
A = cat(3,P,P,P);
A(center(2),:,1) = 255;
S = floor(min(H,W)/2);
for k = -S+1:S-1,
    A(k+center(2),k+center(1),2) = 255;
    A(-k+center(2),k+center(1),2) = 255;
end
A(:,center(1),3) = 255;
imwrite(uint8(A),['Calibration.bmp'],'bmp');
figure('name','Calibration'); imshow(uint8(A))



%% A) Horizontal Stripes plus Vertical Stripes.
close all
P            = zeros(H,W); % pattern inicialization; 
strW         = 4;  % set stripes width (in pixels),
Level        = 250; % maximum level 
% horizontal stripes,
Ph = P; % inicialization
numOfStripes = floor( H / (2*strW) ); 
for i = 0:numOfStripes-1,
    Ph(i*strW*2+[1:strW],:) = Level;
end

% vertical stripes,
Pv = P; % inicializatino 
numOfStripes = floor( W / (2*strW) ); 
for i = 0:numOfStripes-1,
    Pv(:,i*strW*2+[1:strW]) = Level;
end

%P = 1/2*(Ph+Pv);
%imwrite(uint8(P),['Stripes_Sum_strW_' num2str(strW) '_Level' num2str(Level) '.bmp'],'bmp');
%figure('name','Stripes [sum]'); imshow(uint8(P))
P = sqrt(Pv.*Ph); 
imwrite(uint8(P),['Stripes_Pro_strW_' num2str(strW) '_Level' num2str(Level) '.bmp'],'bmp');
figure('name','Stripes [prod]'); imshow(uint8(P))

%imwrite(uint8(cat(3,Ph,zeros(H,W),Pv)),['Stripes_Color_strW_' num2str(strW) '_Level' num2str(Level) '.bmp'],'bmp');
%figure('name','Stripes color'); imshow(uint8(cat(3,Ph,zeros(H,W),Pv)))


clear numOfStripes Ph Pv;

%% B) Binary codded stripes 
if 0,
N = 8; % set num of bits

P = ones(H,1)*mod([1:W],2^N);
figure('name','aux'); imshow(uint8(P)); 

Pn = P; % inicialization 
for n = 1:N;
    for i = 1:H;
        for j = 1:W;
            aux = dec2bin(P(i,j),N);
            Pn(i,j) = 255*str2double(aux(n));
        end
    end 
    imwrite(uint8(Pn),['Stripes_BinCoded_N' num2str(N) '_' num2str(n) '.bmp'],'bmp');
    clf; imshow(uint8(Pn)); drawnow
end
end
%% C) PWM
if 1,
% define and save the sinusoidal images (for comparision)
[X,Y] = meshgrid([1:W],[1:H]);
f0 = 1/50;
S1 = 255/2*( 1 + sin(2*pi*f0*X-2*pi/3) );
S2 = 255/2*( 1 + sin(2*pi*f0*X)        );
S3 = 255/2*( 1 + sin(2*pi*f0*X+2*pi/3) );
imshow(uint8(S1))
imwrite(uint8(S1),['Sin1_T' num2str(uint8(1/f0)) '.bmp'],'bmp');
imwrite(uint8(S2),['Sin2_T' num2str(uint8(1/f0)) '.bmp'],'bmp');
imwrite(uint8(S3),['Sin3_T' num2str(uint8(1/f0)) '.bmp'],'bmp');

fc  = 7*f0;
Tc  = round(1/fc); % period of the triangles function 
aux = [[1:1:round(Tc/2)] [round(Tc/2)-1:-1:1]];
ind = 1+mod(X,length(aux));
Comp = aux(ind);
Comp = 2*( ( Comp-min(Comp(:)) ) / (max(Comp(:))-min(Comp(:))) ) - 1;

figure, plot(Comp(1,:));

S1 = 255 * ( sin(2*pi*f0*X-2*pi/3) > Comp );
S2 = 255 * ( sin(2*pi*f0*X)        > Comp );
S3 = 255 * ( sin(2*pi*f0*X+2*pi/3) > Comp );
imwrite(uint8(S1),['PWMSin1_T' num2str(uint8(1/f0)) '.bmp'],'bmp');
imwrite(uint8(S2),['PWMSin2_T' num2str(uint8(1/f0)) '.bmp'],'bmp');
imwrite(uint8(S3),['PWMSin3_T' num2str(uint8(1/f0)) '.bmp'],'bmp');
end

%% D) Takeda
P  = zeros(H,W);
% horizontal stripes,
Ph = P; % inicialization
numOfStripes = floor( H / (2*strW) ); 
for i = 0:numOfStripes-1,
    Ph(i*strW*2+[1:strW],:) = 255;
end

P = Ph;
imwrite(uint8(P),['Stripes_Tak_strW_' num2str(strW) '.bmp'],'bmp');
figure('name','Stripes [Tak]'); imshow(uint8(P))

clear numOfStripes Ph;


