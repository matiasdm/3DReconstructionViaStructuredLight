%% Evaluate D3D with syntetic and real simulated data.
% 
% -------------------------------------------------------------------------
% matias di martino, matiasdm@fing.edu.uy                              2014
% -------------------------------------------------------------------------

set(0,'defaultfigureposition',[1370 1039 334 501]);

close all
clear all
home 

addpath tools/2dintegration tools/3dreconstruction 
addpath tools/misc tools/visualization 

% /////////////////////////////////////////////////////////////////////////
%% Set Experiment Parameters 
% /////////////////////////////////////////////////////////////////////////
verbose = 2; % 1-display additional text and 2-also some graphics, 

% load texture images (usefull for normalization)
id            = '_1';  % set id = '' for ExpFolder = 'Data/Ramp/
ExpFolder     = 'Data/Faces/';
Texture       = mean(double(imread([ExpFolder 'white.bmp']))  ,3);
Texture_p     = mean(double(imread([ExpFolder 'white_p.bmp'])),3);

par_display.type          = 'gray'; 
par_display.C             = Texture; 
par_display.Normalization = 1;

% Preproces (crop, resize, etc)

% /////////////////////////////////////////////////////////////////////////
%% A) Binary Coded Stripes
% /////////////////////////////////////////////////////////////////////////
fprintf('Binary Coded Stripes> ')
N = 8; % number of bit,
verbose = 2;
for n = 1:N; % load images, 
    P{n} = mean(double(imread(...
        [ExpFolder 'BinaryCodedStripes/Stripes_BinCoded_N8_' ...
        num2str(n) '_p.bmp'])),3);
    if verbose>1;
       imshow(P{n},[]); drawnow, set(gcf,'name',['BinCoded - ' num2str(n)])
       pause(.2)
    end
    
    I{n} = mean(double(imread(...
        [ExpFolder 'BinaryCodedStripes/Stripes_BinCoded_N8_' ...
        num2str(n) id '.bmp'])),3);
    if verbose>1;
       imshow(I{n},[]); drawnow, set(gcf,'name',['BinCoded - ' num2str(n)])
       pause(.2)
    end
    
    % Perform normalization 
    P{n} = P{n}./mean(Texture_p,3);
    I{n} = I{n}./mean(Texture,3);
end
% /////////////////////////////////////////////////////////////////////////
% // perform reconstruction, //////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
parameters.verbose = 2;
[D_BCS] = BinaryCodedStripes(I,P,parameters);
D_BCS   = -D_BCS;
D_BCS = medfilt2(D_BCS,[4 4]);
dmin  = -1; dmax = 26;
D_BCS(D_BCS(:)<dmin | D_BCS(:)>dmax) = NaN;
% Visualization -------------------------------------------------------
display3D(D_BCS,par_display)
[X,Y] = meshgrid([1:size(D_BCS,2)],[1:size(D_BCS,1)]);
fprintf('          Ok \n')
clear I N P n p parameters dmin dmax;

% /////////////////////////////////////////////////////////////////////////
%% C) Sin
% /////////////////////////////////////////////////////////////////////////
%fprintf('Binary Coded Stripes> ')
fprintf('PSI>                  ')
clear P I

verbose = 2;
for n = 1:3; % load images, 
    P{n} = mean(double(imread(...
        [ExpFolder 'Sin/Sin' num2str(n) '_T100_p.bmp'])),3);
    if verbose>1;
       imshow(P{n},[]); drawnow, set(gcf,'name',['PWM - ' num2str(n)])
       pause(.2)
    end
    
    I{n} = mean(double(imread(...
        [ExpFolder 'Sin/Sin' num2str(n) '_T100' id '.bmp'])),3);
    if verbose>1;
       imshow(I{n},[]); drawnow, set(gcf,'name',['Sin - ' num2str(n)])
       pause(.2)
    end
    
    % Perform normalization 
    P{n} = P{n}./mean(Texture_p,3);
    I{n} = I{n}./mean(Texture,3);
end

% /////////////////////////////////////////////////////////////////////////
% // perform reconstruction, //////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
parameters.verbose   = 0;
parameters.Normalize = 0;
[D_SIN] = PSI(I,P,parameters); 
% To use the actual PSI algorithm we need to compile the unwrapping tools. 
print('PSI algorithm requires compiling the unwrapping tools')
%D_SIN = zeros(size(I{1}));

D_SIN = -D_SIN;
D_SIN = medfilt2(D_SIN,[4 4]);
dmin  = -1; dmax = 1;
D_SIN(D_SIN(:)<dmin | D_SIN(:)>dmax) = NaN;
% Visualization -------------------------------------------------------
display3D(D_SIN,par_display)

fprintf('          Ok \n')
clear I P n p parameters;

% /////////////////////////////////////////////////////////////////////////
%% D) Takeda
% /////////////////////////////////////////////////////////////////////////
%fprintf('Binary Coded Stripes> ')
fprintf('Takeda>               ')
clear P I
strW = 4;

verbose = 1;
P = mean(double(imread([ExpFolder 'Takeda/Stripes_Tak_strW_' ...
    num2str(strW) '_p.bmp'])),3);
I = mean(double(imread([ExpFolder 'Takeda/Stripes_Tak_strW_' ...
    num2str(strW) id '.bmp'])),3);
if verbose>1;
   figure, imshow(P,[]); drawnow, set(gcf,'name',['Takeda (P)'])
   figure, imshow(I,[]); drawnow, set(gcf,'name',['Takeda (I)'])
end

P = P./mean(Texture_p,3);
I = I./mean(Texture,3);

% /////////////////////////////////////////////////////////////////////////
% // perform reconstruction, //////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
parameters.verbose = 2;
[D_TAK] = Takeda(I',P',parameters); %
% To use takeda we need to compile unwrapping tools. 
print('Takeda algorithm requires compiling the unwrapping tools')
%D_TAK = zeros(size(I)); % 

D_TAK = D_TAK';
D_TAK = medfilt2(D_TAK,[4 4]);
dmin  = -4; dmax = 5;
D_TAK(D_TAK(:)<dmin | D_TAK(:)>dmax) = NaN;
% Visualization -------------------------------------------------------
display3D(D_TAK,par_display)

fprintf('          Ok \n')
clear I P parameters;

% /////////////////////////////////////////////////////////////////////////
%% E) D3D
% /////////////////////////////////////////////////////////////////////////
fprintf('D3D>                  ')
clear P I
strW = 4;

verbose = 2;
I = double(imread([ExpFolder 'D3D/Stripes_Color_strW_' num2str(strW) id '.bmp']));

if verbose>1;
   figure, imshow(I,[]); drawnow, set(gcf,'name',['D3D (I)'])
end

% /////////////////////////////////////////////////////////////////////////
% // perform reconstruction, //////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
clear parameters;
parameters.verbose           = 2;
parameters.CodeType          = 'color'; %{sum, prod, color}
parameters.FiltWidth         = .1;
parameters.IntegrationMethod = 'AD'; %{'LS', 'WLS', 'AD', 'FC', 'L1', 'M'}
parameters.TextureIm         = Texture_p;
parameters.Theta             = pi/4;%1.42
[D_D3D,~,Dx,Dy,Qx,Qy]        = mt_D3D(I,parameters);
figure, imagesc(D_D3D), 

% Visualization -------------------------------------------------------
display3D(D_D3D,par_display)

fprintf('          Ok \n')
clear parameters I 

% /////////////////////////////////////////////////////////////////////////
%% Compare 
% /////////////////////////////////////////////////////////////////////////.
K             = 0.1;
GroundTruth   = K*mean(double(imread([ExpFolder 'GroundTruth' id '.png'])),3);
D_reg         = GroundTruth;
Mask          = GroundTruth>0;

MethodsName = {'BCS','SIN','TAK','D3D'};
D = {D_BCS, D_SIN, D_TAK, D_D3D}; 

% Evaluate error.
min_v = min(D_reg(:));
max_v = max(D_reg(:));

figure('Position',[1163 568 757 406]),
for k = 1:length(D);
    Daux          = D{k};
    Daux          = mt_Normalize(Daux,[min_v max_v]);
    %Daux(~Mask)   = 0;

    [~, mse(k)]   = SNR(Daux,D_reg); 
    rmse  = sqrt(mse)/(max_v-min_v); % root mean squared error nomalized
    
    subplot(2,3,k); imagesc(Daux); colormap gray; 
    title([MethodsName{k} ' (e = ' num2str(rmse(k)*100,'%4.1f') '%)']); colorbar;
end
subplot(2,3,k+1); imagesc(D_reg); colormap gray; colorbar; title('GroundTruth');

% Show results. 
fprintf('/////////////////////////////////////////////// \n')
fprintf(' |')
for k = 1:length(MethodsName),
        fprintf([' ' MethodsName{k} '    | '])
end, fprintf('\n |')
for k = 1:length(MethodsName),
        fprintf(' %5.2f%% | ',rmse(k)*100);
end
fprintf('\n/////////////////////////////////////////////// \n')
