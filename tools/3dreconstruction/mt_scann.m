%% Aux function,

function [X,Y,Z,C,M] = mt_scann(Iin,CalibrationImage);

addpath ../../matuls/ ../../otherstoolbox/ ../images ../adquiridas/ ../calibration/
tic
verbose = 2;

%% Read Test Image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
fprintf('Loading image ... ');

m = size(Iin,1); n = size(Iin,2);

% show input image,
pos  = [256        1007         640         512];
pos2 = pos;
if verbose>1,
    h = figure('name','Input','NumberTitle','off','Position',pos); imshow(Iin,'InitialMagnification',100);
    h2 = figure('name',' ','NumberTitle','off','Position',pos2); 
end

% release memory ----------------------------------------------
% clear ; % clear auxiliary variables,
% -------------------------------------------------------------

fprintf('  Ok \n');

%% Preprocesing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Preprocesing and image alignment ...')
% using Calibration image, fix posible misalignment
% perform fft, 
FT = fftshift(fft2(mean(CalibrationImage,3)));

% Remove zero frec,
[aux,i0] = max(abs(FT));
[~,j0] = max(aux'); i0 = i0(j0);
FT(i0+[-40:40],:) = 0; 

% find max (correspond to fringes)
[aux,i1] = max(abs(FT));
[~,j1] = max(aux'); i1 = i1(j1);
di = abs(i1-i0); 
dj = j1-j0;
Cangle = pi/2 - atan2(di,dj);
Cangle = Cangle * 180 / pi; % convert to degrees,

Iin = imrotate(Iin,-Cangle,'bilinear');
CalibrationImage = imrotate(CalibrationImage,-Cangle,'bilinear');

% /////////////////////////////////////////////////////////////////////////
% ///// Estimation of fringe width ////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
fw = round( 1/di * size(CalibrationImage,1) ); % 

% /////////////////////////////////////////////////////////////////////////
% ///// Convert to gray and normalize /////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
In = mt_Normalize(double(Iin),[0 1]);

% release memory ----------------------------------------------
clear aux j0 i0 i1 j1 di dj Cangle FT ; % clear auxiliary variables,
% -------------------------------------------------------------

fprintf('    Ok \n')

%% Segment Face %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Face segmentation ... ');

% /////////////////////////////////////////////////////////////////////////
% /////////// Segment face, ///////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
% set smoothing kernel ---------
KSize = ceil(fw*[3 3]); KSigma = .5*fw;
K = fspecial('gaussian',KSize,KSigma);
% ------------------------------

IGray = mean(In,3);
IGray = conv2(IGray,K,'same');
Mask = IGray>3*median(IGray(:));

% Fill holes in the mask,
Mask = imfill(Mask,'holes');
% Keep just the central part,
Mask = bwselect(Mask,round(m/2),round(n/2));
% Errode to remove borders (that introduce noise)

se = strel('disk',ceil(4*fw));  
Mask = imdilate(Mask,se);
se = strel('disk',ceil(4*fw));  
Mask = imerode(Mask,se);

if verbose>1, 
    figure(h),set(h,'name','segmentation mask'), imshow(Mask,[]), drawnow
end

% Just keep usefull range of pixels, --------------------------------------
%   Remove the columns and rows that have all NaNs,
ValidRows = sum(Mask,2)>0; % is valid if Mask intersect that Row
In = In(ValidRows,:,:);
Mask = Mask(ValidRows,:);
ValidColumns = sum(Mask,1)>0; % is valiz if Mask intersect that Row
In = In(:,ValidColumns,:);
Mask = Mask(:,ValidColumns);
% -------------------------------------------------------------------------

% release memory ----------------------------------------------
clear K KSize KSigma KCol KSizeCol KSigmaCol; % clear auxiliary variables,
clear ValidRows ValidColumns se IGray 
% -------------------------------------------------------------

fprintf('  Ok \n');

%% 3D retrieval by our method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('3D retrieval by our ... ');
parameters.Type = 'CosineExpansion'; % {Gauss, Fourier,'SineExpansion'}
parameters.FilterType = 5; %{1,2,3,4 or 5}
parameters.WinWidth = [1 1];
parameters.verbose = verbose;
[DispOur,Texture, Disp_x, Disp_y ] = mt_3DRetrieveBoth(In,parameters);


% release memory ----------------------------------------------
clear parameters ; % clear auxiliary variables,
% -------------------------------------------------------------

fprintf('  Ok \n');

%% Filtering, postprocesing and visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Filtering, postprocesing and visualization ... '); 

Z = DispOur;
Z = medfilt2(Z,[15 15]);

% set smoothing kernel and filter Z ---------
% KSize = [5 5]; KSigma = 1;
% K = fspecial('gaussian',KSize,KSigma);
% Z = conv2(Z,K,'same');
% -------------------------------------------


% Apply mask obtain, 
Z(Mask==0) = NaN; 

% Define mesh colors, (as the original smoothed image), 
C = Texture;

% /////////////////////////////////////////////////////////////////////////
% //// Display results ////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////

% Define the X,Y grid, 
[X,Y] = meshgrid([1:size(Z,2)],[1:size(Z,1)]);

% Normalize Data ----------------------------------------------------------
X = mt_Normalize(X,[0 1]);
Y = mt_Normalize(-Y,[0 1]); % and flip vertical axe 
Z = mt_Normalize(Z,[0 .8]);
C = mt_Normalize(C,[0 1]);
M = Mask;
% -------------------------------------------------------------------------

end

