%% [DeltaPhi] = mt_Takeda3DProfilometry(I,I0,parameters) 
% function that implements takeda's aproach for 3D profilometry [1]. 
%
% Inputs: 
%   - I (mxnx1) gray image with horizontal fringe prattern, 
%   - I0 (mxnx1) gray image of the pattern without object (calibration image)
%                if we don't have a calibration image we can use I0 = 0 and
%                a calibration image is estimated using fringes period.
%   - parameters: struct that may contain, 
%       .verbose {def 0} [1 - display some text, 2 also some graphics]
%
% Outputs:
%   - DeltaPhi (mxnx1) wrapped 3d profile
%
% Refs:
%   [1] Mitsou Takeda and Kazuhiro Mutoh. "Fourier transforma profilometry 
%     for the automatic measurement of 3D object shapes."
%
% -------------------------------------------------------------------------
% matias di martino, matiasdm@fing.edu.uy                            (2013)
% -------------------------------------------------------------------------

function [DeltaPhi] = mt_Takeda3DProfilometry(I,I0,parameters), 

% check inputs and parameters, 
I = mean(I,3); % convert to gray if is colored image, 

if I0==0,
    [I0] = EstimateCalibrationImage(I);
else %use calibration image provided by user,
    I0 = mean(I0,3); % convert to gray if is colored image, 
end

if isfield(parameters,'verbose'),
    verbose = parameters.verbose;
else % set default value
    verbose = 0;
end

%% 
if verbose>1, %display intput images, 
    figure('Position',[267   539   767   224]); 
    subplot(1,2,1), imagesc(I), axis image, axis off, colormap gray, colorbar, title('I');
    subplot(1,2,2), imagesc(I0), axis image, axis off, colormap gray, colorbar, title('I0');
end
    
%% Define Fourier domains 

[m n] = size(I);

ft2 = @(u) fftshift(fft2(u));
ft = @(u) fftshift(fft(u));
ift2 = @(u) ifft2(ifftshift(u));
ift = @(u) ifft(ifftshift(u));

f0x = round(n/2)+1; f0y = round(m/2)+1; % center of spectrum

for x = 1:n, % for each column, perform fft, 
    G(:,x) = ft(I(:,x));
    G0(:,x) = ft(I0(:,x));
end

if verbose>1,
    figure,
    subplot(1,2,1), imagesc(log(abs(G0))), title('G0'), colormap jet;
    subplot(1,2,2), imagesc(log(abs(G))), title('G'), colormap jet;
end

%% Filter to keep just Q1

% find secondary peak:
AuxMask = ones(size(G0));
AuxMask(f0y+[-10:10],:) = 0; % remove DC  

% find secondary peak (higher after DC value)
[~,inds] = max(abs(G0.*AuxMask));

% to make this more robust, keep the mode
fy1 = mode(inds( inds(:)<f0y ));

% filter;
width = abs( (f0y-fy1) ); % filter width 
AuxMask = zeros(m,n);
AuxMask(f0y+round(width/2)+[1:width],:) = 1;

G0 = G0.*AuxMask;
G = G.*AuxMask;

if verbose>1,
    figure,
    subplot(1,2,1), imagesc(log(abs(G0))), colormap jet;
    subplot(1,2,2), imagesc(log(abs(G))), colormap jet;
end
%% Inversion 

for x = 1:n, % for each column, perform fft, 
    g(:,x) = ift(G(:,x));
    g0(:,x) = ift(G0(:,x));
end

DeltaPhi = imag( log( g.*conj(g0) ) );

figure, imagesc(DeltaPhi), colormap jet, colorbar, 

end


%% Auxiliary functions, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% /////////////////////////////////////////////////////////////////////////
function [I0] = EstimateCalibrationImage(I); % ////////////////////////////
% /////////////////////////////////////////////////////////////////////////
verbose = 0;

% Define some shortcuts for image transformations -----
ft2 = @(u) fftshift(fft2(u));
ft = @(u) fftshift(fft(u));
ift2 = @(u) ifft2(ifftshift(u));
ift = @(u) ifft(ifftshift(u));
% -----------------------------------------------------

[m n] = size(I);
FI    = ft2(I); % fourier transform of gray image

% /////////////////////////////////////////////////////////////////////////
% /// Split "DC and AC" components ////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////

% zero frec 
i0  =  ceil((m+1)/2);
j0  =  ceil((n+1)/2);

% define the filter (hanning, rectangular or any other shape of filter
%                   can be used)
FiltWidthx   = 110; 
FiltWidthy   = round(FiltWidthx*m/n); % respect aspect ratio when filtering
Filt         = hanning(FiltWidthy)*hanning(FiltWidthx)';

DCFilt       = zeros(m,n); 
DCFilt(i0-round(FiltWidthy/2)+[1:FiltWidthy],...
       j0-round(FiltWidthx/2)+[1:FiltWidthx]) = Filt;

FIDC         = FI.*DCFilt;
FAC          = FI.*(1-DCFilt);
    
if verbose>1, % show filters,  
    figure('name','AC-DC'), imagesc(log(abs(FI))), colormap jet; 
    AuxLPFImage = cat(3,zeros(m,n),zeros(m,n),ones(m,n));
    AuxBPFImage = cat(3,ones(m,n),zeros(m,n),zeros(m,n));
    hold on, imagesc(AuxLPFImage,'AlphaData',DCFilt*.5), 
    hold on, imagesc(AuxBPFImage,'AlphaData',(1-DCFilt)*.5),
end

% /////////////////////////////////////////////////////////////////////////
% /// From "AC" image, extraxt vertical and horizontal fringes position ///
% /// on fourier domain                                                 ///
% /////////////////////////////////////////////////////////////////////////    
% Find vertical maximum  --------------------------------------------------
EnergyOnYAxe = mean(abs(FAC),2); 
% remove the efect of the orthogonal fringes (which add enercy near Y=0)
EnergyOnYAxe(i0-round(FiltWidthy/4):end) = 0;
[~,i1] = max(EnergyOnYAxe);
% -------------------------------------------------------------------------
% Find horizontal maximum  ------------------------------------------------
EnergyOnXAxe = mean(abs(FAC),1); 
% remove the efect of the orthogonal fringes (which add enercy near Y=0)
EnergyOnXAxe(j0-round(FiltWidthx/4):end) = 0;
[~,j1] = max(EnergyOnXAxe);
% -------------------------------------------------------------------------
   
if verbose>1,
    figure, imagesc(log(abs(FI))), colormap jet; hold on
    scatter(j1,i0,30,'b','filled');
    scatter(j0,i1,30,'b','filled');
end

% We are assuming horizontal fringes! 
fw = round( 1/(i0-i1) * m ); 

if verbose>0, disp(['Estimated fringes width :  ' num2str(fw) '(pixels)']), end

[X Y] = meshgrid([1:n],[1:m]);
I0 = cos(Y*2*pi/fw);

end
