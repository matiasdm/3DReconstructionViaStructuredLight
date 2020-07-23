%% [Z,[Zx]] = mt_3DRetrieve(I,parameters) --------------------------------------
% function that performs differential 3D retrival (I = image adquiered with
% projected fringe patterns)
%
% Inputs:
%   - I (mxnxc) gray or color image, 
%   - paremeters, struct that may contain:
%       .verbose [def 0] (1 some text is display, 2 also some graphics.)
%       .Type [def Gauss]
%           - Gauss: use gauss-seidel interation (on integrateCPP),
%           - Fourier: integration on fourier domain,
%       .MaxIter [def 2e4] (Max iteration on iterative algorithms)
%       .Tol [def 1e-4] (stop criteria on iterative algorithms)
%       .WinWidth [def [2 1]] (width of the windows in median filter (in
%                              terms of the fringes width)
%
% Outputs:
%   - Z (mxnx1) 
%   - Zx (mxnx1) optional output
% Refs:
% 
% -------------------------------------------------------------------------
% matias di martino (c) - matiasdm@fing.edu.uy - 2013
% -------------------------------------------------------------------------


function [Z,varargout] = mt_3DRetrieve(I,parameters);

%% Load and set general parameters, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(parameters,'verbose')
    verbose = parameters.verbose;
else % Set default value, 
    verbose = 0;
end

if isfield(parameters,'Type'),
    Type = parameters.Type;
else % set default value, 
    Type = 'Gauss';
end

if isfield(parameters,'MaxIter')
    MaxIter = parameters.MaxIter;
else % Set default value, 
    MaxIter = 2e4;
end

if isfield(parameters,'Tol')
    Tol = parameters.Tol;
else % Set default value, 
    Tol = 1e-4;
end

if isfield(parameters,'WinWidth'),
    w = parameters.WinWidth;
else % set default value, 
    w = [2 1];
end


% check input, ------------
I = double(I); I = mean(I,3); % convert to gray if is color image, 
I = mt_Normalize(I,[0 1]); % normalization, 
% -------------------------

[m n] = size(I);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute x-partial derivative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /////////////////////////////////////////////////////////////////////////
% //// A) split high and low frecuency components /////////////////////////
% /////////////////////////////////////////////////////////////////////////
% A0) Estimate fringes width --------------------
FT = fftshift(fft2(I));
% Remove zero frec,
% Remove zero frec,
[aux,i0] = max(abs(FT));
[~,j0] = max(aux'); i0 = i0(j0);
FT(i0+[-10:10],j0+[-10:10]) = 0; 

% find max (correspond to fringes)
[aux,i1] = max(abs(FT));
[~,j1] = max(aux'); i1 = i1(j1);

di = abs(i1(1)-i0(1)); 

fw = round( 1/di * m ); 

if verbose>0, disp(['Estimated fringes width :  ' num2str(fw) '(pixels)']), end
% -----------------------------------------------

% A1) Filter on Fourier domain --------------------------------------------
% Define some shortcuts for image transformations -----
ft2 = @(u) fftshift(fft2(u));
ft = @(u) fftshift(fft(u));
ift2 = @(u) ifft2(ifftshift(u));
ift = @(u) ifft(ifftshift(u));
% -----------------------------------------------------

f0x = round(n/2)+1; f0y = round(m/2)+1; % center of spectrum
FTI = ft2(I); 
% Filter in fourier domain ----------------------------------
% so we have fringe on one side and low frecuencies (asociated with object
% reflectance) on the other.

LPF = zeros(m,n); % LowPassFilter 
LPF(f0y+round([-m/fw:m/fw]/2),:) = 1;

BPF = zeros(m,n); % BandPassFilter 
BPF(f0y+round(m/fw)+round([-m/fw:m/fw]/2),:) = 1;
BPF(f0y-round(m/fw)+round([-m/fw:m/fw]/2),:) = 1;
%BPF(1:f0y+round(-m/fw/2),:) = 1;
%BPF(f0y+round(+m/fw/2):end,:) = 1;

if verbose>1, % show filters,  
    figure, imagesc(log(abs(FTI))), colormap jet; 
    AuxLPFImage = cat(3,zeros(m,n),zeros(m,n),ones(m,n));
    AuxBPFImage = cat(3,ones(m,n),zeros(m,n),zeros(m,n));
    hold on, imagesc(AuxLPFImage,'AlphaData',LPF*.5), 
    hold on, imagesc(AuxBPFImage,'AlphaData',BPF*.5),
end

LowFrecI = real(ift2(LPF.*FTI));
FringeI = real(ift2(BPF.*FTI));

if verbose>1,
    figure, imshow(FringeI,[]); drawnow,
end

% A2) Compute x-derivative ------------------------------------------------
% Compute Image derivatives
clear p
p.verbose = 0;
p.Type = 'Centered'; %{'Centered','Forward','FluxLimiter','RotationInvariant'...
                    % 'Sobel','Fourier'};
[Ix,Iy] = mt_ImDerivate(FringeI,p); 

Zx = -Ix./(Iy+eps);

if verbose>1, 
    figure('name','sign(Disp_x).*log(abs(Disp_x)+1)')
    imagesc(sign(Zx).*log(abs(Zx)+1)), colormap jet, colorbar, drawnow,
end

% -------------------------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remove outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /////////////////////////////////////////////////////////////////////////
% //// B) Integrate  //////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
Zx = medfilt2(Zx,[w(1)*fw w(2)*fw]);

if verbose>1,
     figure('name','Zx (with med filt)')
     imagesc(Zx), colormap jet, colorbar, drawnow,
end
varargout{1} = Zx;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Integrate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /////////////////////////////////////////////////////////////////////////
% //// C) Integrate  //////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
switch Type,
    case 'Gauss',
        addpath ../../matulsC++ % add path of cpp tools
        tic
        [Z,flag] = mt_IntegrateCPP(Zx,MaxIter,Tol);
        t = toc; if verbose>0, mt_printtime(t), end
    otherwise,
        p = parameters; % pass all input parametes,   
        Ix = Zx; Iy = 0*Zx; 
        [Z] = mt_ImIntegrate2(Ix,Iy,p);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%