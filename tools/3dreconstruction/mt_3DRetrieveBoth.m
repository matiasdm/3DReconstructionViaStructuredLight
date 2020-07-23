%% [Z,[Texture],[Zx],[Zy]] = mt_3DRetrieve(I,parameters) --------------------------------------
% function that performs differential 3D retrival (I = image adquiered with
% projected fringe patterns)
%
% Inputs:
%   - I1 (mxnxc) gray or color image with horizontal plus vertical fringes
%   - paremeters, struct that may contain:
%       .verbose [def 0] (1 some text is display, 2 also some graphics.)
%       .Type [def Gauss]
%           - Gauss: use gauss-seidel interation (on integrateCPP),
%           - Fourier: integration on fourier domain,
%           - CosineExpansion: using cosine transf props (on CPP),
%       .FilterType {1,2,3,4 or 5}
%       .MaxIter [def 2e4] (Max iteration on iterative algorithms)
%       .Tol [def 1e-4] (stop criteria on iterative algorithms)
%       .WinWidth [def [2 2]] (width of the windows in median filter (in
%                              terms of the fringes width)
%
% Outputs:
%   - Z (mxnx1) 
%   - Texture (mxnxc); optional output
%   - Zy (mxnx1) optional output
%   - Zx (mxnx1) optional output
% Refs:
% 
% -------------------------------------------------------------------------
% matias di martino (c) - matiasdm@fing.edu.uy - 2013
% -------------------------------------------------------------------------


function [Z,varargout] = mt_3DRetrieveBoth(I1,parameters);

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
    w = [2 2];
end

if isfield(parameters,'FilterType');
    FilterType = parameters.FilterType;
else % set default value,
    FilterType = 3;
end


% check input, ------------
Icolor = I1;
I1 = double(I1); I1 = mean(I1,3); % convert to gray if is color image, 
I1 = mt_Normalize(I1,[0 1]); % normalization, 
% -------------------------

[m n col] = size(Icolor);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute x-partial derivative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /////////////////////////////////////////////////////////////////////////
% //// A) split high and low frecuency components /////////////////////////
% /////////////////////////////////////////////////////////////////////////
I = I1;
% A0) Estimate fringes width --------------------
FT = fftshift(fft2(I));
if nargout>1, % if texture image is required,
    % inicialization, 
    FTcol = zeros([m n col]);
    for c = 1:col
        FTcol(:,:,c) = fftshift(fft2(Icolor(:,:,c)));
    end
end

% Remove zero frec,
% Remove zero frec,
waux = 50;
[aux,i0] = max(abs(FT));
[~,j0] = max(aux'); i0 = i0(j0);
FT(i0+[-waux:waux],:) = 0; 

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
LPF(f0y+round([-m/fw:m/fw]/2),f0x+round([-n/fw:n/fw]/2)) = 1;

BPF = zeros(m,n); % BandPassFilter 
switch FilterType,
    case 1,
        BPF(f0y+round(m/fw)+round([-m/fw:m/fw]/2),:) = 1;
        BPF(f0y-round(m/fw)+round([-m/fw:m/fw]/2),:) = 1;
    case 2,
        BPF(1:f0y+round(-m/fw/2),:) = 1;
        BPF(f0y+round(+m/fw/2):end,:) = 1;
    case 3,
        BPF(f0y+round(m/fw)+round([-m/fw:m/fw]/2),f0x+round([-n/fw:n/fw]/2)) = 1;
        BPF(f0y-round(m/fw)+round([-m/fw:m/fw]/2),f0x+round([-n/fw:n/fw]/2)) = 1;
    case 4,
        BPF(f0y+round(m/fw)+round([-m/fw:m/fw]/2),f0x+round([-n/fw:n/fw]/2)) = 1;
    case 5,
        W = 40;
        BPF(f0y+round(m/fw)+round([-m/fw+W:m/fw-W]/2),f0x+round([-n/fw+W:n/fw-W]/2)) = 1;
        BPF(f0y-round(m/fw)+round([-m/fw+W:m/fw-W]/2),f0x+round([-n/fw+W:n/fw-W]/2)) = 1;
end

if verbose>1, % show filters,  
    figure, imagesc(log(abs(FTI))), colormap jet; 
    AuxLPFImage = cat(3,zeros(m,n),zeros(m,n),ones(m,n));
    AuxBPFImage = cat(3,ones(m,n),zeros(m,n),zeros(m,n));
    hold on, imagesc(AuxLPFImage,'AlphaData',LPF*.5), 
    hold on, imagesc(AuxBPFImage,'AlphaData',BPF*.5),
end

FringeI = real(ift2(BPF.*FTI));

% inicialization,
if nargout>1, % if texture image is required,
    Texture = zeros([m n col]);
    for c = 1:col
        Texture(:,:,c) = real(ift2(LPF.*FTcol(:,:,c)));
    end
    varargout{1} = Texture;
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove Texture %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AuxTexture = mt_Normalize( mean(Texture,3),[0 1]);
% FringeI = FringeI./(AuxTexture+.1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    figure('name','[mt_3DRetrieve]')
    subplot(1,2,1),imagesc(Ix), title('Ix'), colormap jet, colorbar, drawnow,
    subplot(1,2,2),imagesc(Iy), title('Iy'), colormap jet, colorbar, drawnow,
end

% -------------------------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute y-partial derivative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /////////////////////////////////////////////////////////////////////////
% //// A) split high and low frecuency components /////////////////////////
% /////////////////////////////////////////////////////////////////////////
% A0) Estimate fringes width --------------------
I = I1; % do the same we did but with the trasposicion
[m,n] = size(I);
FT = fftshift(fft2(I));
% Remove zero frec,
% Remove zero frec,
[aux,i0] = max(abs(FT));
[~,j0] = max(aux'); i0 = i0(j0);
FT(:,j0+[-waux:waux]) = 0; 

% find max (correspond to fringes)
[aux,i1] = max(abs(FT));
[~,j1] = max(aux'); i1 = i1(j1);

di = abs(j1(1)-j0(1)); 

fw = round( 1/di * n ); 

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
LPF(f0y+round([-m/fw:m/fw]/2),f0x+round([-n/fw:n/fw]/2)) = 1;

BPF = zeros(m,n); % BandPassFilter 
switch FilterType,
    case 1,
        BPF(:,f0x+round(n/fw)+round([-n/fw:n/fw]/2)) = 1;
        BPF(:,f0x-round(n/fw)+round([-n/fw:n/fw]/2)) = 1;
    case 2,
        BPF(:,1:f0x+round(-n/fw/2)) = 1;
        BPF(:,f0x+round(n/fw/2):end) = 1;
    case 3,
        BPF(f0y+round([-m/fw:m/fw]/2),f0x+round(n/fw)+round([-n/fw:n/fw]/2)) = 1;
        BPF(f0y+round([-m/fw:m/fw]/2),f0x-round(n/fw)+round([-n/fw:n/fw]/2)) = 1;
    case 4,
        BPF(f0y+round([-m/fw:m/fw]/2),f0x+round(n/fw)+round([-n/fw:n/fw]/2)) = 1;
    case 5,
        BPF(f0y+round([-m/fw+W:m/fw-W]/2),f0x+round(n/fw)+round([-n/fw+W:n/fw-W]/2)) = 1;
        BPF(f0y+round([-m/fw+W:m/fw-W]/2),f0x-round(n/fw)+round([-n/fw+W:n/fw-W]/2)) = 1;
end

if verbose>1, % show filters,  
    figure, imagesc(log(abs(FTI))), colormap jet; 
    AuxLPFImage = cat(3,zeros(m,n),zeros(m,n),ones(m,n));
    AuxBPFImage = cat(3,ones(m,n),zeros(m,n),zeros(m,n));
    hold on, imagesc(AuxLPFImage,'AlphaData',LPF*.5), 
    hold on, imagesc(AuxBPFImage,'AlphaData',BPF*.5),
end

LowFrecI2 = real(ift2(LPF.*FTI));
FringeI = real(ift2(BPF.*FTI));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove Texture %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AuxTexture = mt_Normalize( mean(Texture,3),[0 1]);
% FringeI = FringeI./(AuxTexture+.1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

Zy = Iy./(Ix+eps);

% -------------------------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remove outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /////////////////////////////////////////////////////////////////////////
% //// B) Integrate  //////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
[m,n] = size(Zx);
if size(Zy,1)~=m || size(Zy,2)~=n, error('[dimension mismach]'), end

Zx = medfilt2(Zx,ceil([w(1)*fw w(2)*fw]));
Zy = medfilt2(Zy,ceil([w(1)*fw w(2)*fw]));

if verbose>1,
     figure('name','Zx (with med filt)')
     imshow(Zx,[-1 1]), colormap jet, colorbar, drawnow,
     figure('name','Zy (with med filt)')
     imshow(Zy,[-1 1]), colormap jet, colorbar, drawnow,
end
varargout{2} = Zx;
varargout{3} = Zy;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Integrate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /////////////////////////////////////////////////////////////////////////
% //// C) Integrate  //////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
switch Type,
    case 'Gauss',
        addpath ../../matulsC++ % add path of cpp tools
        tic
        [Z,flag] = mt_IntegrateCPP2(Zx,Zy,MaxIter,Tol);
        t = toc; if verbose>0, mt_printtime(t), end
    otherwise,
        p = parameters; % pass all input parametes,   
        [Z] = mt_ImIntegrate2(Zx,Zy,p);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
