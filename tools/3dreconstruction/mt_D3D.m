%% [D,[Texture],[Dx],[Dy],[Qx],[Qy]] = mt_D3D(I,parameters) 
% function that performs differential 3D retrival (I = image adquiered with
% projected fringe patterns)
%
% Inputs:
%   - I (mxnxc) gray or color image with horizontal and vertical fringes.
%               Fringes can be multiplied, summed of coded in the blue and
%               red channels of the projector. Indicate the kind of
%               codification used in parameters.CodeType. By the default we
%               assume that CodeType = sum.
%       
%   - paremeters, struct that may contain:
%       .CodeType [def sum] {'sum','prod','color'}; in the first two, the
%                           fringes are overlaped in all the channels, in
%                           the case of color we assume that the VERTICAL
%                           fringes are coded in the blue channel and the
%                           HORIZONTAL ones in the red channel.
%       .verbose [def 0] (1 some text is display, 2 also some graphics.)
%       .IntegrationMethod 
%           - 'LS'  : Least Squares Integration
%           - 'WLS' : Weighted Least Squares Integration (-> Qx,Qy)
%           - 'FC'  : Frankot-Chellapa approach (Fourier)
%           - 'AD'  : Anisotropic diffusion integration approach,
%           - 'L1'  : L1-error minimization,
%           - 'M'   : M-estimator integration approach,
%       .FiltWidth [def [2 2]] (width of the windows in median filter (in
%                              terms of the fringes width), this filter is 
%                              applied to the gradient field before integr-
%                              ation)
%       .TextureIm (mxnxc): texture image (used for normalization, if it is
%                           not provided, we estimate it from the low
%                           frecuencies. 
%       .Theta: [def pi/4] angle in which we translate the projector.
%
% Outputs:
%   - D (mxnx1) 
%   - Texture (mxnxc); [optional] 
%   - Dx (mxnx1) [optional] 
%   - Dy (mxnx1) [optional] 
%   - Qx (mxnx1) [optional] Dx quality map 
%   - Qy (mxnx1) [optional] Dy quality map 
%
% Refs:
%   [1] Matias Di Martino, Gaston Ayubi, Alicia Fernandez and Jose Ferrari.
%       "Differential 3d shape retrieval", Opticas and Lasers in
%       Engineering, 58C:114-118, 2014.
%
% see also mt_ImIntegrate2_v2
%
% -------------------------------------------------------------------------
% matias di martino (c) - matiasdm@fing.edu.uy - 2014                v.0.5 
% -------------------------------------------------------------------------
% notes:
% v.0.1 - la que use en inge de muestra y el articulo, esta basicamente en
% la funcion mt_3DRetrieveBoth, mantengo esa funcion para que los programas
% que tengo viejos no den erro.
% v.0.2 - ahora incluyo la posibilidad de que las franjas esten codificadas
% en color y de tener una imagen de referencia para la textura. 
% v.0.3 - Escribiendo la tesis y la parte de relacion geometrica encontre
% la expresion exacta y ahora la uso. Se corrigen un par de sqrt(2) que
% aprecen tambien por correr el proj 45 respecto a la cam. Tambien cambio
% imIntegrate, ahora uso la v.0.2. Se agrega tambien un modo para detectar
% las zonas de discontinuidades.  
% v.0,4. Con esta version hice las pruebas sinteticas.
% v.0.5. Estoy haciendo las pruebas con giros.mat con esta version, Agrego
% la funcionalidad de calcular una maskara y anular los gradientes en
% dichos puntos, Asi cuando viene data que tiene un background ruidoso
function [D,varargout] = mt_D3D(Iin,parameters)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load and set general parameters, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list = {'CodeType' ,'verbose'  ,'IntegrationMethod',...
        'FiltWidth','TextureIm','Theta'            ,'DarkBackground'};
defval = {'sum'    ,0          ,'fourier'          ,...
          1        ,[]         ,pi/4               , 0};
[CodeType,verbose,IntegrationMethod,FiltWidth,TextureIm,Theta,DB] = ...
    GetParameters(parameters,list,defval);

Iin = double(Iin);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Mask (usefull for dark background %%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DB, % if DarkBackground
     Mask = CalculateMask(Iin,verbose,'manual'); % 'manual' or 'automatic'
else Mask = ones(size(Iin,1),size(Iin,2));
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Fringe Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fw] = EstimateFringeWidth(Iin,verbose); % estimate the width of the fringes

% if we don't have a texture image estimate it from the low frecuencies.
if isempty(TextureIm);
    [TextureIm] = EstimateTextureImage(Iin,fw,verbose);
    if strcmpi(CodeType,'color'); % then we just have two channels of
        % the texture image, so we will keep just a gray TextureIm,
        TextureIm = mean(TextureIm,3);
    end        
end
varargout{1} = TextureIm;

% Normalization of the input image (to compensate the effects of surface
% reflectance and obtain Ih and Iv (images with horizontal and vertical 
% fringes)
switch lower(CodeType),
    case 'color',
        Ih   = Iin(:,:,1)./TextureIm(:,:,1);
        if size(TextureIm,3)>1, % this if is becacause we can have an 
                                % acquired texture image (in which case we
                                % know the three channels of the
                                % reflectance. Or we can estimate the
                                % texture image in which case we have a
                                % unique "average" estimation.
            Iv   = Iin(:,:,3)./TextureIm(:,:,3);
        else 
            Iv   = Iin(:,:,3)./TextureIm(:,:);
        end
        k   = fspecial('gaussian',round(.5*[fw fw]),fw/3);
        Ih  = conv2(Ih,k,'same'); Iv = conv2(Iv,k,'same'); 
    otherwise % if the image has the fringes overlapped: filter in the 
              % fourier domain,
        I = Iin./double(TextureIm); I = mean(I,3);
        [Ih,Iv] = ExtractVerticalAndHorizontalFringes(I,fw,verbose);
end

if verbose>1, % show Ih and Iv;
    figure('name','[mt_D3D] Ih and Iv'); 
    subplot(121), imshow(Ih,[]), title('Ih');
    subplot(122), imshow(Iv,[]), title('Iv');
    drawnow,
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Depth Gradient Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Dx, Dy, Qx, Qy] = ComputeGradientField(Ih,Iv,Theta,FiltWidth,fw,verbose); 
Dx = Dx.*Mask; Dy = Dy.*Mask; % if it is necessary remove background
varargout{4} = Qx; varargout{5} = Qy;
varargout{2} = Dx; varargout{3} = Dy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Integrate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = integration(Dx,Dy,IntegrationMethod,Qx,Qy);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end




%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /////////////////////////////////////////////////////////////////////////
function [varargout] = GetParameters(parameters,list,defval) % ////////////
% /////////////////////////////////////////////////////////////////////////

% first check that all the parameters are in the list.
names = fieldnames(parameters);
for par = 1:length(names);
    f        = names{par};
    isinlist = 0;
    for j = 1:length(list);
        if strcmpi(f,list{j});
            isinlist = 1;
        end
    end
    if isinlist == 0;
        warning(['[mt_D3D] input parameter ' f ' in unknown '])
    end
end

% now load input parameters or assign default values,
for par = 1:length(list);
    if isfield(parameters,list{par})
        varargout{par} = getfield(parameters,list{par}); %#ok<GFLD>
    else % Set default value,
        varargout{par} = defval{par};
    end
end

if isfield(parameters,'verbose');
    if parameters.verbose>0;
        % show parameters values, 
        display('[mt_D3D] Input parameters: ');
        for par = 1:length(list);
            if length(varargout{par})>20;
                auxstr = 'ToLargeToBeDisplayed';
            else 
                auxstr = num2str(varargout{par});
            end   
            fprintf(['\t' list{par} ' > ' auxstr '\n']);
        end
    end
end
end

% /////////////////////////////////////////////////////////////////////////
function Mask = CalculateMask(I,verbose,so); % ////////////////////////////
% /////////////////////////////////////////////////////////////////////////
I     = mean(I,3);
I     = mt_Normalize(I,[0 1]);
[m,n] = size(I);
switch lower(so),
    case 'manual',
        figure('name','Iin'); imshow(I);
        fprintf('\n Click on the center of the nose \n ');
        [x0, y0] = ginput(1);
        
        fprintf('\n Click on the right border of the face \n ');
        [x1, y1] = ginput(1);
        
        R     = abs(x1-x0); % radius
        [X,Y] = meshgrid(1:n,1:m);
        Mask  = sqrt( (X-x0).^2 + .4*(Y-y0).^2 ) < R;
    case 'automatic',
        th = graythresh(I); % use otsu's threshold
        Mask = I > th;
        if verbose>1,
            figure('name','Ig>th');
            imshow(Mask);
        end
        
        % Erode and dilatate to remove holes and compute border points
        se = strel('disk',ceil(20));
        Mask = imdilate(Mask,se);
        se = strel('disk',ceil(20));
        Mask = imerode(Mask,se);
end

if verbose>1,
    figure('name','Mask');
    imshow(Mask,[]);
end

end

% /////////////////////////////////////////////////////////////////////////
function [fw] = EstimateFringeWidth(I,verbose)  % /////////////////////////
% /////////////////////////////////////////////////////////////////////////
[m,n] = size(I);
% Define some shortcuts for image transformations -----
ft2  = @(u) fftshift(fft2(u));
ift2 = @(u) real(ifft2(ifftshift(u)));
% -----------------------------------------------------

% A0) Estimate fringes width --------------------
FT = ft2(mean(I,3)); % fourier transform of the gray image
% Remove zero frec,
waux = 7;
[aux,i0] = max(abs(FT));
[~,j0] = max(aux'); i0 = i0(j0);
FT(i0+[-waux:waux],:) = 0; 

% find max (correspond to fringes)
[aux,i1] = max(abs(FT));
[~,j1] = max(aux'); i1 = i1(j1);
di = abs(i1(1)-i0(1)); 
fw = round( 1/di * m ); 

if verbose>0, disp(['[mt_D3D] fringes width :  ' num2str(fw) '(pixels)']), end
% -----------------------------------------------
end

% /////////////////////////////////////////////////////////////////////////
function [T] = EstimateTextureImage(I,fw,verbose) % ///////////////////////
% /////////////////////////////////////////////////////////////////////////
[m,n,col] = size(I);

% Define some shortcuts for image transformations -----
ft2  = @(u) fftshift(fft2(u));
ift2 = @(u) real(ifft2(ifftshift(u)));
% -----------------------------------------------------

% A1) Filter on Fourier domain --------------------------------------------
f0x = round(n/2)+1; f0y = round(m/2)+1; % center of spectrum
% Filter in fourier domain ----------------------------------
% so we have fringe on one side and low frecuencies (asociated with object
% reflectance) on the other.

% Gaussian low pass filter -----------------------
LPF          = zeros(m,n); % LowPassFilter 
LPF(f0y,f0x) = 1;
KSIZE        = [round(m/fw) round(n/fw)];
KSIGMA       = .5 * (round(m/fw)+round(n/fw)) / 2 ;
K            = fspecial('gaussian',KSIZE,KSIGMA);
LPF          = conv2(LPF,K,'same');
% ------------------------------------------------

T = I; % inicialization;
for c = 1:col;
    T(:,:,c) = ift2(LPF.*ft2(I(:,:,c)));
end

T = uint8(mt_Normalize(T,[1 255]));

if verbose>1, % show spectrum, filter and Texture image obtained.
    figure('name','[mt_D3D] ft(Iin) + LPF'), 
    imagesc(log(abs(mean(ft2(I),3)))), colormap jet; % show spectrum
    AuxLPFImage = cat(3,zeros(m,n),zeros(m,n),ones(m,n));
    hold on, imagesc(AuxLPFImage,'AlphaData',LPF/max(LPF(:))*.5), hold off,
    figure('name','[mt_D3D] Estimated TextureIm'), 
    imshow(T);
end

end

% /////////////////////////////////////////////////////////////////////////
function [Ih,Iv] = ExtractVerticalAndHorizontalFringes(I,fw,verbose) % ////
% /////////////////////////////////////////////////////////////////////////
I = mean(I,3); % convert to gray
[m,n] = size(I);

% Define some shortcuts for image transformations -----
ft2  = @(u) fftshift(fft2(u));
ift2 = @(u) real(ifft2(ifftshift(u)));
% -----------------------------------------------------

% A1) Filter on Fourier domain --------------------------------------------
f0x = round(n/2)+1; f0y = round(m/2)+1; % center of spectrum
% Filter in fourier domain ----------------------------------

% Gaussian Band pass filter -----------------------
BPFv         = zeros(m,n); % Band Pass Filter (in the horizontal direction)
BPFh         = zeros(m,n); % Band Pass Filter (in the vertical direction) 
BPFv(f0y,f0x+round(n/fw)) = 1; BPFv(f0y,f0x-round(n/fw)) = 1;
BPFh(f0y+round(m/fw),f0x) = 1; BPFh(f0y-round(m/fw),f0x) = 1;
%KSIZE        = [round(m/fw) round(n/fw)];
%KSIGMA       = .5 * ( round(m/fw)+round(n/fw) ) / 2;
%K            = fspecial('gaussian',KSIZE,KSIGMA);
K            = hanning(round(m/fw))*hanning(round(n/fw))';
BPFv         = conv2(BPFv,K,'same');
BPFh         = conv2(BPFh,K,'same');
% -------------------------------------------------

Ih = ift2(ft2(I).*BPFh); Iv = ift2(ft2(I).*BPFv);
if verbose>1, % show spectrum and filters,
    figure('name','[mt_D3D] ft(Iin) + BPFs'), 
    imagesc(log(abs(mean(ft2(I),3)))), colormap jet; % show spectrum
    AuxBPFhImage = cat(3,ones(m,n),zeros(m,n),zeros(m,n));
    AuxBPFvImage = cat(3,zeros(m,n),ones(m,n),zeros(m,n));
    hold on, imagesc(AuxBPFhImage,'AlphaData',BPFh/max(BPFh(:))*.5), 
             imagesc(AuxBPFvImage,'AlphaData',BPFv/max(BPFv(:))*.5), 
    hold off,
end
end

% /////////////////////////////////////////////////////////////////////////
function [Dx,Dy,varargout] = ComputeGradientField...
                                      (Ih,Iv,Theta,FiltWidth,fw,verbose)
% /////////////////////////////////////////////////////////////////////////
%[Dx,Dy,[Qx],[Qy]] = ComputeGradientField(Ih,Iv,Theta,FiltWidth,fw)  % 
% Qx and Qy are Quality maps,

% define some shortcuts
dx = @(U) 1/2*[0*U(:,1)  U(:,3:end)-U(:,1:end-2)  0*U(:,1)]; 
dy = @(U) 1/2*[0*U(1,:); U(3:end,:)-U(1:end-2,:); 0*U(1,:)]; 

%ta1 = dx(Ih)./dy(Ih);  ta2 = dy(Iv)./dx(Iv);
%Dx = ta1 .* ( 1 + ta2*tan(Theta) ) ./  ( sin(Theta) * (ta1.*ta2 - 1) );
%Dy = ta2 .* ( 1 + ta1/tan(Theta) ) ./  ( cos(Theta) * (ta2.*ta1 - 1) );

% instead of the previous, we compute an equivalent expression that has
% less divisions and so less NaN introduces by 0/0;
Dx = ( dx(Ih) .* ( dx(Iv) + dy(Iv)*tan(Theta) ) ) ./ ...
     ( sin(Theta) * ( dx(Ih).*dy(Iv) - dy(Ih).*dx(Iv) ) );

Dy = ( dy(Iv) .* ( dy(Ih) + dx(Ih)/tan(Theta) ) ) ./ ...
     ( cos(Theta) * ( dx(Ih).*dy(Iv) - dy(Ih).*dx(Iv) ) ); 


% Remove NaN due to 0/0 division;
Dx(isnan(Dx(:))|isinf(Dx(:))) = 0; Dy(isnan(Dy(:))|isinf(Dy(:))) = 0;

if FiltWidth>0,
    %Dx = medfilt2(Dx,max(ceil(FiltWidth*fw),3)*[1 1]);
    %Dy = medfilt2(Dy,max(ceil(FiltWidth*fw),3)*[1 1]);
    Dx = medfilt2(Dx,[7 7]); % solo para hacer las pruebas variando fw
    Dy = medfilt2(Dy,[7 7]); % despues sacar esto...
end

% Find the discontinuty map 
[Qx,Qy] = ComputeQualityMap(Dx,Dy,fw,verbose);
varargout{1} = Qx;
varargout{2} = Qy;

if verbose>1,
     figure('name','[mt_D3D] atan(Dx) (with med filt)')
     imshow(atan(Dx),[]), colormap jet, % colorbar, drawnow, %#ok<DUALC>
     figure('name','[mt_D3D] atan(Dy) (with med filt)')
     imshow(atan(Dy),[]), colormap jet, %colorbar, drawnow, %#ok<DUALC>
     figure('name','[mt_D3D] Qx')
     imshow(Qx,[]), colormap jet, % colorbar, drawnow, %#ok<DUALC>
     figure('name','[mt_D3D] Qy')
     imshow(Qy,[]), colormap jet, %colorbar, drawnow, %#ok<DUALC>
     drawnow,
end

end

% /////////////////////////////////////////////////////////////////////////
function [Qx,Qy] = ComputeQualityMap(Dx,Dy,fw,verbose)  % /////////////////
% /////////////////////////////////////////////////////////////////////////
[m,n] = size(Dx);

Qx = ones(m,n);
Qy = ones(m,n);

Qx(isnan(Dx) | isinf(Dx)) = 0;
Qy(isnan(Dy) | isinf(Dy)) = 0;

Dx2 = Dx;
Dy2 = Dy;

Dx2(isnan(Dx) | isinf(Dx)) = 0;
Dy2(isnan(Dy) | isinf(Dy)) = 0;


auxX = abs(Dx2 - medfilt2(Dx2, round([1 fw*2])));
auxY = abs(Dy2 - medfilt2(Dy2, round([fw*2 1])));

th2  = 2e-1;
w1   = 3; 
w2   = 30;% 30 is a reasonable value
w3   = w2-5;

% Mask for Qx,
AuxMask = auxX>th2;
se      = strel('square',w1);
AuxMask = imerode(AuxMask,se);  %imagesc(AuxMask)
se      = strel('square',w2); 
AuxMask = imdilate(AuxMask,se); %imagesc(AuxMask)
se      = strel('square',w3);
AuxMask = imerode(AuxMask,se);  %imagesc(AuxMask)
Qx(AuxMask) = 0;

% Mask for Qy
AuxMask = auxY>th2;
se      = strel('square',w1);
AuxMask = imerode(AuxMask,se);
se      = strel('square',w2);
AuxMask = imdilate(AuxMask,se);
se      = strel('square',w3);
AuxMask = imerode(AuxMask,se);
Qy(AuxMask) = 0;

if verbose>1,
    figure('name','[mt_D3D] Qmap'), 
    subplot(1,2,1), imagesc(Qx);
    subplot(1,2,2), imagesc(Qy);
end

end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
