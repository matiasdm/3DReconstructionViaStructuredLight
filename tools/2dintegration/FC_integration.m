% Frankot-Chellapa
% 
% Inputs: 
%   - gx (mxnx1): estimated x-partial derivative,
%   - gy (mxnx1): estimated y-partial derivative,
%
% Outputs:
%   - Z (mxnx1) retrieved surface
%
% Refs:
%  [1] Fankot and Chellappa. "A method for enforcing 
%      integrability in shape from shading algorithms"
%      TPAMI, 1988.
%  [2] Agrawal et. al. "What is the range of surface 
%      reconstructions from a gradient field? ". 
%      ECCV 2006.
% ----------------------------------------------------
% Matias Di Martino (c)                           2014
%                                 matiasdm@fing.edu.uy
% ----------------------------------------------------

function Z = FC_integration(gx,gy)

[m,n] = size(gx); 

% Define frecuency domain, 
[wx,wy] = meshgrid(0:n-1,0:m-1);
wx = 2*pi*wx/n;
wy = 2*pi*wy/m;

% Gradient fourier transform, 
Gx = fft2(gx);
Gy = fft2(gy);

% FT(Z) = (ax* Gx + ay* Gy) ./ ( |ax|^2 + |ay|^2 ); % [1]

% Asuming we approximate the derivatives by finite
% central differences [1]: 
% ax = (1/2) e^{j wx} - (1/2) e^{-j wx} = j sin(wx)
% ay = (1/2) e^{j wy} - (1/2) e^{-j wy} = j sin(wy)
j  = sqrt(-1);
ax = j * sin(wx); ay = j * sin(wy);

FTZ      = (conj(ax).*Gx + conj(ay).*Gy ) ...
           ./ ( abs(ax).^2 + abs(ay).^2 );

% Set zero component (which is undifined due to 
% 0/0 division) , 
FTZ(1,1) = 0;

Z = real(ifft2(FTZ));

end