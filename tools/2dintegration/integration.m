% function Z = integration(gx,gy,method,[Qx],[Qy]);
% 
% methos: 'LS', 'WLS', 'AD', 'FC', 'L1', 'M'.

function Z = integration(gx,gy,method,varargin)

switch method,
    case 'LS',
        Z = LS_integration(gx,gy);
    case 'WLS',
        Z = WLS_integration(gx,gy,varargin{1},varargin{2});
    case 'FC',
        Z = FC_integration(gx,gy);
    case 'AD',
        Z = AD_integration(gx,gy);
    case 'L1',
        % as it is to conputationally expensive we must reduce gradient
        % field resolution, 
        %[m,n] = size(gx); 
        %gx = medfilt2(gx,[10 10]);
        %gy = medfilt2(gy,[10 10]);
        %gx    = imresize(gx,[128 128]);
        %gy    = imresize(gy,[128 128]);
        Z     = L1_integrationCVx(gx,gy);
        %Z     = imresize(Z,[m n]);
    case 'M',
        Z = M_integration(gx,gy);
    
end
