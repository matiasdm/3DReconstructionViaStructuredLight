% L1-error minimization,
% 
% Inputs: 
%   - gx (mxnx1): estimated x-partial derivative,
%   - gy (mxnx1): estimated y-partial derivative,
%
% Outputs:
%   - Z (mxnx1) retrieved surface
%
% Refs:
%  [1] Zhouyu Du et al. "Robust surface reconstruction 
%      from gradients field using the L1 norm". 
%      Digital Image Computing Techiques and
%      Applications, IEEE, 2007.
% ----------------------------------------------------
% Matias Di Martino (c)                           2014
%                                 matiasdm@fing.edu.uy
% ----------------------------------------------------

function Z = L1_integrationCVx(gx,gy)


% Define Dx and Dy operators, 
[H,W] = size(gx);
N     = (H+2)*(W+2);
mask  = zeros(H+2,W+2);

mask(2:end-1,2:end-1) = 1;
idx                   = find(mask==1);

%Dx = 1/2 * ( sparse(idx,idx+(H+2),1,N,N) ...
%          - sparse(idx,idx-(H+2),1,N,N) );
%Dy = 1/2 * ( sparse(idx,idx+1    ,1,N,N) ...
%          - sparse(idx,idx-1    ,1,N,N) );
 Dx = ( sparse(idx,idx+(H+2),1,N,N) ...
            - sparse(idx,idx,1,N,N) );
 Dy = ( sparse(idx,idx+1    ,1,N,N) ...
            - sparse(idx,idx,1,N,N) );

Dx = Dx(idx,idx); Dy = Dy(idx,idx); 


A = [Dx; Dy];
b = [gx(:); gy(:)];
n = H*W;
cvx_begin
    variable x(n)
    minimize( norm(A*x-b,1) )
cvx_end

Z = reshape(x,[H W]);


end