% M-estimator integration approach,
% 
% this code is based on the source code provided 
% by Amit Agrawal, 2006
% http://www.umiacs.umd.edu/~aagrawal/
% Permitted for personal use and research purpose only
% Refer to the following citations:
%
%   1.  A. Agrawal, R. Raskar and R. Chellappa, "What 
%       is the Range of Surface Reconstructions 
%       from a Gradient Field? European Conference on
%       Computer Vision (ECCV) 2006
%   2.  A. Agrawal, R. Chellappa and R. Raskar, 
%       "An Algebraic approach to surface reconstructions 
%        from gradient fields? Intenational Conference 
%        on Computer Vision (ICCV) 2006
%
%=========================================================
% Inputs: 
%   - gx (mxnx1): estimated x-partial derivative,
%   - gy (mxnx1): estimated y-partial derivative,
%
% Outputs:
%   - Z (mxnx1) retrieved surface
%
% ----------------------------------------------------
% Matias Di Martino                               2014
%                                 matiasdm@fing.edu.uy
% ----------------------------------------------------

function Z1 = M_integration(gx,gy)

[m,n] = size(gx); mn = m*n;

% Use LS solution as seed,
Z0   = WLS_integration(gx,gy,ones(m,n),ones(m,n));
Z0   = Z0 - mean(Z0(:));

% /////////////////////////////////////////////
% chose weight function ///////////////////////
% /////////////////////////////////////////////
option = 'Huber'; %{'Huber', 'Fair', 'Cauchy'}
switch option,
    case 'Huber',
        k   = 1.345;
        rho = @(e) (abs(e)<=k) .* 1 + ...
                   (abs(e)>k)  .* (k ./ (abs(e)+eps) );
    case 'Fair',
        c   = 1.3998;
        rho = @(e) 1 ./ ( 1 + abs(e)/c );
    case 'Cauchy',
        c   = 1.4;
        rho = @(e) 1 ./ (1 + (e/c).^2 ); 
end

% Inicialization, 
verbose  = 0;
tol      = 1e-6;
diff     = 2*tol+1;
iter     = 0;
max_iter = 50;
dx       = @(U) [U(:,2:end)-U(:,1:end-1)  0*U(:,1)];
dy       = @(U) [U(2:end,:)-U(1:end-1,:); 0*U(1,:)];
wx       = rho(gx-dx(Z0));
wy       = rho(gy-dy(Z0));

% ////////////////////////////////////////
% // Begin iteration /////////////////////
% ////////////////////////////////////////
if verbose, 
    h  = figure('Position',[1464 623 456 351]);
    h2 = figure('Position',[1464 161 456 351]); 
    h3 = figure('Position',[1008 623 456 351]);
end

while diff>tol && iter<max_iter;
    Z1       = WLS_integration(gx,gy,wx,wy);
    Z1       = Z1 - mean(Z1(:));
    wx       = rho(gx-dx(Z1));
    wy       = rho(gy-dy(Z1));
    diff     = mean(abs(Z1(:)-Z0(:)));
    Z0       = Z1;
    iter     = iter+1;
    if verbose,
        fprintf('Iter %i | diff %4.2e \n', iter, diff)
        set(h,'name',num2str(iter)), figure(h),
        imagesc(reshape(Z0,[m n])), drawnow, pause(.1)
        set(h2,'name',num2str(iter)), figure(h2),
        imagesc(wx), drawnow, pause(.1)
        set(h3,'name',num2str(iter)), figure(h3),
        imagesc(wy), drawnow, pause(.1)

    end
end

end

% Aux functions, 

% //////////////////////////////////////////////////////
% // Weighted Least Squares Approach ///////////////////
function Z = WLS_integration(gx,gy,wx,wy) % ////////////
% //////////////////////////////////////////////////////

[m,n] = size(gx); mn = m*n;

A = laplacian_matrix_tensor(m,n,wx,zeros(m,n),zeros(m,n),wy);
A = -A;
b = calculate_f_tensor(gx,gy,wx,zeros(m,n),zeros(m,n),wy);
b = b(:);
  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally, solve the linear (sparse) problem
x = A(:,2:end)\b; 
x = [0;x];
Z = reshape(x,[m n]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%

end



function [A] = laplacian_matrix_tensor(H,W,D11,D12,D21,D22)


if(exist('D11','var') & exist('D12','var') & exist('D21','var') & exist('D22','var'))
    %disp('Weighted Poisson Solver');
else
    D11 = ones(H,W);
    D12 = ones(H,W);
    D21 = ones(H,W);
    D22 = ones(H,W);
    disp('All weights are one in diffusion tensor');
end


D21 = D12;


D11 = padarray(D11,[1 1],0,'both');
D12 = padarray(D12,[1 1],0,'both');
D21 = padarray(D21,[1 1],0,'both');
D22 = padarray(D22,[1 1],0,'both');


N = (H+2)*(W+2);
mask = zeros(H+2,W+2);
mask(2:end-1,2:end-1) = 1;
idx = find(mask==1);

A = sparse(idx,idx+1,-D22(idx),N,N);
A = A + sparse(idx,idx+H+2,-D11(idx),N,N);
A = A + sparse(idx,idx-1,-D22(idx-1),N,N);
A = A + sparse(idx,idx-H-2,-D11(idx-H-2),N,N);

A = A + sparse(idx,idx+1,-D12(idx),N,N);
A = A + sparse(idx,idx-H-2,-D12(idx-H-2),N,N);
A = A + sparse(idx,idx-H-2+1,D12(idx-H-2),N,N);
A = A + sparse(idx,idx+H+2,-D21(idx),N,N);
A = A + sparse(idx,idx-1,-D21(idx-1),N,N);
A = A + sparse(idx,idx-1+H+2,D21(idx-1),N,N);


A = A(idx,idx);
N = size(A,1);
dd = sum(A,2);
idx = [1:N]';
A = A + sparse(idx,idx,-dd,N,N);

end


% kernel should be symmteric
% Calculate weighted divergence (uD in paper)

function [f] = calculate_f_tensor(gx,gy,d11,d12,d21,d22)


[H,W] = size(gx);


gx(:,end) = 0;
gy(end,:) = 0;

if(~(exist('d11','var') & exist('d12','var') & exist('d21','var') & exist('d22','var')))
    disp('Weights are all zeros')
    d11 = ones(H,W);
    d21 = d11;
    d12 = d11;
    d22 = d11;
end

d21 = d12;


gx1 = gx.*d11;
gy1 = gy.*d22;
gx1 = padarray(gx1,[1 1],0,'both');
gy1 = padarray(gy1,[1 1],0,'both');
gxx = zeros(size(gx1)); gyy = gxx;
j = 1:H+1;
k = 1:W+1;
% Laplacian
gyy(j+1,k) = gy1(j+1,k) - gy1(j,k);
gxx(j,k+1) = gx1(j,k+1) - gx1(j,k);
f = gxx + gyy;
f = f(2:end-1,2:end-1);
clear gx1 gy1 gxx gyy


gx1 = gx.*d12;
gy1 = gy.*d21;

gx2 = gy.*d12;
gy2 = gx.*d21;

gx2(end,:) = gx1(end,:);
gy2(:,end) = gy1(:,end);


gx2(:,end) = 0;
gy2(end,:) = 0;


gx2 = padarray(gx2,[1 1],0,'both');
gy2 = padarray(gy2,[1 1],0,'both');
gxx = zeros(size(gx2)); gyy = gxx;
j = 1:H+1;
k = 1:W+1;
% Laplacian
gyy(j+1,k) = gy2(j+1,k) - gy2(j,k);
gxx(j,k+1) = gx2(j,k+1) - gx2(j,k);
f2 = gxx + gyy;
f2 = f2(2:end-1,2:end-1);

f = f + f2;
end

