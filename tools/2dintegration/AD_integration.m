% Anisotropic diffusion integration approach,
% 
% Inputs: 
%   - gx (mxnx1): estimated x-partial derivative,
%   - gy (mxnx1): estimated y-partial derivative,
%
% Outputs:
%   - Z (mxnx1) retrieved surface
%
% Refs:
%   [1] Gilles Aubert, Pierre Kornprobst, 
%       "Mathematical Problems In Image Proccesing". 
%       Springer, sec. edition
%   [2] A. Agrawal. "Scene analysis under variable 
%       illumination using gradient domain methods". 
%       PhD Thesis, University of Mariland 2006.
%   [3] Joachim Weickert. "Anisotropic Diffusion in 
%       Image Processing". Published by B. G. 
%       Teubner Stuttgart. 2008. 
%   [4] Agrawal et. al. "What is the range of surface 
%       reconstructions from a gradient field? ". 
%       ECCV 2006.
% ----------------------------------------------------
% Matias Di Martino (c)                           2014
%                                 matiasdm@fing.edu.uy
% ----------------------------------------------------

function Z = AD_integration(gx,gy)

% Set some parametes, 
verbose     = 0; 
KernelSigma = 1/2; % Sigma of the gaussian kernel
                   % used to compute the smoothed 
                   % gradient field from which 
                   % the tensor is estimated.

[m,n] = size(gx);
mn    = m*n; 

% Calculate diffusion tensor
[d11, d12, d22] = ...
    CalculateDiffusionTensor(gx,gy,verbose,KernelSigma);

A = laplacian_matrix_tensor(m,n,d11,d12,d12,d22);
b = calculate_f_tensor(gx,gy,d11,d12,d12,d22);

b = b(:); % convert to column vector

x = A(:,2:end)\b; % remove one variable to make 
                  % the solution unique (we have
                  % a family of functions difering 
                  % on a constant
x = [0;x]; % arbitrary set as 0 the fist pixel.

Z = reshape(x,[m n]); % reshape the vector solution 
                      % to the original shape.

end

% ///////////////////////////////////////////
function [D11, D12, D22] = ...
    CalculateDiffusionTensor(gx,gy,verbose,sigma);
% ///////////////////////////////////////////
if verbose>0, 
    tic, fprintf('\t Computing diffusion tensors '),
end
    
[m,n] = size(gx);

% Set tensor parameters -------------------------
KernelSigma = sigma;     %0.5 is a typical value;
KernelWidth = max(3,round(2*KernelSigma)); 
Beta        = 0.02; % [2] 
% -----------------------------------------------

GaussianKernel = fspecial('gaussian',KernelWidth,KernelSigma);

% Memory prelocation 
H11 = zeros(m,n); H12 = zeros(m,n); H22 = zeros(m,n);

H11(:,:)    = conv2(gx(:,:).^2,GaussianKernel,'same');
H12(:,:)    = conv2(gx(:,:).*gy(:,:),GaussianKernel,'same');
H22(:,:)    = conv2(gy(:,:).^2,GaussianKernel,'same'); 

% Memory prelocation 
D11 = zeros(m,n); D12 = zeros(m,n); D22 = zeros(m,n);

% compute largest eigen value,
mu1 = 1/2*( H11+H22 + sqrt( (H11-H22).^2 + 4 * H12.^2 ) );
th  = max(mu1(:))/100; % this small threshold is 
                       % computed as the tenth of the mean 
                       % of H11+H22 which is in the order 
                       % of magnitud of the eigen values 
                       % of matrix H.
                              
for i = 1:m,
    for j = 1:n,
        AuxH = [H11(i,j) H12(i,j); H12(i,j) H22(i,j)];
        [V,D] = eigs(AuxH);
        % return de eigenvalus in the diagonal of D
        % in descendent order of magnitud and they
        % corresponding eigenvectors as columns in V.
                                      
        %SANITY CHECK ~~~~~~~~~~~~~~~~~~~~~~
        if D(2,2)>D(1,1),
            error('there is an error in eigenvalues')
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               
        % modify the diagonal acording to ([2],pag 55)
        D(2,2) = 1; % the smallest eigen value (mu2) 
                    % is modify by 1 the largest eigen 
                    % value (mu1) is set to 1 if it is small
        if D(1,1) < th,
            D(1,1) = 1;
        else % if it is ">0" is modified according [1] pag 115:
            D(1,1) = Beta + 1 - exp( -3.315 / (D(1,1))^4 );
        end
                                      
        auxD = V*D*V';
        D11(i,j) = auxD(1,1);
        D12(i,j) = auxD(1,2);
        D22(i,j) = auxD(2,2);
    end
end

if verbose>0, 
    t = toc; fprintf('   Ok \n')
    fprintf('\t\tCalculating diffusion tensor took: \n')
    fprintf('\t\t'), mt_printtime(t)
end

% release memory
clear auxD auxH D V
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
A = -A;

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
