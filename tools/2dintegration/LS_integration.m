% Least Squares Integration method
% 
% Inputs: 
%   - gx (mxnx1): estimated x-partial derivative,
%   - gy (mxnx1): estimated y-partial derivative,
%
% Outputs:
%   - Z (mxnx1) retrieved surface
%
% ----------------------------------------------------
% Matias Di Martino (c)                           2014
%                                 matiasdm@fing.edu.uy
% ----------------------------------------------------

function Z = LS_integration(gx,gy)

[m,n] = size(gx); mn = m*n;

% Define Dx and Dy matrixs such that 
% {d(U)/di}(:) = Di*U(:), i = x,y,
[H,W] = size(gx); p = H*W;

% Define Dx and Dy matrixs such that 
% {d(U)/di}(:) = Di*U(:), i = x,y,
N                              = (H+2)*(W+2);
mask                           = zeros(H+2,W+2);
mask(2:end-1,2:end-1)          = 1;
idx                            = find(mask==1);

Dx = 1/2 * (sparse(idx,idx+(H+2),1,N,N) - sparse(idx,idx-(H+2),1,N,N));
Dy = 1/2 * (sparse(idx,idx+1    ,1,N,N) - sparse(idx,idx-1    ,1,N,N));

L  = sparse(idx,idx,-4,N,N)  ...
   + sparse(idx,idx+1,1,N,N) ...
   + sparse(idx,idx-1,1,N,N) ...
   + sparse(idx,idx+(H+2),1,N,N) ...
   + sparse(idx,idx-(H+2),1,N,N); 

Dx = Dx(idx,idx); Dy = Dy(idx,idx); L  = L(idx,idx); 

% Correct some errors in the border differenctiation definition, 
Dx       = Dx       - sparse(1:mn,1:mn,sum(Dx,2)      ,p,p);
Dy       = Dy       - sparse(1:mn,1:mn,sum(Dy,2)      ,p,p);
L        = L        - sparse(1:mn,1:mn,sum(L ,2)      ,p,p);

% ////////////////////////////
% // Poisson Eq //////////////
% ////////////////////////////
% Eqs Zxx + Zyy = div(gx,gy)
A = L;
b = Dx*gx(:) + Dy*gy(:);

% ///////////////////////////
% // Border Conditions, /////
% ///////////////////////////
% Abord  = []; b_bord = [];
% 
% % Eq Zx-gx = 0 on x = 1, x = n, 
% Abord  = [Abord; Dx([1:m],:)]; 
% b_bord = [b_bord; gx(transpose(1:m))];
% Abord  = [Abord; Dx(m*(n-1) + [1:m],:)]; 
% b_bord = [b_bord; gx(m*(n-1) + transpose(1:m) )];
% 
% % Eqs Zy-gy = 0 on y = 1, y = m,
% Abord  = [Abord; Dy( ([1:n] - 1) * m+1,:)];   
% b_bord = [b_bord; gy((transpose(1:n)-1)*m+1)];
% Abord  = [Abord; Dy( [1:n]*m,:)];      
% b_bord = [b_bord; gy( transpose(1:n)*m )];

% Add border conditions to the linear sistem of equations, 
%A = [A; Abord]; b = [b; b_bord];


% ///////////////////////////
% // Additional constraints /
% ///////////////////////////
% the integration process is an Ill-Posed problem, 
% we have a set of surfaces that differ by a constant 
% and are solution of the system. This causes that the
% rank(A) = mn-1, we can solve this isue by adding an 
% additional constrain, e.g. looking for the solution 
% with zero mean or arbitrary setting the value of some 
% pixel. We found that besides this look reasonable in 
% theory, in practice the numberical method used to 
% solve the linear system is robust and do not need this 
% additional constraint.
%
%A = [A; ones(1,mn) ]; b = [b; 0]; 


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally, solve the linear (sparse) problem
x = A(:,2:end)\b; 
x = [0; x];
Z = reshape(x,[m n]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%

end