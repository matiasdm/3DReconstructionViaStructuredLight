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

function Z = L1_integration(gx,gy)

[m,n]       = size(gx); 
mn          = m*n;

% we want to minimize Sum fi * xi;
% the vector xi contains in the first mn
% entries Z, then u and v. To impose
% min {Sum (ui + vi)} f must be defined as:
f           = zeros(3*mn,1);
f(mn+1:end) = 1;

% Define Dx and Dy operators, 
[H,W] = size(gx);
N     = (H+2)*(W+2);
mask  = zeros(H+2,W+2);

mask(2:end-1,2:end-1) = 1;
idx                   = find(mask==1);

Dx = 1/2 * ( sparse(idx,idx+(H+2),1,N,N) ...
           - sparse(idx,idx-(H+2),1,N,N) );
Dy = 1/2 * ( sparse(idx,idx+1    ,1,N,N) ...
           - sparse(idx,idx-1    ,1,N,N) );

Dx = Dx(idx,idx); Dy = Dy(idx,idx); 

% Alternative circular definition (this impose symetric boundary contd.), 
%se = speye(mn);
%Dx = 1/2 * ( circshift(se,[0 m]) - circshift(se,[0 -m]) );
%Dy = 1/2 * ( circshift(se,[0 1]) - circshift(se,[0 -1]) );

% Now we set the set of inequalities, as
% A*x < b.
A = speye(6*mn,3*mn);
b = zeros(6*mn,1);

for i = 1:mn,
    % Z_x - p < u 
    % ==> 1/2 (Zij+1 - Zij-1) - uij < pij
    i2 = i;
    A(i2,1:mn)    = Dx(i,:);
    A(i2,i+mn)    = -1;
    b(i2)         =  gx(i);
   
    % Z_x - p > -u 
    % ==> -1/2 (Zij+1 - Zij-1) - uij < -pij
    i2 = i+mn; % Stack after the previous mn eqs.
    A(i2,1:mn)    =  -Dx(i,:);
    A(i2,i+mn)    =  -1;
    b(i2)         =  -gx(i);
     
    % Z_y - q < v 
    % ==> 1/2 (Zi+1j - Zi-1j) - vij < qij  
    i2 = i+2*mn; % Stack after the previous 2mn eqs. 
    A(i2,1:mn)    = Dy(i,:);
    A(i2,i+2*mn)  = -1;
    b(i2)         =  gy(i);
    
    % Z_y - q > -v 
    % ==> -1/2 (Zi+1j - Zi-1j) - vij < -qij  
    i2 = i+3*mn; % Stack after the previous 3mn eqs. 
    A(i2,1:mn)    =  -Dy(i,:);
    A(i2,i+2*mn)  =  -1;
    b(i2)         =  -gy(i);
        
    % u > 0 ==> - uij < 0  
    i2 = i+4*mn; % Stack after the previous 4mn eqs. 
    A(i2,i+mn)    = -1;  
    b(i2)         = 0;
    
    % v > 0 ==> - vij < 0  
    i2 = i+5*mn; % Stack after the previous 5mn eqs. 
    A(i2,i+2*mn)  = -1;  
    b(i2)         = 0;
end

x = linprog(f,A,b); % minimize sum(fi*xi) st. Ax<b

% keep the first mn elements, and reshape it.
Z     = reshape(x(1:mn),[m n]);

end