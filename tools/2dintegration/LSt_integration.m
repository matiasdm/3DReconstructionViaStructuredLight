%% [Z,[vx],[vy]] = LSt_integration(gx,gy,[parameters])
% 
% Least Squares Integration method using temporal
% regularization, 
% 
% Inputs: 
%   - gx (HxWxT): estimated x-partial derivative,
%   - gy (HxWxT): estimated y-partial derivative,
%   - parameters: (OPTIONAL) struct that may contain
%       .lambda: weight of tempral reg
%       .mu: smothness of vx and vy fields
%       .algOp:        
%           0 - process each frame indep.
%           1 - uses input scalars vx and vy                
%           2 - estimates scalar vx and vy 
%           3 - estimates fields vx and vy
%       .vx: (scalar) velocity along x of the objects
%       .vy: (scalar) velocity along x of the objects
%            vx and vy must be provided just for 
%            algOp = 1.
%
% Outputs:
%   - Z (HxWxT) retrieved surface
%   - vx (optional)
%   - vy (optional)
%
% Refs:
%   [1] M. Di Martino, A. Fernandez and J. Ferrari
%       "Gradient domanin methods with application to
%        4D scene reconstruction". Opticas and Lasers 
%        in Engineering, 2014
%   [2] Agrawal et. al. "What is the range of surface 
%       reconstructions from a gradient field? ". 
%       ECCV 2006.
% 
% ----------------------------------------------------
% Matias Di Martino (c)                           2014
%                                 matiasdm@fing.edu.uy
% ----------------------------------------------------

function [Z,varargout] = LSt_integration(gx,gy,varargin)

% ///////////////////////////////////////////////
% // Load and set general parameters ////////////
% ///////////////////////////////////////////////
if nargin>2, par = varargin{1}; else p = []; end, 
% list of input par.
list   = {'lambda','mu','algOp','vx','vy'}; 
% list of def. values.
defval = {.1      ,5   ,2      ,0   , 0}; 
[lambda,mu,algOp,vx,vy] = GetParameters(par,list,defval);
[H,W,T] = size(gx);
% set some derivation operators as global 
% as they are use for various algorithms.
global Dx Dy Dt L L_lambda; 
% -----------------------------------------------

% Now perform integration following
% the method indicated by algOp,
switch algOp,
    % ///////////////////////////////////////////////
    case 0, % Proces frame by frame indep.        ///
    % ///////////////////////////////////////////////    
    % this methods follows standard LS approach
    % frame by frame,     
        Z = IntegrateEachFrame(gx,gy);          
    
    % //////////////////////////////////////////////
    case 1, % impose temporal constraint with    ///
            % known velocity.                    ///
    % //////////////////////////////////////////////
    % In this case, Z is the dinamyc surface that 
    % minimizes 
    %
    % E1[u] = int{ (u_x-gx)^2 + (u_y-gy)^2 ...
    %            + lambda ( u_t + (vx*gx + vy*gy))^2   
    %            }dxdydt
    % whose euler lagrange eq. is:
    %    u_xx + u_yy + lambda u_tt = ...
    %    g_xx + g_yy - lambda (vx*gx+vy*gy)_t
    %
    % if we use vx = 0 and vy = 0 this approaches 
    % is a particular case of the techniques presented
    % in ref [1].             
         DefineLapAndDiOperators(H,W,T,lambda);
         Z = IntegrateGiven_vxvy(gx,gy,lambda,vx,vy);
    
    % //////////////////////////////////////////////
    case 2, % impose temporal constraint with    ///
    % unknown (but uniform) velocity.            ///
    % //////////////////////////////////////////////
    % In this case, Z, vx and vy are those who 
    % minimizes 
    %
    % E1[u,vx,vy] = int{ (u_x-gx)^2 + (u_y-gy)^2 ...
    %             + lambda ( u_t + (vx*gx + vy*gy))^2   
    %                  }dxdydt
    % whose euler lagrange eqs. are:
    % (ai) u_xx + u_yy + lambda u_tt ...
    %      + lambda (vx*gx+vy*gy)_t = g_xx + g_yy; 
    % (b)  gx'*u_t + gx'*gx * vx + gx'*gy * vy = 0
    % (c)  gy'*u_t + gy'*gx * vx + gy'*gy * vy = 0    
        DefineLapAndDiOperators(H,W,T,lambda);   
        [Z,vx,vy] = IntegrateAndEtimate_vxvy(gx,gy,lambda);
    
    % //////////////////////////////////////////////
    case 3, % impose temporal constraint with    ///
    % a unknown velocity field.                  ///
    % //////////////////////////////////////////////
    % In this case, Z, vx and vy are those who 
    % minimizes 
    %
    % E1[u,vx,vy] = int{ (u_x-gx)^2 + (u_y-gy)^2 ...
    %             + lambda ( u_t + (vx*gx + vy*gy))^2   
    %             + mu * [|grad(vx)|^2 + |grad(vy)|^2]    }dxdydt
    % whose euler lagrange eqs. are:
    % (ai) u_xx + u_yy + lambda u_tt ...
    %      + lambda (vx.*gx+vy.*gy)_t = g_xx + g_yy; 
    % (bi)  lambda * gx ( u_t + gx.*vx + gy.*vy ) - mu lap(vx) = 0
    % (ci)  lambda * gy ( u_t + gx.*vx + gy.*vy ) - mu lap(vy) = 0
        DefineLapAndDiOperators(H,W,T,lambda);
        [Z,vx,vy] = IntegrateAndEtimate_VXVY(gx,gy,lambda,mu);        
    otherwise
        error('algOp value is invalid.')
end

% set additional optional outputs.
varargout{1} = vx; varargout{2} = vy;

end

%% Auxiliary functions, %%%%%%%%%%%%%%
% /////////////////////////////////////////////////
function [varargout] = ...                % ///////
    GetParameters(parameters,list,defval) % ///////
% /////////////////////////////////////////////////

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
        warning(['input parameter ' f ' in unknown '])
    end
end

% now load input parameters or assign default values,
for par = 1:length(list);
    if isfield(parameters,list{par})
        varargout{par} = getfield(parameters,list{par}); 
    else % Set default value,
        varargout{par} = defval{par};
    end
end

if isfield(parameters,'verbose');
    if parameters.verbose>0;
        % show parameters values, 
        display('Input parameters: ');
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

% ////////////////////////////////////////////////
function DefineLapAndDiOperators(H,W,T,lambda)% //
% ////////////////////////////////////////////////
p = H*W*T;
global Dx Dy Dt L L_lambda; 

% Define Dx and Dy matrixs such that 
% {d(U)/di}(:) = Di*U(:), i = x,y,
N                              = (H+2)*(W+2)*(T+2);
mask                           = zeros(H+2,W+2,T+2);
mask(2:end-1,2:end-1,2:end-1)  = 1;
idx                            = find(mask==1);

Dx = 1/2 * ( sparse(idx,idx+(H+2),1,N,N) ...
           - sparse(idx,idx-(H+2),1,N,N) );
Dy = 1/2 * ( sparse(idx,idx+1    ,1,N,N) ...
           - sparse(idx,idx-1    ,1,N,N) );
Dt = 1/2 * ( sparse(idx,idx+(H+2)*(W+2),1,N,N)...
           - sparse(idx,idx-(H+2)*(W+2),1,N,N) );
                 
L  = sparse(idx,idx      ,-4,N,N)  ...
   + sparse(idx,idx+1    ,1,N,N) ...
   + sparse(idx,idx-1    ,1,N,N) ...
   + sparse(idx,idx+(H+2),1,N,N) ...
   + sparse(idx,idx-(H+2),1,N,N); 
 
L_lambda = L + lambda * ...
             ( sparse(idx,idx,-2,N,N) ... 
             + sparse(idx,idx+((H+2)*(W+2)),1,N,N) ...
             + sparse(idx,idx-((H+2)*(W+2)),1,N,N) );

Dx = Dx(idx,idx); Dy = Dy(idx,idx); Dt = Dt(idx,idx); 
L  = L(idx,idx); L_lambda = L_lambda(idx,idx); 
clear mask idx 

% sanity check
if p~=size(L,1), error('mismatch dimensions'), end 

% Correct some errors in the border definition, 
Dx       = Dx       - sparse(1:p,1:p,sum(Dx,2)      ,p,p);
Dy       = Dy       - sparse(1:p,1:p,sum(Dy,2)      ,p,p);
Dt       = Dt       - sparse(1:p,1:p,sum(Dt,2)      ,p,p);
L        = L        - sparse(1:p,1:p,sum(L ,2)      ,p,p);
L_lambda = L_lambda - sparse(1:p,1:p,sum(L_lambda,2),p,p);
end

% ////////////////////////////////////////////////
function Z0 = IntegrateEachFrame(gx,gy)  % ///////
% ////////////////////////////////////////////////
% Case 0: each frame independently
[H,W,T] = size(gx);

Z0 = zeros(H,W,T);
tic,
for t = 1:T;
    [Z0(:,:,t)] = LS_integration(gx(:,:,t),gy(:,:,t));
    Z0(:,:,t)   = Z0(:,:,t) - mean(mean(Z0(:,:,t)));
end
fprintf('Case 0, took: \n')
mt_printtime(toc)
end

% ////////////////////////////////////////////////
function Z1 = ...                             % //
    IntegrateGiven_vxvy(gx,gy,lambda,vx0,vy0) % //
% ////////////////////////////////////////////////
% // Case 1: vx and vy are known scalar //////////
% solve,
% Eq Z_xx + Z_yy + lambda Z_tt = ...
% gx_x + gy_y + lambda (vx*gx +vy*gy)_t
%  => A = L_lambda ;
%  => b = Dx*gx+Dy*gy+lambda*Dt*(vx*gx + vy*gy);
[H,W,T] = size(gx);
global Dx Dy Dt L L_lambda; 

A  = L_lambda;
b  = Dx*gx(:) + Dy*gy(:) - lambda * Dt*( vx0*gx(:) + vy0*gy(:) );
fprintf('Solving first case ... \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally, solve the linear (sparse) problem
tic
x1 = A(:,2:end)\b;
x1 = [0; x1];
fprintf('Case 1, A\b took: \n')
mt_printtime(toc)
Z1 = reshape(x1,[H W T]); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% /////////////////////////////////////////////////////
function [Z,vx,vy] = ...                      % ///////
    IntegrateAndEtimate_vxvy(gx,gy,lambda)    % ///////
% /////////////////////////////////////////////////////
% // Case 2: vx and vy are unknown scalar /////////////
% Euler lagrange eqs. are:
% (ai) u_xx + u_yy + lambda u_tt ...
%      + lambda (vx*gx+vy*gy)_t = g_xx + g_yy; 
% (b)  gx'*u_t + gx'*gx * vx + gx'*gy * vy = 0
% (c)  gy'*u_t + gy'*gx * vx + gy'*gy * vy = 0% solve Eqs:
%
%  => A = [L_lambda | lambda*Dt*gx | lambda*Dt*gy;
%          gxTDt    | gx'*gx       | gx'*gy      ;
%          gyTDt    | gy'*gx       | gy'*gy      ];
%  => b = [Dx*gx+Dy*gy;
%          0;
%          0]
[H,W,T]   = size(gx); p = H*W*T;
global Dx Dy Dt L L_lambda; 
aux_gxTDt = zeros(1,p);
aux_gyTDt = zeros(1,p);
for i = 1:p
    aux_gxTDt = aux_gxTDt + gx(i)*Dt(i,:);    
    aux_gyTDt = aux_gyTDt + gy(i)*Dt(i,:);
end

A  = [L_lambda, lambda*Dt*gx(:), lambda*Dt*gy(:);
      lambda*([aux_gxTDt, gx(:)'*gx(:), gx(:)'*gy(:)]);
      lambda*([aux_gyTDt, gy(:)'*gx(:), gy(:)'*gy(:)]) 
     ];
  
b  = [Dx*gx(:) + Dy*gy(:); 0; 0];

fprintf('Solving second case ... \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally, solve the linear (sparse) problem
tic
x = A(:,2:end)\b;
x = [0;x];
fprintf('Case 2, A\b took: \n')
mt_printtime(toc)
Z = reshape(x(1:p),[H W T]);
vx = x(p+1);
vy = x(p+2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%   
end
 
% /////////////////////////////////////////////////////
function [Z,vx,vy] = ...                      % ///////
    IntegrateAndEtimate_VXVY(gx,gy,lambda,mu) % ///////
% /////////////////////////////////////////////////////
% // Case 3: vx and vy are unknown scalar /////////////
% Euler lagrange eqs. are:   
% (ai) u_xx + u_yy + lambda u_tt ...
%      + lambda (vx.*gx+vy.*gy)_t = g_xx + g_yy; 
% (bi)  lambda * gx ( u_t + gx.*vx + gy.*vy ) - mu lap(vx) = 0
% (ci)  lambda * gy ( u_t + gx.*vx + gy.*vy ) - mu lap(vy) = 0
[H,W,T] = size(gx); p = H*W*T;
global Dx Dy Dt L L_lambda; 
% Define some additional operators,
Gx  = sparse(1:p,1:p,gx(:),p,p);
Gy  = sparse(1:p,1:p,gy(:),p,p);
Oo  = 0*speye(p,p); %null sparse matrix (pxp)

A  = [L_lambda, lambda*Dt*Gx, lambda*Dt*Gy;
      lambda*(gx(:)*ones(1,3*p)).*[Dt, Gx, Gy] - mu*[Oo, L, Oo];
      lambda*(gy(:)*ones(1,3*p)).*[Dt, Gx, Gy] - mu*[Oo, Oo, L]
     ];
  
b  = [Dx*gx(:) + Dy*gy(:); zeros(p,1); zeros(p,1)];

fprintf('Solving third case ... \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally, solve the linear (sparse) problemtic
x = A(:,2:end)\b;
x = [0; x];
fprintf('Case 3, A\b took: \n')
mt_printtime(toc)
Z  = reshape(x(1:p)      ,[H W T]);
vx = reshape(x(p+[1:p])  ,[H W T]);
vy = reshape(x(2*p+[1:p]),[H W T]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%     
end















%