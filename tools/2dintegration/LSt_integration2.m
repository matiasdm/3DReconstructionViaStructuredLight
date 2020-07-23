%% [Z,[vx],[vy]] = LSt_integration2(gx,gy,[parameters])
% performs LSt_integration for times [-n ... 0 ... n] to integrate the 
% frame 0
% 
% ----------------------------------------------------
% Matias Di Martino (c)                           2014
%                                 matiasdm@fing.edu.uy
% ----------------------------------------------------

function Z = LSt_integration2(gx,gy,par)

[H,W,T] = size(gx);
Z       = zeros(H,W,T);
n       = 2;

gx = padarray(gx,[0 0 n],0,'both');
gy = padarray(gy,[0 0 n],0,'both');

for t = 1:T;
    auxZ     = LSt_integration(gx(:,:,t+n+[-n:n]),gy(:,:,t+n+[-n:n]),par);
    Z(:,:,t) = auxZ(:,:,n+1);
end
    
end
%