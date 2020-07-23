%% [Disparity] = BinaryCodedStripes(Iin,Plane,parameters);
% function that compute the disparity of the binary coded stripes. 
%
% Inputs:
%   - Iin: 1xN cell where Iin{n} code the n bit of the stripes projected
%          over the surface.
%   - Plane: 1xN cell where P{n} code the n bit of the stripes projected
%            over the plan.
%   - parameters: struct that may contain,
%       .verbose {def 0} 1-display some text, 2- also some graphics
%       .???
%
% Outputs:
%   - Disparity: is propo to surf. height, matrix of the same size of
%                Iin{i}
%
% -------------------------------------------------------------------------
% matias di martino, matiasdm@fing.edu.uy                              2014
% -------------------------------------------------------------------------

function [D] = BinaryCodedStripes(Ic,Pc,parameters);

% load input parameters, 
if isfield(parameters,'verbose');
    verbose = parameters.verbose;
else % set default value, 
    verbose = 0;
end



%% Recover coded stripes
N      = length(Ic);
[m n]  = size(Ic{1});
P      = zeros(m,n);
I      = zeros(m,n);

% Binarize input images;
if verbose>1, figure, end;

for k = 1:N;
    Ic{k} = mt_Normalize(Ic{k},[0 1]);
    th    = graythresh(Ic{k});
    Ic{k} = Ic{k}>th;
    Pc{k} = mt_Normalize(Pc{k},[0 1]);
    th    = graythresh(Pc{k});
    Pc{k} = Pc{k}>th; 
    if verbose>1, 
        clf, subplot(1,2,1), imshow(Ic{k},[]); title('Binarized Ic') 
             subplot(1,2,2), imshow(Pc{k},[]); title('Binarized P')
        drawnow, pause(.5);
    end
end

for k = 1:N;
    P = P + 2^(N-k) * Pc{k};
    I = I + 2^(N-k) * Ic{k};
end 

if verbose>1,
    clf, subplot(1,2,1), imshow(uint8(I)); title('decoded Ic')
    subplot(1,2,2), imshow(uint8(P)); title('decoded P')
    drawnow, pause(.5);
end

%% Find the correspondences 
max_d = 70; % half of the width of the "looking window"
% add borders;
I = [zeros(m,max_d) I zeros(m,max_d)];
P = [zeros(m,max_d) P zeros(m,max_d)];
D = zeros(m,n); % inicialization;
for i = 1:m;
    for j = 1:n;
        j2     = j+max_d;
        aux    = P(i,j2+[-max_d:max_d]);
        [~,d]  = min(abs(I(i,j2)-aux));
        D(i,j) = d-max_d;       
    end
end

if verbose>1, 
    figure('name','D'); imagesc(D); colormap jet;
end

    





