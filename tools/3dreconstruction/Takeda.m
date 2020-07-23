%% [Disparity] = Takeda(Iin,Plane,parameters);
% function that compute the disparity of the binary coded stripes. 
%
% Inputs:
%   - Iin: 1x3 cell where Iin{n} code the n image of the PSI sequence
%   - Plane: 1x3 cell where P{n} 
%   - parameters: struct that may contain,
%       .verbose {def 0} 1-display some text, 2- also some graphics
%       .
%
% Outputs:
%   - Disparity: is propo to surf. height, matrix of the same size of
%                Iin{i}
%
% -------------------------------------------------------------------------
% matias di martino, matiasdm@fing.edu.uy                              2014
% -------------------------------------------------------------------------

function [D] = Takeda(I,P,parameters);

% load input parameters, 
if isfield(parameters,'verbose');
    verbose = parameters.verbose;
else % set default value, 
    verbose = 0;
end

%%
parameters.verbose = 2;
[DeltaPhi] = mt_Takeda3DProfilometry(I,P,parameters);

if verbose>1,
    figure, imagesc(DeltaPhi), colormap jet;
end

% Unwrapping tools, 
path = 'tools/UnwrappingGhiglia/code/';
[m n]  = size(DeltaPhi);
% write DeltaPhi as byte image 
fid = fopen('tmp/DeltaPhi','w');
DeltaPhiAux = mt_Normalize(DeltaPhi,[0 255]);
%
DeltaPhiAux = imresize(DeltaPhiAux,[1025 1025]); % when using unwt
%
xsize = size(DeltaPhiAux,2);
ysize = size(DeltaPhiAux,1);
DeltaPhiAux = uint8(DeltaPhiAux);
fwrite(fid,DeltaPhiAux','uint8');
fclose(fid);

fopen('tmp/Phi','w'); fclose('all');

algorithm = 'gold'; % gold, diff, flynn, fmg, lpno, mcut, pcg, qual, unwt
auxstr = [path algorithm ...
         ' -input tmp/DeltaPhi -format byte -output tmp/Phi -xsize ' ...
         num2str(xsize) ' -ysize ' num2str(ysize)];

[status,result]=system(auxstr);

% Read output from disk
fid = fopen('tmp/Phi');
phi = fread(fid,[xsize ysize],'float',0,'n')';
fclose(fid);

fopen('tmp/Phi','w'); fopen('tmp/DeltaPhi','w'); fclose('all'); % eras files.
% back to the original size;
D = imresize(-phi,[m n]);


if verbose>1, 
    figure('name','D'); imagesc(D); colormap jet;
end

    
