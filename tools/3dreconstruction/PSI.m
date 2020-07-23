%% [Disparity] = PSI(Iin,Plane,parameters);
% function that compute the disparity of the binary coded stripes. 
%
% Inputs:
%   - Iin: 1x3 cell where Iin{n} code the n image of the PSI sequence
%   - Plane: 1x3 cell where P{n} 
%   - parameters: struct that may contain,
%       .verbose {def 0} 1-display some text, 2- also some graphics
%       .Normalize {def 0}, if 1, move fringes to zero mean and unitary
%                           variance. 
%
% Outputs:
%   - Disparity: is propo to surf. height, matrix of the same size of
%                Iin{i}
%
% -------------------------------------------------------------------------
% matias di martino, matiasdm@fing.edu.uy                              2014
% -------------------------------------------------------------------------

function [D] = PSI(I,P,parameters);

% load input parameters, 
if isfield(parameters,'verbose');
    verbose = parameters.verbose;
else % set default value, 
   verbose = 0;
end

if isfield(parameters,'Normalize');
    Normalize = parameters.Normalize;
else % set default value, 
    Normalize = 0;
end

%% Normalize inputs, 
% /////////////////////////////////////////////////////////////////////////
if Normalize, % normalize input data //////////////////////////////////////
    % /////////////////////////////////////////////////////////////////////////
    for k = 1:3;
        mu = mean(I{k}(:));
        sg = std(I{k}(:));  
        I{k} = (I{k}-mu)/sg;
 
        mu = mean(P{k}(:));
        sg = std(P{k}(:));        
        P{k} = (P{k}-mu)/sg;        
    end
    if verbose>1,
        aux = round(size(P{1},1));
        figure('name','PSI Ps slice'), plot(P{1}(aux,:),'r'), hold on, 
        plot(P{2}(aux,:),'g'), plot(P{3}(aux,:),'b'), 
        figure('name','PSI Is slice'), plot(I{1}(aux,:),'r'), hold on, 
        plot(I{2}(aux,:),'g'), plot(I{3}(aux,:),'b'), 
    
    end
end
    
%%
% ///////////////////////////////////////////////////////////////////
% // Mask ///////////////////////////////////////////////////////////
% ///////////////////////////////////////////////////////////////////

DeltaPhi   = atan( sqrt(3)*(I{1}-I{3}) ./ ( 2*I{2}-I{1}-I{3} ) );
[m n]   = size(DeltaPhi);

if verbose>1,
    figure('name','[PSI] DeltaPhi'), imagesc(DeltaPhi), colormap jet;
end

% Unwrapping, 
path = 'tools/UnwrappingGhiglia/code/';

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

[status,result]=system(auxstr)

% Read output from disk
fid = fopen('tmp/Phi');
phi = fread(fid,[xsize ysize],'float',0,'n')';
fclose(fid);

fopen('tmp/Phi','w'); fopen('tmp/DeltaPhi','w'); fclose('all'); % eras files.
D = phi;

%%
% ///////////////////////////////////////////////////////////////////
% // Plane ///////////////////////////////////////////////////////////
% ///////////////////////////////////////////////////////////////////

DeltaPhi   = atan( sqrt(3)*(P{1}-P{3}) ./ ( 2*P{2}-P{1}-P{3} ) );
if verbose>1,
    figure('name','[PSI] DeltaPhi_p'), imagesc(DeltaPhi), colormap jet;
end
% Unwrapping, 
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
auxstr = [path algorithm ...
         ' -input tmp/DeltaPhi -format byte -output tmp/Phi -xsize ' ...
         num2str(xsize) ' -ysize ' num2str(ysize)];

[status,result]=system(auxstr)

% Read output from disk
fid = fopen('tmp/Phi');
phi = fread(fid,[xsize ysize],'float',0,'n')';
fclose(fid);
fopen('tmp/Phi','w'); fopen('tmp/DeltaPhi','w'); fclose('all'); % remove files.

D = D-phi;

% back to the original size;
D = imresize(D,[m n]);

if verbose>1, 
    figure('name','D'); imagesc(D); colormap jet;
end

    





