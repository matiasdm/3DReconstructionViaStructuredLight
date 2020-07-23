% funcion [Inom]=mt_Normalize(I,[range],[verbose]),
% 
% % Input:
% I: imagen original (si la imagen es color lo hace por componente)
% [rango]: establece entre que valores queremos la normalizacion
%        min=rango(1), max=rango(2).
%        default=[0 1];
% verbose:[default 0]
% Output:
% Inom: double image with range in [rando(1) rangoo(2)]
%
% -------------------------------------------------------------------------
% matiasdm@fing.edu.uy, 6/6/2013
% -------------------------------------------------------------------------


function [Inorm]=mt_Normalize(I,varargin)

if nargin<2,
    rango=[0 1];
else
    rango=varargin{1};
end

if nargin<3,
    verbose=0;
else
    verbose=varargin{2};
end

I = double(I);


for k=1:size(I,3),
    aux=I(:,:,k);
    % sanity check,
    if min(aux(:))==max(aux(:)),
        % do nothing,
        disp('[mt_normalize]WARNING min==max ')
        Inorm(:,:,k) = I(:,:,k);
    else
        aux_nom=(rango(2)-rango(1)) * ( (aux-min(I(:)))/(max(I(:))-min(I(:))) ) ...
            + rango(1);
        Inorm(:,:,k)=aux_nom;
    end
end

if verbose,
    disp(['Imagen normalizada entre ' num2str(rango(1)) ' y ' num2str(rango(2))])
end
