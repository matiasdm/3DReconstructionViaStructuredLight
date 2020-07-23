function [] = mydisplay(Z,view_dir)

[H,W] = size(Z);
global maxZ;

if(isempty(maxZ))
    disp('maxZ is not defined');
    maxZ = max(Z(:));
end


ss = [-60 21]';  %ramp peaks


if(exist('view_dir','var'))
    ss = view_dir;
end


if(exist('ss','var'))
    figure;surfl(Z,ss);axis ij;shading interp;colormap gray;axis off;view(ss(1),ss(2));
else
    figure;surfl(Z);axis ij;shading interp;colormap gray;axis off;
end
axis off


