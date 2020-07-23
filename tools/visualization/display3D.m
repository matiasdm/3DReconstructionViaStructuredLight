function display3D(Z,p)

if isfield(p,'type');
    type = p.type;
else 
    type = 'gray';
end

if isfield(p,'X');
    X = p.X;
    Y = p.Y;
else 
    [X,Y] = meshgrid([1:size(Z,2)],[1:size(Z,1)]);
end 

if isfield(p,'Normalization')
    if p.Normalization == 1,
        X = -mt_Normalize(X,[0 1]);
        Y = mt_Normalize(Y,[0 1]);
        Z = mt_Normalize(Z,[0 1]);
    end
end

% Open figura and set plot propieties -----------------
set(gcf,'color',[1 1 1]);
clf
switch lower(type), 
    case 'color', 
        C = double(mt_Normalize(p.C,[0 1]));
        h = surface(X,Y,Z,C,'EdgeColor','none');
        shading interp, colormap gray    
        %light
    case 'gray'
        h = surface(X,Y,Z,'EdgeColor','none','FaceColor',[.5 .5 .5]);
        light;
    case 'mesh'
        set(h3,'color',[0 0 0]);
        h = surface(X(1:5:end,1:5:end),Y(1:5:end,1:5:end),Z(1:5:end,1:5:end),...
            'EdgeColor',[0 1 0],'FaceColor','none');
    case 'jet'        
        h = surface(X,Y,Z,'EdgeColor','none','FaceColor','interp',...
            'CData',Z); colormap jet, 
end
camorbit(180,0,'data',[0 0 1])
camorbit(40,0,'data',[0 1 0])

drawnow, axis vis3d, axis off, 
xlabel('x'); ylabel('y');
%daspect([1 1 2.5]); 

%lighting phong

% Set light propieties
%material metal 
%material metal
material dull
%set(h,'SpecularStrength',.1); % in [0 1]

%set(h,'SpecularExponent',10); % in [1 500]
%set(h,'SpecularColorReflectance', 0.5); % in [0 1]

%camlight
% -----------------------------------------------------
