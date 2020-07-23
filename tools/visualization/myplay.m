function myplay(Z,p);
%% myplay(Z,parameters)
% parameters is a struct that may contain
%   .FileName
%   .FigName
%   .Position
%   .VideoName: give a name if you want the video to be saved
%   .PrintFrames: [1 2 3 etc] give the number of the frames to be printed

if isfield(p,'FileName');
    FileName = p.FileName;
else 
    FileName = 'default';
end

h = figure; 

if isfield(p,'FigName')    
    set(h,'name',p.FigName);
end

if isfield(p,'Position')
    set(h,'Position',p.Position);
end

num_frames = size(Z,3);

if isfield(p,'VideoName');
    VideoName = p.VideoName;
    MakeVideo = 1;
else 
    MakeVideo = 0;
end

if isfield(p,'PrintFrames');
    Frames = p.PrintFrames;
    PrintFrames = 1;
    if isfield(p,'FramesName'),
        frname = p.FramesName;
    else 
        frname = [];
    end
else 
    PrintFrames = 0;
end

if MakeVideo, 
    frame=0;
    clear F
end

for k = 1:num_frames,
    switch FileName,         
        % //////////////////////////////////////////////    
        case 'visor1';
            imshow(Z(:,:,k),[],'InitialMagnification',500); 
            colormap gray;
            drawnow; pause(.15);

        % //////////////////////////////////////////////    
        case 'RampPeaks';
            zmin = min(Z(:)); zmax = max(Z(:));        
            surfl(Z(:,:,k)), zlim([zmin zmax]);
            axis vis3d, shading interp;colormap gray;axis off;
            camzoom(1.5); camorbit(-40,0);
            drawnow, pause(.3)           
            
        % //////////////////////////////////////////////    
        case 'Gaussian',
            imagesc(Z(:,:,k)), 
            %surf(Z(:,:,k),'FaceColor',[.5 .5 .5],'EdgeColor','none'), 
            %shading interp, axis equal, axis vis3d, axis off, 
            %material dull, light, colormap summer            
            drawnow; pause(.3);   
            
        % //////////////////////////////////////////////    
        otherwise
            imshow(Z(:,:,k),[-.1 .1]), axis equal, axis off, 
            colormap gray; drawnow; pause(.1);
    end
    if MakeVideo,
        frame = frame+1;
        F(frame)=getframe(h);
    end
    if PrintFrames, 
        if sum((Frames(:)-k)==0)>=1,
            print('-depsc',[frname '_frame' num2str(k) '.eps'])
        end
    end
end


% Genero video de salida,
if MakeVideo,
    writerObj = VideoWriter([VideoName '.avi']);
    writerObj.FrameRate=5;
    writerObj.Quality=50;
    open(writerObj);
    writeVideo(writerObj,F) 
    close(writerObj);
end
