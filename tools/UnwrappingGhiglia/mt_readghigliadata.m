%    Phase Unwrapping Data From
% Two-dimensional Phase Unwrapping: Theory, Algorithms, and Software
%               by Dennis C. Ghiglia and Mark D. Pritt
% 
%               Published by John Wiley and Sons, Inc.
%                Copyright 1998, all rights reserved
% 
% 
%    These are the datasets used for the experiments in the book.
%    For detailed information regarding these data, we refer the
%    user to chapter 3.
%    
%    The file naming convention is as follows:
%             name.NNNxNNN.ext
%    where "name" is the dataset name (e.g., head), "NNNxNNN" is
%    the file dimension (e.g., 256x256), and "ext" is the file
%    extension (e.g., phase).
%    
%    All files are in simple "raster" format, and the dimensions
%    are part of the file name.  For example, the mask values for
%    the Longs Peak example are contained in longs.152x458.mask.
%    There are 458 rows and 152 values in each row.  The file
%    extension "phase" indicates a file of phase values, "corr"
%    indicates a file of correlation values for weighting, "mask"
%    indicates a file that contains 0's where the phase should
%    be masked out or ignored, and "surf" indicates a file that
%    contains the elevation values for the actual (true) surface.
%    
%    For the actual data sets (ifsar, head and knee), only the
%    phase data is available.  For the "artificial" data sets
%    (shear, spiral and noise) only the phase data and the actual
%    surface data used to generate the phase data are available.
%    The simulated IFSAR data sets (longs and isola) include the
%    phase data, mask data, correlation data and actual surface
%    data.  
%    
%    All of the files consist of 1-byte values, except for the
%    surface (surf) files, which consist of 4-byte floating-point
%    values.  The 1-byte phase values denote the phase scaled and
%    quantized to the range 0-255.  The 1-byte correlation and
%    mask values denote the weights 0-1 scaled and quantized to
%    the range 0-255.
%    
%    The simulated IFSAR data sets (longs and isola) were created
%    by Dr. Michael Roth of the Johns Hopkins University Applied
%    Physics Laboratory in Laurel, Maryland.  The MRI data sets
%    (head and knee) are courtesy of Dr. Gang Zhu and come from
%    the ESTEEM clinical MRI scanner manufactured by Elscint MR
%    Inc., Ft. Collins, Colorado.  We thank them for allowing us
%    to post these files at this ftp site.
%    
%    The following table lists the files, along with a short
%    description and the book figure that shows the data.
%    
%    
%                                                          Book
%         Filename                 Description            Figure 
%         --------                 -----------            ------ 
%    
%    ifsar.512x512.phase  Phase data for IFSAR example     3.2a
%                               
%    head.256x256.phase   Phase data for MRI head example  3.3a 
%    
%    knee.256x256.phase   Phase data for MRI knee example  3.4a 
%    
%    longs.152x458.phase  Phase data for Longs Peak        3.7a
%    longs.152x458.mask   Mask data                         -
%    longs.152x458.corr   Correlation data                 3.11a
%    longs.152x458.surf   True surface data                3.6a 
%    
%    isola.157x458.phase  Phase data for Isolation Peak    3.8a
%    isola.157x458.mask   Mask data                         - 
%    isola.157x458.corr   Correlation data                 3.12a
%    isola.157x458.surf   True surface data                3.6b
%    
%    shear.257x257.phase  Phase data for shear example     3.9b
%    shear.257x257.surf   True surface data                3.9a
%    
%    spiral.257x257.phase Phase data for spiral shear      3.10b
%    spiral.257x257.surf  True surface data                3.10a
%    
%    noise.257x257.phase  Phase data for noisy shear       4.15a
   


clear all
close all

fid = fopen('knee.256x256.phase');

dx = 256;
dy = 256;

size = dx*dy;

M = fread(fid,size,'uint8');

M = reshape(M,[dy dx]);

figure, imshow(M,[])
fclose('all');


%%
fid = fopen('salida');

M = fread(fid,size,'float');

M = reshape(M,[dy dx]);

figure, imshow(M,[])
fclose('all');

