function displayD(D,names);



figure,
subplot(2,3,1); imagesc(D{1}); colormap jet; title(names{1}); colorbar;
subplot(2,3,2); imagesc(D{2}); colormap jet; title(names{2}); colorbar;
subplot(2,3,3); imagesc(D{3}); colormap jet; title(names{3}); colorbar;
subplot(2,3,4); imagesc(D{4}); colormap jet; title(names{4}); colorbar;
subplot(2,3,5); imagesc(D{5}); colormap jet; title(names{5}); colorbar;
subplot(2,3,6); imagesc(D{6}); colormap jet; title(names{6}); colorbar;
 
    