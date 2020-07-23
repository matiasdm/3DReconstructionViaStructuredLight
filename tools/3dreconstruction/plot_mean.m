%% plot with mean and variance, 
% [[m],[v]] = plot_mean(Data,[xvalues]);
% Data should be>
% [y1   y2   y3   y4   ... 
%  y1'  y2'  y3'  y4'  ...
%  y1'' y2'' y3'' y4'' ...
% is Data has a third dimension each frame is plot in the same chart with
% different colors, (in this implementation size(Data,3) should be < 8)

function [varargout] = plot_mean(D,varargin);

if nargin>1,
    xvalues = varargin{1};
else 
    xvalues = 1:size(D,2);
end

colors  = {'r','g','b','y','c','k','m'};
markers = {'s','d','^','v','o','>','<'};
figure, 

for k = 1:size(D,3),
    m = mean(D(:,:,k));
    s = std(D(:,:,k));
    plot(xvalues,m,['-' colors{k} markers{k}],'MarkerSize',10,'LineWidth',2);
    hold on, 
end
% this is done in two steps so we can set the legend of the first elements,
for k = 1:size(D,3),
    m = mean(D(:,:,k));
    s = std(D(:,:,k));
    highlow((m+s/2)',(m-s/2)',(m-s/2)',(m+s/2)',colors{k},xvalues');
    highlow((m+s/2)',(m-s/2)',(m+s/2)',(m-s/2)',colors{k},xvalues');
end
set(gca,'xTick',xvalues);
set(gca,'xTickLabel',xvalues);

