function [p] = verpendiente;


[x,y] = ginput(2);

p = (y(2)-y(1)) / (x(2)-x(1));