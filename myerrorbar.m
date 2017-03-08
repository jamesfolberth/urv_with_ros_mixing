function [] = myerrorbar(x,y,l,u,width,color)
% MYERRORBAR  A dumb function to draw error bars (only) with a specified width and color

vlines_x = [x(:), x(:)].';
vlines_y = [l(:), u(:)].';

hlines_x = [x(:)-width/2, x(:)+width/2].';
hlines_y1 = [l(:), l(:)].';
hlines_y2 = [u(:), u(:)].';

plot(vlines_x, vlines_y, color, hlines_x, hlines_y1, color, hlines_x, hlines_y2, color);

end
