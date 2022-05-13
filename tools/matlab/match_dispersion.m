function [x1,y1]=match_dispersion(x,y)
 % y delta for dispersion
 
x = x - mean(x);
y = y - mean(y);
x = x - mean(x.*y)/mean(y.^2).*y;
x1 = x;
y1 =y;
end
