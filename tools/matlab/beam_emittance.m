%% beam_emittance

function [emit,beta,alpha]=beam_emittance(x,y)

% fit_coefs=polyfit(x,delta,1);
% x= x- fit_coefs(1)*delta;
% 
% fit_coefs=polyfit(y,delta,1);
% y= y- fit_coefs(1)*delta;


x=x-mean(x);
y=y-mean(y);
emit=sqrt(mean(x.^2)*mean(y.^2)-(mean(x.*y))^2);
beta = mean(x.^2)/emit;
alpha = -mean(x.*y)/emit;

end

