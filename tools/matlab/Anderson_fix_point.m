
clear all
input_parameters.wrf = 2*pi*166e6;
input_parameters.alphac = 2.491432e-05;
input_parameters.E0 = 3.5e9;
input_parameters.T0 = 972/3e8;
input_parameters.U0 = 9.010539e-01*1e6;
input_parameters.energy_spread = 1.2e-3;
input_parameters.n2 = 2;

U0 = 9.010539e-01*1e6;
alphac = 2.491432e-05;
h = 540;
E0 = 3.5e9;


%%  Anderson acceleration 

alpha = -0.6; % self-adaption? from large abs(alpha) to small abs(alpha)

x1 = (0.0391)^2/(U0/(pi*h*alphac*E0));  % 4.196e-2

err_bar = 1e-9;
v0  = -1.05;
err = 1;
jj = 1;
while err > err_bar
 jj = jj +1;
 f0= asin(1./v0);
 fb = -pi-2*f0;
   g0 = x1 - (fb+v0.*(cos(fb+f0)-cos(f0)))+ v0;
v1 = alpha*g0 + (1-alpha)*v0;
 
 err = abs(v1-v0);
v0 = v1;
end

vout = v1

jj






[x,iter,res_hist] = AndAcc(@g_find_v,v0)














