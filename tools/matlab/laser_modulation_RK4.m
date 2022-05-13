
function [PS]=laser_modulation_RK4(phase,Lw,lu,K,lL,w0,Power)

%% Laser parameter

ku=2*pi/lu;
kr=2*pi/(lL);
C=3e8;




%% mean 
PL=Power;
x=phase(:,1);
x_=phase(:,2);
y=phase(:,3);
y_=phase(:,4);
s=phase(:,5);
gamma=phase(:,6);
s=-s;
t=s/C;
t=t-min(t);
Ps0=[x_,y_,x,y,t,gamma];
Ps=Ps0;
zc=Lw/2;
 z=0;
theta =0 ;
nstep=200;
dz=Lw/nstep;
for I=1:nstep
    
k1=dz*tilt_laser_eq(Ps,z,ku,kr,K,PL,w0,zc,theta);
k2=dz*tilt_laser_eq(Ps,z,ku,kr,K,PL,w0,zc,theta);
k3=dz*tilt_laser_eq(Ps,z,ku,kr,K,PL,w0,zc,theta);
k4=dz*tilt_laser_eq(Ps,z,ku,kr,K,PL,w0,zc,theta);
z=z+dz;
Ps=Ps+(k1+2*k2+2*k3+k4)/6;
PS=Ps;

end
PS(:,5)=PS(:,5);
PS(:,5)=-PS(:,5); % head-tail to tail head (bunch head at larger time)
PS(:,5)=(PS(:,5)-min(PS(:,5)))*C; 



 

