
function [U]=tilt_laser_eq(Ps,z,ku,kr,K,PL,w0,zc,theta)

P0=8.7e9;
C=3e8;
z0=kr*w0^2/2;
ux=Ps(:,1); %betax
uy=Ps(:,2);
x=Ps(:,3);
y=Ps(:,4);
t=Ps(:,5);
gamma=Ps(:,6);

% undulator field (leading terms)
By=K*ku*(1+(ku*y*1).^2).*cos(ku*z);
Bz=-K*ku*ku.*y.*sin(ku*z);
% laser filed


%% tilt y-z theta (y,z) to (y0,z0)

%%
z1=z-zc;

%fxy=1./(1+(z1/z0).^2).^(1/2).*(exp(-(x.^2+y.^2)./(w0^2*(1+(z1/z0).^2))));
fxy=1;
Eamp=2*sqrt(4*PL/P0)/w0;
fi=kr*(z*cos(theta)+y*sin(theta)-C*(t));
Ex=Eamp*sin(fi).*fxy;
%% 

% diff.eq.
uz=sqrt(1-1./gamma.^2-ux.^2-uy.^2);
uxp=(Ex.*(1-uz*cos(theta)))./gamma./uz+By./gamma+uy.*Bz./gamma./uz;  % d(batax*gamma)/dz = -e(E+VxB)/mc;
uyp=-(Ex.*ux*sin(theta))./(gamma)./uz-ux.*Bz./gamma./uz;
xp=ux; % dx/dz= (batax*gamma)/(bataz*gamma)
yp=uy;
tp=1./(uz*C)*1; % t;
gp=-Ex.*ux;% gamma;
U=[uxp,uyp,xp,yp,tp,gp];
end



