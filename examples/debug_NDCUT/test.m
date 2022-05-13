clear all
clc

data0 = textread('wake_CU.beam','','headerlines',4,'commentstyle','shell');
data = textread('lcls.dist','','headerlines',5,'commentstyle','shell');

%%
ndcut = 1;

t00  = data0(:,1)/3e8;
cur0 = data0(:,2);
t = data(:,5);
charge = 9.8e-10;
tdmin = min(t);
tdmax = max(t);

npart = 16384;
nbins = 4;

%mpart =npart/nbins;
dtd = (tdmax-tdmin)/ndcut;
dcharge = charge/(length(t));
zsep = 32;
nslice = 14002;
ntail = 30;
xlamds = 1.5e-10;
tn =[];
for ii =1:nslice

t0 = tdmax - (ntail+ii-1)*zsep*xlamds/3e8-0.5*dtd;
t1 = t0 +dtd;
indx = find(t>=t0 & t<= t1);
mget = length(indx);
xcur(ii)=dcharge*mget/dtd;
tn = [tn,t0];
end

%%
figure(1)
% plot(t00,cur0)
% hold on
plot(tn-min(tn),xcur)



