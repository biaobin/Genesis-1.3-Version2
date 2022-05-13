%% postprocess
clear all

filename='./rad.out';% output/0.4m_100as_1e6/
%filename='../../benchmark/radpslice260nm.out';
fig=1;
lambdar=0.8;
%[out]=read_genesis_output(filename);
[data0,datas,current,xlamd,delts,nslice,datainfo]=readflie(filename,fig);

% max_power along position
z=data0(1:1,1);
for I=1:length(z)
% power
power(I)=mean(datas(I,1,:));

% bunch 

bunch(I)=mean(datas(I,8,:));

% beam energy 


%beam szie
beamx(I)=mean(datas(I,9,:));

beamy(I)=mean(datas(I,10,:));
 
end


% output radiation pulse
s=delts*(1:nslice);
powers=datas(length(z),1,:);
powers=reshape(powers, [1 ,length(s)]);
M=find(powers>=0.5*max(powers));
FWHM0=(s(M(end))-s(M(1)))/3e8*1e15
sm = s(find(powers==max(powers)));
mean_s=sum(s.*powers)/sum(powers);
std_s=sqrt(sum((s-mean_s).^2.*powers)/sum(powers))
%
% spectrum 
powers_mid = datas(length(z),3,:);
powers_mid=reshape(powers_mid, [1 ,length(s)]);

bunchs = datas(1,8,:);
bunchs=reshape(bunchs, [1 ,length(s)]);

energys=datas(length(z),7,:);
energys=reshape(energys, [1 ,length(s)]);

phase=datas(length(z),4,:);
phase=reshape(phase, [1 ,length(s)]);
[l,plunit]=spectrim(s,powers_mid,phase,nslice,lambdar);
mean_l=sum(l.*plunit)/sum(plunit);
bw=sqrt(sum((l-mean_l).^2.*plunit)/sum(plunit))/lambdar


M=find(plunit>=0.5);
FWHM=(l(M(end))-l(M(1)))/lambdar




%
figure(1)
hold on
plot(data0(:,1),power*1,'linewidth',1)
xlabel('z (m)')
ylabel('power (W)')
set(gca,'FontSize',16,'yscale','log');
box on

figure(2)
hold on
plot((s)/3e8*1e15,(powers))
xlabel('t (fs)')
ylabel('power (W)')
%set(gca,'yscale','log','FontSize',16);
set(gca,'FontSize',16);
%xlim([2,4])
%title(['$ pulse length $ =',num2str(FWHM0),'$ fs $'],'Interpreter','latex','FontSize',20,'FontWeight','bold');
box on
figure(3)
hold on
plot(l,plunit)
xlabel('wavelength (nm)')
ylabel('P(\lambdar) (arb.units)')
set(gca,'FontSize',16);
%xlim([0.795, 0.815])
 %title(['$ BW $ =',num2str(FWHM*100),'$ \% $'],'Interpreter','latex','FontSize',20,'FontWeight','bold');
box on

 figure(4)
hold on
plot(s,bunchs,'linewidth',1)
xlabel('s (m)')
ylabel('bunch')
set(gca,'FontSize',16);
box on




 figure(5)
hold on
plot(z,bunch,'linewidth',1)
xlabel('s (m)')
ylabel('bunch')
set(gca,'FontSize',16);
box on
 






 
 figure(7)

 hold on
plot(z,beamx,z,beamy,'linewidth',2)
set(gca,'FontSize',16);




