
function [outp]=Genesis_outptut_analysis(filename,lambdar,fig)
%filename='../output/radp.out';
%fig=1; time-dep
 %lambdar=13;% nm
 
 [data0,datas,current,xlamd,delts,nslice]=readflie(filename,fig);
% 

 bunch=max(max(datas(:,8,:)));
 current_out=max(current);
 %% std t-duration 
 s=delts*(1:nslice);
powers=reshape(datas(end,1,:), [1 ,length(s)]);
mean_s=sum(s.*powers)/sum(powers);
std_t=sqrt(sum((s-mean_s).^2.*powers)/sum(powers))/3e8;
  power=max(powers);
%% std - bw duraiton 
phase=reshape(datas(end,4,:), [1 ,length(s)]);
[l,plunit]=spectrim(s,powers,phase,nslice,lambdar);
mean_l=sum(l.*plunit)/sum(plunit);
bw=sqrt(sum((l-mean_l).^2.*plunit)/sum(plunit))/(lambdar);
std_e=bw*2*pi*3e8/(lambdar*1e-9)*6.582119569e-16 ; %(eV)%
pluse_e=sum(powers*delts)/3e8;
M=find(powers>=0.5*max(powers));
fwt=(s(M(end))-s(M(1)))/3e8*1e15;
 outp=[power;bunch;current_out;std_t;std_e;pluse_e;bw;fwt];
end