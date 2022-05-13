%% current_for_bunch

function [cur,zn]=bunch_current(z,Q,n_slice,q)
[ct]=physic_constants;

if q ==0
dz = (max(z)-min(z))/n_slice;
 bins0=min(z)+dz/2:dz:max(z)-dz/2;
[ntot,zn0]=hist(z,bins0);
cur0=ntot/(1*dz)*ct.c*Q/sum(ntot);
dz1 = (max(z)-min(z))/n_slice;
zn  = min(z)+dz1/2:dz1:max(z)-dz1/2;
cur = interp1(zn0,cur0,zn,'cubic');
else
 dz = (max(z)-min(z))/n_slice;
 bins0=min(z)+dz/2:dz:max(z)-dz/2;
for ii = 1:length(bins0)-1
    ind = find(z>=bins0(ii)&z<=bins0(ii+1));
    ntot(ii) =sum(q(ind));
end
 cur0=ntot/(1*dz)*ct.c;
 zn0=bins0(1:end-1);
 dz1 = (max(z)-min(z))/n_slice;
zn  = min(z)+dz1/2:dz1:max(z)-dz1/2;  
cur = interp1(zn0,cur0,zn,'cubic');
end




end



