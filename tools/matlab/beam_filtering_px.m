%% beam filtering (tr_mometum)

function [new]=beam_filtering_px(old,Q,down_px,up_px)

px = old(:,2);
py = old(:,4);
np_old=length(px);
ind = find(px>=down_px & px <= up_px& py>=down_px & py <= up_px);
x = old(ind,1);
px = old(ind,2);
y = old(ind,3);
py = old(ind,4);
z  = old(ind,5);
gamma = old(ind,6);
Q = Q/np_old*length(gamma);


new = struct('x',x,'px',px,'y',y,'py',py,'z',z,'gamma',gamma,'Q',Q);
end