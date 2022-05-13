%% slice analysis

function [std_delta,emitx,emity]=slice_enegy_spread(z,delta,x,px,y,py,n_slice)

dz = (max(z)-min(z))/n_slice;
 bins0=min(z):dz:max(z)-dz;
for i = 1:length(bins0)-1
    
  ind = find(z>=bins0(i) & z<= bins0(i+1));
  std_delta(i) = std(delta(ind)); 
    emitx(i)=beam_emittance(x(ind),px(ind));
    emity(i)=beam_emittance(y(ind),py(ind));
    
    
end





end

