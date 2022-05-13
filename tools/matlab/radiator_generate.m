function radiator_generate(tdp,nslice,ntail,beam_param,wavelength,nham,lambdu,Nw,np,gajitter)



gamma =beam_param.Ene/0.511e6;  % Gamma, corresponding to beam energy 500 MeV
gamma = gamma*(1+gajitter);
rdegama=beam_param.delE;
dgamma = rdegama*gamma;  % Energy spread, in unit of gamma

lambdau2 = lambdu;  % undulator period, unit: m
lambdar2 = wavelength; % radiation wavelength, unit: m
aw3 = sqrt(lambdar2*2*gamma^2/lambdau2-1); % undulator strength parameter, dimenionless
nwig3=Nw;
aw3 = 2.69261;
peakcurrent=beam_param.Ipk;  % peak current, unit: A
bl = beam_param.bl; % bunch length, unit: m

% twiss 
betax  = beam_param.betax;
betay  = beam_param.betay;
alphax = beam_param.alphax;
alphay = beam_param.alphay;
emitx=beam_param.emit_x;
emity=beam_param.emit_y;
ex = emitx*gamma;
ey= emity*gamma;% normalized emittance, unit: m.rad


beamx = sqrt(emitx*betax);
beamy = sqrt(emity*betay);


 %% write rad file
 fid = fopen('../genesis/rad.in','w');
 fprintf(fid,'  $NEWRUN \r\n');
 fprintf(fid,'  partfile = ''../output/modu_out.dpa'' \r\n');
 fprintf(fid,'  trama   = 0 \r\n');
 fprintf(fid,'  itram16   = %g \r\n',0);
 fprintf(fid,'convharm= %g,\r\n ',nham);
 fprintf(fid,'multconv=0,\r\n ');
 fprintf(fid,' AW0=  %g,\r\n ',aw3);
 fprintf(fid,' AWD=  %g,\r\n ',aw3);
 fprintf(fid,' XKX= 0.0,\r\n ');
 fprintf(fid,' XKY= 1,\r\n ');
 fprintf(fid,' XLAMD=  %g \r\n ',lambdau2);
 fprintf(fid,'  NPART= %g \r\n ',np);
 fprintf(fid,' GAMMA0=  %g,\r\n ',gamma);
 fprintf(fid,' DELGAM=  %g,\r\n ',dgamma);
 fprintf(fid,' RXBEAM=  %g ,\r\n ',beamx);
 fprintf(fid,' RYBEAM=  %g ,\r\n ',beamy);
 fprintf(fid,' ALPHAX=  %g,\r\n ',alphax);
 fprintf(fid,' ALPHAY=  %g,\r\n ',alphay);
 fprintf(fid,' XBEAM=  0.,\r\n ');
 fprintf(fid,' PXBEAM=  0.,\r\n ');
 fprintf(fid,' YBEAM=  0.,\r\n ');
 fprintf(fid,' PYBEAM=  0., \r\n ');
 fprintf(fid,' EMITX=  %g,\r\n ',ex);
 fprintf(fid,' EMITY=  %g,\r\n ',ey);
 fprintf(fid,' XLAMDS=  %g,\r\n ',lambdar2);
 fprintf(fid,' PRAD0=   %g,\r\n ',0);
 fprintf(fid,' ZRAYL=  10 \r\n ');
 fprintf(fid,' ZWAIST= 0 \r\n ');
 fprintf(fid,' NCAR= 151,\r\n ');
 fprintf(fid,' rmax0= 9.,\r\n ');
 fprintf(fid,' DELZ= 1,\r\n ');
 fprintf(fid,' ZSEP=%g,\r\n ',nham);
 fprintf(fid,' NSEC= 1,\r\n ');
   fprintf(fid,'itgaus= 1,\r\n ');
   fprintf(fid,' NBINS= 4 \r\n ');
  fprintf(fid,'igamgaus = 1 \r\n ');
  fprintf(fid,' LOUT= 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 1,\r\n ');
  fprintf(fid,'iotail= 1,\r\n ');
    fprintf(fid,' NHARM= 1,\r\n ');
    fprintf(fid,'iallharm = 1,\r\n ');
     fprintf(fid,'iharmsc = 1,\r\n ');
   fprintf(fid,' CURPEAK=  %g,\r\n ',peakcurrent);
 fprintf(fid,' CURLEN= %g,\r\n ',bl);
 fprintf(fid,' ntail= %g,\r\n ',ntail);
   fprintf(fid,' SHOTNOISE= 1 ,\r\n ');
  fprintf(fid,' ITDP= %g,\r\n ',tdp);
   fprintf(fid,'nslice=  %g,\r\n ',nslice);
  fprintf(fid,' ippart= 0,\r\n ');
   fprintf(fid,' IDMPPAR= 1,\r\n ');
 fprintf(fid,' IDMPFLD= 1,\r\n ');
 fprintf(fid,' filetype = ''ORIGINAL'' ,\r\n ');
 fprintf(fid,' OUTPUTFILE= ''../output/radp.out'' ,\r\n ');
 fprintf(fid,' MAGINFILE= ''../lat/rad1.lat'' ,\r\n ');
 fprintf(fid,' BEAMFILE= ''../output/curr.beam'',\r\n');
 fprintf(fid,' IPHSTY=1,\r\n ');
 fprintf(fid,' IWITYP= 0,\r\n ');
 fprintf(fid,'  quadf = 0 ,\r\n ');
  fprintf(fid,'  quadd = 0 ,\r\n ');
 fprintf(fid,' $end \r\n');
 fclose(fid);
 
%aw3 = 2.72156;

 % write lat file   
fid2 = fopen('../lat/rad1.lat','w');

 fprintf(fid2,'?VERSION= 1.0 \n');
 fprintf(fid2,'?UNITLENGTH=  %g  \n ',lambdau2);
 fprintf(fid2,'aw %g %g %g \n',aw3,nwig3,0);
 fclose(fid2);    
    
pause(1);   
unix('echo ../genesis/rad.in | ./genesis_c');









end