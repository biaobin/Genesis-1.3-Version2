function beam_generate(tdp,nslice,ntail,beam_param,lambda,np,ildpsi,gajitter)


gamma =beam_param.Ene/0.511e6;  % Gamma, corresponding to beam energy 500 MeV
gamma = gamma*(1+gajitter);
rdegama=beam_param.delE;
dgamma = rdegama*gamma;  % Energy spread, in unit of gamma
lambdau1 = 1e-5;  % undulator period, unit: m
lambdar1 = lambda; % radiation wavelength, unit: m
%aw1 = sqrt(lambdar1*2*gamma^2/lambdau1-1); % undulator strength parameter, dimenionless
aw1=1;

aw0=1e-7*aw1;
nwig0=1;



peakcurrent=beam_param.Ipk*1e-8;  % peak current, unit: A
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

xbeam = beam_param.xbeam;
ybeam = beam_param.ybeam;
 px    = beam_param.px;
 py    = beam_param.py;




NPART= np;



  %% write 1.in file / 
 fid = fopen('../genesis/1.in','w');
 fprintf(fid,'  $NEWRUN \r\n');
 fprintf(fid,' AW0=  %g,\r\n ',aw0);
 fprintf(fid,' AWD=  %g,\r\n ',aw0);
 fprintf(fid,' XKX= 0.0,\r\n ');
 fprintf(fid,' XKY= 1,\r\n ');
  fprintf(fid,' IALL=0,\n ');
    fprintf(fid,' ILDPSI=%g \r\n ',ildpsi(1));
    fprintf(fid,' ildgam=%g \r\n ',ildpsi(2));
    fprintf(fid,' ildx=%g \r\n ',ildpsi(3));
        fprintf(fid,' ildy=%g \r\n ',ildpsi(4));
        fprintf(fid,' ildpx=%g \r\n ',ildpsi(5));
        fprintf(fid,' ildpy=%g \r\n ',ildpsi(6));
    
  fprintf(fid,' XLAMD=  %g \r\n ',lambdau1);
 fprintf(fid,'  NPART= %g \r\n ',NPART);
 fprintf(fid,' GAMMA0=  %g,\r\n ',gamma);
 fprintf(fid,' DELGAM=  %g,\r\n ', dgamma);
 fprintf(fid,' RXBEAM=  %g ,\r\n ',beamx);
 fprintf(fid,' RYBEAM=  %g ,\r\n ',beamy);
 fprintf(fid,' ALPHAX=  %g,\r\n ',alphax);
 fprintf(fid,' ALPHAY=  %g,\r\n ',alphay);
 fprintf(fid,' XBEAM=  %g,\r\n ',xbeam);
 fprintf(fid,' PXBEAM=  %g,\r\n ',px);
 fprintf(fid,' YBEAM=   %g,\r\n ',ybeam);
 fprintf(fid,' PYBEAM=  %g, \r\n ',py);
 fprintf(fid,' EMITX=  %g,\r\n ', ex);
 fprintf(fid,' EMITY=  %g,\r\n ', ey);

 fprintf(fid,' XLAMDS=  %g,\r\n ',lambdar1);
 fprintf(fid,' PRAD0=   %g,\r\n ',0);
 fprintf(fid,' ZRAYL=  1\r\n ');
 fprintf(fid,' ZWAIST= 0.15 \r\n ');
 fprintf(fid,' NCAR= 151,\r\n ');
 fprintf(fid,' rmax0= 9.,\r\n ');
 fprintf(fid,' DELZ= 1,\r\n ');
 fprintf(fid,' ZSEP= 1,\r\n ');
 fprintf(fid,' NSEC= 1,\r\n ');
   fprintf(fid,'itgaus= 1,\r\n ');
   fprintf(fid,' NBINS= 4 \r\n ');
  fprintf(fid,'igamgaus = 1 \r\n ');
  fprintf(fid,' LOUT= 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 1,\r\n ');
    fprintf(fid,'iotail= 1,\r\n ');
    fprintf(fid,' NHARM= 1,\r\n ');
    fprintf(fid,'iallharm = 0,\r\n ');
     fprintf(fid,'iharmsc = 0,\r\n ');
   fprintf(fid,' CURPEAK=  %g,\r\n ',peakcurrent);
 fprintf(fid,' CURLEN= %g,\r\n ',bl);
  fprintf(fid,' ntail= %g,\r\n ',ntail);
   fprintf(fid,' SHOTNOISE= 0 ,\r\n ');
  fprintf(fid,' ITDP= %g,\r\n ',tdp);
   fprintf(fid,'nslice=  %g,\r\n ',nslice);
  fprintf(fid,' ippart= 0,\r\n ');
   fprintf(fid,' IDMPPAR= 1,\r\n ');
 fprintf(fid,' IDMPFLD= 0,\r\n ');
 fprintf(fid,' filetype = ''ORIGINAL'' ,\r\n ');
 fprintf(fid,' OUTPUTFILE= ''../output/1.out'' ,\r\n ');
 fprintf(fid,' MAGINFILE= ''../lat/1.lat'' ,\r\n ');
 fprintf(fid,' IPHSTY=1,\r\n ');
 fprintf(fid,' IWITYP= 0,\r\n ');
 fprintf(fid,'  quadf = 0 ,\r\n ');
  fprintf(fid,'  quadd = 0 ,\r\n ');
  fprintf(fid,'  fl = 0 ,\r\n ');
   fprintf(fid,'  SHOTNOISE= 0 ,\n ');
 fprintf(fid,' $end \r\n');
 fclose(fid);

 % write lat file
 
  fid2 = fopen('../lat/1.lat','w');

 fprintf(fid2,'?VERSION= 1.0 \n');
 fprintf(fid2,'?UNITLENGTH=  %g  \n ',lambdau1);
  fprintf(fid2,'aw %g %g %g \n',aw0,nwig0,0);
  
  fclose(fid2);

unix('echo ../genesis/1.in | ./genesis_c');
 pause(1);
 




end
