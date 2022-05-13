function [data0,datas,current,xlamd,delts,nslice,datainfo]=readflie(filename,fig)
%for time-depending fig=1,
%read outflie
%Author: Frank.liu


fidin = fopen(filename,'r');
%read xlamd
    line=fgetl(fidin); %use this function one read  line one/use again read line 2 ...so on
   temp=split(line);
   while ~strcmp(temp{1},'xlamd')% find key word xlamd
       line=fgetl(fidin);
      temp=split(line);
   end
   xlamd=str2num(temp{3});
 %read nslice
  while ~strcmp(temp{1},'nslice=')
       line=fgetl(fidin);
      temp=split(line);
   end
   nslice=str2num(temp{2});
   if fig ~=1
        nslice=1;
   end
%read seperation of output slices 
   while ~strcmp(temp{length(temp)},'slices')
       line=fgetl(fidin);
      temp=split(line);
   end
   delts=str2num(temp{1});

% read z[m]=data0(1),aw,qfld
   while ~strcmp(temp{1},'z[m]')
       line=fgetl(fidin);
      temp=split(line);
   end
   line=fgetl(fidin);
   temp=split(line);
   data0=[];
      while ~strcmp(temp{1},'')
          t=str2num(line);
          data0=[data0;t];
          line=fgetl(fidin);
          temp=split(line);
         
   end
%read current and data
 
  datas=[];
current=[]; 
islice=1;
if fig==1
while islice <= nslice-1
      while ~strcmp(temp{length(temp)},'current')
      line=fgetl(fidin);
      temp=split(line);
      end
      current=[current;str2num(temp{1})];
    %read data
      
       while ~strcmp(temp{1},'power')
      line=fgetl(fidin);
      temp=split(line);
       end
        datainfo = temp;
      data1=[];
    while ~strcmp(temp{1},'')
          t1=str2num(line);
         
          data1=[data1;t1];
          line=fgetl(fidin);
          temp=split(line);
    end
    datas(:,:,islice)=data1;
    islice=islice+1;
    
end
end

 while ~strcmp(temp{length(temp)},'current')
      line=fgetl(fidin);
      temp=split(line);
      end
      current=[current;str2num(temp{1})];
      
  while ~strcmp(temp{1},'power')
      line=fgetl(fidin);
      temp=split(line);
  end
      datainfo = temp;
data1=[];      
    while ~feof(fidin)
        line=fgetl(fidin);
          t1=str2num(line);
          
          data1=[data1;t1];    
    end
    datas(:,:,nslice)=data1;  
fclose(fidin);
   
   
     