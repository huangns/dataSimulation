%const value
%Twc
clear;
clc;
format long;
fidwi=fopen('pwi.txt','w');fclose(fidwi);
fidwi = fopen('pwi.txt','a+');
fidam=fopen('am.txt','w');fclose(fidam);
fidam = fopen('am.txt','a+');
fidwm=fopen('wm.txt','w');fclose(fidwm);
fidwm = fopen('wm.txt','a+');
fidpcm=fopen('pcm.txt','w');fclose(fidpcm);
fidpcm = fopen('pcm.txt','a+');
fidqcm=fopen('qcm.txt','w');fclose(fidqcm);
fidqcm = fopen('qcm.txt','a+');

fidba=fopen('ba.txt','w');fclose(fidba);
fidba = fopen('ba.txt','a+');
fidbw=fopen('bw.txt','w');fclose(fidbw);
fidbw = fopen('bw.txt','a+');
fidqiw=fopen('qiw.txt','w');fclose(fidqiw);
fidqiw = fopen('qiw.txt','a+');
fidvwi=fopen('vwi.txt','w');fclose(fidvwi);
fidvwi = fopen('vwi.txt','a+');

qwc=[0.1527 0.3481 -0.8499 0.3649];
qwc = quatnormalize(qwc);
Rwc=quat2dcm(qwc);
twc=[0 0 0]';
Twc=[Rwc twc;0 0 0 1];
%Tmi
qmi=[0.7 0 0 -0.7];
qmi = quatnormalize(qmi);
%Rmi=[1 0 0;0 1 0;0 0 1];
Rmi=quat2dcm(qmi);
tmi=[0.05 0.05 0.01]';
Tmi=[Rmi tmi;0 0 0 1];

%IMU

bw{1}=[0.00006 -0.00008 0.00001]';
ba{1}=[-0.18 0.15 0.3]';
sig_na=0.0025*[1 1 1]';
sig_nba=sig_na./10;
sig_nw=0.000025*[1 1 1]';
sig_nbw=sig_nw./10;



g=[0 0 9.81]';

%init value
qcm_init=[-0.0193613 0.932408 0.244601 -0.182484];
qcm_init=quatnormalize(qcm_init);
Rcm_init=quat2dcm(qcm_init);
tcm_init=[-0.00985595 0.0310215 0.396001]';
Tcm_init=[Rcm_init tcm_init;0 0 0 1];

Twi_init=Twc*Tcm_init*Tmi;

idt=0.00001;
idtt=0.01;
vdt=1/30;
w{1}=[0 0 0]';
a{1}=Twi_init(1:3,1:3)'*g;

wm{1}=w{1}+bw{1}+sig_nw*randn;
am{1}=a{1}+ba{1}+sig_na*randn;

qiw{1}=(dcm2quat(Twi_init(1:3,1:3)'))';
pwi{1}=Twi_init(1:3,4);
vwi{1}=[0,0,0]';
(dcm2quat(Twi_init(1:3,1:3)'))';

fprintf(fidwi,'%f %f %f %f \r\n',0.0,pwi{1}(1,1),pwi{1}(2,1),pwi{1}(3,1));
%fprintf(fidqiw,'%f %f %f %f %f \r\n',0.0,qiw{1}(1,1),qiw{1}(2,1),qiw{1}(3,1),qiw{1}(4,1));

fprintf(fidam,'%f %f %f %f \r\n',0.0,am{1}(1,1),am{1}(2,1),am{1}(3,1));
fprintf(fidwm,'%f %f %f %f \r\n',0.0,wm{1}(1,1),wm{1}(2,1),wm{1}(3,1));
fprintf(fidba,'%f %f %f %f \r\n',0.0,ba{1}(1,1),ba{1}(2,1),ba{1}(3,1));
fprintf(fidbw,'%f %f %f %f \r\n',0.0,bw{1}(1,1),bw{1}(2,1),bw{1}(3,1));
fprintf(fidvwi,'%f %f %f %f \r\n',0.0,vwi{1}(1,1),vwi{1}(2,1),vwi{1}(3,1));



%init static state
Rcm{1}=Rcm_init;
ci=2;
cii=2;
vi=2;
count_v=0;
timesum=idt;
timesumdt=idt;

%vision part
pcm{1}=tcm_init;
qcm{1}=qcm_init';
sig_pcm=0.001*[1 1 1]';
sig_qcm=0.001*[1 1 1]';
pcm_m{1}=pcm{1}(:,:)+sig_pcm*randn;
deltaq=quatFromSmallAngle(sig_qcm*randn);
qcm_m{1}=dcm2quat(quat2dcm(deltaq)*quat2dcm(qcm{1}(:,:)'))';
%qcm{1}(:,:)
%qcm_m{1}(:,:)
%pcm_m{1}(:,:)
%pcm{1}(:,:)
DV{1}=[0 0 0]';

fprintf(fidpcm,'%f %f %f %f \r\n',0.0,pcm_m{1}(1,1),pcm_m{1}(2,1),pcm_m{1}(3,1));
fprintf(fidqcm,'%f %f %f %f %f\r\n',0.0,qcm_m{1}(1,1),qcm_m{1}(2,1),qcm_m{1}(3,1),qcm_m{1}(4,1));
fprintf(fidqiw,'%f %f %f %f %f\r\n',0.0,qiw{1}(1,1),qiw{1}(2,1),qiw{1}(3,1),qiw{1}(4,1));

%fprintf(fidwm,'%f %f %f %f \r\n',0.0,wm{1}(1,1),wm{1}(2,1),wm{1}(3,1));
cilast=1;
for t=idt:idt:5-idt
%for t=idt:idt:0.02
    if(ci>1000)
        clast=1000;
        ci=1;
    else
        clast=ci-1;
    end
    a{ci}=a{clast};
    
    w{ci}=w{clast};
    
    a{ci}(:,:)';
    %atemp =a{ci}(1:3,1)'
    %w{ci}(1:3,1)' 
    ba{ci}=ba{clast}+sig_nba*randn*idt;
    bw{ci}=bw{clast}+sig_nbw*randn*idt;
    am{ci}=a{ci}+ba{ci}+sig_na*randn;
    wm{ci}=w{ci}+bw{ci}+sig_nw*randn;
    
    ew=w{ci}(:,:);
    ea=a{ci}(:,:);
    ewold=w{clast}(:,:);
    eaold=a{clast}(:,:);
    Omega=omegaMatJPL(ew);
    OmegaOld=omegaMatJPL(ewold);
    OmegaMean=omegaMatJPL((ew+ewold)/2);
    
    div=1;
    MatExp=eye(4);
    OmegaMean=OmegaMean*0.5*idt;
    for ki=1:1:4
        div=div*ki;
        MatExp=MatExp+OmegaMean./div;
        OmegaMean=OmegaMean*OmegaMean;
    end
    
    quat_int=MatExp+1.0/48.0*(Omega*OmegaOld-OmegaOld*Omega)*idt*idt;
    qcoeffs1=[qiw{clast}(2:4,1)',qiw{clast}(1,1)]';
    qcoeffs2=quat_int*qcoeffs1;
    qiw{ci}=( quatnormalize([qcoeffs2(4,1),qcoeffs2(1:3,1)']) )';
    
    dv=(quat2dcm(qiw{ci}(:,:)')' * ea +  quat2dcm(qiw{clast}(:,:)')'*eaold)/2;
    DV{ci}=dv-g;
    vwi{ci}=vwi{clast}(:,:)+(dv-g)*idt;
    pwi{ci}=pwi{clast}(:,:)+(vwi{ci}(:,:)+vwi{clast}(:,:))/2.0*idt;
    
    DV{ci}=( vwi{ci}(:,:) + vwi{clast}(:,:) )/2.0*idt;
    

   
    timesum=timesum+idt;
 if(timesum>0.03)
        %Twc*Tcm_init*Tmi
        tempTwi=[ (quat2dcm(qiw{ci}(:,:)'))' pwi{ci}(:,:);0 0 0 1];
        temptcm= inv(Twc)*tempTwi*inv(Tmi);
        pcm{vi}=temptcm(1:3,4);
        qcm{vi}=dcm2quat(temptcm(1:3,1:3))';
        pcm_m{vi}=pcm{vi}(:,:)+sig_pcm*randn;
        deltaq_=quatFromSmallAngle(sig_qcm*randn);
        qcm_m{vi}=dcm2quat(quat2dcm(deltaq_)*quat2dcm(qcm{vi}(:,:)'))';
        
        %tempTwi
        timesum=idt;
       % pcm_m{vi}(:,:)
       % qcm_m{vi}(:,:)
       
       
        fprintf(fidpcm,'%f %f %f %f \r\n',t,pcm_m{vi}(1,1),pcm_m{vi}(2,1),pcm_m{vi}(3,1));
        fprintf(fidqcm,'%f %f %f %f %f \r\n',t,qcm_m{vi}(1,1),qcm_m{vi}(2,1),qcm_m{vi}(3,1), qcm_m{vi}(4,1));
       
       
        %vi=vi+1;
 end
    timesumdt=timesumdt+idt;
 if(timesumdt>0.01)
    fprintf(fidwi,'%f %f %f %f \r\n',t,pwi{ci}(1,1),pwi{ci}(2,1),pwi{ci}(3,1));
    fprintf(fidam,'%f %f %f %f \r\n',t,am{ci}(1,1),am{ci}(2,1),am{ci}(3,1));
    fprintf(fidwm,'%f %f %f %f \r\n',t,wm{ci}(1,1),wm{ci}(2,1),wm{ci}(3,1));
    
    fprintf(fidba,'%f %f %f %f \r\n',t,ba{ci}(1,1),ba{ci}(2,1),ba{ci}(3,1));
    fprintf(fidbw,'%f %f %f %f \r\n',t,bw{ci}(1,1),bw{ci}(2,1),bw{ci}(3,1));
    fprintf(fidvwi,'%f %f %f %f \r\n',t,vwi{ci}(1,1),vwi{ci}(2,1),vwi{ci}(3,1));
    fprintf(fidqiw,'%f %f %f %f %f\r\n',t,qiw{ci}(1,1),qiw{ci}(2,1),qiw{ci}(3,1),qiw{ci}(4,1));
    
    
    timesumdt=idt;
    %cii=cii+1;
   end
 %fprintf(fidam,'%f %f %f %f \r\n',t,am{ci}(1,1),am{ci}(2,1),am{ci}(3,1));
    ci=ci+1;
end

initv=(5);
lastax=0.1;
lastay=0.13;
lastaz=0.2;
for t=5:idt:9-idt
     if(ci>1000)
        clast=1000;
        ci=1;
    else
        clast=ci-1;
    end
   % a{ci}=a{ci-1}(:,:)+[randn*0.01 randn*0.001 randn*0.1]'*0.0000001;
    %a{ci}=[0.5*cos(5*t+1),0.5*sin(5*t-1),-0.5*sin(5*t-2)]'+quat2dcm(qiw{ci-1}(:,:)')*g;
    %(quat2dcm(qiw{ci-1}(:,:)')*g)';
   % disp('aci');
  %  a{ci}(:,:);
  %  disp('wci');
    %w{ci}=[randn*0.5*sin(2.5*(t-initv)) randn*0.6*sin(1.5*(t-initv)) randn*0.6*sin(3*(t-initv))]';
    w{ci}=[0.95*sin(2.5*(t-initv)) -0.9*sin(1.5*(t-initv)) 0.62*sin(3*(t-initv))]';
    %w{ci}=w{ci-1};
    w{ci}(:,:);
    % w{ci}=w{ci-1}+[randn*0.1 randn*0.1 randn]'*0.000001;
   % a{ci}=a{ci-1}(:,:)+randn*[0.002*sin(4*t+1),0.001*sin(10*t+3),0.001*sin(-6*t+2)]'/100.0;
    %w{ci}=w{ci-1}+[randn*0.25*sin(t+1),randn*0.3*sin(1.5*t+3),randn*0.5*cos(-2*t+10)]';
    % a{ci}=a{ci-1}(:,:)+[randn*0.05,randn*0.05,randn*0.1]';
    %w{ci}=w{ci-1}+[randn*0.25*sin(t+1),randn*0.3*sin(1.5*t+3),randn*0.5*cos(-2*t+10)]';
   % a{ci}(:,:)';
    %atemp =a{ci}(1:3,1)'
    %w{ci}(1:3,1)' 
    ba{ci}=ba{clast}+sig_nba*randn*idt;
    bw{ci}=bw{clast}+sig_nbw*randn*idt;
   
    wm{ci}=w{ci}+bw{ci}+sig_nw*randn;
    ew=w{ci}(:,:);
    ewold=w{clast}(:,:);
    Omega=omegaMatJPL(ew);
    OmegaOld=omegaMatJPL(ewold);
    OmegaMean=omegaMatJPL((ew+ewold)/2);
    
    div=1;
    MatExp=eye(4);
    OmegaMean=OmegaMean*0.5*idt;
    for ki=1:1:4
        div=div*ki;
        MatExp=MatExp+OmegaMean./div;
        OmegaMean=OmegaMean*OmegaMean;
    end
    
    quat_int=MatExp+1.0/48.0*(Omega*OmegaOld-OmegaOld*Omega)*idt*idt;
    qcoeffs1=[qiw{clast}(2:4,1)',qiw{clast}(1,1)]';
    qcoeffs2=quat_int*qcoeffs1;
    qiw{ci}=( quatnormalize([qcoeffs2(4,1),qcoeffs2(1:3,1)']) )';
    %disp('dv');
   a{ci}=quat2dcm(qiw{ci}(:,:)')*[0,0.1,0]'+quat2dcm(qiw{ci}(:,:)')*g;
   % a{ci}=quat2dcm(qiw{ci}(:,:)')*[0.4*sin(pi*(t-initv)+asin(lastax))-0.4*lastax,0.4*sin(pi*(t-initv)+asin(lastay))-0.4*lastay,0.85*sin(pi*(t-initv)+asin(lastaz))-0.85*lastaz]'+quat2dcm(qiw{ci}(:,:)')*g;
   % lastax=0.4*cos(pi*(t-initv));
  %  lastay=0.4*cos(pi*(t-initv));
  %  lastaz=0.85*cos(pi*(t-initv));
    %disp('awi')
   % [0.4*cos(3*(t-initv)),0.4*cos(3*(t-initv)),0.85*cos(4*(t-initv))]
    am{ci}=a{ci}+ba{ci}+sig_na*randn;
    ea=a{ci}(:,:);
    eaold=a{clast}(:,:);
     
    dv=(quat2dcm(qiw{ci}(:,:)')' * ea +  quat2dcm(qiw{clast}(:,:)')'*eaold)/2;
    dv';
    %DV{ci}=dv-g;
    vwi{ci}=vwi{clast}(:,:)+(dv-g)*idt;
    %disp('dv-g:')
    %(dv-g)'
    %disp('vwi');
   % vwi{ci-1}(:,:)';
    pwi{ci}=pwi{clast}(:,:)+( vwi{ci}(:,:) + vwi{clast}(:,:) )/2.0*idt;
    DV{ci}=( vwi{ci}(:,:) + vwi{clast}(:,:) )/2.0*idt;
%fprintf(fid,'%f %f %f %f \r\n',t,pwi{ci}(1,1),pwi{ci}(2,1),pwi{ci}(3,1));
    %disp('pwi=');
    timesum=timesum+idt;
    if(timesum>0.03)

        %Twc*Tcm_init*Tmi
        tempTwi=[ (quat2dcm(qiw{ci}(:,:)'))' pwi{ci}(:,:);0 0 0 1];
        temptcm= inv(Twc)*tempTwi*inv(Tmi);
        pcm{vi}=temptcm(1:3,4);
        qcm{vi}=dcm2quat(temptcm(1:3,1:3))';
        pcm_m{vi}=pcm{vi}(:,:)+sig_pcm*randn;
        deltaq_=quatFromSmallAngle(sig_qcm*randn);
        qcm_m{vi}=dcm2quat(quat2dcm(deltaq_)*quat2dcm(qcm{vi}(:,:)'))';
        % timesum=timesum-0.03;
        %tempTwi
       timesum=idt;
       % pcm_m{vi}(:,:)
       % qcm_m{vi}(:,:)
       
        fprintf(fidpcm,'%f %f %f %f \r\n',t,pcm_m{vi}(1,1),pcm_m{vi}(2,1),pcm_m{vi}(3,1));
        fprintf(fidqcm,'%f %f %f %f %f \r\n',t,qcm_m{vi}(1,1),qcm_m{vi}(2,1),qcm_m{vi}(3,1), qcm_m{vi}(4,1));
       
      %  vi=vi+1;
    end
    
    %fprintf(fidwi,'%f %f %f %f \r\n',t,pwi{ci}(1,1),pwi{ci}(2,1),pwi{ci}(3,1));
   % fprintf(fidam,'%f %f %f %f \r\n',t,am{ci}(1,1),am{ci}(2,1),am{ci}(3,1));
    %fprintf(fidwm,'%f %f %f %f \r\n',t,wm{ci}(1,1),wm{ci}(2,1),wm{ci}(3,1));
    %    timesumdt=timesumdt+idt;
    timesumdt=timesumdt+idt;
    if(timesumdt>0.01)
    fprintf(fidwi,'%f %f %f %f \r\n',t,pwi{ci}(1,1),pwi{ci}(2,1),pwi{ci}(3,1));
    fprintf(fidam,'%f %f %f %f \r\n',t,am{ci}(1,1),am{ci}(2,1),am{ci}(3,1));
    fprintf(fidwm,'%f %f %f %f \r\n',t,wm{ci}(1,1),wm{ci}(2,1),wm{ci}(3,1));
    
    fprintf(fidba,'%f %f %f %f \r\n',t,ba{ci}(1,1),ba{ci}(2,1),ba{ci}(3,1));
    fprintf(fidbw,'%f %f %f %f \r\n',t,bw{ci}(1,1),bw{ci}(2,1),bw{ci}(3,1));
    fprintf(fidvwi,'%f %f %f %f \r\n',t,vwi{ci}(1,1),vwi{ci}(2,1),vwi{ci}(3,1));
    fprintf(fidqiw,'%f %f %f %f %f\r\n',t,qiw{ci}(1,1),qiw{ci}(2,1),qiw{ci}(3,1),qiw{ci}(4,1));
    
    
   timesumdt=idt;
%    cii=cii+1;
    end
    ci=ci+1;
end





    lastax=0.4*cos(pi*(t-initv));
    lastay=0.4*cos(pi*(t-initv));
    lastaz=0.85*cos(pi*(t-initv));

    lastwx=0.95*sin(2.5*(t-initv));
    lastwy=-0.9*sin(1.5*(t-initv));
    lastwz=0.62*sin(3*(t-initv));
initv=(9);    
    
    for t=9:idt:15
     if(ci>1000)
        clast=1000;
        ci=1;
    else
        clast=ci-1;
    end
   % a{ci}=a{ci-1}(:,:)+[randn*0.01 randn*0.001 randn*0.1]'*0.0000001;
    %a{ci}=[0.5*cos(5*t+1),0.5*sin(5*t-1),-0.5*sin(5*t-2)]'+quat2dcm(qiw{ci-1}(:,:)')*g;
    %(quat2dcm(qiw{ci-1}(:,:)')*g)';
   % disp('aci');
  %  a{ci}(:,:);
  %  disp('wci');
    %w{ci}=[randn*0.5*sin(2.5*(t-initv)) randn*0.6*sin(1.5*(t-initv)) randn*0.6*sin(3*(t-initv))]';
   % w{ci}=[0.5*sin(2.5*(t-initv)) -0.6*sin(1.5*(t-initv)) 0.6*sin(3*(t-initv))]';
   
    w{ci}=[0.95*sin(2.5*(t-initv)+asin(lastwx)), -0.9*sin(1.5*(t-initv)+asin(lastwy)), 0.62*sin(3*(t-initv)+asin(lastwz))]';
    %w{ci}=w{ci-1};
    w{ci}(:,:);
    % w{ci}=w{ci-1}+[randn*0.1 randn*0.1 randn]'*0.000001;
   % a{ci}=a{ci-1}(:,:)+randn*[0.002*sin(4*t+1),0.001*sin(10*t+3),0.001*sin(-6*t+2)]'/100.0;
    %w{ci}=w{ci-1}+[randn*0.25*sin(t+1),randn*0.3*sin(1.5*t+3),randn*0.5*cos(-2*t+10)]';
    % a{ci}=a{ci-1}(:,:)+[randn*0.05,randn*0.05,randn*0.1]';
    %w{ci}=w{ci-1}+[randn*0.25*sin(t+1),randn*0.3*sin(1.5*t+3),randn*0.5*cos(-2*t+10)]';
   % a{ci}(:,:)';
    %atemp =a{ci}(1:3,1)'
    %w{ci}(1:3,1)' 
    ba{ci}=ba{clast}+sig_nba*randn*idt;
    bw{ci}=bw{clast}+sig_nbw*randn*idt;
   
    wm{ci}=w{ci}+bw{ci}+sig_nw*randn;
    ew=w{ci}(:,:);
    ewold=w{clast}(:,:);
    Omega=omegaMatJPL(ew);
    OmegaOld=omegaMatJPL(ewold);
    OmegaMean=omegaMatJPL((ew+ewold)/2);
    
    div=1;
    MatExp=eye(4);
    OmegaMean=OmegaMean*0.5*idt;
    for ki=1:1:4
        div=div*ki;
        MatExp=MatExp+OmegaMean./div;
        OmegaMean=OmegaMean*OmegaMean;
    end
    
    quat_int=MatExp+1.0/48.0*(Omega*OmegaOld-OmegaOld*Omega)*idt*idt;
    qcoeffs1=[qiw{clast}(2:4,1)',qiw{clast}(1,1)]';
    qcoeffs2=quat_int*qcoeffs1;
    qiw{ci}=( quatnormalize([qcoeffs2(4,1),qcoeffs2(1:3,1)']) )';
    %disp('dv');
   a{ci}=quat2dcm(qiw{ci}(:,:)')*[-0.4*cos((t-initv)),-0.4*sin(0.2*pi*(t-initv)),0.85*cos(0.5*pi*(t-initv))]'+quat2dcm(qiw{ci}(:,:)')*g;
   % a{ci}=quat2dcm(qiw{ci}(:,:)')*[0.4*sin(pi*(t-initv)+asin(lastax))-0.4*lastax,0.4*sin(pi*(t-initv)+asin(lastay))-0.4*lastay,0.85*sin(pi*(t-initv)+asin(lastaz))-0.85*lastaz]'+quat2dcm(qiw{ci}(:,:)')*g;
   % lastax=0.4*cos(pi*(t-initv));
  %  lastay=0.4*cos(pi*(t-initv));
  %  lastaz=0.85*cos(pi*(t-initv));
    %disp('awi')
   % [0.4*cos(3*(t-initv)),0.4*cos(3*(t-initv)),0.85*cos(4*(t-initv))]
    am{ci}=a{ci}+ba{ci}+sig_na*randn;
    ea=a{ci}(:,:);
    eaold=a{clast}(:,:);
     
    dv=(quat2dcm(qiw{ci}(:,:)')' * ea +  quat2dcm(qiw{clast}(:,:)')'*eaold)/2;
    dv';
    %DV{ci}=dv-g;
    vwi{ci}=vwi{clast}(:,:)+(dv-g)*idt;
    %disp('dv-g:')
    %(dv-g)'
    %disp('vwi');
   % vwi{ci-1}(:,:)';
    pwi{ci}=pwi{clast}(:,:)+( vwi{ci}(:,:) + vwi{clast}(:,:) )/2.0*idt;
    DV{ci}=( vwi{ci}(:,:) + vwi{clast}(:,:) )/2.0*idt;
%fprintf(fid,'%f %f %f %f \r\n',t,pwi{ci}(1,1),pwi{ci}(2,1),pwi{ci}(3,1));
    %disp('pwi=');
    timesum=timesum+idt;
    if(timesum>0.03)

        %Twc*Tcm_init*Tmi
        tempTwi=[ (quat2dcm(qiw{ci}(:,:)'))' pwi{ci}(:,:);0 0 0 1];
        temptcm= inv(Twc)*tempTwi*inv(Tmi);
        pcm{vi}=temptcm(1:3,4);
        qcm{vi}=dcm2quat(temptcm(1:3,1:3))';
        pcm_m{vi}=pcm{vi}(:,:)+sig_pcm*randn;
        deltaq_=quatFromSmallAngle(sig_qcm*randn);
        qcm_m{vi}=dcm2quat(quat2dcm(deltaq_)*quat2dcm(qcm{vi}(:,:)'))';
        % timesum=timesum-0.03;
        %tempTwi
       timesum=idt;
       % pcm_m{vi}(:,:)
       % qcm_m{vi}(:,:)
       
        fprintf(fidpcm,'%f %f %f %f \r\n',t,pcm_m{vi}(1,1),pcm_m{vi}(2,1),pcm_m{vi}(3,1));
        fprintf(fidqcm,'%f %f %f %f %f \r\n',t,qcm_m{vi}(1,1),qcm_m{vi}(2,1),qcm_m{vi}(3,1), qcm_m{vi}(4,1));
       
      %  vi=vi+1;
    end
    
    %fprintf(fidwi,'%f %f %f %f \r\n',t,pwi{ci}(1,1),pwi{ci}(2,1),pwi{ci}(3,1));
   % fprintf(fidam,'%f %f %f %f \r\n',t,am{ci}(1,1),am{ci}(2,1),am{ci}(3,1));
    %fprintf(fidwm,'%f %f %f %f \r\n',t,wm{ci}(1,1),wm{ci}(2,1),wm{ci}(3,1));
    %    timesumdt=timesumdt+idt;
    timesumdt=timesumdt+idt;
    if(timesumdt>0.01)
    fprintf(fidwi,'%f %f %f %f \r\n',t,pwi{ci}(1,1),pwi{ci}(2,1),pwi{ci}(3,1));
    fprintf(fidam,'%f %f %f %f \r\n',t,am{ci}(1,1),am{ci}(2,1),am{ci}(3,1));
    fprintf(fidwm,'%f %f %f %f \r\n',t,wm{ci}(1,1),wm{ci}(2,1),wm{ci}(3,1));
    
    fprintf(fidba,'%f %f %f %f \r\n',t,ba{ci}(1,1),ba{ci}(2,1),ba{ci}(3,1));
    fprintf(fidbw,'%f %f %f %f \r\n',t,bw{ci}(1,1),bw{ci}(2,1),bw{ci}(3,1));
    fprintf(fidvwi,'%f %f %f %f \r\n',t,vwi{ci}(1,1),vwi{ci}(2,1),vwi{ci}(3,1));
    fprintf(fidqiw,'%f %f %f %f %f\r\n',t,qiw{ci}(1,1),qiw{ci}(2,1),qiw{ci}(3,1),qiw{ci}(4,1));
    
    
   timesumdt=idt;
%    cii=cii+1;
    end
    ci=ci+1;
end






fclose(fidwi);
fclose(fidam);
fclose(fidwm);
fclose(fidpcm);
fclose(fidqcm);



fclose(fidba);
fclose(fidbw);
fclose(fidqiw);
fclose(fidvwi);



