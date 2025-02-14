%This code is used to analyze the influence of rotation axis inclination Angle on MDU
clc;
clear all;

% % %    z/|\  /y
% % %      |  /
% % %      | / 
% % %       ---------->x
% % %      /
% % %     /  The ray source is located on the -y axis
% % %    /   

u0 = 0;
v0 = 0;

Sx = u0;
Sy = -1000;  %SDD=1000��mm��
SDD=abs(Sy);
Sz = v0;
zuobiao_S = [Sx,Sy,Sz];

Ox = u0;
Oy = -500;  %SOD=500(mm)
Oz = v0;
zuobiao_O = [Ox,Oy,Oz];

objVoxel = 0.05;
detVoxel = 0.1;

tm=zeros(500,500,5);
  a =1;
  b =1;
  c =3;
tm(a,b,c)=10000;

[obj_lenth  obj_width obj_height] = size(tm);
obj_lenth_center = (obj_lenth+1)/2;
obj_width_center = (obj_width+1)/2;
obj_height_center = (obj_height+1)/2;  
d2r = pi/180; 

Rsitemax=46;
nummax=50;

s1site=zeros(Rsitemax,nummax); 
s1fai=zeros(Rsitemax,nummax); 
s1eta=zeros(Rsitemax,nummax); 
s1Dx=zeros(Rsitemax,nummax); 
s1Dz=zeros(Rsitemax,nummax); 
s1Dy=zeros(Rsitemax,nummax); 
s1Ry=zeros(Rsitemax,nummax); 
sinu=zeros(1,360);     
sinv=zeros(1,360);     

for Rsitea=0:1:Rsitemax-1
     Rsite=Rsitea*pi/180;
for num=1:nummax
    site=0.01*num*pi/180;
    fai=0.01*num*pi/180;
    eta=0.01*num*pi/180;
    dx=0.1*num*detVoxel;
    dz=0.1*num*detVoxel;
    dy=0.1*num*detVoxel;
    ry=0.1*num*detVoxel;
  
for sita = 0:1:359 
    
                l = (a-obj_lenth_center)*objVoxel;  %z
                w = (b-obj_width_center)*objVoxel;  %x
                h = (c-obj_height_center)*objVoxel; %y
                
                x0 = w*cos(sita*d2r)+h*sin(sita*d2r);  
                y0 = -w*sin(sita*d2r)+h*cos(sita*d2r);
                objXYZ = [x0;y0;l];  
              
                Rx = [1 0 0; 0 cos(Rsite) sin(Rsite); 0 -sin(Rsite) cos(Rsite)];
                wRotateXYZ = Rx*objXYZ; 
                
                x1 = wRotateXYZ(1)+Ox;
                y1 = wRotateXYZ(2)+Oy;
                z1 = wRotateXYZ(3)+Oz;
                
                dy0 = -Sy/(y1-Sy);  
                du = dy0*(x1-Sx)+Sx; 
                dv = dy0*(z1-Sz)+Sz; 
  
                sinu(sita+1)=du;     
                sinv(sita+1)=dv;     
               
                du1=du*(1-dv*sin(site)/(SDD*cos(site)+dv*sin(site)));
                dv1=dv*SDD/(SDD*cos(site)+dv*sin(site));
       
                sinu1(sita+1)=du1;     
                sinv1(sita+1)=dv1;    
 
                du2=du*SDD/(SDD*cos(fai)+du*sin(fai));
                dv2=dv*(1-du*sin(fai)/(SDD*cos(fai)+du*sin(fai)));
              
                sinu2(sita+1)=du2;     
                sinv2(sita+1)=dv2;     
            
                du3=du*cos(eta)+dv*sin(eta);
                dv3=-du*sin(eta)+dv*cos(eta);
              
                sinu3(sita+1)=du3;    
                sinv3(sita+1)=dv3;    
               
                du4 = du-dx; 
                dv4 = dv; 
                sinu4(sita+1)=du4;     
                sinv4(sita+1)=dv4; 
           
                du5 = du; 
                dv5 = dv-dz; 
                sinu5(sita+1)=du5;     
                sinv5(sita+1)=dv5; 
            
                dy6 = (-Sy+dy)/(y1-Sy);  
                du6 = dy6*(x1-Sx)+Sx; 
                dv6 = dy6*(z1-Sz)+Sz; 
                sinu6(sita+1)=du6;     
                sinv6(sita+1)=dv6; 
           
                dy7 = -Sy/(y1-Sy+ry);   
                du7 = dy7*(x1-Sx)+Sx; 
                dv7 = dy7*(z1-Sz)+Sz; 
                sinu7(sita+1)=du7;     
                sinv7(sita+1)=dv7; 
              
end
  siteu= sinu1-sinu;
  sitev= sinv1-sinv;
  msiteu=max(abs(siteu));
  msitev=max(abs(sitev));
  s1site(Rsitea+1,num)=max(msiteu,msitev);
 
  
  faiu= sinu2-sinu;
  faiv= sinv2-sinv;
  mfaiu=max(abs(faiu));
  mfaiv=max(abs(faiv));
  s1fai(Rsitea+1,num)=max(mfaiu,mfaiv);
  
  etau= sinu3-sinu;
  etav= sinv3-sinv;
  metau=max(abs(etau));
  metav=max(abs(etav));
  s1eta(Rsitea+1,num)=max(metau,metav);
  
  dxu= sinu4-sinu;
  dxv= sinv4-sinv;
  mdxu=max(abs(dxu));
  mdxv=max(abs(dxv));
  s1Dx(Rsitea+1,num)=max(mdxu,mdxv);
  
  dzu= sinu5-sinu;
  dzv= sinv5-sinv;
  mdzu=max(abs(dzu));
  mdzv=max(abs(dzv));
  s1Dz(Rsitea+1,num)=max(mdzu,mdzv);
  
  dyu= sinu6-sinu;
  dyv= sinv6-sinv;
  mdyu=max(abs(dyu));
  mdyv=max(abs(dyv));
  s1Dy(Rsitea+1,num)=max(mdyu,mdyv);
  
  ryu= sinu7-sinu;
  ryv= sinv7-sinv;
  mryu=max(abs(ryu));
  mryv=max(abs(ryv));
  s1Ry(Rsitea+1,num)=max(mryu,mryv);
end
end






