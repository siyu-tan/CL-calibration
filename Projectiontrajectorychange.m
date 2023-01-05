%This code is used to analyze the change of the projection trajectory caused by the deviation of 7 parameters.
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
Sy = -1000;  %SDD=1000£¨mm£©
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

nummax=30;

sinu1=zeros(7,360); 
sinv1=zeros(7,360); 

sinu2=zeros(7,360); 
sinv2=zeros(7,360); 

sinu3=zeros(7,360); 
sinv3=zeros(7,360); 

sinu4=zeros(7,360); 
sinv4=zeros(7,360); 

sinu5=zeros(7,360); 
sinv5=zeros(7,360); 

sinu6=zeros(7,360); 
sinv6=zeros(7,360); 

sinu7=zeros(7,360); 
sinv7=zeros(7,360); 

sinu=zeros(1,360);     
sinv=zeros(1,360);    

for num=0:5:nummax
    
    site=num*pi/180;
    fai=num*pi/180;
    eta=num*pi/180;
    dx=num*detVoxel;
    dz=num*detVoxel;
    dy=num*detVoxel;
    ry=num*detVoxel;
  
for sita = 1:1:360 
    
                l = (a-obj_lenth_center)*objVoxel;  %z
                w = (b-obj_width_center)*objVoxel;  %x
                h = (c-obj_height_center)*objVoxel; %y
          
                x0 = w*cos(sita*d2r)+h*sin(sita*d2r);  
                y0 = -w*sin(sita*d2r)+h*cos(sita*d2r);
                objXYZ = [x0;y0;l];  
              
                x1 = objXYZ(1)+Ox;
                y1 = objXYZ(2)+Oy;
                z1 = objXYZ(3)+Oz;
                
                dy0 = -Sy/(y1-Sy);   %magnification,SDD£ºSy/SOD:y1-Sy
              
                du = dy0*(x1-Sx)+Sx; 
                dv = dy0*(z1-Sz)+Sz; 
               
                sinu(sita+1)=du;      
                sinv(sita+1)=dv;     
              
                du1=du*(1-dv*sin(site)/(SDD*cos(site)+dv*sin(site)));
                dv1=dv*SDD/(SDD*cos(site)+dv*sin(site));
       
                sinu1(num/5+1,sita)=du1;     
                sinv1(num/5+1,sita)=dv1;       
     
                du2=du*SDD/(SDD*cos(fai)+du*sin(fai));
                dv2=dv*(1-du*sin(fai)/(SDD*cos(fai)+du*sin(fai)));
              
                sinu2(num/5+1,sita)=du2;      
                sinv2(num/5+1,sita)=dv2;      
   
               
                du3=du*cos(eta)+dv*sin(eta);
                dv3=-du*sin(eta)+dv*cos(eta);
              
                sinu3(num/5+1,sita)=du3;     
                sinv3(num/5+1,sita)=dv3;     
            
                du4 = du-dx; 
                dv4 = dv; 
                sinu4(num/5+1,sita)=du4;     
                sinv4(num/5+1,sita)=dv4; 

                du5 = du; 
                dv5 = dv-dz; 
                sinu5(num/5+1,sita)=du5;     
                sinv5(num/5+1,sita)=dv5; 
               
                dy6 = (-Sy+dy)/(y1-Sy);  
                du6 = dy6*(x1-Sx)+Sx; 
                dv6 = dy6*(z1-Sz)+Sz; 
                sinu6(num/5+1,sita)=du6;     
                sinv6(num/5+1,sita)=dv6; 
              
                dy7 = -Sy/(y1-Sy+ry);  
                du7 = dy7*(x1-Sx)+Sx; 
                dv7 = dy7*(z1-Sz)+Sz; 
                sinu7(num/5+1,sita)=du7;     
                sinv7(num/5+1,sita)=dv7;       
end 
end






