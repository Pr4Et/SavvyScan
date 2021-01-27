%Calibration of OPAL detector, producing montage of segments and calculating mean and STD of active detection regions
%Mat files should be acquired without sample, in direct mode (not diffraction mode) with the beam scanning on detector area.
%Written by Shahar Seifer, Weizmann Institute of Science, 2021
[filename,path] = uigetfile('E:\Shadowdata\intensity_calib_no_haadf\CH*.mat','Fetch savvyscan file CH1.mat');
display([path filename]);
keyload_flag=false;
%load scan file to header0, series0,1,2...
load([path 'CH1.mat'],'-mat');
img1=-savvyscan_decipher_mat(series0,header0,keyload_flag);
load([path 'CH2.mat'],'-mat');
img2=-savvyscan_decipher_mat(series0,header0,keyload_flag);
load([path 'CH3.mat'],'-mat');
img3=-savvyscan_decipher_mat(series0,header0,keyload_flag);
load([path 'CH4.mat'],'-mat');
img4=-savvyscan_decipher_mat(series0,header0,keyload_flag);
load([path 'CH5.mat'],'-mat');
img5=-savvyscan_decipher_mat(series0,header0,keyload_flag);
load([path 'CH6.mat'],'-mat');
img6=-savvyscan_decipher_mat(series0,header0,keyload_flag);
%x-left to right, y-bottom up
R1=16;
R2=55;
R3=87;
R4=118;
[yN,xN]=size(img1);
xc=xN/2+5;
yc=yN/2-1;
ymesh=ones(yN,1)*(yN-0.5:-1:0.5);
xmesh=(0.5:xN-0.5)'*ones(1,xN);
rmesh=sqrt((xmesh-xc).^2+(ymesh-yc).^2);
del_alpha=pi/2;
alpha_mesh=atan((ymesh-yc)./(xmesh-xc+0.000001))+pi*((xmesh-xc)<0);
alpha_mesh(alpha_mesh>0)=alpha_mesh(alpha_mesh>0)-2*pi;
alpha0=-2*pi/180;
img=img1; %use first channel as background. CH1 is at 4st quadrature, 
img(rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-0*del_alpha & alpha_mesh>alpha0-del_alpha)=img1(rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-0*del_alpha & alpha_mesh>alpha0-1*del_alpha);
img(rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-1*del_alpha & alpha_mesh>alpha0-2*del_alpha)=img4(rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-1*del_alpha & alpha_mesh>alpha0-2*del_alpha);
img(rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-2*del_alpha & alpha_mesh>alpha0-3*del_alpha)=img3(rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-2*del_alpha & alpha_mesh>alpha0-3*del_alpha);
img(rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-3*del_alpha & alpha_mesh>alpha0-4*del_alpha)=img2(rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-3*del_alpha & alpha_mesh>alpha0-4*del_alpha);
img(rmesh>R2 & rmesh<R3)=img5(rmesh>R2 & rmesh<R3);
img(rmesh>R3 & rmesh<R4)=img6(rmesh>R3 & rmesh<R4);

Nshades=1024;
mapvector=linspace(0,1,Nshades)';
cmap=zeros(Nshades,3);
for loop=1:3
    cmap(:,loop)=mapvector;
end

imgcrop=img(390:644,388:638);
showpic=balance_pic(imgcrop,Nshades);
figure(1);
imshow(showpic',cmap);
%line ([x1,x2], [y1,y2])
xcm=xc-390+1.5;
ycm=yc-390+1.5;
line([xcm+R1*cos(alpha0-0*del_alpha) xcm+R2*cos(alpha0-0*del_alpha)],[ycm-R1*sin(alpha0-0*del_alpha) ycm-R2*sin(alpha0-0*del_alpha)],'Color','red','LineWidth',2);
line([xcm+R1*cos(alpha0-1*del_alpha) xcm+R2*cos(alpha0-1*del_alpha)],[ycm-R1*sin(alpha0-1*del_alpha) ycm-R2*sin(alpha0-1*del_alpha)],'Color','red','LineWidth',2);
line([xcm+R1*cos(alpha0-2*del_alpha) xcm+R2*cos(alpha0-2*del_alpha)],[ycm-R1*sin(alpha0-2*del_alpha) ycm-R2*sin(alpha0-2*del_alpha)],'Color','red','LineWidth',2);
line([xcm+R1*cos(alpha0-3*del_alpha) xcm+R2*cos(alpha0-3*del_alpha)],[ycm-R1*sin(alpha0-3*del_alpha) ycm-R2*sin(alpha0-3*del_alpha)],'Color','red','LineWidth',2);
viscircles([xcm,ycm],R1,'Color','red','LineWidth',3);
viscircles([xcm,ycm],R2,'Color','red','LineWidth',3);
viscircles([xcm,ycm],R3,'Color','red','LineWidth',3);
viscircles([xcm,ycm],R3,'Color','red','LineWidth',3);
viscircles([xcm,ycm],R4,'Color','red','LineWidth',3);

figure(2)
threshold =0;
maxs=0;
series=img1(img1>threshold & rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-0*del_alpha & alpha_mesh>alpha0-1*del_alpha);
vect_temp=0:50:(max(series)*1.2);
result=histogram(series,vect_temp);
histvalues1=[0 result.Values];
series=img2(img2>threshold & rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-3*del_alpha & alpha_mesh>alpha0-4*del_alpha);
result=histogram(series,vect_temp);
histvalues2=[0 result.Values];
series=img3(img3>threshold & rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-2*del_alpha & alpha_mesh>alpha0-3*del_alpha);
result=histogram(series,vect_temp);
histvalues3=[0 result.Values];
series=img4(img4>threshold & rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-1*del_alpha & alpha_mesh>alpha0-2*del_alpha);
result=histogram(series,vect_temp);
histvalues4=[0 result.Values];
series=img5(img5>threshold & rmesh>R2 & rmesh<R3);
result=histogram(series,vect_temp);
histvalues5=[0 result.Values];
series=img6(img6>threshold & rmesh>R3 & rmesh<R4);
result=histogram(series,vect_temp);
histvalues6=[0 result.Values];
plot(vect_temp,histvalues1/sum(histvalues1),vect_temp,histvalues2/sum(histvalues2),vect_temp,histvalues3/sum(histvalues3),vect_temp,histvalues4/sum(histvalues4),vect_temp,histvalues5/sum(histvalues5),vect_temp,histvalues6/sum(histvalues6));
xlabel('Intensity');
ylabel('Normalized Count');
legend('Ch1','Ch2','Ch3','Ch4','Ch5','Ch6');

figure(3)
series=reshape(imgcrop,[length(imgcrop(1,:))*length(imgcrop(:,1)) 1]);
thresholdH=0.01*max(series);
vect_temp=0:200:max(series);
result=histogram(series,vect_temp);
histvalues=[0 result.Values];
p0=mean(vect_temp(histvalues==max(histvalues(vect_temp>0.1*max(series))) & vect_temp>0.1*max(series)));
p1=max(vect_temp(histvalues<thresholdH & vect_temp<p0));
p2=min(vect_temp(histvalues<thresholdH & vect_temp>p0));
result2=histogram(series,p1:50:p2);
threshold=p1;
ch1_mean=mean(mean(img1(img1>threshold & rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-0*del_alpha & alpha_mesh>alpha0-1*del_alpha)));
ch2_mean=mean(mean(img2(img2>threshold & rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-3*del_alpha & alpha_mesh>alpha0-4*del_alpha)));
ch3_mean=mean(mean(img3(img3>threshold & rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-2*del_alpha & alpha_mesh>alpha0-3*del_alpha)));
ch4_mean=mean(mean(img4(img4>threshold & rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-1*del_alpha & alpha_mesh>alpha0-2*del_alpha)));
ch5_mean=mean(mean(img5(img5>threshold & rmesh>R2 & rmesh<R3)));
ch6_mean=mean(mean(img6(img6>threshold & rmesh>R3 & rmesh<R4)));
disp(sprintf('ch1=%d, ch2=%d, ch3=%d, ch4=%d, ch5=%d, ch6=%d',ch1_mean,ch2_mean,ch3_mean,ch4_mean,ch5_mean,ch6_mean));

ch1_std=std((img1(img1>threshold & rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-0*del_alpha & alpha_mesh>alpha0-1*del_alpha)));
ch2_std=std((img2(img2>threshold & rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-3*del_alpha & alpha_mesh>alpha0-4*del_alpha)));
ch3_std=std((img3(img3>threshold & rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-2*del_alpha & alpha_mesh>alpha0-3*del_alpha)));
ch4_std=std((img4(img4>threshold & rmesh>R1 & rmesh<R2 & alpha_mesh<=alpha0-1*del_alpha & alpha_mesh>alpha0-2*del_alpha)));
ch5_std=std((img5(img5>threshold & rmesh>R2 & rmesh<R3)));
ch6_std=std((img6(img6>threshold & rmesh>R3 & rmesh<R4)));
disp(sprintf('Dch1=%d, Dch2=%d, Dch3=%d, Dch4=%d, Dch5=%d, Dch6=%d',ch1_std,ch2_std,ch3_std,ch4_std,ch5_std,ch6_std));


