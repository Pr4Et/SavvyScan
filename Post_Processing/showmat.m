[filename,path] = uigetfile('CH*.mat','Fetch savvyscan file CH*.mat');
%load scan file to header0, series0,1,2...
load([path filename],'-mat');
display([path filename]);
keyload_flag=false;
if exist('header','var') %check if key loaded
    if header.scan_type==header0.scan_type && header.netXsize==header0.netXsize
      keyload_flag=true;
    end 
end
if keyload_flag==false
   img=savvyscan_decipher_mat(series0,header0,keyload_flag);
else
   img=savvyscan_decipher_mat(series0,header0,keyload_flag,locationsX,locationsY);
end

Nshades=1024;
mapvector=linspace(0,1,Nshades)';
cmap=zeros(Nshades,3);
for loop=1:3
    cmap(:,loop)=mapvector;
end
showpic2=balance_pic(img,Nshades);
figure(1)
imshow(showpic2',cmap);