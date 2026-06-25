%Load Savvyscan MAT file and generate full MRC file with metadata
%Written by Shahar Seifer, Weizmann Institute of Science, 2021
%Requires @MRCImage library from MatTomo, PEET project: https://bio3d.colorado.edu/imod/matlab.html
clear;
%load([path key]); %Here you can load a custom key with header included according to scan type
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
savename=['F' filename];
savename(length(savename)-2:length(savename))='mrc';
newFilename=[path savename];

mRCImage = MRCImage;%Instentiate MRCImage object
mRCImage.filename=newFilename;
nY=header0.netYsize;
nX=header0.netXsize;
pixelx_Angstrom=header0.pixelXnm*10;
pixely_Angstrom=header0.pixelYnm*10;
for index=0:1000
    name=sprintf('series%d',index);
    if ~exist(name, 'var')
        break;
    end %if exist
end %for index
nZ=index;

vol=zeros(nX,nY,nZ,'int16');

mRCImage.header.nBytesExtended=0;
mRCImage.header.machineStamp=[68 68 0 0]; %for FEI machine, change otherwise 
mRCImage.header.nLabels=0; %no labels included
mRCImage.header.serialEMType=0;
mRCImage.header.nBytesPerSection=0; %extended per section

%min_density=0;
%max_density=0;
%sum_density=0;
for index=0:nZ-1
    name=sprintf('series%d',index);
    series=eval(name);
    if keyload_flag==false
        img=savvyscan_decipher_mat(series,header0,keyload_flag);
    else
        img=savvyscan_decipher_mat(series,header0,keyload_flag,locationsX,locationsY);
    end
    vol(:,:, index+1) = img;
end %for index

mRCImage = setVolume(mRCImage, vol); %enter to mRCImage, do statistics, and fill many details to the header

mRCImage.header.cellDimensionX = nX * pixelx_Angstrom;
mRCImage.header.cellDimensionY = nY * pixely_Angstrom;
mRCImage.header.cellDimensionZ = nZ * pixelx_Angstrom;

oversampling=header0.samples_per_pixel;
pixeltime_us=header0.time_resolution_us*oversampling;
net_time=header0.netXsize*header0.netYsize*pixeltime_us/1000000;
iLabel=1;
mRCImage.header.nLabels=iLabel;
ScanLabel=sprintf('scan time=%g sec, scan type:%s',net_time,header0.scan_type);
mRCImage.header.labels(1,1:length(ScanLabel)) =ScanLabel;

save(mRCImage, newFilename)
close(mRCImage);
disp(sprintf('Saved to file: %s',newFilename));

