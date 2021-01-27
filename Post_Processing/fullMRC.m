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

mRCImage = MRCImage;%Initiate MRCImage object
mRCImage.filename=newFilename;
nY=header0.netYsize;
nX=header0.netXsize;
pixelx_Angstrom=header0.pixelXnm*10;
pixely_Angstrom=header0.pixelYnm*10;
mRCImage.header.mX=1;
mRCImage.header.mY=1;
mRCImage.header.mZ=1;
mRCImage.header.nX=nX;
mRCImage.header.nY=nY;
mRCImage.header.mode=1; %signed 16 bits
for index=0:1000
    name=sprintf('series%d',index);
    if ~exist(name, 'var')
        break;
    end %if exist
end %for index
nZ=index;
mRCImage.header.nZ=nZ;
mRCImage.header.cellDimensionX = nX * pixelx_Angstrom;
mRCImage.header.cellDimensionY = nY * pixely_Angstrom;
mRCImage.header.cellDimensionZ = nZ * pixelx_Angstrom;
mRCImage.flgVolume=true;
mRCImage.header.nXStart=0;
mRCImage.header.nYStart=0;
mRCImage.header.nZStart=0;
mRCImage.header.nBytesExtended=0;
mRCImage.header.cellAngleX=90;
mRCImage.header.cellAngleY=90;
mRCImage.header.cellAngleZ=90;

mRCImage.header.densityRMS=-1;
mRCImage.header.spaceGroup=0;
mRCImage.header.nSymmetryBytes=0;
mRCImage.header.nBytesExtended=0;
mRCImage.header.creatorID=0;

mRCImage.header.xOrigin=0;
mRCImage.header.yOrigin=0;
mRCImage.header.zOrigin=0;
mRCImage.header.imodStamp = 0;%1146047817
mRCImage.header.imodFlags=0;%13
mRCImage.header.map='MAP';
mRCImage.header.machineStamp=[68 68 0 0]; %for FEI machine, change otherwise 
mRCImage.header.nLabels=0; %no labels included
mRCImage.header.serialEMType=0;
mRCImage.header.nBytesPerSection=0; %extended per section

min_density=0;
max_density=0;
sum_density=0;
for index=0:nZ-1
    name=sprintf('series%d',index);
    series=eval(name);
    if keyload_flag==false
        img=savvyscan_decipher_mat(series,header0,keyload_flag);
    else
        img=savvyscan_decipher_mat(series,header0,keyload_flag,locationsX,locationsY);
    end
    mRCImage.volume(:,:, index+1) = img';
    maxd=max(max(img));
    mind=min(min(img));
    if mind<min_density || index==0
        min_density=mind;
    end
    if maxd>max_density || index==0
        max_density=maxd;
    end
    sum_density=sum_density+mean(mean(img));
    %[nX,nY] = size(img);
    %putImage(mRCImage, img, index+1) % NOGO: Only suitable for replacing exisiting file
end %for index
mRCImage.header.meanDensity=sum_density/nZ;
mRCImage.header.minDensity=min_density;
mRCImage.header.maxDensity=max_density;

oversampling=header0.samples_per_pixel;
pixeltime_us=header0.time_resolution_us*oversampling;
net_time=header0.netXsize*header0.netYsize*pixeltime_us/1000000;
iLabel=1;
mRCImage.header.nLabels=iLabel;
ScanLabel=sprintf('scan time=%g sec, scan type:%s',net_time,header0.scan_type);
mRCImage.header.labels(1,1:length(ScanLabel)) =ScanLabel;

save(mRCImage, newFilename)
close(mRCImage);
disp(sprintf('Saved to file: %s',newFilename);

