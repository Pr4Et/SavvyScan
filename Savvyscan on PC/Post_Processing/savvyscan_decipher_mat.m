function showpic=savvyscan_decipher_mat(series,header0, keyload_flag,locationsX,locationsY)
%Deciphering Savvyscan mat files according to recognized scan patterns
%Written by Shahar Seifer, Weizmann Institute of Science

correction=1;%with correction for delay
show_margins=0;
extra_correct=1;

len=length(series);
height=header0.netYsize;
width=header0.netXsize;
%oversampling=10;
oversampling=header0.samples_per_pixel;
%pixeltime_us=1000000*integration_time/(width*height);%declared and also actual (sampling time per pixel)
pixeltime_us=header0.time_resolution_us*oversampling;
%scXsize=floor(width*sqrt(2));
%scYsize=floor(height*sqrt(2));
scXsize=header0.fullXsize;
scYsize=header0.fullYsize;
net_time=header0.netXsize*header0.netYsize*pixeltime_us/1000000;
%disp(sprintf('net time=%g sec',net_time));


if keyload_flag==false

    locationsX=zeros(len,1);
    locationsY=zeros(len,1);
    locationsX_rot=zeros(len,1);
    locationsY_rot=zeros(len,1);
    hide=zeros(len,1);
    rnbs=floor((scXsize/2));
    marginx=floor(0.5*(scXsize-width));
    marginy=floor(0.5*(scYsize-height));
    

    %binning=input("binning=");%4
    delay_lowfreq=0.000265;%in sec 
    chi=0;%0.025*delay_lowfreq^2;
    delay=delay_lowfreq;
    dt=0.000001*pixeltime_us/oversampling;
    flyback_us=40.0;
    angle_rot=0*pi/180;
    out_amplitude=2000;
    
    if strcmp(header0.scan_type,'savvy-spiral1')
        if show_margins
            marginx=0; %!!! for demo of showing margins
            marginy=0;
        end
        locationsX0_rot=(scXsize-marginx);%assume the beam starts from parking position
        locationsY0_rot=(scYsize-marginy);
        counter=1;
        for scanrvar=rnbs:-1:0  %radius of circle, enclosing a square widthXwidth
            if scanrvar~=0
                phnmbs=(floor(2*pi*scanrvar))*oversampling; %number of phase points: so the arc steps are 1pixel/oversampling long
            else
                phnmbs=1*oversampling;
            end
            for scanphvar=0:phnmbs-1
                scanph=scanphvar*2*pi/phnmbs;
                xkey_request=(scanrvar*sin(scanph)+rnbs-marginx);
                ykey_request=(scanrvar*cos(scanph)+rnbs-marginy);
                xkey_request_rot=cos(angle_rot)*xkey_request+sin(angle_rot)*ykey_request;
                ykey_request_rot=cos(angle_rot)*ykey_request-sin(angle_rot)*xkey_request;
                delayx=delay_lowfreq*(1-extra_correct*0.17-0.00000010*(out_amplitude/width)/dt);%something to do with trajectory velocity
                delayy=delay_lowfreq*(1-0.00000010*(out_amplitude/width)/dt);%something to do with trajectory velocity
                if counter>1 && correction==1 
                    locationsX_rot(counter)=(xkey_request_rot+(delayx/dt)*locationsX_rot(counter-1))/(1+delayx/dt);
                    locationsY_rot(counter)=(ykey_request_rot+(delayy/dt)*locationsY_rot(counter-1))/(1+delayy/dt);
                elseif counter==1 && correction==1
                    locationsX_rot(counter)=(xkey_request_rot+(delayx/dt)*locationsX0_rot)/(1+delayx/dt);%assume the beam starts from parking position
                    locationsY_rot(counter)=(ykey_request_rot+(delayy/dt)*locationsY0_rot)/(1+delayy/dt);
                else
                    locationsX_rot(counter)=xkey_request_rot;
                    locationsY_rot(counter)=ykey_request_rot;
                end
                locationsX(counter)=cos(-angle_rot)*locationsX_rot(counter)+sin(-angle_rot)*locationsY_rot(counter);
                locationsY(counter)=cos(-angle_rot)*locationsY_rot(counter)-sin(-angle_rot)*locationsX_rot(counter);
              
                counter=counter+1;
            end
        end %for scanvar
        if show_margins==1
            width=scXsize;% !!!  only after! for demo of showin margins
            height=scYsize;% !!!  only after! for demo of showin margins
        end
    elseif strcmp(header0.scan_type,'savvy-raster1')
       counter=1;
       %for scanyvar=scYsize-1:-1:0
        if show_margins
            marginx=0; %!!! for demo of showing margins
            marginy=0;
        end
        locationsX0_rot=(scXsize-marginx);%assume the beam starts from parking position
        locationsY0_rot=(scYsize-marginy);
        for scanyvar=0:scYsize-1
 			for scanxvar=0:(scXsize*oversampling)-1
				xkey_request=((scanxvar/oversampling)-marginx);
				ykey_request=scanyvar-marginy;
                xkey_request_rot=cos(angle_rot)*xkey_request+sin(angle_rot)*ykey_request;
                ykey_request_rot=cos(angle_rot)*ykey_request-sin(angle_rot)*xkey_request;
                delayx=delay_lowfreq*(1-extra_correct*0.17-0.00000010*(out_amplitude/width)/dt);%something to do with trajectory velocity
                delayy=delay_lowfreq*(1-0.00000010*(out_amplitude/width)/dt);%something to do with trajectory velocity
                if counter>1 && correction==1 
                    locationsX_rot(counter)=(xkey_request_rot+(delayx/dt)*locationsX_rot(counter-1))/(1+delayx/dt);
                    locationsY_rot(counter)=(ykey_request_rot+(delayy/dt)*locationsY_rot(counter-1))/(1+delayy/dt);
                elseif counter==1 && correction==1
                    locationsX_rot(counter)=(xkey_request_rot+(delayx/dt)*locationsX0_rot)/(1+delayx/dt);%assume the beam starts from parking position
                    locationsY_rot(counter)=(ykey_request_rot+(delayy/dt)*locationsY0_rot)/(1+delayy/dt);
                else
                    locationsX_rot(counter)=xkey_request_rot;
                    locationsY_rot(counter)=ykey_request_rot;
                end
                locationsX(counter)=cos(-angle_rot)*locationsX_rot(counter)+sin(-angle_rot)*locationsY_rot(counter);
                locationsY(counter)=cos(-angle_rot)*locationsY_rot(counter)-sin(-angle_rot)*locationsX_rot(counter);
              
                counter=counter+1;
            end % for scanxvar
            
            for scanxvar=floor(flyback_us/pixeltime_us)*oversampling-1:-1:0
				xkey_request=scXsize*scanxvar/(floor(flyback_us/pixeltime_us)*oversampling)-marginx;
				ykey_request=scanyvar-marginy;
                locationsX(counter)=xkey_request;
                locationsY(counter)=ykey_request;
                hide(counter)=1;
                counter=counter+1;
            end
            if show_margins==1
                 width=scXsize;% !!!  only after! for demo of showin margins
                 height=scYsize;% !!!  only after! for demo of showin margins
            end

         end % for scanyvar
         
    elseif strcmp(header0.scan_type,'savvy-mandala1')% case of 2/3 AspectRatio
        delay_lowfreq=0.000265*1.51;
        if show_margins
            marginx=0; %!!! for demo of showing margins
        end
        locationsX0_rot=(scXsize-marginx);%assume the beam starts from parking position
        locationsY0_rot=(scYsize-marginy);
        rnbs=width;%//number of circling times,case of 2/3 AspectRatio
        scanr=floor(height/2); %//radius of the circle
        phnmbs=floor(2.0*pi*scanr)*oversampling; %//number of phase points in one round: so the arc steps are 1pixel/oversampling long
        counter=1;
        for scanphvar=0:phnmbs*rnbs-1
         	scanph=scanphvar*2.0*pi/phnmbs; 
			scanphcor=scanph - floor(scanph/(2.0*pi))*2.0*pi; %//modulus 2PI
            xdrift=-width+scanph/pi; %//drift smoothly, case of 2/3 AspectRatio
            xkey_request=(-scanr*sin(scanphcor)-xdrift+scXsize/2.0-marginx);
			ykey_request=(scanr*cos(scanphcor)+scYsize/2.0-marginy);
            xkey_request_rot=cos(angle_rot)*xkey_request+sin(angle_rot)*ykey_request;
            ykey_request_rot=cos(angle_rot)*ykey_request-sin(angle_rot)*xkey_request;
            delayx=delay_lowfreq*(1-extra_correct*0.17-0.00000010*(out_amplitude/width)/dt);%something to do with trajectory velocity
            delayy=delay_lowfreq*(1-0.00000010*(out_amplitude/width)/dt);%something to do with trajectory velocity
            if counter>1 && correction==1 
                locationsX_rot(counter)=(xkey_request_rot+(delayx/dt)*locationsX_rot(counter-1))/(1+delayx/dt);
                locationsY_rot(counter)=(ykey_request_rot+(delayy/dt)*locationsY_rot(counter-1))/(1+delayy/dt);
            elseif counter==1 && correction==1
                locationsX_rot(counter)=(xkey_request_rot+(delayx/dt)*locationsX0_rot)/(1+delayx/dt);%assume the beam starts from parking position
                locationsY_rot(counter)=(ykey_request_rot+(delayy/dt)*locationsY0_rot)/(1+delayy/dt);
            else
                locationsX_rot(counter)=xkey_request_rot;
                locationsY_rot(counter)=ykey_request_rot;
            end
            locationsX(counter)=cos(-angle_rot)*locationsX_rot(counter)+sin(-angle_rot)*locationsY_rot(counter);
            locationsY(counter)=cos(-angle_rot)*locationsY_rot(counter)-sin(-angle_rot)*locationsX_rot(counter);
              
            counter=counter+1;
        end %for scanphvar
        if show_margins==1
            width=scXsize;% !!!  only after! for demo of showin margins
        end
    elseif strcmp(header0.scan_type,'savvy-mandala2') %Normal case
        if show_margins
            marginx=0; %!!! for demo of showing margins
        end
        locationsX0_rot=(scXsize-marginx);%assume the beam starts from parking position
        locationsY0_rot=(scYsize-marginy);
        rnbs=floor(width/2);%//number of circling times, so the right hand of the circle interlaces with left hand, jumping 2 pixels per cycle, over 2 width distance
        scanr=floor(height/2); %//radius of the circle
        phnmbs=floor(2.0*pi*scanr)*oversampling; %//number of phase points in one round: so the arc steps are 1pixel/oversampling long
        counter=1;
        for scanphvar=0:phnmbs*rnbs-1
         	scanph=scanphvar*2.0*pi/phnmbs; 
			scanphcor=scanph - floor(scanph/(2.0*pi))*2.0*pi; %//modulus 2PI
            xdrift=-floor(width/2)+scanph/pi; %//drift smoothly one pixel per half a round
            xkey_request=(-scanr*sin(scanphcor)-xdrift+scXsize/2.0-marginx);
			ykey_request=(scanr*cos(scanphcor)+scYsize/2.0-marginy);
            xkey_request_rot=cos(angle_rot)*xkey_request+sin(angle_rot)*ykey_request;
            ykey_request_rot=cos(angle_rot)*ykey_request-sin(angle_rot)*xkey_request;
            delayx=delay_lowfreq*(1-extra_correct*0.17-0.00000010*(out_amplitude/width)/dt);%something to do with trajectory velocity
            delayy=delay_lowfreq*(1-0.00000010*(out_amplitude/width)/dt);%something to do with trajectory velocity
            if counter>1 && correction==1 
                locationsX_rot(counter)=(xkey_request_rot+(delayx/dt)*locationsX_rot(counter-1))/(1+delayx/dt);
                locationsY_rot(counter)=(ykey_request_rot+(delayy/dt)*locationsY_rot(counter-1))/(1+delayy/dt);
            elseif counter==1 && correction==1
                locationsX_rot(counter)=(xkey_request_rot+(delayx/dt)*locationsX0_rot)/(1+delayx/dt);%assume the beam starts from parking position
                locationsY_rot(counter)=(ykey_request_rot+(delayy/dt)*locationsY0_rot)/(1+delayy/dt);
            else
                locationsX_rot(counter)=xkey_request_rot;
                locationsY_rot(counter)=ykey_request_rot;
            end
            locationsX(counter)=cos(-angle_rot)*locationsX_rot(counter)+sin(-angle_rot)*locationsY_rot(counter);
            locationsY(counter)=cos(-angle_rot)*locationsY_rot(counter)-sin(-angle_rot)*locationsX_rot(counter);
              
            counter=counter+1;
        end %for scanphvar
        if show_margins==1
            width=scXsize;% !!!  only after! for demo of showin margins
        end
    else
        disp('Scan type key not supplied');
    end
    
    
end %of if keyload_flag==false
    
    
fillresults=zeros(height*width,1);
fillcounter=zeros(height*width,1);
		
for j = 1:len
    for delx=0:1
        for dely=0:1
            if hide(j)==0
                xk=floor(locationsX(j))+delx;
                yk=floor(locationsY(j))+dely;
                factorx=1-abs(locationsX(j)-xk)/1;
                factory=1-abs(locationsY(j)-yk)/1;
                if xk>=1 && xk<=width && yk>=1 && yk<=height
                            index=xk+(yk-1)*width;
                            fillresults(index)=fillresults(index)+series(j)*factorx*factory;
                            fillcounter(index)=fillcounter(index)+1*factorx*factory;
                end
            end
        end
    end
end
fillresults=round(fillresults./(fillcounter+0.1*(fillcounter==0)));
meanval=median(fillresults(fillresults>0));
imtable=zeros(width,height);
for j = 1:height*width
            x=mod(j-1,width)+1;
            y=floor((j-1)/width)+1;
			if fillcounter(j)>0
			   imtable(x,y)=fillresults(j); %floor(fillresults(j)/fillcounter(j));
            else
               imtable(x,y)=meanval;
            end
end
Nshades=1024;
mapvector=linspace(0,1,Nshades)';
cmap=zeros(Nshades,3);
for loop=1:3
    cmap(:,loop)=mapvector;
end
imtable(isnan(imtable))=0;
showpic=imtable; %RESULT RETURNED TO CALLER, NOT NORMALLIZED
%showpic=balance_pic(imtable,Nshades);% RETURN NORMALIZED RESULTS BETWEEN 0 AND 1024
%Remove remarks to show images: 
%showpic2=balance_pic(imtable,Nshades); 
%figure(1)
%imshow(showpic2',cmap);
