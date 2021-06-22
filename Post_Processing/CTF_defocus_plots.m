%calculate CTF for iDPC and ADF
%Written by Shahar Seifer, based on formulas found in Lazic & Bosch, Analytical Review of Direct Stem
%Imaging Techniques for Thin Samples, 2017
global N dr dk x y kx ky mx my nx ny;
lambda=0.0025;%nm
defocus_vect=-2000:25:2000;
CS=2000000; %our spherical aberration nm
%CS=2.7*10^6; %in lazic's example, nm
C2=30;%um
NA=0.5*C2/4000;
kBF=NA/lambda; %1/nm
L=1500; %camera length, mm
OPAL_4Q_inner_radius=3.075; %mm ## Change to zero for the case without hole in the DPC detector ##
OPAL_4Q_outer_radius=6.15; %mm
OPAL_E_outer_radius=10; %mm
OPAL_F_outer_radius=14; %mm
k_DPC_max=atan(OPAL_4Q_outer_radius/L)/lambda;
k_DPC_min=atan(OPAL_4Q_inner_radius/L)/lambda;
N=512;%1024;
VN=0:N-1;
Nc=N/2; % Not mean(VN) since FFT is defined with /N
%Calculating Fourier transform from FFT is according to defintions:
%x=n_x/(dk*N), y=n_y/(dk*N), kx=m_x*dk, ky=m_y*dk. So kx*dk_x=1/N.
%FFT2=sum(exp(-i*2*pi*(n_x*m_x+n_y*m_y)/N))
%IFFT2=sum(exp(i*2*pi*(n_x*m_x+n_y*m_y)/N)/N^2)
[my,mx]=meshgrid(VN,VN);
[ny,nx]=meshgrid(VN,VN);
[ky,kx]=meshgrid((VN-Nc)*2*kBF/Nc,(VN-Nc)*2*kBF/Nc); %in 1/nm
%Here the span of k is 4*kBF in total.
dk=2*kBF/Nc;
dr=1/(dk*N);
[y,x]=meshgrid((VN-Nc)/(dk*N),(VN-Nc)/(dk*N)); %in nm
knorm=sqrt(kx.^2+ky.^2);
kangle=(atan(abs(ky)./(abs(kx)+0.00001)).*(-2*(kx<0)+1)+pi*(kx<0)).*(-2*(ky<0)+1);
theta0=pi/4;
if theta0==0
W_DPC_kx=(pi*kBF/sqrt(8))*(knorm>=k_DPC_min & knorm<=k_DPC_max).*(+1*(kangle<pi/4 & kangle>-pi/4)-1*(kangle>3*pi/4 | kangle<-3*pi/4)); 
W_DPC_ky=(pi*kBF/sqrt(8))*(knorm>=k_DPC_min & knorm<=k_DPC_max).*(+1*(kangle<3*pi/4 & kangle>pi/4)-1*(kangle>-3*pi/4 & kangle<-pi/4)); 
elseif theta0==pi/4
W_DPC_kx=(pi*kBF/2)*(knorm>=k_DPC_min & knorm<=k_DPC_max).*(+1*(kx>0)-1*(kx<0)); 
W_DPC_ky=(pi*kBF/2)*(knorm>=k_DPC_min & knorm<=k_DPC_max).*(+1*(ky>0)-1*(ky<0)); 
end

for focus_ind=1:length(defocus_vect)
    defocus=defocus_vect(focus_ind);

    psi_k=(knorm<=kBF).*exp(-1i*(-pi*lambda*defocus*knorm.^2+(pi/2)*CS*lambda^3*knorm.^4));
    psi_in_r=F_r(psi_k);
    %psi_in_k_accurate=F_k_fast(psi_in_r);
    %disp(mean(mean(abs(psi_in_k_accurate-psi_in_k)))); %check accuracy of fast fourier transform; delicate!
    psi_in_k=F_k(psi_in_r); 
    psi_r_gradx=(1/dr)*0.5*([zeros(1,N); psi_in_r(2:N,:)-psi_in_r(1:N-1,:)]+[psi_in_r(2:N,:)-psi_in_r(1:N-1,:); zeros(1,N)]);
    psi_r_grady=(1/dr)*0.5*([zeros(N,1) psi_in_r(:,2:N)-psi_in_r(:,1:N-1)]+[psi_in_r(:,2:N)-psi_in_r(:,1:N-1) zeros(N,1)]);

    CTFSx_sh=-1i*(conj(F_k(psi_in_r.*(F_r(W_DPC_kx.*(invF_k(conj(psi_in_r)))))))...
    -conj(F_k(conj(psi_in_r).*(invF_r(W_DPC_kx.*(F_k(psi_in_r)))))));
    CTFSy_sh=-1i*(conj(F_k(psi_in_r.*(F_r(W_DPC_ky.*(invF_k(conj(psi_in_r)))))))...
    -conj(F_k(conj(psi_in_r).*(invF_r(W_DPC_ky.*(F_k(psi_in_r)))))));
    CTFCx_sh=-1*(conj(F_k(psi_in_r.*(F_r(W_DPC_kx.*(invF_k(conj(psi_in_r)))))))...
    +conj(F_k(conj(psi_in_r).*(invF_r(W_DPC_kx.*(F_k(psi_in_r)))))));
    CTFCy_sh=-1*(conj(F_k(psi_in_r.*(F_r(W_DPC_ky.*(invF_k(conj(psi_in_r)))))))...
    +conj(F_k(conj(psi_in_r).*(invF_r(W_DPC_ky.*(F_k(psi_in_r)))))));

    CTFiS=(kx.*CTFSx_sh+ky.*CTFSy_sh)./(2*pi*1i*knorm.^2);
    CTFiC=(kx.*CTFCx_sh+ky.*CTFCy_sh)./(2*pi*1i*knorm.^2);

    CTFicomCx=conj(F_k(psi_in_r.*conj(psi_r_gradx)-conj(psi_in_r).*psi_r_gradx));
    CTFicomCy=conj(F_k(psi_in_r.*conj(psi_r_grady)-conj(psi_in_r).*psi_r_grady));
    CTFphi2=CTFiC-(kx.*CTFicomCx+ky.*CTFicomCy)./(4*pi^2*knorm.^2);

    CTFicomS=(1/(2*pi))*conj(F_k(psi_in_r.*conj(psi_in_r)));
    CTFphi3=CTFicomS-CTFiS;

    CTFiS_profile=radialAverage((CTFiS),Nc,Nc,Nc);
    CTFphi2_profile=radialAverage((CTFphi2),Nc,Nc,Nc);
    CTFphi3_profile=radialAverage((CTFphi3),Nc,Nc,Nc);
    k_profile=(dk*(0:Nc-1))/kBF;
    figure(1);
    plot([k_profile' k_profile' k_profile'],[real(CTFiS_profile') real(CTFphi2_profile') real(CTFphi3_profile')],'-');
    xlabel('k / k_{BF}');
    ylabel('CTF');
    legend('CTFiS','CTF\phi^{2}','CTF\phi^{3}');

    pointsv=[2 4 6 10];%[7 28 28*4];
    kprof=k_profile(pointsv);
    if focus_ind==1
        CTF1=real(CTFiS_profile(pointsv)');
        CTF2=real(CTFphi2_profile(pointsv)');
        CTF3=real(CTFphi3_profile(pointsv)');
    else
        CTF1=[CTF1 real(CTFiS_profile(pointsv)')];
        CTF2=[CTF2 real(CTFphi2_profile(pointsv)')];
        CTF3=[CTF3 real(CTFphi3_profile(pointsv)')];
    end
    
    %Calculate CTF of any ADF detector (without constant prefactor)
    psi_k_ADF=(knorm<=kBF).*exp(-1i*(-pi*lambda*defocus*knorm.^2+(pi/2)*CS*lambda^3*knorm.^4));
    psi_in_r_ADF=F_r(psi_k_ADF);
    Fk_psi_r_square=F_k(psi_in_r_ADF.*conj(psi_in_r_ADF));
    CTF_ADF_profile=radialAverage(real(Fk_psi_r_square),Nc,Nc,Nc);
    if focus_ind==1
        CTF_ADF=real(CTF_ADF_profile(pointsv)');
    else
        CTF_ADF=[CTF_ADF real(CTF_ADF_profile(pointsv)')];
    end
    
    
end %for focus_ind
text{1}=sprintf('k/k_{BF}=%g',0.001*round(1000*kprof(1)/kBF));
text{2}=sprintf('k/k_{BF}=%g',0.001*round(1000*kprof(2)/kBF));
text{3}=sprintf('k/k_{BF}=%g',0.001*round(1000*kprof(3)/kBF));
close all;

figure(2)
plot(defocus_vect,CTF1,'-');
xlabel('defocus [nm]');
ylabel('CTFiS');
legend(text);

figure(3)
plot(defocus_vect,CTF2,'-');
xlabel('defocus [nm]');
ylabel('CTF\phi^{2}');
legend(text);

figure(4)
plot(defocus_vect,CTF3,'-');
xlabel('defocus [nm]');
ylabel('CTF\phi^{3}');
legend(text);

figure(4)
plot(defocus_vect,CTF_ADF,'-');
xlabel('defocus [nm]');
ylabel('CTF-ADF');
legend(text);


function result=F_k_accurate(array)
%for testing by comparing F_k to this result
    global N dr x y kx ky;
    result=zeros(N,N);
    for ind1=1:N
        for ind2=1:N
            result(ind1,ind2)=dr^2*sum(sum(array.*exp(-1i*2*pi*(x*kx(ind1,ind2)+y*ky(ind1,ind2)))));
        end
    end
end
function result=F_k(array)
    global dr;
    result=dr^2*ifftshift(fft2(fftshift(array)));
    %Note that using fftshift is equivalent to :
    %result=dr^2*(exp(1i*pi*(mx+my-N/2)).*fft2(exp(1i*pi*(nx+ny)).*(array)));
end
function result=F_r(array)
    global dk;
    result=dk^2*ifftshift(fft2(fftshift(array)));
end
function result=invF_k(array)
    global dr N;
    result=dr^2*N^2*ifftshift(ifft2(fftshift(array)));
    %note that ifft2 is defnied with excessive 1/N^2 that we replace with integration differential 
end
function result=invF_r(array)
    global dk N;
    result=dk^2*N^2*ifftshift(ifft2(fftshift(array)));
end

function profile = radialAverage(IMG, cx, cy, w)
    % computes the radial average of the image IMG around the cx,cy point
    % w is the size of vector of radii starting from zero
    [a,b] = size(IMG);
    [Y, X] = meshgrid( (1:a)-cx, (1:b)-cy);
    R = sqrt(X.^2 + Y.^2);
    profile = [];
    for i = 1:w % radius of the circle
        mask = (i-1<R & R<i+1); % smooth 1 px around the radius
        values = (1-abs(R(mask)-i)) .* double(IMG(mask)); % smooth based on distance to ring
        % values = IMG(mask); % without smooth
        profile(end+1) = mean( values(:) );
    end
end
