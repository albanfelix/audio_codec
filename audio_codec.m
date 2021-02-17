% Alban-Félix Barreteau, M1 CORO SIP
% Mail : alban-felix.barreteau@eleves.ec-nantes.fr
% Matlab R2020a Update 5 (9.8.0.1451342), Student license
% Signal processing, Subject 2 (A basic codec)
%--------------------------------------------------------%

%% Preliminary settings

Nf=4096; %Number of frequencies for FT computations
fmin=70;
fmax=3500;

originalfile='musique.wav'; codedfile = 'handel.wav';
[y, fs] = audioread(originalfile);
%sound(y, fs);
aInfo = audioinfo(originalfile);
nbits = aInfo.BitsPerSample; %nbits=16
npt = aInfo.TotalSamples; %npt=length(y)=total nb of samples=300000
T = aInfo.Duration; %T=6,8027s



%% Original signal analysis


tmin=0; Ts=1/fs; tmax=6.8027; %tmax=T doesn't work while T=6.8027
t=tmin:Ts:tmax;
figure (10); plot(t,y,'b'); 
legend('y_{original}(t)',"location",'southeast')
xlabel('Time t'); ylabel('Signal x(t)');
set(figure(10),'name','Spectrum of the original wav music file')


% %Perform and blot (in blue) a FT on the original signal
% [f, tfx]=FT(y, Nf, fs);
% figure(20); plot(f, abs(tfx),'b')
%     legend('|tfx_{original}(f)|',"location",'southeast')
%     xlabel('f=k/Nf*fs'); ylabel('Magnitude of the FT');
%     set(figure(20),'name','Magnitude spectrum of the FT of y_{original}(t)')
    
    
 %Perform and plot (in blue) a zccsFT on the original signal
[fshift, tfxshift]=zccsFT(y, Nf, fs);
figure(30); plot(fshift, abs(tfxshift),'b')
        legend('|tfxshift_{original}(f)|',"location",'southeast')
        xlabel('f=k/Nf*fs'); ylabel('Magnitude of the zero centered circular shifted FT of y_{original}(t)');
        set(figure(30),'name','Magnitude spectrum of the zccs FT of y_{original}(t)')





%% Codec signal analysis


%Code the signal in the codedfile
[npt, scale]=coder(y, fs, nbits, fmin, fmax, codedfile);

%Decode and plot (in blue grey) the codedfile and sound the codec signal
[ycodec, fscodec, nbits_codec]=decoder(codedfile, fmin, fmax, npt, scale, t); 
%sound(ycodec,fscodec)
figure (40); plot(t,ycodec);
    legend('y_{codec}(t)',"location",'southeast')
    xlabel('Time t'); ylabel('Signal codec y(t)');
    set(figure(40),'name','Spectrum of the wav codec music file')
    
% %Perform and blot (in blue) a FT on the codec signal
% [fshift, tfxshift]=FT(ycodec, Nf, fscodec);
% figure(50); plot(fshift, abs(tfxshift))
%     legend('|tfx_{codec}(f)|',"location",'southeast')
%     xlabel('f=k/Nf*fs'); ylabel('Magnitude of the FT');
%     set(figure(50),'name','Magnitude spectrum of the FT of y_{codec}(t)')

%Perform and plot (in blue grey) a zccsFT on the codec signal
[fshift, tfxshift]=zccsFT(ycodec, Nf, fscodec);
figure(60); plot(fshift, abs(tfxshift))
        legend('|tfxshift_{codec}(f)|',"location",'southeast')
        xlabel('f=k/Nf*fs'); ylabel('Magnitude of the zccs FT of y_{codec}(t)');
        set(figure(60),'name','Magnitude spectrum of the zccs FT of y_{codec}(t)')
        
        
       
%% Functions
        
%Compute a Fourier Transform
function [f, tfx]=FT(y, Nf, fs)
    YNk=fft(y); %fft function computed the Nt values of the signal
    f=(0:length(YNk)-1)*fs/length(YNk); %Vector that corresponds to the array of frequencies for which the transform is calculated
    tfx=1/fs*YNk; %Calculation of the approximation of the FT of x(t)
end


%Compute a zero-centered circular shifted (zccs) Fourier Transform
function [fshift, tfxshift]=zccsFT(y, Nf, fs)
    YNk=fft(y);
    YNkshift=fftshift(YNk);
    fshift=(-length(YNk)/2:length(YNk)/2-1)*(fs/length(YNk));
    tfxshift=1/fs*YNkshift; %Performs a zero-centered circular shift on the transform
end 


%Create a smaller coded file from the original one
function [npt, scale] = coder(y, fs, nbits, fmin, fmax, codedfile)
tfy = fft(y); %Tfd computation,
tfy = tfy(:); %Then ensures column vector

%From frequency mask to matlab index, and mask
npt = length(y);
kmin = round(npt*fmin/fs) + 1;
kmax = round(npt*fmax/fs) + 1;
tfymask = tfy(kmin:kmax);
tfymask = [real(tfymask) imag(tfymask)]; %Real and imaginary part are stored as both channels of a stereo signal

scale = max(max(abs(tfymask)))*1.01;
tfymask = tfymask/scale; % Ensures maximum value lower than 1

audiowrite(codedfile,tfymask,fs,'BitsPerSample',nbits); % Saves with same quantization
end 


%Decode the coded file
function [y, fs, nbits]=decoder(codedfile, fmin, fmax, npt, scale, t) 
[tfymask, fs] = audioread(codedfile); %We read tfymask and the sample frequency fs
aInfo = audioinfo(codedfile); %Read the info of the file
nbits = aInfo.BitsPerSample; %Read the quantization

tfymask=scale*tfymask; %We de-standardize
tfymask=tfymask(:,1)+1i*tfymask(:,2); %Switch back to a mono signal

%From frequency mask to matlab index, and mask
kmin = round(npt*fmin/fs) + 1;
kmax = round(npt*fmax/fs) + 1; 
tfy = zeros(npt,1);
tfy(kmin:kmax)=tfymask;

y=ifft(tfy,'symmetric'); %Inverse discrete Fourier transform
end



%% End
