 %%discard unwanted samples in audio file OFDM QPSK
 close all
 x=txout;
 %x=x(:,1);
 if sum(x(1,:))==1
 x=x';
 end
[xr,locr]=findpeaks(x(1:end));
meanp = mean(xr);
[row,colmn]=find(xr>(meanp));
%x=x(locr(colmn(1)):locr(colmn(end)));
x=x(locr(colmn(1)):end);
x=x(8001:end);
plot(x)


 %%

 xdet=[1 1];
 Nchar = 1024;  % make power of 2
bits = Nchar * 8;
 Nguard=512;
N2 = (bits+1)/3;
Nfft = N2 + 1024; %for qpsk
Nfft = (2*Nfft)  - 1 ; 

%% x = x(:,1);
%x = x.';

fswav=fs;

if fswav ~= fs
    n = gcd(fs,fswav);
    p = fs / n;
    q = fswav / n;
    x = resample(x,p,q);   % resample file to expected sample rate
end


sigtime = Nguard + 2*Nfft + 2*Nfft + 2*Nguard;
 plotcrap=1;
 
[maxv, maxi] = max(abs(xdet));
xstart = maxi ;
xend = xstart + sigtime ;

if plotcrap
    figure(101)
    xplot = abs(x) / max(abs(x));
    xplot = 20 * log10(xplot);
    yplot = length(xdet);
    yplot = zeros(1,yplot);
    yplot(maxi) = -40;
    plot(xplot)
    hold on
    plot(yplot)
    hold off
end

%%
x = x(xstart:xend);
a = Nguard  + round(0.8*Nfft) + 1 ;
b = a + Nfft   - 1 ;
xp = x(a:b);

a1 =  Nguard + 2*Nfft + round(0.8*Nfft) + 1 ;
b2 = a1 + Nfft  - 1 ;
xd = x(a1:b2);







Xp = fft(xp);
Xp = Xp( 1:floor(end/2) );
Xd = fft(xd);
Xd = Xd( 1:floor(end/2) );

if plotcrap
    n = length(Xd);
    fvec = linspace(0,fs/2,n);
    figure(103)
    plot(fvec,20*log10(abs(Xp)))
    hold on
    plot(fvec,20*log10(abs(Xd)))
    legend('Pilot FFT', 'Data FFT');
    hold off
end










detect = angle(Xp ./ Xd);          % ratio to get delta phase diff
%detect = Xp ./ Xd;          % ratio to get delta phase diff
%detect = abs(detect);    %only for bpsk

thres = pi/2;

% detect = ( detect <= thres) .* 1   +  (detect > thres) .* 0;
detect = detect(513:end);
detect = detect(1:(bits+1)/3);   % bits for char  

if plotcrap
    figure(104)
    plot(detect)
    title('demodulate phase vector')
end



%detect(0<=detect <= thres) = 0;
%detect(detect > thres) = 1;

%detect(detect <= thres) = 1;
%detect(detect > thres) = 0;

%detect=pskdemod(detect,4); matlab function
det2=detect;
det3=[];
 for rr=1:2731
if abs(det2(rr))<pi/8
det3(rr,:)=[0 0 0];
elseif abs(pi/4-det2(rr))<pi/8
det3(rr,:)=[0 0 1];
elseif abs(pi/2-det2(rr))<pi/8
det3(rr,:)=[0 1 0];
elseif abs(3*pi/4-det2(rr))<pi/8
det3(rr,:)=[0 1 1];
elseif abs(det2(rr))>3
det3(rr,:)=[1 0 0];
elseif abs(-3*pi/4-det2(rr))<pi/8
det3(rr,:)=[1 0 1];
elseif abs(-pi/2-det2(rr))<pi/8
det3(rr,:)=[1 1 0];
else % abs(3*pi/2-det2(rr))<=pi/4
det3(rr,:)=[1 1 1];
end
end
detect=det3;

%%


%%
if plotcrap
    figure(105)
    plot(detect)
    title('decoded bits vector')
end
pause on

%detect=dec2bin(detect)
pause()
%detect = reshape(detect,8,[]);
%det1=[512,8];
%for z=0:511
   % z4=z*4+1;
    %zz4=z4+3;
%det1(z,:)=detect(z*4+1:z*4+3);
%end
detect = detect';
detect=detect(:,1:2728);
detect=reshape(detect,8,1023);
detect=detect';

charmes = zeros(1,Nchar);
w = [7:-1:0];
w = 2.^w;


for k = 1:Nchar-1   %for8psk
    x = detect(k,:);
    x = x .* w;
    x = sum(x);
    charmes(k) = x;
end



    clc

   xstr =  char(charmes)

   msgbox(xstr,'replace')
