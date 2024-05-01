close all
clear all
clc

fs=8000
Nchar=1024   %for psk
guardN=512  

%%
%N = 6000;  %number of samples for preamble
t = [0:7999]/fs;
bits=Nchar*8;
Nfft = (bits+1)/3 + 1024;   %for 8psk

pream=sin(2*pi*256*t).*exp(t);
pream = pream / max(pream);

%%get message from txt file
fid = fopen('D:\Users\Student\Desktop\0_readme.important.txt');
mtxt = fread(fid);
fclose(fid);
mtxt = mtxt';
if length(mtxt) >= Nchar
    mtxt = mtxt(1:Nchar);
else
    z = length(mtxt);
    z = Nchar - z;
    mtxt = [mtxt zeros(1,z)];
end
mbits = dec2bin(mtxt,8);
mbits = reshape(mbits',1,[]);
mbits = mbits(1:bits);
mbits=[mbits 0];     %for 8psk

%%
mbits=mbits';
mbits1=mbits(1:3:8191);
mbits2=mbits(2:3:8192);
mbits3=mbits(3:3:8193);
mbits=[mbits1,mbits2,mbits3];
%mbits=reshape(mbits,2048,2);
 mbits=bin2dec(mbits);
txbits = zeros(1,bits/2); %for qpsk

%%
for k = 1:(bits+1)/3 %for 8psk
    if mbits(k) == '1'
        txbits(k) = 1;
    elseif mbits(k) == '2'
         txbits(k) = 2;
    elseif mbits(k) == '3'
         txbits(k) = 3;
          elseif mbits(k) == '4'
         txbits(k) = 4;
          elseif mbits(k) == '5'
         txbits(k) = 5;
          elseif mbits(k) == '6'
         txbits(k) = 6;
          elseif mbits(k) == '7'
         txbits(k) = 7;
    else
        txbits(k) = 0;
    end
end
% Above section has no effect

txbits=mbits;



% Perform psk mod
psk_modulated_data = pskmod(txbits, 8);

psk_modulated_data=psk_modulated_data';
psk_modulated_data= [zeros(1,512)   psk_modulated_data    zeros(1,512) ];  %512 guard on each side
%%

% make pilot vector
xp = randn(1,Nfft);
xp = xp / max(xp);
xp = xp * 12;


% pilot vector with 512 guard on each side
xp = exp(1i*xp);
xp(1:512) = 0;
xp(end-511:end) = 0;

% relate data vector to pilot vector
psk_modulated_data = xp .* psk_modulated_data;  % e^a * e^b =  e^(a + b)

xp1 = conj( xp(end:-1:2));
xp = [xp1 xp];

psk_modulated_data1 = conj(psk_modulated_data(end:-1:2));
psk_modulated_data = [psk_modulated_data1 psk_modulated_data];



td = 2;
tdsamples = round(td * fs);
tdvec = zeros(1,tdsamples);

%%Generate random data
     data_in=randi([0 1],2048,1);       
   % data_in=reshape(data_in,[2048,8]);
    data_in=reshape(data_in,[256,8]);
    data_in(:,8)=0; %CP
    data_in=reshape(data_in,[2048,1]);


%%
figure(1),stem(txbits); grid on; xlabel('Data Points'); ylabel('Amplitude')
title('Original Data ')

%%

% pilot time vector
xp = ifftshift(xp);
xp = ifft(xp);
xp = xp / max(abs(xp));


%% IFFT 256 subcarriers

psk_modulated_data=ifftshift(psk_modulated_data);
 txout = ifft(psk_modulated_data);
 txout = txout / max(abs(txout));
 %txout=txout';
xz = zeros(1,guardN);
txout=[tdvec pream xz xp xp txout txout xz tdvec];
%sigtime = Nguard + 2*Nfft + 2*Nfft + 2*Nguard



