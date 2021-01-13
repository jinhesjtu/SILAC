clear
close all
%% Radar application for INCM-SILAC Beamformer

%%%%%%%%%  data generation
lambda = 0.03; N = 256; fs = 25600;
Sensor = 10; Sp = [0:Sensor-1]*0.5;

thetat = 5;   fd = 0.0911;  At = exp(-j*2*pi*Sp'*sind(thetat));  vt = lambda*fd*fs/2;
thetaj = [-17, 25]; Aj = exp(-j*2*pi*Sp'*sind(thetaj));
SnrdB = 20; InrdB = 50;
Snr=sqrt(10.^(SnrdB/10)); Inr = sqrt(10.^(InrdB/10));
st = exp(j*2*pi*fd*(0:N-1));
sj = exp(j*2*rand(2,N)*pi);
Noise = sqrt(0.5)* (randn(Sensor,N) + j* randn(Sensor,N));
g = [1,1,1.0139 + 0.0323i,0.9992 - 0.0178i,0.9956 + 0.0179i,1.0068 + 0.0040i,1.0096 - 0.0132i,1.0008 + 0.0098i,0.9938 - 0.0164i,1];
G = diag(g);
X = [Snr*G*At, Inr*G*Aj] * [st;sj] + Noise;
Xs = Snr*G*At*st;
%%%%%%%%%  data generation

%%%%
bt = At; yo = bt'*X; Y = abs(fft(yo)); yos = bt'*Xs; Ys = abs(fft(yos));
% Y = sum(abs(fft(X.').'));
Y = 0.5*db(Y);
% Ys = sum(abs(fft(Xs.').'));
Ys = 0.5*db(Ys);
x = linspace(-192,192,N);
figure; set(gcf,'DefaultLineLineWidth',1.5)
plot(x,fftshift(Y)); hold on; grid on;
plot(x,fftshift(Ys), '-o')
xlabel('Velocity, in m/s')
ylabel('Dopple map, in dB');
legend('Received signals', 'Target signals',  'Location', 'NorthWest')
title('(a) Doppler map of detection beam, without SILAC')
%%%

sref = 10;SigAll = 3;Scal = [1,2,10];
C0 = cum4mtx(X(1,:),X(1,:),X,X); C1 = cum4mtx(X(1,:),X(sref,:),X,X);
CC = [C0;C1]; [UU SS VV] = svd(CC);
Es = UU(:,1:SigAll - 1); E1 = Es(1:Sensor,:); E2 = Es(Sensor+1:2*Sensor,:);
[T, Phi] = eig(pinv(E1)*E2);  A1 = E1*T + E2*T*inv(Phi); phi = diag(Phi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Jamming Localization
for num = 1:SigAll - 1
    n1 = -round(Sp(sref)):round(Sp(sref));
    pall = angle(phi(num)) + 2*pi*n1 ;
    uamall = pall/(2*pi*Sp(sref));
    uall = uamall(find(abs(uamall) < 1));
    for nc = 1:length(uall)
        a(:,nc) = exp(-j*2*pi*Sp'*uall(nc));
        p(:,nc) =  pinv(diag(a(:,nc)))*A1(:,num);
    end
    q = p(Scal,:);
    vall = mean(abs(diff(q)));
    [bb aa] = min(vall);
    uhat(num) = uall(aa);
    Aje((num-1)*Sensor+1: num*Sensor, :) = diag(a(:,aa));
    A1lg((num-1)*Sensor+1: num*Sensor, :) = A1(:,num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Jamming Localization
ghat = pinv(Aje) * A1lg;    ghat = ghat/ghat(1);  %% Array error
Ghat = diag(ghat);

%%%%% INC Reconstruction
R = X*X'/N;   [U S V] = svd(R);  Es = U(:,1:SigAll); Ds = S(1:SigAll, 1:SigAll);
Evall = diag(S);  evs = Evall(1:SigAll); evn = Evall(SigAll + 1:end);
sn = mean(evn);
for num = 1:SigAll - 1
    a = exp(-j*2*pi*Sp'*uhat(num));
    b = Ghat*a;
    Ajhat(:,num) = b;
    sigma(num) = inv(b'*pinv(Es*(Ds - sn*eye(SigAll))*Es')*b);
end
Rine = Ajhat * diag(sigma) * (Ajhat)' + sn*eye(Sensor);
%%%%% INC Reconstruction
bsoi = Ghat*At;
w = inv(Rine)*bsoi/ (bsoi'*inv(Rine)*bsoi);

yo = w'*X; Y = abs(fft(yo));
Y = 0.5*db(Y) - 20;
x = linspace(-192,192,N);
figure; set(gcf,'DefaultLineLineWidth',1.5)
plot(x,fftshift(Y),'-ro'); hold on; grid on;
xlabel('Velocity, in m/s')
ylabel('Dopple map, in dB');
legend('Detection SNR \approx 15 dB', 'Location', 'NorthWest')
title('(b) Doppler map of INCM-SILAC beamforming')

