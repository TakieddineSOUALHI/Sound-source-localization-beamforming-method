close all
clear all

%% Etude du diagramme d'antenne
% -------------------------------------------------------------------------
ANTENNE.N = 2; % Nombre de microphones
ANTENNE.C = 340; % Vitesse du son
ANTENNE.D = 0.15; % Distance entre les microphones
for i=1:ANTENNE.N % Position des microphones
    ANTENNE.Pos(i) = (i-(ANTENNE.N + 1)/2)*ANTENNE.D;
end

% Source
load('source1.mat');

% COMPLETER LE CODE ...

Fs = 44100;
dt = 1/Fs;
t = 0:dt:0.05;

F = 100:500:8000;

reponse_antenne = zeros(length(F),length(0:180));
for i = 1:length(F)
    SOURCE.signal = sin(2*pi*F(i)*t);
    e = sum(SOURCE.signal.*conj(SOURCE.signal));
%     reponse_antenne = [];
    for theta = 0:180
        y = diagramme(SOURCE,theta,ANTENNE);
        reponse_antenne(i,theta+1) = y'*y/e;
    end
%     figure();
%     polarplot(deg2rad(0:180),reponse_antenne)
%     thetalim([0 180])
end

surf(0:180,F,reponse_antenne)
shading interp

%% formation de voix
load('data1.mat')
theta_0 = 0:180;

energie = zeros(1,length(theta_0));
e = sum(SOURCE.signal.*conj(SOURCE.signal));
for i = 1:length(theta_0)
    y = beamforming(MICROS.Signal,MICROS.fe,ANTENNE,theta_0(i));
    energie(i) = y'*y/e;
end

polarplot(deg2rad(0:180),energie)

close all;

load('data2.mat')

K = 512;

nb_trames = floor(length(MICROS.Signal)/K);

energie = zeros(nb_trames,length(theta_0));
% e = sum(SOURCE.signal.*conj(SOURCE.signal));
e=1;
for trame = 1:nb_trames
    for i = 1:length(theta_0)
        y = beamforming(MICROS.Signal(K*(trame-1)+1:K*trame,:),MICROS.fe,ANTENNE,theta_0(i));
        energie(trame,i) = y'*y/e;
    end
%     polarplot(deg2rad(0:180),energie(trame,:))
%     hold on;
%     pause(0.1)
end
 
dt = MICROS.t(2) - MICROS.t(1);
dK = K*dt;
trames2time = dK.*(1:nb_trames);
figure();
surf(theta_0,trames2time,energie);
shading interp;
colorbar();

idxes = [];
for i = 1:size(energie,1)
    [val,idx]=max(energie(i,:));
    idxes = [idxes, idx-1];
end
figure();
plot(trames2time,idxes);
