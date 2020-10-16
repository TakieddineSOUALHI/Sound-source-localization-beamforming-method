function [out,y] = beamforming(x_mic,fe,ANTENNE,theta_0)
% BEAMFORMING Compute the output of a classical beamformer
%
% BEAMFORMING(MICRO,ARRAY, theta_0) computes the output of a
% classical beamformer polarized in the azimuth THETA_0. MICRO is a 
% structure containing the microphones outputs and their sampling frequency
% ARRAY is a structure which must at leat contain :
%     - ARRAY.N : microphone number,
%     - ARRAY.C : speed of sound value,
%     - ARRAY.D : microphone interspace (in meter),
%     - ARRAY.Pos : microphones position (linear vector)
%
% [y,out] = BEAMFORMING(x_mic,fs,ARRAY, theta_0) also computes the
% beamformer filter outputs in the OUT vector.

N = length(x_mic);
X_source = fft(x_mic);
f = (0:fe/N:(N-1)/N*fe).';

tau_0 = -ANTENNE.Pos./ANTENNE.C.*cos(theta_0*pi/180);
delay = exp(j*2*pi*kron(f,tau_0));
Y = X_source.*delay;

if (mod(N,2)==0) %N est pair
    Y(N/2+2:end,:) = conj(Y(N/2:-1:2,:));
else % N est impair
    Y(round(N/2)+1:end,:) = conj(Y(round(N/2):-1:2,:));
end

y = real(ifft(Y));
out = sum(y,2)/ANTENNE.N;