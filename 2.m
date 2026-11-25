clear; clc; close all;

N = 1000; trials = 1e5;
snr_dB = -10:2:15;
Pfa = 1e-3;
thr = sqrt(-2*log(Pfa));
sigma = 10.^(-snr_dB/20);
L = numel(snr_dB);


maxG = @(r,c) sqrt(2)*erfcinv(-2*expm1(log(rand(r,c))/N));

Mx = maxG(trials, L);                 
Mn = maxG(trials, L);                 

Pd_sim  = mean(1 + sigma.*Mx > thr*sigma);
Pfa_sim = mean(sigma.*Mn > thr*sigma);

Pd_th = marcumq(sqrt(2*10.^(snr_dB/10)), thr);

figure;
plot(snr_dB, Pd_sim, 'ro-', snr_dB, Pd_th, 'b--', 'LineWidth', 1.5);
grid on; xlabel('SNR (dB)'); ylabel('P_d');
title(['P_d vs SNR (P_{fa} = ', num2str(Pfa), ')']);
legend('Simulated', 'Theoretical');

pfa = logspace(-4, 0, 100);
[P, S] = meshgrid(pfa, sqrt(2*10.^([0 5 10]/10)));
pd = marcumq(S, sqrt(-2*log(P)));

figure;
semilogx(pfa, pd', 'LineWidth', 1.5);
grid on; xlabel('P_{fa}'); ylabel('P_d');
title('ROC Curves'); legend('0 dB','5 dB','10 dB');
