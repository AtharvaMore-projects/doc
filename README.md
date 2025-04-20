Link Power 

clc;  
clear all; 
close all; 
pin=input('Enter the mean input optical power launched into the fiber in dBm'); 
po=input('Enter the mean incident optical power required at receiver in dBm');  
af =input('Enter the in the attenuation factor in the fiber');  
Lsp=input('Enter the in the splice loss per km');  
Lc=input('Enter the in the total connector losses');  
M=input('Enter the required safety margin');  
L=(pin-po-Lc-M)/( af +Lsp);  
disp('The maximum link length is');  
disp(L);


pin=input('Enter the mean input optical power launched into the fiber in dBm'); 
po=input('Enter the mean incident optical power required at receiver in dBm');  
af =input('Enter the in the attenuation factor in the fiber');  
Lsp=input('Enter the in the splice loss per km');  
Lc=input('Enter the in the total connector losses');  
M=input('Enter the required safety margin');  
L=(pin-po-Lc-M)/( af +Lsp);  
disp('The maximum link length is');  
disp(L);

pin=input('Enter the mean input optical power launched into the fiber in dBm'); 
po=input('Enter the mean incident optical power required at receiver in dBm');  
af =input('Enter the in the attenuation factor in the fiber');  
Lsp=input('Enter the in the splice loss per km');  
Lc=input('Enter the in the total connector losses');  
M=input('Enter the required safety margin'); 
DEP=input('Enter the required DEP');
L=(pin-po-Lc-M-DEP)/( af +Lsp);  
disp('The maximum link length is');  
disp(L);




















Rise Time
clc;
clear all; 
close all; 
td=6;
ts=8;
l=input("fiber lenght: "); 
inter=5; 
tinter=(l*inter); 
intra=1; 
tintra =(l*intra);
tsys=  sqrt (ts^2+td^2+tinter^2+tintra^2); 
fprintf('Total system Rise time:'); 
disp (tsys) 
nrz=(0.7/tsys*1000); 
fprintf('bit rate for NRZ format' ); 
disp (nrz) 
rz=(0.35/tsys*1000); 
fprintf('Bit rate for RZ format: '); 
disp(rz)








SI & GI
clc;
close all;
clear all;
% Parameters
n1 = 1.5; delta = 0.01; a = 30;
n2 = n1 * sqrt(1 - 2 * delta);

% Step-Index Data
r_full = linspace(0, 50, 500);
n_step_full = n1 * ones(size(r_full));
n_step_full(r_full > a) = n2;

% Graded-Index Data
r = linspace(0, a, 500);
m = r / a;
n_step = n1 * ones(size(r));
n_tri = n1 .* sqrt(1 - 2 * delta * m);       % Triangular (α = 1)
n_para = n1 .* sqrt(1 - 2 * delta * m.^2);   % Parabolic (α = 2)

% Plotting
figure;

% Subplot 1: Step-Index
subplot(1,2,1);
plot(r_full, n_step_full, 'b', 'LineWidth', 1.5); hold on;
xline(a, 'r--', 'Core Boundary');
xlabel('Radial Distance (\mum)');
ylabel('Refractive Index');
title('Step-Index Fiber');
legend('Step Index'); grid on;

% Subplot 2: Graded-Index
subplot(1,2,2);
plot(r, n_step, 'k--', r, n_tri, 'g', r, n_para, 'm', 'LineWidth', 1.5); hold on;
xline(a, 'm--', 'Core Boundary');
xlabel('Radial Distance (mm)');
ylabel('Refractive Index');
title('Graded-Index Profiles');
legend('Step (\alpha=\infty)', 'Triangular (\alpha=1)', 'Parabolic (\alpha=2)'); grid on;
