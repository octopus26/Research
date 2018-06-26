close all;clear all;
%% Parameters
Buoy_R=3; %[m]
Buoy_Draft_o=20; %[m]
    
gravity = 9.80665; 
rho=1027;
Cd = 0.7;

z = linspace(0,(-1)*Buoy_Draft_o,500);
x = 0;
%% Incident Wave
Hs=10; %[m]
Tp=14; %[s]
gamma = 1.7;
domega=0.01; %[rad/s]
time_record = linspace(0,2*pi/domega,1200);
omega_max = sqrt(2*gravity/Hs);
omega_min=0.1;
N_omega=1+ceil(omega_max-omega_min)/domega;
omega_s = linspace(omega_min,omega_max,N_omega);
% Johnswap spectrum
S_zeta = 155*Hs^2./(Tp^4*omega_s.^5).*exp(-944./(Tp^4.*omega_s.^4))*3.3^gamma;

% A(omega) & k(omega)
amp_s =sqrt(2.*S_zeta*domega);
%k_w_s =omega_s.^2./gravity;  %deep water
dispersion = @(k_w_s) omega_s.^2 - k_w_s*gravity.*tanh(k_w_s*Buoy_Draft_o);
k_w_s = fsolve(dispersion,ones(1,length(omega_s)));

fileID1 = fopen('amp.txt','at');
fprintf(fileID1,'\n%e',amp_s');
fileID2 = fopen('omega.txt','at');
fprintf(fileID2,'\n%e',omega_s');
%Incident wave potential 

phase = 2*pi*rand(1,length(omega_s));
fileID = fopen('phase.txt','at');
fprintf(fileID,'\n%e',phase');
Phio = zeros(length(time_record),length(z));
for count_z = 1:length(z)
    for count_t = 1:length(time_record)
        Phio(count_t,count_z) = real(sum(1i*gravity*amp_s./omega_s.*exp(k_w_s.*z(count_z)-1i*k_w_s.*x + 1i*omega_s.*time_record(count_t)+ 1i*phase)));
    end 
end 
%%
%Incident wave elevation 
zetao = zeros(1,length(time_record));
for count_t = 1:length(time_record)
zetao(count_t) = real(sum(amp_s.*exp(-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
end
%%
%Incident wave velocity
u1 = zeros(length(time_record),length(z));
u1dot = zeros(length(time_record),length(z));
u1x = zeros(length(time_record),length(z));
u1z = zeros(length(time_record),length(z));
u3 = zeros(length(time_record),length(z));
for count_z = 1:length(z)
    for count_t = 1:length(time_record)
        u1(count_t,count_z) = real(sum((-1i)*k_w_s.*1i*gravity.*amp_s./omega_s.*exp(k_w_s.*z(count_z)-1i*k_w_s.*x + 1i*omega_s.*time_record(count_t)+ 1i*phase)));
        u1dot(count_t,count_z) = real(sum((1i*omega_s).*((-1i)*k_w_s).*1i*gravity.*amp_s./omega_s.*exp(k_w_s*z(count_z)-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
        u1x(count_t,count_z) = real(sum(((-1i)*k_w_s).^2.*1i*gravity.*amp_s./omega_s.*exp(k_w_s*z(count_z)-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
        u1z(count_t,count_z) = real(sum(k_w_s.*((-1i)*k_w_s).*1i*gravity.*amp_s./omega_s.*exp(k_w_s*z(count_z)-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
        u3(count_t,count_z) = real(sum(k_w_s.*1i*gravity.*amp_s./omega_s.*exp(k_w_s*z(count_z)-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
       
    end 
end 

%Derivatives of incident wave elevation
zetaot = zeros(1,length(count_t));
zetaox = zeros(1,length(count_t));
for count_t = 1:length(time_record)
zetaot(count_t) = real(sum(1i*omega_s.*amp_s.*exp(-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
zetaox(count_t) = real(sum(-1i*k_w_s.*amp_s.*exp(-1i*k_w_s*x + 1i*omega_s*time_record(count_t)+ 1i*phase)));
end

%% Morison Force
F_Morison = 2*rho*pi*Buoy_R^2*(trapz(u1dot,2)');  %Classic Morison 
F_body_nonlinear = 2*rho*pi*Buoy_R^2*trapz(u1.*u1x+u3.*u1z,2)';
F_wave_ele_nonl = 2*rho*pi*Buoy_R^2*zetao.*(u1dot(:,1)+u1(:,1).*u1x(:,1)+u3(:,1).*u1z(:,1))';

F_crest = rho*pi*Buoy_R^2*u1(:,1)'.*(zetaot-u1(:,1)'.*zetaox);

F_viscous = 0.5*rho*2*Buoy_R*Cd*trapz(u1.*abs(u1),2)'...
                                    +0.5*rho*2*Buoy_R*Cd*zetao.*(u1(:,1).*abs(u1(:,1)))';
%% Plots
figure;
plot(time_record,F_Morison,'r');
hold on;
plot(time_record,F_Morison_cross,'b');
hold on;
plot(time_record,F_Morison_wave,'g');
hold on;
plot(time_record,F_Morison_viscous,'y');
legend('Classic Morison','Nonlinear Morison','Wave Effect Morison','Viscous Effect Morison');

