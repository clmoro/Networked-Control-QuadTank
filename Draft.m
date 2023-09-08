%  General Parameters
A1 = 28;    % [cm^2]
A2 = 32;    % [cm^2]
A3 = A1;    % [cm^2]
A4 = A2;    % [cm^2]
a1 = 0.071; % [cm^2]
a2 = 0.057; % [cm^2]
a3 = a1;    % [cm^2]
a4 = a2;    % [cm^2]
kc = 0.50;  % [V/cm]
g = 981;    % [cm/s^2]

% % Parameters P- MINIMUM PHASE
% h10 = 12.4; % [cm]
% h20 = 12.7; % [cm]
% h30 = 1.8;  % [cm]
% h40 = 1.4;  % [cm]
% v10 = 3.00; % [V]
% v20 = 3.00; % [V]
% k1 = 3.33;  % [cm^3/V*s]
% k2 = 3.35;  % [cm^3/V*s]
% gamma1 = 0.7;
% gamma2 = 0.6;

% Parameters P+ NON-MINIMUM PHASE
h10 = 12.6; % [cm]
h20 = 13; % [cm]
h30 = 4.8;  % [cm]
h40 = 4.9;  % [cm]
v10 = 3.15; % [V]
v20 = 3.15; % [V]
k1 = 3.14;  % [cm^3/V*s]
k2 = 3.19;  % [cm^3/V*s]
gamma1 = 0.43;
gamma2 = 0.34;

% Parameters that depends on the operating condition
T1 = A1/a1*sqrt(2*h10/g);
T2 = A2/a2*sqrt(2*h20/g);
T3 = A3/a3*sqrt(2*h30/g);
T4 = A4/a4*sqrt(2*h40/g);
c1 = T1*k1*kc/A1;
c2 = T2*k2*kc/A2;

% Transfer function
s = tf('s');
G11 = gamma1*c1/(1+s*T1);
G12 = (1-gamma2)*c1/((1+s*T3)*(1+s*T1));
G21 = (1-gamma1)*c2/((1+s*T4)*(1+s*T2));
G22 = gamma2*c2/(1+s*T2);
G = [G11 G12;G21 G22];

%% Zero location
% condition: 1-(T3+T4)^2/(4*T3*T4)<=k<=1 
% if k>1 --> s>0 
% if k<-0.0334 --> s are imaginaries (Re<0) 
% if k = 0 --> s are -1/T3 and -1/T4
k = 0;  
syms x;
double(solve((1+x*T3)*(1+x*T4)-k==0));

%% State-space matrices
A = [-1/T1 0 A3/(A1*T3) 0; 0 -1/T2 0 A4/(A2*T4); 0 0 -1/T3 0; 0 0 0 -1/T4];
B = [(gamma1*k1)/A1 0; 0 (gamma2*k2)/A2; 0 ((1-gamma2)*k2)/A3; ((1-gamma1)*k1)/A4 0];
C = [kc 0 0 0; 0 kc 0 0];
D = 0;
% The system is input-decoupled since B is block diagonal after switching row 2 and row 4
% A = [-1/T1 0 A3/(A1*T3) 0; 0 0 0 -1/T4; 0 0 -1/T3 0; 0 -1/T2 0 A4/(A2*T4)];
% B1 = [(gamma1*k1)/A1 0; ((1-gamma1)*k1)/A4 0];
% B2 = [0 (gamma2*k2)/A2; 0 ((1-gamma2)*k2)/A3];
% B = [B1; B2];
% C1 = [kc 0 0 0];
% C2 = [0 kc 0 0];
% C = [C1; C2];

%% C-T system analysis
systemCT = ss(A,B,C,D);
% Step response
figure;
step(systemCT);
% Impulse response
figure;
impulse(systemCT);
% Eigenvalues
eig(A);
% Spectral abscissa
rho = max(real(eig(A)));
% Eigenvalues on the diagonal of D and eigenvectors on columns of V
[Vec, Dc] = eig(A);
% Columns of V are the generalized eigenvectors
[Vjc, Jc] = jordan(A);

%% D-T system analysis
% sampling time TS chose considering Ttransient = 5/rho
TS = 1.0;
systemDT = c2d(systemCT, TS);
[F,G,H,L,Ts]=ssdata(systemDT);

% Step response
figure;
step(systemDT);
% Impulse response
figure;
impulse(systemDT);
% Eigenvalues
eig(F);
% Spectral abscissa
rho = max(abs(eig(F)));
% Eigenvalues on the diagonal of D and eigenvectors on columns of V
[Ved, Dd] = eig(F);
% Columns of V are the generalized eigenvectors
[Vjd, Jd] = jordan(F);

%%




