%% Initialization
clc
clear
close all

%% Parameters
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
h10 = 12.4; % [cm]
h20 = 12.7; % [cm]
h30 = 1.8;  % [cm]
h40 = 1.4;  % [cm]
v10 = 3.00; % [V]
v20 = 3.00; % [V]
k1 = 3.33;  % [cm^3/V*s]
k2 = 3.35;  % [cm^3/V*s]
gamma1 = 0.7;
gamma2 = 0.6;

% Parameters P+ NON-MINIMUM PHASE
% h10 = 12.6; % [cm]
% h20 = 13; % [cm]
% h30 = 4.8;  % [cm]
% h40 = 4.9;  % [cm]
% v10 = 3.15; % [V]
% v20 = 3.15; % [V]
% k1 = 3.14;  % [cm^3/V*s]
% k2 = 3.19;  % [cm^3/V*s]
% gamma1 = 0.43;
% gamma2 = 0.34;

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
% k = 0;  
% syms x;
% double(solve((1+x*T3)*(1+x*T4)-k==0));

%% State-space matrices
% The system is input-decoupled
Atot0 =[-1/T1 0 A3/A1/T3 0
    0 -1/T2 0 A4/A2/T4
    0 0 -1/T3 0
    0 0 0 -1/T4];         
Bdec0{1}=[gamma1*k1/A1 0 0 (1-gamma1)*k1/A4]';
Bdec0{2}=[0 gamma2*k2/A2 (1-gamma2)*k2/A3 0]';
% C_real = [kc 0 0 0; 0 kc 0 0];
Ctot=eye(4);
Cdec0{1}=Ctot([1,4],:);
Cdec0{2}=Ctot(2:3,:);
D = 0;

% use suitable change of coordinates to make the state partition 
%  xnew=[xnew1',xnew2']', where xnew1 includes h1 and h4, while
% x2 includes h2 and h3

T=[1 0 0 0;
    0 0 0 1
    0 1 0 0
    0 0 1 0];
Atot=T*Atot0/T;
Bdec{1}=T*Bdec0{1};
Bdec{2}=T*Bdec0{2};
Cdec{1}=Cdec0{1}/T;
Cdec{2}=Cdec0{2}/T;

B(:,1)=Bdec{1};
B(:,2)=Bdec{2};
C(1:2,:)=Cdec{1};
C(3:4,:)=Cdec{2};

N = 2;
rounding_n = 3;

%% C-T system analysis
systemCT = ss(Atot,B,C,D);
% Step response
figure;
step(systemCT);
% Impulse response
figure;
impulse(systemCT);
% Eigenvalues
eig(Atot);
% Spectral abscissa
rho = max(real(eig(Atot)));
% Eigenvalues on the diagonal of D and eigenvectors on columns of V
[Vec, Dc] = eig(Atot);
% Columns of V are the generalized eigenvectors
[Vjc, Jc] = jordan(Atot);

%% D-T system analysis
% sampling time TS chose considering Ttransient = 5/rho
TS = 1.0;
systemDT = c2d(systemCT, TS);
[Ftot,G,H,L,Ts]=ssdata(systemDT);
Gdec{1} = G(:,1);
Gdec{2} = G(:,2);
Hdec{1} = H(1:2,:);
Hdec{2} = H(3:4,:);

% Step response
figure;
step(systemDT);
% Impulse response
figure;
impulse(systemDT);
% Eigenvalues
eig(Ftot);
% Spectral abscissa
rho = max(abs(eig(Ftot)));
% Eigenvalues on the diagonal of D and eigenvectors on columns of V
[Ved, Dd] = eig(Ftot);
% Columns of V are the generalized eigenvectors
[Vjd, Jd] = jordan(Ftot);

%% Centralized Control
ContStruc = ones(N,N);
[CFM_CT]=di_fixed_modes(Atot,Bdec,Cdec,N,ContStruc,rounding_n);
[CFM_DT]=di_fixed_modes(Ftot,Gdec,Hdec,N,ContStruc,rounding_n);
[K_c_CT, rho_c_CT, feas_c_CT] = LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc);
[K_c_DT, rho_c_DT, feas_c_DT] = LMI_DT_DeDicont(Ftot,Gdec,Hdec,N,ContStruc);

%% Decentralized Control
ContStruc = diag(ones(N,1));

[DFM_CT]=di_fixed_modes(Atot,Bdec,Cdec,N,ContStruc,rounding_n);
[DFM_DT]=di_fixed_modes(Ftot,Gdec,Hdec,N,ContStruc,rounding_n);
[K_dec_CT, rho_dec_CT, feas_dec_CT] = LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc);
[K_dec_DT, rho_dec_DT, feas_dec_DT] = LMI_DT_DeDicont(Ftot,Gdec,Hdec,N,ContStruc);

%% Distributed Control
% ContStruc = [1 1
%               1 1];





