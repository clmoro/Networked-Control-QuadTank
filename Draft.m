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
Gtf = [G11 G12;G21 G22];

%% Zero location
% k is the eta parameter of the paper, where k =
% (1-gamma1)*(1-gamma2)/(gamma1*gamma2)
% condition: 1-(T3+T4)^2/(4*T3*T4)<=k<=1 
% if k>1 --> s>0 
% if k<-0.0334 --> s are imaginaries (Re<0) 
% if k = 0 --> s are -1/T3 and -1/T4
k = 0;  
syms x;
TheoreticalZero = double(solve((1+x*T3)*(1+x*T4)-k==0));
% disp(['k = ',num2str(k),'     zero =',num2str(TheoreticalZero(1)),'   ',num2str(TheoreticalZero(2))])

k = (1-gamma1)*(1-gamma2)/(gamma1*gamma2);
ActualZero = double(solve((1+x*T3)*(1+x*T4)-k==0));
disp(['Actual Zeros are in ',num2str(ActualZero(1)),' and ',num2str(ActualZero(2)),'.'])

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
eig_OL_CT=eig(Atot);
% Spectral abscissa
rho_CT = max(real(eig(Atot)));
% Eigenvalues on the diagonal of D and eigenvectors on columns of V
[Vec, Dc] = eig(Atot);
% Columns of V are the generalized eigenvectors
[Vjc, Jc] = jordan(Atot);

% Eigenvalues plot
figure
plot(real(eig_OL_CT),imag(eig_OL_CT),'*b','MarkerSize',8,'LineWidth',2);
grid;title('CT Open loop eigenvalue Positions','FontSize',20);legend('Open loop','FontSize',18)
xlabel('Real Axis','FontSize',18);ylabel('Imaginary Axis','FontSize',18);

%% D-T system analysis
% sampling time TS chose considering Ttransient = 5/rho
TS = 1;
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
eig_OL_DT=eig(Ftot);
% Spectral abscissa
rho_DT = max(abs(eig(Ftot)));
% Eigenvalues on the diagonal of D and eigenvectors on columns of V
[Ved, Dd] = eig(Ftot);
% Columns of V are the generalized eigenvectors
[Vjd, Jd] = jordan(Ftot);

% Eigenvalues plot
figure
x = [-1:0.01:1];
y = sqrt(1-x.^2);
hold on
plot(zeros(size(x)),x,':k',x,zeros(size(x)),':k',x,y,':k',x,-y,':k',0,0,'.b')
plot(real(eig_OL_DT),imag(eig_OL_DT),'*r','MarkerSize',8,'LineWidth',2);legend('','','','','','Open loop','FontSize',18)
title('DT Open loop eigenvalue positions','FontSize',20);
xlabel('Real Axis','FontSize',18);ylabel('Imaginary Axis','FontSize',18);
hold off

%% STABLE CONTROLLER
NumCont = [2 4];    % [The number of type controller implemented, for each type how many]
Gain{1} = [];

% Centralized
ContStruc = ones(N,N);
[CFM_CT]=di_fixed_modes(Atot,Bdec,Cdec,N,ContStruc,rounding_n);
[CFM_DT]=di_fixed_modes(Ftot,Gdec,Hdec,N,ContStruc,rounding_n);
[K_c_CT, rho_c_CT, feas_c_CT] = LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc);
[K_c_DT, rho_c_DT, feas_c_DT] = LMI_DT_DeDicont(Ftot,Gdec,Hdec,N,ContStruc);
Gain{1} = [K_c_CT,K_c_DT];

% Decentralized
ContStruc = diag(ones(N,1));
[DFM_CT]=di_fixed_modes(Atot,Bdec,Cdec,N,ContStruc,rounding_n);
[DFM_DT]=di_fixed_modes(Ftot,Gdec,Hdec,N,ContStruc,rounding_n);
[K_dec_CT, rho_dec_CT, feas_dec_CT] = LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc);
[K_dec_DT, rho_dec_DT, feas_dec_DT] = LMI_DT_DeDicont(Ftot,Gdec,Hdec,N,ContStruc);
Gain{1} = [Gain;
        K_dec_CT,K_dec_DT];

% Distributed 1
ContStruc = [1 0
              1 1];
[Dist1FM_CT]=di_fixed_modes(Atot,Bdec,Cdec,N,ContStruc,rounding_n);
[Dist1FM_DT]=di_fixed_modes(Ftot,Gdec,Hdec,N,ContStruc,rounding_n);
[K_dist1_CT, rho_dist1_CT, feas_dist1_CT] = LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc);
[K_dist1_DT, rho_dist1_DT, feas_dist1_DT] = LMI_DT_DeDicont(Ftot,Gdec,Hdec,N,ContStruc);
Gain{1} = [Gain;
         K_dist1_CT,K_dist1_DT];

% Distributed 2
ContStruc = [1 1
              0 1];
[Dist2FM_CT]=di_fixed_modes(Atot,Bdec,Cdec,N,ContStruc,rounding_n);
[Dist2FM_DT]=di_fixed_modes(Ftot,Gdec,Hdec,N,ContStruc,rounding_n);
[K_dist2_CT, rho_dist2_CT, feas_dist2_CT] = LMI_CT_DeDicont(Atot,Bdec,Cdec,N,ContStruc);
[K_dist2_DT, rho_dist2_DT, feas_dist2_DT] = LMI_DT_DeDicont(Ftot,Gdec,Hdec,N,ContStruc);
Gain{1} = [Gain;
         K_dist2_CT,K_dist2_DT];

%% REDUCTION OF THE CONTROL EFFORT CONTROLLER
Gain{2} = [];
% Centralized
ContStruc = ones(N,N);
[K_c_CT_2, rho_c_CT_2, feas_c_CT_2] = LMI_CT_ContrEffort(Atot,Bdec,Cdec,N,ContStruc);
[K_c_DT_2, rho_c_DT_2, feas_c_DT_2] = LMI_DT_ContrEffort(Ftot,Gdec,Hdec,N,ContStruc);
% Decentralized
ContStruc = diag(ones(N,1));
[K_dec_CT_2, rho_dec_CT_2, feas_dec_CT_2] = LMI_CT_ContrEffort(Atot,Bdec,Cdec,N,ContStruc);
[K_dec_DT_2, rho_dec_DT_2, feas_dec_DT_2] = LMI_DT_ContrEffort(Ftot,Gdec,Hdec,N,ContStruc);
% Distributed 1
ContStruc = [1 0
              1 1];
[K_dist1_CT_2, rho_dist1_CT_2, feas_dist1_CT_2] = LMI_CT_ContrEffort(Atot,Bdec,Cdec,N,ContStruc);
[K_dist1_DT_2, rho_dist1_DT_2, feas_dist1_DT_2] = LMI_DT_ContrEffort(Ftot,Gdec,Hdec,N,ContStruc);
% Distributed 2
ContStruc = [1 1
              0 1];
[K_dist2_CT_2, rho_dist2_CT_2, feas_dist2_CT_2] = LMI_CT_ContrEffort(Atot,Bdec,Cdec,N,ContStruc);
[K_dist2_DT_2, rho_dist2_DT_2, feas_dist2_DT_2] = LMI_DT_ContrEffort(Ftot,Gdec,Hdec,N,ContStruc);

%% LIMITATION ON THE SPECTRAL ABSCISSA CONTROLLER

%% LIMITATION ON THE SPECTRAL RADIUS CONTROLLER

%% LIMITATION IN A DISK CONTROLLER

%% LIMITATION IN A REGION CONTROLLER

%% Results
% Gain = table(Gain{1},Gain{2},'VariableNames',{'Stable controller','Reduct of control effort'});
% STABLE CONTROLLER
% Continuos time
clc
disp('Results STABLE CONTROLLER (Continuous-time):')
disp(['-  Centralized: Feasibility=',num2str(feas_c_CT),', rho=',num2str(rho_c_CT),', FM=',num2str(CFM_CT),'.'])
disp(['-  Decentralized: Feasibility=',num2str(feas_dec_CT),', rho=',num2str(rho_dec_CT),', FM=',num2str(DFM_CT),'.'])
disp(['-  Distributed1 (u2-x1): Feasibility=',num2str(feas_dist1_CT),', rho=',num2str(rho_dist1_CT),', FM=',num2str(Dist1FM_CT),'.'])
disp(['-  Distributed2 (u1-x2): Feasibility=',num2str(feas_dist2_CT),', rho=',num2str(rho_dist2_CT),', FM=',num2str(Dist2FM_CT),'.'])
% Discrete time
disp('Results STABLE CONTROLLER (Discrete-time):')
disp(['-  Centralized: Feasibility=',num2str(feas_c_DT),', rho=',num2str(rho_c_DT),', FM=',num2str(CFM_DT),'.'])
disp(['-  Decentralized: Feasibility=',num2str(feas_dec_DT),', rho=',num2str(rho_dec_DT),', FM=',num2str(DFM_DT),'.'])
disp(['-  Distributed (u2-x1): Feasibility=',num2str(feas_dist1_DT),', rho=',num2str(rho_dist1_DT),', FM=',num2str(Dist1FM_DT),'.'])
disp(['-  Distributed (u1-x2): Feasibility=',num2str(feas_dist2_DT),', rho=',num2str(rho_dist2_DT),', FM=',num2str(Dist2FM_DT),'.'])

% REDUCTION OF THE CONTROL EFFORT
% Continuos time
disp('Results REDUCTION OF THE CONTROL EFFORT (Continuous-time):')
disp(['-  Centralized: Feasibility=',num2str(feas_c_CT_2),', rho=',num2str(rho_c_CT_2),', FM=',num2str(CFM_CT),'.'])
disp(['-  Decentralized: Feasibility=',num2str(feas_dec_CT_2),', rho=',num2str(rho_dec_CT_2),', FM=',num2str(DFM_CT),'.'])
disp(['-  Distributed (u2-x1): Feasibility=',num2str(feas_dist1_CT_2),', rho=',num2str(rho_dist1_CT_2),', FM=',num2str(Dist1FM_CT),'.'])
disp(['-  Distributed (u1-x2): Feasibility=',num2str(feas_dist2_CT_2),', rho=',num2str(rho_dist2_CT_2),', FM=',num2str(Dist2FM_CT),'.'])
% % Discrete time
disp('Results REDUCTION OF THE CONTROL EFFORT (Discrete-time):')
disp(['-  Centralized: Feasibility=',num2str(feas_c_DT_2),', rho=',num2str(rho_c_DT_2),', FM=',num2str(CFM_DT),'.'])
disp(['-  Decentralized: Feasibility=',num2str(feas_dec_DT_2),', rho=',num2str(rho_dec_DT_2),', FM=',num2str(DFM_DT),'.'])
disp(['-  Distributed (u2-x1): Feasibility=',num2str(feas_dist1_DT_2),', rho=',num2str(rho_dist1_DT_2),', FM=',num2str(Dist1FM_DT),'.'])
disp(['-  Distributed (u1-x2): Feasibility=',num2str(feas_dist2_DT_2),', rho=',num2str(rho_dist2_DT_2),', FM=',num2str(Dist2FM_DT),'.'])

%% Analysis closed-loop stability STABLE CONTROLLER

% CT closed-loop stability
% Eigenvalues CENTRALIZED 
A_c_cl=Atot+B*K_c_CT;
eig_c_CL_CT=eig(A_c_cl);
% Eigenvalues DECENTRALIZED
A_dec_cl=Atot+B*K_dec_CT;
eig_dec_CL_CT=eig(A_dec_cl);
% Eigenvalues DISTRIBUTED 1
A_dist1_cl=Atot+B*K_dist1_CT;
eig_dist1_CL_CT=eig(A_dist1_cl);
% Eigenvalues DISTRIBUTED 2
A_dist2_cl=Atot+B*K_dist2_CT;
eig_dist2_CL_CT=eig(A_dist2_cl);

% PLOT Continuos time closed loop STABLE CONTROLLER
figure
plot(real(eig_c_CL_CT),imag(eig_c_CL_CT),'xk',real(eig_dec_CL_CT),imag(eig_dec_CL_CT),'xr',real(eig_dist1_CL_CT),imag(eig_dist1_CL_CT),'*g',real(eig_dist2_CL_CT),imag(eig_dist2_CL_CT),'xb','MarkerSize',8,'LineWidth',2);
grid
xlabel('Real Axis','FontSize',18);ylabel('Imaginary Axis','FontSize',18);
title('CT Closed-loop eigenvalue positions BASIC','FontSize',20)
legend('Centralized','Decentralized','Distributed 1','Distributed 2','FontSize',18)

% DT closed-loop stability
% Eigenvalues CENTRALIZED 
F_c_cl=Ftot+G*K_c_DT;
eig_c_CL_DT=eig(F_c_cl);
% Eigenvalues DECENTRALIZED
F_dec_cl=Ftot+G*K_dec_DT;
eig_dec_CL_DT=eig(F_dec_cl);
% Eigenvalues DISTRIBUTED 1
F_dist1_cl=Ftot+G*K_dist1_DT;
eig_dist1_CL_DT=eig(F_dist1_cl);
% Eigenvalues DISTRIBUTED 2
F_dist2_cl=Ftot+G*K_dist2_DT;
eig_dist2_CL_DT=eig(F_dist2_cl);

% PLOT Discrete time closed loop STABLE CONTROLLER
figure
x = [-1:0.01:1];
y = sqrt(1-x.^2);
hold on
plot(zeros(size(x)),x,':k',x,zeros(size(x)),':k',x,y,':k',x,-y,':k',0,0,'.b')
plot(real(eig_c_CL_DT),imag(eig_c_CL_DT),'xk',real(eig_dec_CL_DT),imag(eig_dec_CL_DT),'xr',real(eig_dist1_CL_DT),imag(eig_dist1_CL_DT),'*g',real(eig_dist2_CL_DT),imag(eig_dist2_CL_DT),'xb','MarkerSize',8,'LineWidth',2);
legend('','','','','','Centralized','Decentralized','Distributed 1','Distributed 2','FontSize',18)
title('DT Closed-loop eigenvalue positions BASIC','FontSize',20);
xlabel('Real Axis','FontSize',18);ylabel('Imaginary Axis','FontSize',18);
hold off

%% Analysis closed-loop stability REDUCTION CONTROL EFFORT CONTROLLER

% CT closed-loop stability
% Eigenvalues CENTRALIZED 
A_c_cl_2=Atot+B*K_c_CT_2;
eig_c_CL_CT_2=eig(A_c_cl_2);
% Eigenvalues DECENTRALIZED
A_dec_cl_2=Atot+B*K_dec_CT_2;
eig_dec_CL_CT_2=eig(A_dec_cl_2);
% Eigenvalues DISTRIBUTED 1
A_dist1_cl_2=Atot+B*K_dist1_CT_2;
eig_dist1_CL_CT_2=eig(A_dist1_cl_2);
% Eigenvalues DISTRIBUTED 2
A_dist2_cl_2=Atot+B*K_dist2_CT_2;
eig_dist2_CL_CT_2=eig(A_dist2_cl_2);

% PLOT Continuos time closed loop STABLE CONTROLLER
figure
plot(real(eig_c_CL_CT_2),imag(eig_c_CL_CT_2),'xk',real(eig_dec_CL_CT_2),imag(eig_dec_CL_CT_2),'xr',real(eig_dist1_CL_CT_2),imag(eig_dist1_CL_CT_2),'*g',real(eig_dist2_CL_CT_2),imag(eig_dist2_CL_CT_2),'xb','MarkerSize',8,'LineWidth',2);
grid
xlabel('Real Axis','FontSize',18);ylabel('Imaginary Axis','FontSize',18);
title('CT Closed-loop eigenvalue positions REDUCTION CONTROL EFFORT','FontSize',20)
legend('Centralized','Decentralized','Distributed 1','Distributed 2','FontSize',18)

% DT closed-loop stability
% Eigenvalues CENTRALIZED 
F_c_cl_2=Ftot+G*K_c_DT_2;
eig_c_CL_DT_2=eig(F_c_cl_2);
% Eigenvalues DECENTRALIZED
F_dec_cl_2=Ftot+G*K_dec_DT_2;
eig_dec_CL_DT_2=eig(F_dec_cl_2);
% Eigenvalues DISTRIBUTED 1
F_dist1_cl_2=Ftot+G*K_dist1_DT_2;
eig_dist1_CL_DT_2=eig(F_dist1_cl_2);
% Eigenvalues DISTRIBUTED 2
F_dist2_cl_2=Ftot+G*K_dist2_DT_2;
eig_dist2_CL_DT_2=eig(F_dist2_cl_2);

% PLOT Discrete time closed loop STABLE CONTROLLER
figure
x = [-1:0.01:1];
y = sqrt(1-x.^2);
hold on
plot(zeros(size(x)),x,':k',x,zeros(size(x)),':k',x,y,':k',x,-y,':k',0,0,'.b')
plot(real(eig_c_CL_DT_2),imag(eig_c_CL_DT_2),'xk',real(eig_dec_CL_DT_2),imag(eig_dec_CL_DT_2),'xr',real(eig_dist1_CL_DT_2),imag(eig_dist1_CL_DT_2),'*g',real(eig_dist2_CL_DT_2),imag(eig_dist2_CL_DT_2),'xb','MarkerSize',8,'LineWidth',2);
legend('','','','','','Centralized','Decentralized','Distributed 1','Distributed 2','FontSize',18)
title('DT Closed-loop eigenvalue positions REDUCTION CONTROL EFFORT','FontSize',20);
xlabel('Real Axis','FontSize',18);ylabel('Imaginary Axis','FontSize',18);
hold off

%% Analysis closed-loop stability REDUCTION CONTROL EFFORT CONTROLLER

%% Analysis closed-loop stability LIMITATION ON THE SPECTRAL ABSCISSA CONTROLLER

%% Analysis closed-loop stability LIMITATION ON THE SPECTRAL RADIUS CONTROLLER

%% Analysis closed-loop stability LIMITATION IN A DISK CONTROLLER

%% Analysis closed-loop stability LIMITATION IN A REGION CONTROLLER




%% SIMULATIONs

% Simulation BASIC CONTROLLER

Tfinal=20;
T=[0:0.01:Tfinal];
h=0.01;

x0=[];
x0=[x0;randn(4,1)];

k=0;
for t=T
    k=k+1;
    x_c_CT(:,k)=expm((A_c_cl)*t)*x0;
    x_dec_CT(:,k)=expm((A_dec_cl)*t)*x0;
    x_dist1_CT(:,k)=expm((A_dist1_cl)*t)*x0;
    x_dist2_CT(:,k)=expm((A_dist2_cl)*t)*x0;
end
for k=1:Tfinal/h
    x_c_DT(:,k)=((F_c_cl)^k)*x0;
    x_dec_DT(:,k)=((F_dec_cl)^k)*x0;
    x_dist1_DT(:,k)=((F_dist1_cl)^k)*x0;
    x_dist2_DT(:,k)=((F_dist2_cl)^k)*x0;
end

figure
for i=1:4
    subplot(2,2,i)
    hold on
    grid on
    subtitle(['CT BASIC X_{',num2str(i),'}'])
    plot(T,[x_c_CT((i),:)],'k')
    plot(T,[x_dec_CT((i),:)],'m')
    plot(T,[x_dist1_CT((i),:)],'b')
    plot(T,[x_dist2_CT((i),:)],'r')
    axis([0 T(end) min(x0)-4 max(x0)+4])
end
legend('Centralized','Decentralized','Distributed 1','Distributed 2')

figure
for i=1:4
    subplot(2,2,i)
    hold on
    grid on
    subtitle(['DT BASIC X_{',num2str(i),'}'])
    plot([h:h:Tfinal],[x_c_DT((i),:)],'k.-')
    plot([h:h:Tfinal],[x_dec_DT((i),:)],'m.-')
    plot([h:h:Tfinal],[x_dist1_DT((i),:)],'b.-')
    plot([h:h:Tfinal],[x_dist2_DT((i),:)],'r.-')
    axis([0 T(end) min(x0)-4 max(x0)+4])
end
legend('Centralized','Decentralized','Distributed 1','Distributed 2')

% Simulation REDUCTION CONTROL EFFORT CONTROLLER

k=0;
for t=T
    k=k+1;
    x_c_CT_2(:,k)=expm((A_c_cl_2)*t)*x0;
    x_dec_CT_2(:,k)=expm((A_dec_cl_2)*t)*x0;
    x_dist1_CT_2(:,k)=expm((A_dist1_cl_2)*t)*x0;
    x_dist2_CT_2(:,k)=expm((A_dist2_cl_2)*t)*x0;
end
for k=1:Tfinal/h
    x_c_DT_2(:,k)=((F_c_cl_2)^k)*x0;
    x_dec_DT_2(:,k)=((F_dec_cl_2)^k)*x0;
    x_dist1_DT_2(:,k)=((F_dist1_cl_2)^k)*x0;
    x_dist2_DT_2(:,k)=((F_dist2_cl_2)^k)*x0;
end

figure
for i=1:4
    subplot(2,2,i)
    hold on
    grid on
    subtitle(['CT REDUCED CONTROL EFFORT X_{',num2str(i),'}'])
    plot(T,[x_c_CT_2((i),:)],'k')
    plot(T,[x_dec_CT_2((i),:)],'m')
    plot(T,[x_dist1_CT_2((i),:)],'b')
    plot(T,[x_dist2_CT_2((i),:)],'r')
    axis([0 T(end) min(x0)-4 max(x0)+4])
end
legend('Centralized','Decentralized','Distributed 1','Distributed 2')

figure
for i=1:4
    subplot(2,2,i)
    hold on
    grid on
    subtitle(['DT REDUCED CONTROL EFFORT X_{',num2str(i),'}'])
    plot([h:h:Tfinal],[x_c_DT_2((i),:)],'k.-')
    plot([h:h:Tfinal],[x_dec_DT_2((i),:)],'m.-')
    plot([h:h:Tfinal],[x_dist1_DT_2((i),:)],'b.-')
    plot([h:h:Tfinal],[x_dist2_DT_2((i),:)],'r.-')
    axis([0 T(end) min(x0)-4 max(x0)+4])
end
legend('Centralized','Decentralized','Distributed 1','Distributed 2')














