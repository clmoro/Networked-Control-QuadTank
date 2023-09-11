function [K,rho,feas]=LMI_CT_Disk(A,B,C,N,ContStruc)
% Computes, using LMIs, the distributed "state feedback" control law for the continuous-time system, with reference to the control
% information structure specified by 'ContStruc'.
%
% Inputs:
% - A: system matrix.
% - B: input matrices (i.e., B{1},..., B{N} are the input matrices of the decomposed system, one for each channel).
% - C: output matrices  (i.e., C{1},..., C{N} are the output matrices of the decomposed system, one for each channel, where [Cdec{1}',...,
% Cdec{N}']=I).
% - N: number of subsystems.
% - ContStruc: NxN matrix that specifies the information structure
% constraints (ContStruc(i,j)=1 if communication is allowed between channel
% j to channel i, ContStruc(i,j)=0 otherwise).
%
% Output:
% - K: structured control gain
% - rho: spectral abscissa of matrix (A+B*K) - note that [C{1}',...,
% C{N}']=I
% - feas: feasibility of the LMI problem (=0 if yes)

Btot=[];
for i=1:N
    m(i)=size(B{i},2);
    n(i)=size(C{i},1);
    Btot=[Btot,B{i}];
end
ntot=size(A,1);
mtot=sum(m);

yalmip clear

if ContStruc==ones(N,N)
    % Centralized design
    P=sdpvar(ntot);
    L=sdpvar(mtot,ntot);
else
    % Dentralized/distributed design
    P=[];
    L=sdpvar(mtot,ntot);
    minc=0;
    for i=1:N
        P=blkdiag(P,sdpvar(n(i)));
        ninc=0;
        for j=1:N
            if ContStruc(i,j)==0
                L(minc+1:minc+m(i),ninc+1:ninc+n(j))=zeros(m(i),n(j));
            end
            ninc=ninc+n(j);
        end
        minc=minc+m(i);
    end
end

alpha = 8; % We constraints Re(eigs(A))<alpha
rho = 1; % We constraints all eigs to have radius < rho
% Dimensions of the matrices: Kx(mtot,ntot), L(mtot,ntot), Y(ntot,ntot).
LMIconstr = [[(rho^2-alpha^2)*P-A*P*A'-A*L'*B'-B*L*A'-alpha*(P*A'+A*P+L'*B'+B*L),B*L;
                L'*B',P]>=-1e-2*eye(2*ntot)];
% options=sdpsettings('solver','sedumi');
J=optimize(LMIconstr);%,[],options);
feas=J.problem;
L=double(L);
P=double(P);

K=L/P;
rho=max(real(eig(A+Btot*K)));
