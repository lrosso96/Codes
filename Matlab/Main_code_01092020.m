% ---------------------------------------------- %
%% Two asset Life Cycle Consumption Model %%
% Based on Two Asset Kinked and Life Cycle codes %
% Continuous Income Proces %
% Author: Lucas Rosso %
% Date: 01-09-2020 %
% ---------------------------------------------- %
clear all; close all; clc;
global chi0 chi1 ga

tic

%Parameters
ga = 2;         % CRRA utility with parameter gamma
ra = 0.08;      % Return illiquid asset
rb = 0.03;      % Return liquid asset
rho = 0.06;     % Discount rate
chi0 = 0.03;       % Linear parameter adj. cost
chi1 = 1;       % Convex parameter adj. cost

xi = 0;         % Fraction of income that is automatically deposited

if ra - 1/chi1 > 0
    disp('Warning: ra - 1/chi1 > 0')
end

sig2 = (0.8)^2;  % sigma^2 O-U
Corr = exp(-0.9);  % persistence -log(Corr) O-U
w = 1;

zmean = exp(sig2/2);

T=75;       %maximum age
N=300;      %number of age steps
dt=T/N;


%simulation parameters
maxit  = 100;     %maximum number of iterations in the HJB loop
crit = 10^(-10);  %criterion HJB loop

%ORNSTEIN-UHLENBECK IN LOGS
the = -log(Corr);

%grids
I = 20;
bmin = 0;
bmax = 10;
b = linspace(bmin,bmax,I)';
db = (bmax-bmin)/(I-1);

J= 20;
amin = 0;
amax = 10;
a = linspace(amin,amax,J);
da = (amax-amin)/(J-1);

Nz=10;         % number of z points 
zmin = 0.75;   % Range z
zmax = 2.5;
z = linspace(zmin,zmax,Nz);   % productivity vector
dz = (zmax-zmin)/(J-1);
dz2 = dz^2;

bb = b*ones(1,J);
aa = ones(I,1)*a;
zz = ones(J,1)*z;

bbb = zeros(I,J,Nz); aaa = zeros(I,J,Nz); zzz = zeros(I,J,Nz);
for nz=1:Nz
    bbb(:,:,nz) = bb;
    aaa(:,:,nz) = aa;
    zzz(:,:,nz) = z(nz);
end

mu = -the.*z.*log(z)+sig2/2*z;   %DRIFT (FROM ITO'S LEMMA)
s2 = sig2.*z.^2;                   %VARIANCE (FROM ITO'S LEMMA)

%Finite difference approximation of the partial derivatives
Vaf = zeros(I,J,Nz);             
Vab = zeros(I,J,Nz);
Vbf = zeros(I,J,Nz);             
Vbb = zeros(I,J,Nz);
Vzf = zeros(I,J,Nz);
Vzb = zeros(I,J,Nz);
Vzz = zeros(I,J,Nz);
c = zeros(I,J,Nz);

%% CONSTRUCT MATRIX Zswitch SUMMARIZING EVOLUTION OF z
yy = - s2/dz2 - mu/dz;
chi =  s2/(2*dz2);
zeta = mu/dz + s2/(2*dz2);


%This will be the upperdiagonal of the matrix Zswitch
updiag=zeros(I*J,1); %This is necessary because of the peculiar way spdiags is defined.
for j=1:Nz
    updiag=[updiag;repmat(zeta(j),I*J,1)];
end

%This will be the center diagonal of the matrix Aswitch
centdiag=repmat(chi(1)+yy(1),I*J,1);
for j=2:Nz-1
    centdiag=[centdiag;repmat(yy(j),I*J,1)];
end
centdiag=[centdiag;repmat(yy(Nz)+zeta(Nz),I*J,1)];

%This will be the lower diagonal of the matrix Aswitch
lowdiag=repmat(chi(2),I*J,1);
for j=3:Nz
    lowdiag=[lowdiag;repmat(chi(j),I*J,1)];
end

%Add up the upper, center, and lower diagonal into a sparse matrix
Zswitch=spdiags(centdiag,0,I*J*Nz,I*J*Nz)+spdiags(lowdiag,-I*J,I*J*Nz,I*J*Nz)+spdiags(updiag,I*J,I*J*Nz,I*J*Nz);
% ------------------------------------

%% CONSTRUCT MATRIX Rswitch SUMMARIZING EVOLUTION OF RISKY ASSET

mu_r = ra;   % rb + adjustment for risk (and liquidity)
sig2_r = (0.04)^2;
s2_r = sig2_r.*z.^2;                   %VARIANCE (FROM ITO'S LEMMA)

%% Preallocation
v = zeros(I,J,Nz,N);
gg = cell(N+1,1);
convergence_criterion = 10^(-5);
up_diag = zeros(I*J,Nz);
low_diag = zeros(I*J,Nz);
cent_diag = zeros(I*J,Nz);
% AAi = cell(Nz,1);
% BBi = cell(Nz,1);
AA = zeros(I*J*Nz,I*J*Nz);
BB = zeros(I*J*Nz,I*J*Nz);
% ------------------------------------

%% HJB Loop
% terminal condition on value function: value of death \approx 0
small_number1 = 10^(-8); small_number2 = 10^(-8);
v_terminal = small_number1*(small_number2 + aaa + bbb).^(1-ga)/(1-ga);

V = v_terminal;

    for n=N:-1:1
        disp(['age = ', num2str(n*dt)])
        v(:,:,:,n)=V;
        
        % LIQUID ASSET
        % forward difference 
        Vbf(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db;
        Vbf(I,:,:) = (w*zz + rb.*bmax).^(-ga); %will never be used, but impose state constraint a<=amax just in case
        % backward difference 
        Vbb(2:I,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db;
        Vbb(1,:,:) = (w*zz + rb.*bmin).^(-ga); %state constraint boundary condition
        % --------------------
        
        % ILLIQUID ASSET
        % forward difference
        Vaf(:,1:J-1,:) = (V(:,2:J,:)-V(:,1:J-1,:))/da;
        % backward difference
        Vab(:,2:J,:) = (V(:,2:J,:)-V(:,1:J-1,:))/da;
        % --------------------
        
        % consumption
        cf = Vbf.^(-1/ga);
        cb = Vbb.^(-1/ga);
        c0 = (1-xi)*w*zzz + rb.*bbb;
        
        % deposits
        dBB = two_asset_kinked_FOC(Vab,Vbb,aaa);
        dFB = two_asset_kinked_FOC(Vab,Vbf,aaa);
        dBF = two_asset_kinked_FOC(Vaf,Vbb,aaa);
        dFF = two_asset_kinked_FOC(Vaf,Vbf,aaa);
        
        %UPWIND SCHEME
        d_B = (dBF>0).*dBF + (dBB<0).*dBB;
        %state constraints at amin and amax
        d_B(:,1,:) = (dBF(:,1,:)>10^(-12)).*dBF(:,1,:); %make sure d>=0 at amin, don't use VaB(:,1,:)
        d_B(:,J,:) = (dBB(:,J,:)<-10^(-12)).*dBB(:,J,:); %make sure d<=0 at amax, don't use VaF(:,J,:)
        d_B(1,1,:)= max(d_B(1,1,:),0);
    %split drift of b and upwind separately
        sc_B = (1-xi)*w*zzz + rb.*bbb - cb;
        sd_B = (-d_B - two_asset_kinked_cost(d_B,aaa));
    
        d_F = (dFF>0).*dFF + (dFB<0).*dFB;
    %state constraints at amin and amax
        d_F(:,1,:) = (dFF(:,1,:)>10^(-12)).*dFF(:,1,:); %make sure d>=0 at amin, don't use VaB(:,1,:)
        d_F(:,J,:) = (dFB(:,J,:)<-10^(-12)).*dFB(:,J,:); %make sure d<=0 at amax, don't use VaF(:,J,:)
    
    %split drift of b and upwind separately
        sc_F = (1-xi)*w*zzz + rb.*bbb - cf;
        sd_F = (-d_F - two_asset_kinked_cost(d_F,aaa));
        sd_F(I,:,:) = min(sd_F(I,:,:),0);
        
        Ic_B = (sc_B < -10^(-12));
        Ic_F = (sc_F > 10^(-12)).*(1- Ic_B);
        Ic_0 = 1 - Ic_F - Ic_B;
    
        Id_F = (sd_F > 10^(-12));
        Id_B = (sd_B < -10^(-12)).*(1- Id_F);
        Id_B(1,:,:)=0;
        Id_F(I,:,:) = 0; Id_B(I,:,:) = 1; %don't use VbF at bmax so as not to pick up articial state constraint
        Id_0 = 1 - Id_F - Id_B;
        
        c = cf.*Ic_F + cb.*Ic_B + c0.*Ic_0;
        u = c.^(1-ga)/(1-ga);
        
    %CONSTRUCT MATRIX BB SUMMARING EVOLUTION OF b
    X = -Ic_B.*sc_B/db -Id_B.*sd_B/db;
    Y = (Ic_B.*sc_B - Ic_F.*sc_F)/db + (Id_B.*sd_B - Id_F.*sd_F)/db;
    Z = Ic_F.*sc_F/db + Id_F.*sd_F/db;
    
    for i = 1:Nz
        center_diag(:,i) = reshape(Y(:,:,i),I*J,1);
    end

    low_diag(1:I-1,:) = X(2:I,1,:);
    up_diag(2:I,:) = Z(1:I-1,1,:);
    for j = 2:J
        low_diag(1:j*I,:) = [low_diag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Nz)];
        up_diag(1:j*I,:) = [up_diag(1:(j-1)*I,:);zeros(1,Nz);squeeze(Z(1:I-1,j,:))];
    end
    
    for nz=1:Nz
        BBi{nz}=spdiags(center_diag(:,nz),0,I*J,I*J)+spdiags([up_diag(:,nz);0],1,I*J,I*J)+spdiags([low_diag(:,nz);0],-1,I*J,I*J);
    end
    
      BB = [BBi{1}, sparse(I*J,(Nz-1)*I*J)];
 
    for nz=2:Nz
        BB = [BB; sparse(I*J,(nz-1)*I*J), BBi{nz},  sparse(I*J,(Nz-nz)*I*J)]; 
    end
    
    % -------------------------------------------------
    
    %CONSTRUCT MATRIX AA SUMMARIZING EVOLUTION OF a
    dB = Id_B.*dBB + Id_F.*dFB;
    dF = Id_B.*dBF + Id_F.*dFF;
    d = Id_B.*d_B + Id_F.*d_F; % deposits
    MB = min(dB,0);
    MF = max(dF,0) + xi*w*zzz + ra.*aaa;
    MB(:,J,:) = xi*w*zzz(:,J,:) + dB(:,J,:) + ra.*amax; %this is hopefully negative
    MF(:,J,:) = 0;
    chi = -MB/da;
    yy =  (MB - MF)/da;
    zeta = MF/da;
    
    %MATRIX AAi
    for nz=1:Nz
        %This will be the upperdiagonal of the matrix AAi
        AAupdiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined.
        for j=1:J
            AAupdiag=[AAupdiag;zeta(:,j,nz)];
        end
        
        %This will be the center diagonal of the matrix AAi
        AAcentdiag= yy(:,1,nz);
        for j=2:J-1
            AAcentdiag=[AAcentdiag;yy(:,j,nz)];
        end
        AAcentdiag=[AAcentdiag;yy(:,J,nz)];
        
        %This will be the lower diagonal of the matrix AAi
        AAlowdiag=chi(:,2,nz);
        for j=3:J
            AAlowdiag=[AAlowdiag;chi(:,j,nz)];
        end
        
        %Add up the upper, center, and lower diagonal into a sparse matrix
        AAi{nz} = spdiags(AAcentdiag,0,I*J,I*J)+spdiags(AAlowdiag,-I,I*J,I*J)+spdiags(AAupdiag,I,I*J,I*J);
        
    end
    
      AA = [AAi{1}, sparse(I*J,(Nz-1)*I*J)];
 
    for nz=2:Nz
        AA = [AA; sparse(I*J,(nz-1)*I*J), AAi{nz},  sparse(I*J,(Nz-nz)*I*J)]; 
    end    
% ----------------------------------------

%% Transition Matriz
    A = AA + BB + Zswitch;
    
% Check all rows sum 0
error = full(sum(A,2));
check = abs(error)>1e-9;
check_ = sum(check);

if max(abs(sum(A,2)))>10^(-9)
       disp('Improper Transition Matrix')
       break
end
% -----------------------------------------
           %%Note the syntax for the cell array
        A_t{n} = A;
        B = (1/dt + rho)*speye(I*J*Nz) - A;
        
        u_stacked = reshape(u,I*J*Nz,1);
        V_stacked = reshape(V,I*J*Nz,1);
        
        b = u_stacked + V_stacked/dt;
        V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
        
        V = reshape(V_stacked,I,J,Nz);
        c_t{n} = c;
        d_t{n} = d;
        sb_t{n} = w*zzz + rb.*bbb - d - two_asset_kinked_cost(d,aaa) - c;   
        sa_t{n} = d + xi*w*zzz + ra.*aaa;
        share_t{n} = (aaa + sa_t{n})./(aaa +sa_t{n} + bbb + sb_t{n}); % risky asset share
        part_t{n} = (share_t{n} >0); % participation
        
    end
toc
