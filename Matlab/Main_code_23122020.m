% ---------------------------------------------- %
%% Two asset Portfolio Allocation Model %%
% Based on Liquid-Illiquid LCP code (from Benjamin Moll) %
% Stochastic return on illiquid asset
% Author: Lucas Rosso %
% Date: 27-12-2020 %
% Extremely Preliminar %
% ---------------------------------------------- %

clear all; close all; clc;

tic 

%% Parameters and initial objetcts

UseNoAdjustmentAsGuess = 1;
AssumeConcavity = 0; % this is relevant in the part with adjustment
BilinW = 1; % when all neighboring points of (b',a') are non-adjustment points, use bilinear weights in dividing mass to be rebalanced 
            % across them in the KF algorithm

% Parameters:
s     = 5;         % risk aversion in CRRA utility
rbPos = 0.02;      % return on liquid asset (e.g. checking account)
rbNeg = 0.10;
rho   = 0.042;     % discount rate (beta = 0.96)
zeta  = 1/60;      % death rate  1/180;     % death rate (if time is quarterly, 45 year lifespan)

rho = rho + zeta;   

ka    = 0.2;       % cost of adjustment

% Liquid Asset
I    = 100; 
bmin    = 0;
bmax    = 10;
b = linspace(bmin,bmax,I)'; % column vector
db = b(2)-b(1);

% Illiquid asset 
J   = 100; 
amin   = 0;
amax   = 10; 
a  = linspace(amin,amax,J); % row vector
da = a(2)-a(1);

% Stock returns
N_ra=10;         % number of r_a points 
rmin = -0.06;   % Range z
rmax = 0.18;
r_a = linspace(rmin,rmax,N_ra);   % productivity vector
dr_a = (rmax-rmin)/(N_ra-1);
dr_a2 = dr_a^2;

% Drift and Variance of the stock return
mu = 0.06;  
s2 = 0.18^2;                   

% Poisson shocks:
Nz = 2; % number of exog. states
w = 1;
z      = [0.8*w,1.2*w]; % state values (row vector)
la = [1/15,1/15]; % Poisson intensities

% Create mesh: b in 1st dimension, a in 2nd dimension, z in 3rd dimension,
% r_a in 4th dimension
bb   = b*ones(1,J);   % repeat b column vector J times in the 2nd dimension
aa   = ones(I,1)*a;   % repeat a row vector I times in the 1st dimension
zz   = ones(J,1)*z;   % repeat z row vector J times in the 1st dimension (will then take the 2 columns with 1st and 2nd state values separately)
rr_a = r_a'*ones(1,Nz); % repeat r_a row vector Z times in the 1st dimension

bbb = zeros(I,J,N_ra); aaa = zeros(I,J,N_ra); zzz = zeros(J,N_ra,Nz);

for n_ra = 1:N_ra
    bbb(:,:,n_ra) = bb;
    aaa(:,:,n_ra) = aa;
    zzz(:,n_ra,:) = zz;
end

rrr_a = zeros(J,N_ra,Nz);
for j = 1:J
    rrr_a(j,:,:) = rr_a;
end

zzzz = zeros(I,J,N_ra,Nz);
rrrr_a = zeros(I,J,N_ra,Nz);
for i= 1:I
    zzzz(i,:,:,:) = zzz;
    rrrr_a(i,:,:,:) = rrr_a;
end

aaaa = zeros(I,J,N_ra,Nz); bbbb = zeros(I,J,N_ra,Nz);
for nz = 1:Nz
    aaaa(:,:,:,nz) = aaa;
    bbbb(:,:,:,nz) = bbb;
end
    
%allow for differential borrowing and lending rates (only matters when bmin<0)
Rb = rbPos.*(bbbb>=0) + rbNeg.*(bbbb<0); % mesh of returns on liquid asset to allow for different borrowing rate if borrowing is allowed

L = I*J*N_ra*Nz; % total number of nodes

% Adjustment decision set-up:
Nx      = 600; % number of grid points for cash-in-hand when adjusting
NaP     = 600; % number of grid points for a' when adjusting
xmin      = amin + bmin + ka; % minimum cash-in-hand
xmax      = amax + bmax; % maximum
x     = linspace(xmin,xmax,Nx);

sumGrid_aux   = aa + bb; % mesh of cash-in-hand from (b, a) nodes
sumGrid_aux   = sumGrid_aux(:); % transform to column vector

sumGrid = sumGrid_aux;

% Iteration set-up:
Delta     = 1000; % step size (for time) in the implicit method (as opposed to db and db which are steps in the state dimensions)
maxit     = 1000;
tol       = 1e-6;
iter      = 0;
dist      = zeros(maxit,1); % hold errors during iteration

% Initial guess
for n_ra=1:N_ra
V0(:,:,n_ra,1) = (1-s)^(-1)*(z(1) + rbPos.*bb).^(1-s)/rho; % assuming u(c) = c^(1 - s)/(1 - s). This is similar to value of "staying put" u(z + rb*b)/rho                                                          
V0(:,:,n_ra,2) = (1-s)^(-1)*(z(2) + rbPos.*bb).^(1-s)/rho;
end
v         = V0;  

% tau = 15; % from two_assets_kinked.m, this is justified as: if ra>>rb, impose tax on ra*a at high a, otherwise some households accumulate infinite illiquid wealth
%           % (not needed if ra is close to - or less than - rb). It implies aDrift = ra*a*(1 - (a/(amax * 0.999))^(tau - 1))
% tau0 = ra.*(amax*.999)^(1-tau); % in two_assets_kinked.m the multiplier is 1.33 instead of 0.999. That multiplier makes the tax much less aggressive, as with 
%                                 % 0.999 the drift of a goes down to slightly negative towards amax, while with 1.33 the drift becomes just slightly concave near amax
% T = tau0.*a.^tau;
% aDrift = ra.*a - T;
% plot(a,aDrift,a,zeros(J,1),a,ra.*a)

aDrift = r_a'*a;

%% Build matrix Aswitch summarizing evolution of (a,r_a,z) % i.e. states which are exogenous in the no adjustment case: risky asset, its return and the income process
chi = -min(aDrift,0)/da; % this is analogous to X, Y, Z constructed below for b (but call chi instead of x since x is cash-in-hand). This (and below) has 
                         % similarities with p. 16/17 in Achdou et al. 2017 - Online Appendix
yy  = - max(aDrift,0)/da + min(aDrift,0)/da;
zeta= max(aDrift,0)/da;

%This will be the upperdiagonal of the A_switch
updiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined (since this will be the Ith upper diagonal, the first I elements are dropped)
for n_ra = 1:N_ra
    for j=1:J-1
    updiag=[updiag;repmat(zeta(n_ra,j),I,1)]; % this is analogous to zeta(j) * ones(I, 1) i.e. repeating zeta(j) I times in a vector
    end
    updiag = [updiag; zeros(I,1)]; % at each break I*J (when r_a changes), need I zeros for zeta in order to generate the correct discretization.
end 

%This will be the center diagonal of the A_switch
    centdiag=repmat(chi(1,1)+yy(1,1),I,1);
    for j = 2:J-1    
        centdiag=[centdiag;repmat(yy(1,j),I,1)];
    end
    centdiag=[centdiag;repmat(yy(1,J)+zeta(1,J),I,1)];
        
for n_ra = 2:N_ra
    centdiag=[centdiag;repmat(chi(n_ra,1)+yy(n_ra,1),I,1)];
    for j=2:J-1
    centdiag=[centdiag;repmat(yy(n_ra,j),I,1)];
    end
    centdiag=[centdiag;repmat(yy(n_ra,J)+zeta(n_ra,J),I,1)];
end

lowdiag = repmat(chi(1,2),I,1);
for j = 3:J
    lowdiag = [lowdiag;repmat(chi(1,j),I,1)];
end

for n_ra = 2:N_ra
    lowdiag = [lowdiag; zeros(I,1)];
    for j = 2:J
        lowdiag = [lowdiag; repmat(chi(n_ra,j),I,1)];
    end
end
%Add up the upper, center, and lower diagonal into a sparse matrix
 Aaux=spdiags(centdiag,0,I*J*N_ra,I*J*N_ra)+spdiags(lowdiag,-I,I*J*N_ra,I*J*N_ra)+spdiags(updiag,I,I*J*N_ra,I*J*N_ra);

Aswitch = [Aaux - speye(I*J*N_ra)*la(1), speye(I*J*N_ra)*la(1);speye(I*J*N_ra)*la(2),Aaux - speye(I*J*N_ra)*la(2)]; % this adds the transitions between z1 and z2 (keep in mind they have
                                                                                                                % to balance along the row and sum to 0)
                                                                                                                
%% CONSTRUCT MATRIX Rswitch SUMMARIZING MISSING ELEMENTS OF r_a

yy = - s2/dr_a2 - mu/dr_a;
chi =  s2/(2*dr_a2);
zeta = mu/dr_a + s2/(2*dr_a2);

%This will be the upperdiagonal of the matrix Rswitch
updiag=zeros(I*J,1); %This is necessary because of the peculiar way spdiags is defined. First I*J elements will be dropped.

for n_ra=1:N_ra
    updiag=[updiag;repmat(zeta,I*J,1)];
end

%This will be the center diagonal of the matrix Rswitch
centdiag=repmat(chi+yy,I*J,1);
for n_ra=2:N_ra-1
    centdiag=[centdiag;repmat(yy,I*J,1)];
end
centdiag=[centdiag;repmat(yy+zeta,I*J,1)];

%This will be the lower diagonal of the matrix Rswitch. 
lowdiag=repmat(chi,I*J,1);
for n_ra=3:N_ra
    lowdiag=[lowdiag;repmat(chi,I*J,1)];
end

%Add up the upper, center, and lower diagonal into a sparse matrix (particular case with z two-state poisson process)
Rswitch=[spdiags(centdiag,0,I*J*N_ra,I*J*N_ra)+spdiags(lowdiag,-I*J,I*J*N_ra,I*J*N_ra)+spdiags(updiag,I*J,I*J*N_ra,I*J*N_ra), 0*speye(I*J*N_ra); ...
         0*speye(I*J*N_ra), spdiags(centdiag,0,I*J*N_ra,I*J*N_ra)+spdiags(lowdiag,-I*J,I*J*N_ra,I*J*N_ra)+spdiags(updiag,I*J,I*J*N_ra,I*J*N_ra)];


%% SOLVE WITHOUT ADJUSTMENT
if UseNoAdjustmentAsGuess == 1

for n=1:maxit
    V = v;
    
    % forward difference for the derivative of v wrt b
    Vbf(1:I-1,:,:,:) = (V(2:I,:,:,:)-V(1:I-1,:,:,:))/db; % remember b changes along the 1st dimension (column vector)
    Vbf(I,:,:,:) = (zzzz(I,:,:,:) + Rb(I,:,:,:).*bmax).^(-s); %will never be used, but impose state constraint b<=bmax just in case
                                                           % this is v'b(I,j,k)_F = u'(z(k) + rb(I) * b(I))
    % backward difference
    Vbb(2:I,:,:,:) = (V(2:I,:,:,:)-V(1:I-1,:,:,:))/db;
    Vbb(1,:,:,:) = (zzzz(1,:,:,:) + Rb(1,:,:,:).*bmin).^(-s); %state constraint boundary condition: v'b(1,j,k)_B = u'(z(k) + rb(1) * b(1))

    Vbf = max(Vbf,10^(-6)); % this is to avoid dealing with too small numbers in parts where v is almost flat wrt b
    Vbb = max(Vbb,10^(-6));
    
     %consumption and savings with forward difference
     cf = Vbf.^(-1/s); % from FOC
     sf = zzzz + Rb.*bbbb - cf;
     %consumption and savings with backward difference
     cb = Vbb.^(-1/s);
     sb = zzzz + Rb.*bbbb - cb;
     %consumption and derivative of value function at steady state
     c0 = zzzz + Rb.*bbbb; % this is where s=0 i.e. db = 0
     
     % make a choice of forward or backward differences based on the sign of the drift (upwind scheme)
     Ib = (sb < 0); %negative drift --> backward difference. This is where we assume concavity, because we don't check sf <= 0 as well (see *** reference)
     If = (sf > 0).*(1-Ib); %positive drift --> forward difference. Note that the (1-Ib) makes sure we take only points where also sb>=0, so Ib and If are 
                               % disjoint. Concavity of v wrt b would imply that, because with concavity sf <= sb, but this enforces it in a setting
                               % where there could be some numerical instability
     I0 = (1-If-Ib); %at steady state. These have to be points where sb>=0 & sf<=0. If we didn't impose disjointness in the previous step, this could = -1 where
                        % (sb < 0 & sf > 0) and we would need to split these points (see *** reference). Instead we assume concavity i.e. these points don't arise 

     c = cf.*If + cb.*Ib + c0.*I0;
     u = c.^(1-s)/(1-s);
     
     %CONSTRUCT MATRIX AA (see p. 5 Achdou et al. 2017 - Online Appendix)
     X = -Ib.*sb/db; % lower
     Y = -If.*sf/db + Ib.*sb/db; % center
     Z = If.*sf/db; % upper    

% upper diagonal of AA
updiag = 0; % this element will be dropped when using spdiags
for n_z = 1:Nz
    for n_ra = 1:N_ra
        for j = 1:J
            for i = 1:I-1
                updiag = [updiag; Z(i,j,n_ra,n_z)];
            end
            updiag = [updiag; 0];
        end
    end
end

% center diagonal of AA
centdiag = 0; % only to generate matrix
for n_z = 1:Nz
    for n_ra = 1:N_ra
        for j = 1:J
            centdiag = [centdiag; X(1,j,n_ra,n_z)+Y(1,j,n_ra,n_z)];
            for i = 2:I-1
                centdiag = [centdiag; Y(i,j,n_ra,n_z)];
            end
            centdiag = [centdiag; Y(I,j,n_ra,n_z)+Z(I,j,n_ra,n_z)];
        end
    end
end
centdiag = centdiag(2:length(centdiag));

% lower diagonal of AA
lowdiag = 0;
for n_z = 1:Nz
    for n_ra = 1:N_ra
        for j = 1:J
            lowdiag = [lowdiag; 0];
            for i= 2:I
                lowdiag = [lowdiag; X(i,j,n_ra,n_z)];
            end
        end
    end
end
lowdiag = lowdiag(3:length(lowdiag));

AA = spdiags(centdiag,0,I*J*N_ra*Nz,I*J*N_ra*Nz) + spdiags(updiag,1,I*J*N_ra*Nz,I*J*N_ra*Nz)+spdiags(lowdiag,-1,I*J*N_ra*Nz,I*J*N_ra*Nz); % notice that when filling the upper diagonal
                                                                                         % the first element of the vector is dropped. When filling the lower diagonal
                                                                                         % the last element is dropped. This is consistent with how lowdiag and updiag
                                                                                         % are constructed initially

  A = AA + Aswitch + Rswitch; % adds the transition of "exogenous states" (in the no adjustment case, risky asset evolves endogenously)

 % checking that rows in A sum to zero
error = full(sum(A,2));
check = abs(error)>1e-9;
check_ = sum(check);

if max(abs(sum(A,2)))>10^(-9)
disp('Improper Transition Matrix')
max(abs(sum(A,2)))
how_many = check_ 
break
end
 
  B = (1/Delta + rho)*speye(I*J*N_ra*Nz) - A; % this and the next 3 steps implement eq. 15 p.6 Achdou et al. 2017 - Online Appendix

  u_stacked = [reshape(u(:,:,:,1),I*J*N_ra,1);reshape(u(:,:,:,2),I*J*N_ra,1)]; % stack into vector I*J*N_ra*Nz
  V_stacked = [reshape(V(:,:,:,1),I*J*N_ra,1);reshape(V(:,:,:,2),I*J*N_ra,1)];

  vec = u_stacked + V_stacked/Delta;

  %IMPLICIT UPDATING
  V_stacked_12 = B\vec; %SOLVE SYSTEM OF EQUATIONS

  % bring back into mesh form
  V(:,:,:,1) = reshape(V_stacked_12(1:I*J*N_ra),I,J,N_ra); 
  V(:,:,:,2) = reshape(V_stacked_12(I*J*N_ra+1:I*J*N_ra*Nz),I,J,N_ra);

  Vchange = V - v;
  v = V;
  
  dist(n) = max(max(max(max(abs(Vchange)))));
  disp(['Initial guess, no adjustment, iter ', int2str(n), ', dist ' num2str(dist(n)) ]);
      if dist(n)<tol
         disp(['Value Function Converged, Iteration = ' int2str(n)]);
         break
      end
end
end

%% SOLVE WITH ADJUSTMENT

for n=1:maxit

    % except for Hf, Hb, H0 and the AssumeConcavity == 0, this part (up to adjustment decision for x grid) is the same as above without adjustment
    V = v;
    % forward difference for the derivative of v wrt b
    Vbf(1:I-1,:,:,:) = (V(2:I,:,:,:)-V(1:I-1,:,:,:))/db; % remember b changes along the 1st dimension (column vector)
    Vbf(I,:,:,:) = (zzzz(I,:,:,:) + Rb(I,:,:,:).*bmax).^(-s); %will never be used, but impose state constraint b<=bmax just in case
                                                           % this is v'b(I,j,k)_F = u'(z(k) + rb(I) * b(I))
    % backward difference
    Vbb(2:I,:,:,:) = (V(2:I,:,:,:)-V(1:I-1,:,:,:))/db;
    Vbb(1,:,:,:) = (zzzz(1,:,:,:) + Rb(1,:,:,:).*bmin).^(-s); %state constraint boundary condition: v'b(1,j,k)_B = u'(z(k) + rb(1) * b(1))

    Vbf = max(Vbf,10^(-6)); % this is to avoid dealing with too small numbers in parts where v is almost flat wrt b
    Vbb = max(Vbb,10^(-6));
    
    %consumption and savings with forward difference
    cf = Vbf.^(-1/s);
    sf = zzzz + Rb.*bbbb - cf;
    Hf = cf.^(1-s)/(1-s) +Vbf.*sf;
    
    %consumption and savings with backward difference
    cb = Vbb.^(-1/s);
    sb = zzzz + Rb.*bbbb - cb;
    Hb = cb.^(1-s)/(1-s) +Vbb.*sb;
    
    %consumption and derivative of value function at steady state
    c0 = zzzz + Rb.*bbbb;
    %hamiltonians;
    
    % make a choice of forward or backward differences based on the sign of the drift
    if AssumeConcavity == 1 % same as above
        Ib = (sb < 0); %negative drift --> backward difference
        If = (sf > 0).*(1-Ib); %positive drift --> forward difference
        I0 = (1-If-Ib); %at steady state
    else % (*** reference)
        I0 = (1-(sf>0)) .* (1-(sb<0)); % points where sf <= 0 & sb >= 0
        Iunique = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0); % points where (sf <= 0 & sb < 0) U (sf > 0 & sb >= 0). "Unique" because we know here we can easily take sb
                                                           % over (sf <= 0 & sb < 0) and sf over (sf > 0 & sb >= 0)
        Iboth = (sb<0).*(sf>0);
        Ib = Iunique.*(sb<0) + Iboth.*(Hb>Hf); % split the points where (sb < 0 & sf > 0) based on u(cb) + Vbb * sb </> u(cf) + Vbf * sf: if > then backward. Note
                                               % that the other elements in the HJB are the same between backward and forward
        If = Iunique.*(sf>0) + Iboth.*(Hb<=Hf);
    end
    
    c = cf.*If + cb.*Ib + c0.*I0;
    u = c.^(1-s)/(1-s);
    
    %CONSTRUCT MATRIX
    X = -Ib.*sb/db;
    Y = -If.*sf/db + Ib.*sb/db;
    Z = If.*sf/db;
    
% upper diagonal of AA
updiag = 0; % this element will be dropped when using spdiags
for n_z = 1:Nz
    for n_ra = 1:N_ra
        for j = 1:J
            for i = 1:I-1
                updiag = [updiag; Z(i,j,n_ra,n_z)];
            end
            updiag = [updiag; 0];
        end
    end
end

% center diagonal of AA
centdiag = 0; % only to generate matrix
for n_z = 1:Nz
    for n_ra = 1:N_ra
        for j = 1:J
            centdiag = [centdiag; X(1,j,n_ra,n_z)+Y(1,j,n_ra,n_z)];
            for i = 2:I-1
                centdiag = [centdiag; Y(i,j,n_ra,n_z)];
            end
            centdiag = [centdiag; Y(I,j,n_ra,n_z)+Z(I,j,n_ra,n_z)];
        end
    end
end
centdiag = centdiag(2:length(centdiag));

% lower diagonal of AA
lowdiag = 0;
for n_z = 1:Nz
    for n_ra = 1:N_ra
        for j = 1:J
            lowdiag = [lowdiag; 0];
            for i= 2:I
                lowdiag = [lowdiag; X(i,j,n_ra,n_z)];
            end
        end
    end
end
lowdiag = lowdiag(3:length(lowdiag));

AA = spdiags(centdiag,0,I*J*N_ra*Nz,I*J*N_ra*Nz) + spdiags(updiag,1,I*J*N_ra*Nz,I*J*N_ra*Nz)+spdiags(lowdiag,-1,I*J*N_ra*Nz,I*J*N_ra*Nz); % notice that when filling the upper diagonal
                                                                                         % the first element of the vector is dropped. When filling the lower diagonal
                                                                                         % the last element is dropped. This is consistent with how lowdiag and updiag
                                                                                         % are constructed initially

A = AA + Aswitch + Rswitch; 

 % checking that rows in A sum to zero
error = full(sum(A,2));
check = abs(error)>1e-9;
check_ = sum(check);

if max(abs(sum(A,2)))>10^(-9)
disp('Improper Transition Matrix')
max(abs(sum(A,2)))
how_many = check_ 
break
end

%%%%% Adjustment decision for x grid: this constructs objects on the predetermined grid for x, before interpolating at the x values corresponding to nodes (b, a)
vstarAux = zeros(N_ra,Nz,Nx); %. The Aux elements are those defined on the grid for x
bAdjAux  = zeros(N_ra,Nz,Nx);
aAdjAux  = zeros(N_ra,Nz,Nx);
vAdj = zeros(N_ra,Nz,NaP); % for each x in the grid, search over a grid for a' - backing out b' = x - a' - ka (see ** reference) - to find the maximum. This vadj
                         % holds the values over such grid for a'
aPmin = amin;

for n_ra = 1:N_ra
G1{n_ra} = griddedInterpolant(bb,aa,V(:,:,n_ra,1)); % this is the interpolant for v over the nodes (b, a) at z(1)
G2{n_ra} = griddedInterpolant(bb,aa,V(:,:,n_ra,2));
end

for n_ra = 1:N_ra
    for i = 1:Nx
        aPmax = min(x(i) - ka - bmin,amax); %CONSTRAINT a'<=aMax, taking into account the minimum value for b
        aP = linspace(aPmin,aPmax,NaP);
        aP = max(aP,x(i)-ka-bmax); %CONSTRAINT b'<=bMax
        bP = x(i) - ka - aP; % (** reference)
        vAdj(n_ra,1,:) = G1{n_ra}(bP',aP'); % collect values over the grid for a' (and implied grid for b', given x)
        vAdj(n_ra,2,:) = G2{n_ra}(bP',aP');

        [VstarAux(n_ra,1,i),idx1] = max(vAdj(n_ra,1,:)); % find maximum and argmax, and store maximum
        [VstarAux(n_ra,2,i),idx2] = max(vAdj(n_ra,2,:));
        aAdjAux(n_ra,1,i) = aP(idx1); aAdjAux(n_ra,2,i) = aP(idx2); % store argmax
        bAdjAux(n_ra,1,i) = bP(idx1); bAdjAux(n_ra,2,i) = bP(idx2);
    end
end
    
%%%%% Interpolate for all possible values of a+b:
Vstar = 0;
for n_z = 1:Nz
    for n_ra = 1:N_ra
        Vstar = [Vstar; lininterp1(x',squeeze(VstarAux(n_ra,n_z,:))',sumGrid)];
    end
end
Vstar = Vstar(2:length(Vstar));
      
% SOLVE USING LCP  
    
  B = (1/Delta + rho)*speye(I*J*N_ra*Nz) - A;

  u_stacked = [reshape(u(:,:,:,1),I*J*N_ra,1);reshape(u(:,:,:,2),I*J*N_ra,1)]; % stack into vector I*J*N_ra*Nz
  V_stacked = [reshape(V(:,:,:,1),I*J*N_ra,1);reshape(V(:,:,:,2),I*J*N_ra,1)];

  vec = u_stacked + V_stacked/Delta;
  
  q = -vec + B*Vstar;
    
  %using Yuval Tassa's Newton-based LCP solver, download from http://www.mathworks.com/matlabcentral/fileexchange/20952
  z0 = V_stacked-Vstar; lbnd = zeros(I*J*N_ra*Nz,1); ubnd = Inf*ones(I*J*N_ra*Nz,1); % set initial value z0 where we guess v(n+1)=v(n) and lower/upper bound for solution z
  z = LCP(B,q,lbnd,ubnd,z0,0); % last argument is display of iteration data
   
  LCP_error = max(abs(z.*(B*z + q))); % check the accuracy of the LCP solution (should all be 0)
    if LCP_error > 10^(-6)
        disp('LCP not solved')
        %break
    end
    
  V_stacked = z+Vstar; %calculate value function
    
  % bring back into mesh form (exactly as in no adjustment case)
  V(:,:,:,1) = reshape(V_stacked(1:I*J*N_ra),I,J,N_ra); 
  V(:,:,:,2) = reshape(V_stacked(I*J*N_ra+1:I*J*N_ra*Nz),I,J,N_ra);
  
  Vchange = V - v;
  v = V;

    dist(n) = max(max(max(max(abs(Vchange)))));
    disp(['With adjustment, adjustment, iter ', int2str(n), ', dist ' num2str(dist(n)) ]);
    if dist(n)<tol
        disp(['Value Function Converged, Iteration = ' int2str(n)]);
        break
    end

end 

toc

%plot(dist(n-20:n))
ss = zzzz + Rb.*bbbb - c; % solution for savings in liquid asset b

adj = V_stacked==Vstar; %INDICATOR FOR ADJUSTING
adj = (abs(V_stacked - Vstar)<10^(-6)); % allow for some tolerance in defining solution for adjustment decision

%AJDUSTMENT REGION
for n_ra = 1:N_ra
adjRegion1{n_ra} = reshape(adj(I*J*(n_ra-1)+1:I*J*n_ra),I,J).*ones(I,J); % multiplication by ones(I,J) transforms logicals into floats
adjRegion2{n_ra} = reshape(adj(I*J*N_ra+I*J*(n_ra-1)+1:I*J*N_ra+I*J*n_ra),I,J).*ones(I,J);
adjRegion(:,:,n_ra,1)=adjRegion1{n_ra}; adjRegion(:,:,n_ra,2)=adjRegion2{n_ra}; % this reshapes the adj. decision into (I, J, N_ra, Nz)

%ADJUSTMENT TARGETS
bAdj1{n_ra} = lininterp1(x',squeeze(bAdjAux(n_ra,1,:))',sumGrid); % this defines the solution for (b',a') conditional on adjusting on the original grid (b, a) from the grid on x
bAdj2{n_ra} = lininterp1(x',squeeze(bAdjAux(n_ra,2,:))',sumGrid);
aAdj1{n_ra} = lininterp1(x',squeeze(aAdjAux(n_ra,1,:))',sumGrid);
aAdj2{n_ra} = lininterp1(x',squeeze(aAdjAux(n_ra,2,:))',sumGrid);

bAdj(:,:,n_ra,1) = reshape(bAdj1{n_ra},I,J);
bAdj(:,:,n_ra,2) = reshape(bAdj2{n_ra},I,J);
aAdj(:,:,n_ra,1) = reshape(aAdj1{n_ra},I,J);
aAdj(:,:,n_ra,2) = reshape(aAdj2{n_ra},I,J);
end

%% KF EQUATION

if max(abs(sum(A, 2))) > 1e-10
    warning('Some rows of A do not sum to 0')
end

% get list of points in the adjustment region
pointslist = (1:L)';
adjpoints = pointslist(adj);
noadjpoints = pointslist(~adj);

% BUILD M MATRIX

% transforming bAdj, aAdj into a matrix of dimension I*J x N_ra
bAdj1_ = bAdj1{1};
bAdj2_ = bAdj2{1};
aAdj1_ = aAdj1{1};
aAdj2_ = aAdj2{1};

for n_ra = 2:N_ra
    bAdj1_ = [bAdj1_, bAdj1{n_ra}];
    bAdj2_ = [bAdj2_, bAdj2{n_ra}];
    aAdj1_ = [aAdj1_, aAdj1{n_ra}];
    aAdj2_ = [aAdj2_, aAdj2{n_ra}];
end
  
bAdj1 = bAdj1_;
bAdj2 = bAdj2_;
aAdj1 = aAdj1_;
aAdj2 = aAdj2_;
clear bAdj1_ bAdj2_ aAdj1_ aAdj2_
 % first, put 1's on diagonal for points not in adjustment region
M = sparse(noadjpoints, noadjpoints, ones(length(noadjpoints),1), L, L);
% next, take care of points in adjustment region
pointsmat = reshape(pointslist, I, J, N_ra, Nz);
count_idwpoints = 0; % count number of times inverse distance weights are used to allocate mass to be rebalanced
for mcount = 1:length(adjpoints)
    m = adjpoints(mcount);
    if m <= I*J*N_ra % get (b',a') from point (b,a)
        zind = 1;
        ra_ind = ceil(m/(I*J));
        bprime = bAdj1(m-I*J*(ra_ind-1),ra_ind);
        aprime = aAdj1(m-I*J*(ra_ind-1),ra_ind);
    else
        zind = 2;
        ra_ind = ceil((m - I*J*N_ra)/(I*J));
        bprime = bAdj2(m - I*J*N_ra -I*J*(ra_ind-1),ra_ind);
        aprime = aAdj2(m - I*J*N_ra -I*J*(ra_ind-1),ra_ind);      
    end
    bprime_left = discretize(bprime, b); % this gives the index of the left edge for liquid asset
    aprime_left = discretize(aprime, a); % same for illiquid asset
    bprime_right = bprime_left + 1;
    aprime_right = aprime_left + 1;
    
    % map from grid indexes to points
    point11 = pointsmat(bprime_left, aprime_left, ra_ind, zind);
    point12 = pointsmat(bprime_left, aprime_right, ra_ind, zind);
    point21 = pointsmat(bprime_right, aprime_left, ra_ind, zind);
    point22 = pointsmat(bprime_right, aprime_right, ra_ind, zind);
    neighpoints = [point11 point12;
                    point21 point22];
    % check if each point is not an adjustment point
    pointscheck = reshape(ismember(neighpoints(:), adjpoints), 2, 2); % 0 = not an adjustment point, so keep it
    pointscheck = 1 - pointscheck;
    if sum(sum(pointscheck)) == 4 && BilinW == 1 % use weights proportional to area opposite the point
        totarea = (b(bprime_right) - b(bprime_left)) * (a(aprime_right) - a(aprime_left)); % this is just da * db with equispaced grid 
        weights = totarea^(-1) * [b(bprime_right) - bprime; bprime - b(bprime_left)] * [a(aprime_right) - aprime aprime - a(aprime_left)];
    else
        count_idwpoints = count_idwpoints + 1; % use inverse distance weights
        neigh_bvals = [ones(1,2) * b(bprime_left);
                     ones(1,2) * b(bprime_right)];
        neigh_avals = [ones(2,1) * a(aprime_left) ones(2,1) * a(aprime_right)];
        bprime_mat = ones(2,2) * bprime;
        aprime_mat = ones(2,2) * aprime;
        dist_points = sqrt((neigh_bvals - bprime_mat).^2 + (neigh_avals - aprime_mat).^2); % this is the Euclidean distance
        inverse_dist = (1 ./ dist_points) .* pointscheck; % kill the points which are adjustment points
        if sum(sum(isinf(inverse_dist))) == 1
            weights = isinf(inverse_dist) .* ones(2,2);
        else
            weights = inverse_dist / sum(sum(inverse_dist));
        end
    end
    M = M + sparse(m * ones(4, 1), neighpoints(:), weights(:), L, L);
end

display(['Number of points where IDWs were used = ', num2str(count_idwpoints)])

breaki=1


AMT = transpose(A*M);

% Solve for stationary distribution (time-iteration)
Deltat = 25;
g_stacked = ones(L,1)./(L*da*db*dr_a); % initial distribution
for n = 1:maxit
    g_nplus12 = transpose(M)*g_stacked; % apply M to ensure g has no mass in adjustment region
    g_nplus1 = (speye(L) - Deltat*AMT)\g_nplus12;   % implicit method
    if (sum(abs(g_nplus1-g_stacked))) < tol
        disp(['Distribution Converged, Iteration = ' int2str(n)]);
        break
    end
    g_stacked = g_nplus1;
end

if max(abs(AMT * g_stacked)) > 1e-6
    warning('Solution to KFE not accurate')
end
if sum(g_stacked < -1e-6) ~= 0
    warning('Solution to KFE has negative values')
end

g(:,:,:,1) = reshape(g_stacked(1:I*J*N_ra),I,J,N_ra);
g(:,:,:,2) = reshape(g_stacked(I*J*N_ra+1:I*J*N_ra*Nz),I,J,N_ra);

fprintf('Made it !!')
save('DATA_23122020') % to avoid going through the iteration again.

run figures.m
