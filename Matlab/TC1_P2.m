%% Pregunta 2: Costos de Ajuste No Convexos
clear all
clc

%Parametros
beta = 0.9; 
R = 0.04; 
theta = 1/3; 
rho = 0.85; 
%rho = 0; %parte 4
sigma = 0.05; 
%sigma = 0.15; %parte 6
delta = 0; 
F = 0.03; 
%F = 0.06; %parte 5
P = 0.02;
%P = 0.04; %parte 5
Kmin=0.1; 
Kmax=90; 
nK=1000; 
nZ=7;  
tol = 10^(-6); 
V0=zeros(nK,nZ);
Vnew=zeros(nK,nZ);
Kpolicy=zeros(nK,nZ);
display=1;
maxIter=1000;
rng('default');

% Tauchen
[log_z, Prob] = tauchen(0,rho,sigma,nZ);
z=exp(log_z);

% Grilla K
grilla_k=linspace(Kmin, Kmax, nK);

Pi_a = zeros(nK, nK, nZ); % Caso (a)
Pi_b = zeros(nK, nK, nZ); % Caso (b)

for i=1:nK %capital hoy
    for j=1:nK %capital mañana
        for h=1:nZ %shock hoy
            if i==j
                Pi_a(i,j,h)= z(h)*grilla_k(i)^theta -R*grilla_k(i);
                Pi_b(i,j,h)= Pi_a(i,j,h);
            else
                Pi_a(i,j,h)= z(h)*grilla_k(i)^theta -R*grilla_k(i) -F;
                Pi_b(i,j,h)= (1-P)*z(h)*grilla_k(i)^theta -R*grilla_k(i);
            end
        end
    end
end


%% caso (a)
tic
for iter=1:maxIter
    for j=1:nZ %z hoy
        Vesperado=ones(nK,1)*sum(Prob(j,:).*V0,2)';
        [maxx,indx]=max(Pi_a(:,:,j)+beta*Vesperado,[],2);
        Vnew(:,j)=maxx;
        Kpolicy(:,j)=grilla_k(indx);
    end 
        dif=max(max(abs(Vnew-V0)));
    V0=Vnew;
    if display==1
        fprintf('\n Diferencia: %g, Iteración: %g \n',dif,iter)
    end
    
    if dif<tol
        fprintf('\n Value converge luego de: %g iteraciones \n',iter)
        break;
    end
end
fprintf('\n Value converge en: %g segundos \n',toc)

umbral = zeros(2,7);
for i=1:nZ
    [minn, indd] = min(Kpolicy(:,i));
    [maxi, indi] = max(Kpolicy(:,i));
    umbral(1,i) = grilla_k(indd -1);
    if indi ==1000
        umbral(2,i) = 90;
    else
    umbral(2,i) = grilla_k(indi +1);
    end
end
figure(1)
plot(z, Kpolicy(1,:), '-o') % k*(z)
hold on
plot(z, umbral(1,:),'--o')
plot(z, umbral(2,:), '--o')
title('Funciones de umbral y capital óptimo, caso (a)')
xlabel('z','FontSize',15)
ylabel('Nivel de capital','FontSize',15);
legend('Capital óptimo','Cota inferior','Cota superior')

% Policy function
figure(2)
plot(grilla_k,Kpolicy(:,1))
hold on
plot(grilla_k,Kpolicy(:,4))
plot(grilla_k,Kpolicy(:,7))
title('Policy Function, caso(a)')
xlabel('K hoy','FontSize',15)
ylabel('K mañana','FontSize',15);
legend('z alto','z promedio','z bajo')

%Value function
figure(3)
plot(grilla_k,Vnew(:,1));hold on;
plot(grilla_k,Vnew(:,4));hold on
plot(grilla_k,Vnew(:,7));hold on
title('Value Function, caso (a)')
xlabel('K hoy','FontSize',15)
ylabel('Value function','FontSize',15);
legend('z alto','z promedio','z bajo')

%% (b)
tic
for iter=1:maxIter
    for j=1:nZ %z hoy
        Vesperado=ones(nK,1)*sum(Prob(j,:).*V0,2)';
        [maxx,indx]=max(Pi_b(:,:,j)+beta*Vesperado,[],2);
        Vnew(:,j)=maxx;
        Kpolicy(:,j)=grilla_k(indx);
    end 
        dif=max(max(abs(Vnew-V0)));
    V0=Vnew;
    if display==1
        fprintf('\n Diferencia: %g, Iteración: %g \n',dif,iter)
    end
    
    if dif<tol
        fprintf('\n Value converge luego de: %g iteraciones \n',iter)
        break;
    end
end
fprintf('\n Value converge en: %g segundos \n',toc)

umbral = zeros(2,7);
for i=1:nZ
    [minn, indd] = min(Kpolicy(:,i));
    [maxi, indi] = max(Kpolicy(:,i));
    umbral(1,i) = grilla_k(indd -1);
    if indi == 1000
        umbral(2,i) = 90;
    else
    umbral(2,i) = grilla_k(indi +1);
    end
end

figure(4)
plot(z, Kpolicy(1,:), '-o') % k*(z)
hold on
plot(z, umbral(1,:),'--o')
plot(z, umbral(2,:), '--o')
title('Funciones de umbral y capital óptimo, caso (b)')
xlabel('z','FontSize',15)
ylabel('Nivel de capital','FontSize',15);
legend('Capital óptimo','Cota inferior','Cota superior')

%Policy function
figure(5)
plot(grilla_k,Kpolicy(:,1))
hold on
plot(grilla_k,Kpolicy(:,4))
plot(grilla_k,Kpolicy(:,7))
title('Policy Function, caso(b)')
xlabel('K hoy','FontSize',15)
ylabel('K mañana','FontSize',15);
legend('z alto','z promedio','z bajo')

%Value function
figure(6)
plot(grilla_k,Vnew(:,1));hold on;
plot(grilla_k,Vnew(:,4));hold on
plot(grilla_k,Vnew(:,7));hold on
title('Value Function, caso (b)')
xlabel('K hoy','FontSize',15)
ylabel('Value function','FontSize',15);
legend('z alto','z promedio','z bajo')

%% Simulando k 
T = 100;
z_index = [4; zeros(T-1,1)];
z_sim = [z(z_index(1)); zeros(T-1,1)];
for i = 2:T
    z_index(i) = randsample(7, 1, true, Prob(z_index(i-1),:));
    z_sim(i) = z(z_index(i));
end

ksim = zeros(T,1);
ksim(1) = randsample(grilla_k,1);
for t=2:T
    ksim(t) = Kpolicy(ksim(t-1) == grilla_k, z_index(t-1));
end

figure(7)
plot(1:T, ksim);


            
            
