clear all; close all; clc
tic
%% Rate Constants
% Rate Constants at 300K
kado = 10^3;
kbo = 1;
kbdo = 10^4/2;
kro = 10^5/4;
PB = 10^-4;

% Activation parameters
% unit convert eV/atom to kJ/mol = 0.0103642688 %
kB = 8.617333262*10^-5; % eV/K (Boltzmann)
h = 4.135667696*10^-15; % eV*s (Planck)
Ea = [0.4; 0.7; 0.4; 0.6; 0.5]; % eV (Choice of activation energy)
H = (Ea - kB*300)'; % eV (Convert to activation enthalpy)
Ent = [-124.6358   29.2753 -124.6358   10.4962   -8.2830]; % J/mol (Choice of activation entropy)
T = 300; % K
Temp = T; % Because later T is a time array

new = kB*T/h.*exp(-H./kB./T).*exp(Ent./8.314); % new rate constants

% "Adjusted" rate constants
adj = new(1); kado = new(2); kbo = PB*new(3); kbdo = new(4); kro = new(5);

delta = 10^-4; % for calculation of reaction sensitivities
tspan = [0 10^5]; 
eps = 10^-5; % initial A* and B* coverages

%% Set up arrays
KA = adj*logspace(-6,-4,10); % adj = rate constant and logspace array = pressure of A
DORC = zeros(length(KA),5)'; r = zeros(1,length(KA)); A = r; B = r; sigma_app = r; TRC_A = r; TRC_B = r; AB = B; AO = A; BB = B; BO = B; AA = B; OO = B; O = B;
%% Loop
for w = 1:length(KA)
    disp(w)
    for m = 1:1
        k0 = [KA(w) kado kbo kbdo kro];
        for j = 1:length(k0)
            K = zeros(2,length(k0));
            K(1,:) = k0; K(2,:) = k0;
            K(2,j) = k0(j)*(1+delta);
            for i = 1:2
                %% Parameters
                ka = K(i,1);
                kad = K(i,2);
                kb = K(i,3);
                kbd = K(i,4);
                kr = K(i,5);
                
                %% Define Sites
                a = @(x) x(1); b = @(x) x(2);
                o = @(x) 1 - a(x) - b(x);
                
                %% A adsorption
                A_ads_a = @(x) ka.*o(x);
                A_ads_b = @(x) 0;

                %% A desorption
                A_des_a = @(x) -kad.*a(x);
                A_des_b = @(x) 0;
                
                %% B adsorption
                B_ads_a = @(x) 0;
                B_ads_b = @(x) 4*kb.*o(x).^2;
                
                %% B desorption
                B_des_a = @(x) 0;
                B_des_b = @(x) -4*kbd.*b(x).^2;

                %% Reaction
                R_a = @(x) -4*kr.*a(x).*b(x);
                R_b = @(x) -4*kr.*a(x).*b(x);
                
                %% Collate ODEs
                F = @(t,x) [A_ads_a(x) + A_des_a(x) + B_ads_a(x) + B_des_a(x) + R_a(x);
                    A_ads_b(x) + A_des_b(x) + B_ads_b(x) + B_des_b(x) + R_b(x)];
                
                %% Solve
                C0 = [eps eps];
                [T,C] = ode23s(F,tspan,C0);
                t = T;
                a = C(:,1);
                b = C(:,2);
                o = 1 - a - b;
                ab = a.*b; aa = a.*a; ao = a.*o; bb = b.*b; bo = b.*o; oo = o.*o;
                
                R(i,j) = 4*kr.*mean(ab(end-round(length(ab)*0.20):end));
            end
            DORC(j,w) = (R(2,j) - R(1,j))./(delta*R(1,j));
            r(m,w) = 4*kro.*mean(ab(end-round(length(ab)*0.20):end));
            A(m,w) = mean(a(end-round(length(a)*0.20):end));
            B(m,w) = mean(b(end-round(length(a)*0.20):end));
            AB(m,w) = mean(ab(end-round(length(a)*0.20):end));
            BB(m,w) = mean(bb(end-round(length(a)*0.20):end));
            BO(m,w) = mean(bo(end-round(length(a)*0.20):end));
            AO(m,w) = mean(ao(end-round(length(a)*0.20):end));
            AA(m,w) = mean(aa(end-round(length(a)*0.20):end));
            OO(m,w) = mean(oo(end-round(length(a)*0.20):end));
            O(m,w) = mean(o(end-round(length(a)*0.20):end));
        end
    end
end
%% Extract results
zA = kado.*A(1,:)./(KA.*O(1,:));
zB = kbdo.*BB(1,:)./(kbo.*OO(1,:));
zA_AA = kado*AA(1,:)./(KA.*AO(1,:));
zA_AB = kado*AB(1,:)./(KA.*BO(1,:));
zA_AO = kado*(AA(1,:)-AO(1,:))./(KA.*(AO(1,:)-OO(1,:)));
zA_OO = kado*AO(1,:)./(KA.*OO(1,:));

zB_AB = (kbdo.*BB(1,:).*AB(1,:)./B(1,:))./(kbo.*OO(1,:).*AO(1,:)./O(1,:));

mu_AB = AB(1,:)./A(1,:)./B(1,:);
mu_AA = AA(1,:)./A(1,:)./A(1,:);
mu_BB = BB(1,:)./B(1,:)./B(1,:);
mu_OO = OO(1,:)./O(1,:)./O(1,:);
mu_AO = AO(1,:)./A(1,:)./O(1,:);
mu_BO = BO(1,:)./B(1,:)./O(1,:);

%% plotting
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(log10(KA),sum(DORC,1),'k-','linewidth',1)
ylim([0 2]);
title('DORC check');

%% Plotting
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),DORC(5,:),'k-','linewidth',1);
plot(log10(KA),zA,'b--','linewidth',1);
plot(log10(KA),zB,'r--','linewidth',1);
xlabel('{\it log k_{1}P_{A}} [a.u.]');
ylabel('Sensitivity and Reversibility');
legend('s_{rxn}','z_{A}','z_{B_{2}}');
legend boxoff
ylim([0 1]);
