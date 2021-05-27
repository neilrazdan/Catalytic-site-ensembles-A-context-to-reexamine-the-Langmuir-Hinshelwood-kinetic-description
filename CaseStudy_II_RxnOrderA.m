clear all; close all; clc
tic
%% Rate Constants
% Rate Constants at 300K
kado = 10^3;
kbo = 1;
kbdo = 10^4/2;
kro = 10^5/4;
PB = 10^-1;

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
                aa = @(x) x(1); ao = @(x) x(2); ab = @(x) x(3); bo = @(x) x(4); bb = @(x) x(5); oo = @(x) x(6);
                
                a = @(x) aa(x) + ab(x) + ao(x);
                b = @(x) ab(x) + bb(x) + bo(x);
                o = @(x) ao(x) + bo(x) + oo(x);
                
                %% A adsorption
                A_ads_aa = @(x) 2*ka.*ao(x);
                A_ads_ab = @(x) ka.*bo(x);
                A_ads_ao = @(x) ka.*(oo(x) - ao(x));
                
                A_ads_bo = @(x) -ka.*bo(x);
                A_ads_bb = @(x) 0;
                
                A_ads_oo = @(x) -2*ka.*oo(x);
                
                %% A desorption
                A_des_aa = @(x) -2*kad.*aa(x);
                A_des_ab = @(x) -kad.*ab(x);
                A_des_ao = @(x) kad.*(aa(x) - ao(x));
                
                A_des_bo = @(x) kad.*ab(x);
                A_des_bb = @(x) 0;
                
                A_des_oo = @(x) 2*kad.*ao(x);
                
                %% B adsorption
                B_ads_aa = @(x) 0;
                B_ads_ab = @(x) 3.*kb.*oo(x)./o(x).*ao(x);
                B_ads_ao = @(x) -3.*kb.*oo(x)./o(x).*ao(x);
                
                B_ads_bo = @(x) 3*kb.*oo(x)./o(x).*(oo(x)-bo(x));
                B_ads_bb = @(x) kb.*oo(x)./o(x).*(o(x) + 6*bo(x));
                
                B_ads_oo = @(x) -kb.*oo(x)./o(x).*(o(x) + 6*oo(x));
                
                %% B desorption
                B_des_aa = @(x) 0;
                B_des_ab = @(x) -3.*kbd.*bb(x)./b(x).*ab(x);
                B_des_ao = @(x) 3.*kbd.*bb(x)./b(x).*ab(x);
                
                B_des_bo = @(x) 3*kbd.*bb(x)./b(x).*(bb(x) - bo(x));
                B_des_bb = @(x) -kbd.*bb(x)./b(x).*(b(x) + 6*bb(x));
                
                B_des_oo = @(x) kbd.*bb(x)./b(x).*(b(x) + 6*bo(x));
                
                %% Reaction
                R_aa = @(x) -6*kr.*ab(x).*aa(x)./a(x);
                R_ab = @(x) -kr.*ab(x).*(1+3*ab(x)./b(x)+3*ab(x)./a(x));
                R_ao = @(x) 3*kr.*ab(x).*(ab(x)./b(x) + aa(x)./a(x) - ao(x)./a(x));
                
                R_bo = @(x) 3*kr.*ab(x).*(bb(x)./b(x) + ab(x)./a(x) - bo(x)./b(x));
                R_bb = @(x) -6*kr.*ab(x).*bb(x)./b(x);
                
                R_oo = @(x) kr.*ab(x)*(2 + 6.*ao(x)./a(x) + 6*bo(x)./b(x));
                
                %% Collate ODEs
                F = @(t,x) [A_ads_aa(x) + A_des_aa(x) + B_ads_aa(x) + B_des_aa(x) + R_aa(x);
                    A_ads_ao(x) + A_des_ao(x) + B_ads_ao(x) + B_des_ao(x) + R_ao(x);
                    A_ads_ab(x) + A_des_ab(x) + B_ads_ab(x) + B_des_ab(x) + R_ab(x);
                    A_ads_bo(x) + A_des_bo(x) + B_ads_bo(x) + B_des_bo(x) + R_bo(x);
                    A_ads_bb(x) + A_des_bb(x) + B_ads_bb(x) + B_des_bb(x) + R_bb(x);
                    A_ads_oo(x) + A_des_oo(x) + B_ads_oo(x) + B_des_oo(x) + R_oo(x)];
                
                %% Solve
                C0 = [eps.^2 eps.*(1-eps) eps.^2 eps.*(1-eps) eps.^2 (1-eps).*(1-eps)];
                [T,C] = ode23s(F,tspan,C0);
                t = T;
                aa = C(:,1);
                ao = C(:,2);
                ab = C(:,3);
                bo = C(:,4);
                bb = C(:,5);
                oo = C(:,6);
                
                a =  aa + ab + ao;
                b =  ab + bb + bo;
                o =  ao + bo + oo;
                
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