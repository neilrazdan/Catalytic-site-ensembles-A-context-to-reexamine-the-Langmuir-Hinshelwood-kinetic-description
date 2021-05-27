% clear all; close all; clc
tic
%% Variables and Parameters
KA = logspace(-6,3,50); % k_ads array
kado = 1; % k_des
kro = 10^2/2; % k_r
kreo = 1; % k_ER
tspan = [0 10^10]; % time for ode solver
delta = 10^-4; % differential change for calculation of reaction sensitivities
eps = 10^-10; % infinitesmal amount of a on surface at t = 0
%% Arrays for results
S = zeros(length(KA),4)'; r = zeros(1,length(KA)); A = r; O = r; AO = A; AA = A; OO = A; re = r; Sel = r; SEL = r;
S_ER = zeros(length(KA),4)';
S_Sel = S_ER;
%% Solvers and Nested Loops
for w = 1:length(KA)
    disp(w)
    for m = 1:1
        k0 = [KA(w) kado kro kreo];
        for j = 1:length(k0)
            K = zeros(2,length(k0));
            K(1,:) = k0; K(2,:) = k0;
            K(2,j) = k0(j)*(1+delta);
            for i = 1:2
                %% Rate Constants
                ka = K(i,1);
                kad = K(i,2);
                kr = K(i,3);
                kre = K(i,4);           
                %% Sites
                oo = @(x) x(1); aa = @(x) x(2);
                ao = @(x) 1/2 - (aa(x) + oo(x))/2;
                a = @(x) aa(x) + ao(x);
                o = @(x) oo(x) + ao(x);
                
                % A adsorption
                ads_aa = @(x) 2*ka.*ao(x);
                ads_ao = @(x) ka.*(oo(x)-ao(x));
                ads_oo = @(x) -2*ka.*oo(x);
                
                % A desorption
                des_aa = @(x) -2*kad.*aa(x);
                des_ao = @(x) kad.*(aa(x) - ao(x));
                des_oo = @(x) 2*kad.*ao(x);
                
                % rxn
                r_aa = @(x) -kr.*aa(x).*(1+6.*aa(x)./a(x));
                r_ao = @(x) 3.*kr.*aa(x).*(aa(x)./a(x)-ao(x)./a(x));
                r_oo = @(x) kr.*aa(x).*(1+6.*ao(x)./a(x));
                
                % er rxn
                re_aa = @(x) -2*kre.*aa(x);
                re_ao = @(x) kre.*(aa(x) - ao(x));
                re_oo = @(x) 2*kre.*ao(x);
                
                F = @(t,x) [ads_oo(x) + des_oo(x) + r_oo(x) + re_oo(x); ads_aa(x) + des_aa(x) + r_aa(x) + re_aa(x)];
                
                C0 = [(1-eps)^2 eps^2];
              
                [T,C] = ode23s(F,tspan,C0);
                t = T;
                oo = C(:,1);
                aa = C(:,2);
                ao = 1/2 - (aa+oo)/2;
                a = ao + aa;
                o = ao + oo;
                
                em = round(length(aa)*0.20);
                R(i,j) = 4*kr.*mean(aa((end-em):end));
                RE(i,j) =  kre.*mean(a((end-em):end));
                SEL(i,j) = 100*R(i,j)./(R(i,j)+RE(i,j));
            end
            S(j,w) = (R(2,j) - R(1,j))./(delta*R(1,j));
            S_ER(j,w) = (RE(2,j) - RE(1,j))./(delta*RE(1,j));
            S_Sel(j,w) = (SEL(2,j) - SEL(1,j))./(delta*SEL(1,j));
            r(m,w) = 4*kro.*mean(aa((end-em):end));
            re(m,w) = kreo.*mean(a((end-em):end));
            Sel(m,w) = 100*r(m,w)./(r(m,w)+re(m,w));
            A(m,w) = mean(a((end-em):end)); % a
            AO(m,w) = mean(ao((end-em):end)); % ao
            AA(m,w) = mean(aa((end-em):end)); % aa
            OO(m,w) = mean(oo((end-em):end)); % oo
            O(m,w) = mean(o((end-em):end)); % o
        end
    end
end
%% Extract results
% reversibilities
zads = kado.*A(1,:)./(KA.*O(1,:));
zads_AA = kado.*AA(1,:)./(KA.*AO(1,:));
zads_OO = kado.*AO(1,:)./(KA.*OO(1,:));
zads_AO = -kado.*(AA(1,:)-AO(1,:))./(KA.*(OO(1,:)-AO(1,:)));

% mean-field metrics
mu_AA = AA(1,:)./A(1,:)./A(1,:);
mu_AO = AO(1,:)./A(1,:)./O(1,:);
mu_OO = OO(1,:)./O(1,:)./O(1,:);
%% DoRC sum to unity 
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(log10(KA),sum(S,1),'k-','linewidth',1)
plot(log10(KA),sum(S_ER,1),'ko','linewidth',1)
ylim([0 2]);
legend('\Sigma X_{RC} for A_{2} form.','\Sigma X_{RC} for AB form.');
legend boxoff
title('DORC check');
xlabel('log_{10}k_{ads}');
%% Mean-field metrics
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(log10(KA),mu_AA,'k--','linewidth',1.5)
plot(log10(KA),mu_AO,'b-.','linewidth',1.5)
plot(log10(KA),mu_OO,'r-','linewidth',1.5)
xlabel('log_{10}k_{ads}');
ylabel('\mu_{ij}');
legend('\mu_{aa}','\mu_{ao}','\mu_{oo}')
legend boxoff
legend('location','northwest');
%% Selectivity
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(log10(KA),Sel(1,:),'k-','linewidth',1)
title('selectivity');
xlabel('log_{10}k_{ads}');
%% DORC for A2 formation
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),S(1,:)+S(2,:),'k-','linewidth',1.5);
plot(log10(KA),S(3,:),'r--','linewidth',1.5);
plot(log10(KA),S(4,:),'b-.','linewidth',1.5);
xlabel('log_{10}k_{ads}');
legend('{\it X_{RC,ads}}', '{\it X_{RC,r}}', '{\it X_{RC,r,ER}}');
legend boxoff
title('DoRC for {\it A}_{2(g)} formation');
%% DORC for AB formation
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),S_ER(1,:)+S_ER(2,:),'k-','linewidth',1.5);
plot(log10(KA),S_ER(3,:),'r--','linewidth',1.5);
plot(log10(KA),S_ER(4,:),'b-.','linewidth',1.5);
xlabel('log_{10}k_{ads}');
legend('{\it X_{RC,ads}}', '{\it X_{RC,r}}', '{\it X_{RC,r,ER}}');
legend boxoff
title('DoRC for {\it AB}_{(g)} formation');
%% Degree of Selectivity Control
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),S_Sel(1,:)+S_Sel(2,:),'k-','linewidth',1.5);
plot(log10(KA),S_Sel(3,:),'r--','linewidth',1.5);
plot(log10(KA),S_Sel(4,:),'b-.','linewidth',1.5);
xlabel('log_{10}k_{ads}');
ylabel('Degree of Selectivity Control');
legend('{\it X_{SC,ads}}', '{\it X_{SC,r}}', '{\it X_{SC,ER}}');
legend boxoff
ylim([-1.5 1.5])
xlim([-6 3])
title('1^{st}-order')
toc