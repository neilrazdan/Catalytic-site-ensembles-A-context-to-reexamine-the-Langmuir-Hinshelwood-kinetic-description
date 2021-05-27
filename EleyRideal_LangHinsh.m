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
S = zeros(length(KA),4)'; r = zeros(1,length(KA)); A = r; O = r; AO = A; AA = A; OO = A; re = r; Sel = r;
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
                a = @(x) x(1); 
                o = @(x) 1 - a(x);
                
                % A adsorption
                ra_a = @(x) ka.*o(x) - kad*a(x);
                
                % lh rxn
                rr_a = @(x) -4*kr.*a(x).^2;
                
                % eley-rideal rxn
                re_a = @(x) -kre.*a(x);
                
                F = @(t,x) [ra_a(x) + rr_a(x) + re_a(x)];
                
                C0 = [eps];
              
                [T,C] = ode23s(F,tspan,C0);
                t = T;
                a = C(:,1);
                o = 1 - a;
                aa = a.^2;
                ao = a.*o;
                oo = o.^2;
                
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
title('Langmuir-Hinshelwood')
toc