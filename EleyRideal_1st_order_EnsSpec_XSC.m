clear all; close all; clc
tic
%% Variables and Parameters
KA = logspace(-6,3,50); % k_ads array
kad = 1; % k_des
kr = 10^2/2; % k_r
kre = 1; % k_ER
tspan = [0 10^10]; % time for ode solver
delta = 10^-4; % differential change for calculation of reaction sensitivities
deldel = 10^-5;
eps = 10^-10; % infinitesmal amount of a on surface at t = 0
%% Arrays for results
S = zeros(length(KA),4)'; r = zeros(1,length(KA)); A = r; O = r; AO = A; AA = A; OO = A; re = r; Sel = r; SEL = r;
S_ER = zeros(length(KA),4)';
S_Sel = S_ER;
%% Solvers and Nested Loops
for w = 1:length(KA)
    disp(w)
    ka = KA(w);
    for m = 1:9
        %% Sites
        ao = @(x) x(1); aa = @(x) x(2);
        a = @(x) aa(x) + ao(x);
        o = @(x) 1-a(x);
        oo = @(x) o(x) - ao(x);
        
        % A adsorption
        if m == 2
            ads_aa = @(x) 2*ka*(1+deldel).*ao(x);
        else
            ads_aa = @(x) 2*ka.*ao(x);
        end
        
        if m == 4
            ads_ao = @(x) ka*(1+deldel).*(oo(x)-ao(x));
        else
            ads_ao = @(x) ka.*(oo(x)-ao(x));
            
        end
        
        % A desorption
        if m == 3
            des_aa = @(x) -2*kad*(1+deldel).*aa(x);
        else
            des_aa = @(x) -2*kad.*aa(x);
        end
        
        if m == 5
            des_ao = @(x) kad*(1+deldel).*(aa(x)-ao(x));
        else
            des_ao = @(x) kad*(aa(x)-ao(x));
        end
        
        % rxn
        if m == 6
            r_aa = @(x) -kr.*(1+deldel)*aa(x).*(1+6.*aa(x)./a(x));
        else
            r_aa = @(x) -kr.*aa(x).*(1+6.*aa(x)./a(x));
        end
        
        if m == 7
            r_ao = @(x) 3*kr*(1+deldel).*aa(x).*(aa(x)./a(x) - ao(x)./a(x));
        else
            r_ao = @(x) 3*kr.*aa(x).*(aa(x)./a(x) - ao(x)./a(x));
        end
        
        % A desorption
        if m == 8
            re_aa = @(x) -2*kre*(1+deldel).*aa(x);
        else
            re_aa = @(x) -2*kre.*aa(x);
        end
        
        if m == 9
            re_ao = @(x) kre*(1+deldel).*(aa(x)-ao(x));
        else
            re_ao = @(x) kre*(aa(x)-ao(x));
        end
        
        F = @(t,x) [ads_ao(x) + des_ao(x) + r_ao(x) + re_ao(x); ads_aa(x) + des_aa(x) + r_aa(x) + re_aa(x)];
        
        eps = 10^-10;
        C0 = [eps*(1-eps) eps^2];
        
        
        [T,C] = ode23s(F,tspan,C0);
        t = T;
        ao = C(:,1);
        aa = C(:,2);
        a = ao + aa;
        o = 1 - a;
        oo = o - ao;
        
        em = round(length(aa)*0.20);
        
        if m < 6 || m > 7
            r(m,w) = 4*kr.*mean(aa(end-em:end));
        else
            r(m,w) = (ka.*mean(o(end-em:end)) - kad.*mean(a(end-em:end))) - kre.*mean(a(end-em:end));
        end
        
        if m < 8
            re(m,w) = kre.*mean(a(end-em:end));
        else
            re(m,w) = (ka.*mean(o(end-em:end)) - kad.*mean(a(end-em:end))) - 4*kr.*mean(aa(end-em:end));
        end
        
        SEL(m,w) = 100*r(m,w)./(r(m,w)+re(m,w));
        
        A(m,w) = mean(a(end-em:end));
        AO(m,w) = mean(ao(end-em:end));
        AA(m,w) = mean(aa(end-em:end));
        OO(m,w) = mean(oo(end-em:end));
        O(m,w) = mean(o(end-em:end));
    end
end
%% Extract results
% DoRC for A2 formation
sAA_ads = (r(2,:) - r(1,:))./(deldel*r(1,:));
sAA_des = (r(3,:) - r(1,:))./(deldel*r(1,:));

sAO_ads = (r(4,:) - r(1,:))./(deldel*r(1,:));
sAO_des = (r(5,:) - r(1,:))./(deldel*r(1,:));

sAA_rxn = (r(6,:) - r(1,:))./(deldel*r(1,:));
sAO_rxn = (r(7,:) - r(1,:))./(deldel*r(1,:));

sAA_ER = (r(8,:) - r(1,:))./(deldel*r(1,:));
sAO_ER = (r(9,:) - r(1,:))./(deldel*r(1,:));

DORC_AA = sAA_ads + sAA_des + sAA_rxn + sAA_ER;
DORC_AO = sAO_ads + sAO_des + sAO_rxn + sAO_ER;

% DoRC for AB formation
seAA_ads = (re(2,:) - re(1,:))./(deldel*re(1,:));
seAA_des = (re(3,:) - re(1,:))./(deldel*re(1,:));

seAO_ads = (re(4,:) - re(1,:))./(deldel*re(1,:));
seAO_des = (re(5,:) - re(1,:))./(deldel*re(1,:));

seAA_rxn = (re(6,:) - re(1,:))./(deldel*re(1,:));
seAO_rxn = (re(7,:) - re(1,:))./(deldel*re(1,:));

seAA_ER = (re(8,:) - re(1,:))./(deldel*re(1,:));
seAO_ER = (re(9,:) - re(1,:))./(deldel*re(1,:));

% DoSC 
sSelAA_ads = (SEL(2,:) - SEL(1,:))./(deldel*SEL(1,:));
sSelAA_des = (SEL(3,:) - SEL(1,:))./(deldel*SEL(1,:));

sSelAO_ads = (SEL(4,:) - SEL(1,:))./(deldel*SEL(1,:));
sSelAO_des = (SEL(5,:) - SEL(1,:))./(deldel*SEL(1,:));

sSelAA_rxn = (SEL(6,:) - SEL(1,:))./(deldel*SEL(1,:));
sSelAO_rxn = (SEL(7,:) - SEL(1,:))./(deldel*SEL(1,:));

sSelAA_ER = (SEL(8,:) - SEL(1,:))./(deldel*SEL(1,:));
sSelAO_ER = (SEL(9,:) - SEL(1,:))./(deldel*SEL(1,:));

% reversibilities
zads = kad.*A(1,:)./(KA.*O(1,:));
zads_AA = kad.*AA(1,:)./(KA.*AO(1,:));
zads_OO = kad.*AO(1,:)./(KA.*OO(1,:));
zads_AO = -kad.*(AA(1,:)-AO(1,:))./(KA.*(OO(1,:)-AO(1,:)));

% mean-field metrics
mu_AA = AA(1,:)./A(1,:)./A(1,:);
mu_AO = AO(1,:)./A(1,:)./O(1,:);
mu_OO = OO(1,:)./O(1,:)./O(1,:);

%% DoRC adsorption for A2 formation
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),sAO_ads+sAO_des+sAA_ads+sAA_des,'k-','linewidth',1);
plot(log10(KA),sAA_ads+sAA_des,'r--','linewidth',1);
plot(log10(KA),sAO_ads+sAO_des,'b-.','linewidth',1);
xlabel('log_{10}k_{A,ads}');
legend('{\it X_{RC,ads}}','{\it X_{RC,ads,[aa]}}', '{\it X_{RC,ads,[ao]}}');
legend boxoff

%% DoRC surface reaction for A2 formation
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),sAO_rxn+sAA_rxn,'k-','linewidth',1);
plot(log10(KA),sAA_rxn,'r--','linewidth',1);
plot(log10(KA),sAO_rxn,'b-.','linewidth',1);
xlabel('log_{10}k_{A,ads}');
legend('{\it X_{RC,r}}','{\it X_{RC,r,[aa]}}', '{\it X_{RC,r,[ao]}}');
legend boxoff

%% DoRC Eley-Rideal for A2 formation
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),sAO_ER+sAA_ER,'k-','linewidth',1);
plot(log10(KA),sAA_ER,'r--','linewidth',1);
plot(log10(KA),sAO_ER,'b-.','linewidth',1);
xlabel('log_{10}k_{A,ads}');
legend('{\it X_{RC,ER}}','{\it X_{RC,ER,[aa]}}', '{\it X_{RC,ER,[ao]}}');
legend boxoff

%% DoRC adsorption for AB formation
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),seAO_ads+seAO_des+seAA_ads+seAA_des,'k-','linewidth',1);
plot(log10(KA),seAA_ads+seAA_des,'r--','linewidth',1);
plot(log10(KA),seAO_ads+seAO_des,'b-.','linewidth',1);
xlabel('log_{10}k_{A,ads}');
legend('{\it X_{RC,ER,ads}}','{\it X_{RC,ER,ads,[aa]}}', '{\it X_{RC,ER,ads,[ao]}}');
legend boxoff

%% DoRC surface reaction for AB formation
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),seAO_rxn+seAA_rxn,'k-','linewidth',1);
plot(log10(KA),seAA_rxn,'r--','linewidth',1);
plot(log10(KA),seAO_rxn,'b-.','linewidth',1);
xlabel('log_{10}k_{A,ads}');
legend('{\it X_{RC,ER,r}}','{\it X_{RC,ER,r,[aa]}}', '{\it X_{RC,ER,r,[ao]}}');
legend boxoff

%% DoRC Eley-Rideal for AB formation
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),seAO_ER+seAA_ER,'k-','linewidth',1);
plot(log10(KA),seAA_ER,'r--','linewidth',1);
plot(log10(KA),seAO_ER,'b-.','linewidth',1);
xlabel('log_{10}k_{A,ads}');
legend('{\it X_{RC,ER,ER}}','{\it X_{RC,ER,ER,[aa]}}', '{\it X_{RC,ER,ER,[ao]}}');
legend boxoff

%% DoSC adsorption
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),sSelAO_ads+sSelAO_des+sSelAA_ads+sSelAA_des,'k-','linewidth',1);
plot(log10(KA),sSelAA_ads+sSelAA_des,'r--','linewidth',1);
plot(log10(KA),sSelAO_ads+sSelAO_des,'b-.','linewidth',1);
ylabel('{\it X_{SC,ads,[j]}}');
xlabel('log_{10}k_{A,ads}');
legend('{\it X_{SC,ads}}','{\it X_{SC,ads,[aa]}}', '{\it X_{SC,ads,[ao]}}');
legend boxoff
ylim([-0.5 1.5])
xlim([-6 3])

%% DoSC surface reaction
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),sSelAO_rxn+sSelAA_rxn,'k-','linewidth',1);
plot(log10(KA),sSelAA_rxn,'r--','linewidth',1);
plot(log10(KA),sSelAO_rxn,'b-.','linewidth',1);
xlabel('log_{10}k_{A,ads}');
legend('{\it X_{SC,r}}','{\it X_{SC,r,[aa]}}', '{\it X_{SC,r,[ao]}}');
legend boxoff

%% DoSC Eley-Rideal
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),sSelAO_ER+sSelAA_ER,'k-','linewidth',1);
plot(log10(KA),sSelAA_ER,'r--','linewidth',1);
plot(log10(KA),sSelAO_ER,'b-.','linewidth',1);
xlabel('log_{10}k_{A,ads}');
ylabel('{\it X_{SC,ER,[j]}}');
ylim([-1.5 0.5])
xlim([-6 3])
legend('{\it X_{SC,ER}}','{\it X_{SC,ER,[aa]}}', '{\it X_{SC,ER,[ao]}}');
legend boxoff