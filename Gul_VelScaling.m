%% START

clc;
clear;
close all;

set(0,'defaultTextInterpreter','latex');

%% Load Data

genData = load('./ExperimentalData/Gul 2022/generalData.txt');

% P24 to S row 1
% P60 to P24 row 7

% P24 to Smooth
P24toS_Uinf = genData(1,1);
P24toS_utau1 = genData(1,2);
P24toS_delta0 = genData(1,3);

% P60 to P24
P60toP24_Uinf = genData(7,1);
P60toP24_utau1 = genData(7,2);
P60toP24_delta0 = genData(7,3);

P24toS_data = load('./ExperimentalData/Gul 2022/P24toS/P24toS.txt');
P24toS_xhat = P24toS_data(:,1);
P24toS_dibl = P24toS_data(:,2).*P24toS_delta0;
P24toS_utau2 = P24toS_data(:,3).*P24toS_utau1;

myDir = dir('./ExperimentalData/Gul 2022/P24toS/xhat*');

N = length(P24toS_xhat);

for i = 1:N
   
    thisName = myDir(i).name;
    thisFolder = myDir.folder;
    
    thisFile = strcat(thisFolder,'/',thisName);
    thisData = load(thisFile);
    
    P24toS_ydelta{i} = thisData(:,1);
    P24toS_velDefect{i} = thisData(:,2);
    
end

clear this*

P60toP24_data = load('./ExperimentalData/Gul 2022/P60toP24/P60toP24.txt');
P60toP24_xhat = P60toP24_data(2:end,1);
P60toP24_dibl = P60toP24_data(2:end,2).*P60toP24_delta0;
P60toP24_utau2 = P60toP24_data(2:end,3).*P60toP24_utau1;

myDir = dir('./ExperimentalData/Gul 2022/P60toP24/xhat*');

Np60 = length(P60toP24_xhat);


for i = 1:Np60
   
    thisName = myDir(i).name;
    thisFolder = myDir.folder;
    
    thisFile = strcat(thisFolder,'/',thisName);
    thisData = load(thisFile);
    
    P60toP24_ydelta{i} = thisData(:,1);
    P60toP24_velDefect{i} = thisData(:,2);
    
end

clear this* myDir

%% Mean Velocity

for i = 1:N
   
    thisDefect = P24toS_velDefect{i};
    
    P24toS_U{i} = P24toS_Uinf - (thisDefect.*P24toS_utau1);
    P24toS_ydelta0{i} = P24toS_ydelta{i}.*P24toS_delta0;
    
    
end

for i = 1:Np60
   
    thisDefect = P60toP24_velDefect{i};
    
    P60toP24_U{i} = P60toP24_Uinf - (thisDefect.*P60toP24_utau1);
    P60toP24_ydelta0{i} = P60toP24_ydelta{i}.*P60toP24_delta0;
    
    
end

clear this* ii

%% Need to calculate an approximated local BL height

nu = 1.52e-5; %Nu for air at T = 20C (just a guess)
P24toS_dnu = nu./P24toS_utau2;
P60toP24_dnu = nu./P60toP24_utau2;

% using delta/x = 0.16/(Re_x)^1/7

for i = 2:N
    
    P24toS_x(i) = P24toS_xhat(i)*P24toS_delta0;
    P24toS_Rex(i) = P24toS_Uinf*P24toS_x(i)/nu;
    P24toS_deltaloc(i) = (0.16*P24toS_x(i))/(P24toS_Rex(i)^(1/7)) ...
        + P24toS_delta0;
    thisYdelta = P24toS_ydelta{i};
    P24toS_y{i} = thisYdelta.*P24toS_deltaloc(i);
    
    P24toS_yibl{i} = P24toS_y{i}./P24toS_dibl(i);
    
    ii = find(P24toS_yibl{i} > 1.0, 1);
    P24toS_UinftyIBL(i) = P24toS_U{i}(ii);
    
end

for i = 1:Np60
    
    P60toP24_x(i) = P60toP24_xhat(i)*P60toP24_delta0;
    P60toP24_Rex(i) = P60toP24_Uinf*P60toP24_x(i)/nu;
    P60toP24_deltaloc(i) = (0.16*P60toP24_x(i))/(P60toP24_Rex(i)^(1/7)) ...
        + P60toP24_delta0;
    thisYdelta = P60toP24_ydelta{i};
    P60toP24_y{i} = thisYdelta.*P60toP24_deltaloc(i);
    
    P60toP24_yibl{i} = P60toP24_y{i}./P60toP24_dibl(i);
    
    ii = find(P60toP24_yibl{i} > 1.0, 1);
    P60toP24_UinftyIBL(i) = P60toP24_U{i}(ii);
    
end


%%  Mean Velocity Profiles

close all;

figure();
subplot(1,2,1);
for i = 2:N
    semilogx(P24toS_y{i}./P24toS_dnu(i),...
        P24toS_U{i}./P24toS_utau2(i),'^','MarkerSize',8); hold on;
end
set(gca,'FontSize',16);
xlabel('$z^+$','FontSize',18);
ylabel('$U^+$','FontSize',18);
title('Gul and Ganapathisubramani $P24 \rightarrow S$','FontSize',20);

subplot(1,2,2)
for i = 2:N
    semilogx(P24toS_yibl{i},P24toS_U{i}./P24toS_UinftyIBL(i)...
        ,'^','MarkerSize',8); hold on;
end
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$U/U_\infty(z = \delta_i)$','FontSize',18);
title('Gul and Ganapathisubramani $P24 \rightarrow S$','FontSize',20);

%%%% P60 to P24

figure();
subplot(1,2,1);
for i = 1:Np60
    semilogx(P60toP24_y{i}./P60toP24_dnu(i),...
        P60toP24_U{i}./P60toP24_utau2(i),'^','MarkerSize',8); hold on;
end
set(gca,'FontSize',16);
xlabel('$z^+$','FontSize',18);
ylabel('$U^+$','FontSize',18);
title('Inner Scaling','FontSize',20);

subplot(1,2,2)
for i = 1:Np60
    semilogx(P60toP24_yibl{i},P60toP24_U{i}./P60toP24_UinftyIBL(i)...
        ,'^','MarkerSize',8); hold on;
end
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$U/U_i$','FontSize',18);
title('IBL Scaling','FontSize',20);

%% Three Panel

myColors = ["#2800c7",  "#9d0076", "#c70043", "#cc4400", "#958b00", "#31b800"];

close all;

figure();
subplot(1,3,1);
for i = 2:N
    semilogx(P24toS_y{i}./P24toS_dnu(i),...
        P24toS_U{i}./P24toS_utau2(i),'^','MarkerSize',8,...
        'Color',myColors(i-1)); hold on;
end
set(gca,'FontSize',16);
xlabel('$z^+$','FontSize',18);
ylabel('$U^+$','FontSize',18);
title('Viscous Scaling','FontSize',20)
legend('$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','Interpreter','Latex',...
    'Location','NorthWest','NumColumns',2);

subplot(1,3,2)
for i = 2:N
    plot(P24toS_U{i}./P24toS_Uinf,P24toS_ydelta{i}...
        ,'^','MarkerSize',8,...
        'Color',myColors(i-1)); hold on;
end
set(gca,'FontSize',16);
ylabel('$z/\delta$','FontSize',18);
xlabel('$U/U_\infty$','FontSize',18);
title('Outer Scaling','FontSize',20)

subplot(1,3,3)
for i = 2:N
    plot(P24toS_U{i}./P24toS_UinftyIBL(i),P24toS_yibl{i}...
        ,'^','MarkerSize',8,...
        'Color',myColors(i-1)); hold on;
end
% xlim([0 1]);
set(gca,'FontSize',16);
ylabel('$z/\delta_i$','FontSize',18);
xlabel('$U/U_i$','FontSize',18);
title('IBL Scaling','FontSize',20)

%%
%%%%%% P60 to P24

figure();
subplot(1,3,1);

%%%% SEMILOG
% for i = 1:Np60
%     semilogx(P60toP24_y{i}./P60toP24_dnu(i),...
%         P60toP24_U{i}./P60toP24_utau2(i),'o','MarkerSize',8,...
%         'Color',myColors(i)); hold on;
% end
% set(gca,'FontSize',16);
% xlabel('$z^+$','FontSize',18);
% ylabel('$U^+$','FontSize',18);

for i = 1:Np60
    plot(P60toP24_U{i}./P60toP24_utau2(i),P60toP24_y{i}./P60toP24_dnu(i),...
        'o','MarkerSize',8,...
        'Color',myColors(i)); hold on;
end
set(gca,'FontSize',16);
ylabel('$z^+$','FontSize',18);
xlabel('$U^+$','FontSize',18);
title('Viscous Scaling','FontSize',20)
legend('$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','Interpreter','Latex',...
    'Location','NorthWest','NumColumns',2);

subplot(1,3,2)
for i = 1:Np60
    plot(P60toP24_U{i}./P60toP24_Uinf,P60toP24_ydelta{i}...
        ,'o','MarkerSize',8,...
        'Color',myColors(i)); hold on;
end
set(gca,'FontSize',16);
ylabel('$z/\delta$','FontSize',18);
xlabel('$U/U_\infty$','FontSize',18);
title('Outer Scaling','FontSize',20)

subplot(1,3,3)
for i = 1:Np60
    plot(P60toP24_U{i}./P60toP24_UinftyIBL(i),P60toP24_yibl{i}...
        ,'o','MarkerSize',8,...
        'Color',myColors(i)); hold on;
end
% xlim([0 1]);
set(gca,'FontSize',16);
ylabel('$z/\delta_i$','FontSize',18);
xlabel('$U/U_i$','FontSize',18);
title('IBL Scaling','FontSize',20)

%% Velocity Defect Law

smoothPI = 0.55;
roughPI = 0.70;
kappa = 0.39;

P24toS_YD = P24toS_ydelta{end};
P60toP24_YD = P60toP24_ydelta{end};

P24toS_YDi = P24toS_yibl{end};
P60toP24_YDi = P60toP24_yibl{end};

ii1 = find(P24toS_YDi >= 1.2,1);
ii2 = find(P60toP24_YDi >= 1.2,1);

P24toS_YDi = P24toS_YDi(1:ii1);
P60toP24_YDi = P60toP24_YDi(1:ii2);

P24toS_wakeIBL = defectLaw(kappa,smoothPI,P24toS_YDi);
P24toS_wake = defectLaw(kappa,smoothPI,P24toS_YD);

P60toP24_wakeIBL = defectLaw(kappa,roughPI,P60toP24_YDi);
P60toP24_wake = defectLaw(kappa,roughPI,P60toP24_YD);


%% Velocity Defect

color12 = [ "#292f56", "#6c3c6c", "#a94b67", ...
   "#caa73d", "#9ddf62", ...
    "#70fa8e"];

close all;

% Upstream friction velocity
figure();
subplot(1,3,1);
for i = 2:N
    semilogx(P24toS_ydelta{i},...
        P24toS_velDefect{i}.*(P24toS_utau1/P24toS_utau2(i)),...
        '^','MarkerSize',8,'Color',myColors(i-1)); hold on
end
semilogx(P24toS_YD,P24toS_wake,'k--','LineWidth',3);
ylim([0 15]);
set(gca,'FontSize',16);
xlabel('$z/\delta$','FontSize',18);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18);
title('Classic Scaling with Local $u_\tau$','FontSize',20);
legend('$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','Interpreter','Latex',...
    'Location','NorthEast','NumColumns',2);


% Now using local friction velocity
subplot(1,3,2);
for i = 2:N
    semilogx(P24toS_ydelta{i},P24toS_velDefect{i},...
        '^','MarkerSize',8,'Color',myColors(i-1)); hold on
end
semilogx(P24toS_YD,P24toS_wake,'k--','LineWidth',3);
ylim([0 15]);
set(gca,'FontSize',16);
xlabel('$z/\delta$','FontSize',18);
ylabel('$(U_\infty - U)/u_{\tau,0}$','FontSize',18);
title('Classic Scaling with Upstream $u_\tau$','FontSize',20);


% Finally with IBL Scaling
subplot(1,3,3);
for i = 2:N
    
    thisVel = P24toS_U{i};
    thisDefect = (P24toS_UinftyIBL(i) - thisVel)./P24toS_utau2(i);
    
    semilogx(P24toS_yibl{i},thisDefect,...
        '^','MarkerSize',8,'Color',myColors(i-1)); hold on
end
% ylim([0 5]);
semilogx(P24toS_YDi,P24toS_wakeIBL,'k--','LineWidth',3);
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$','FontSize',20);



%%%%%%% P60 to P24


figure();
subplot(1,2,1);
for i = 1:Np60
    semilogx(P60toP24_ydelta{i},...
        P60toP24_velDefect{i}.*(P60toP24_utau1/P60toP24_utau2(i)),...
        'o','MarkerSize',8,'Color',color12(i)); hold on
end
semilogx(P60toP24_YD,P60toP24_wake,'k--','LineWidth',3);
ylim([0 15]);
set(gca,'FontSize',16);
xlabel('$z/\delta$','FontSize',18);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18);
title('Classic Scaling with Local $u_\tau$','FontSize',20);
legend('$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','Interpreter','Latex',...
    'Location','NorthEast','NumColumns',2);

% Now using local friction velocity
% subplot(1,3,2);
% for i = 1:Np60
%     semilogx(P60toP24_ydelta{i},P60toP24_velDefect{i},...
%         '^','MarkerSize',8,'Color',myColors(i)); hold on
% end
% semilogx(P60toP24_YD,P60toP24_wake,'k--','LineWidth',3);
% ylim([0 15]);
% set(gca,'FontSize',16);
% xlabel('$z/\delta$','FontSize',18);
% ylabel('$(U_\infty - U)/u_{\tau,0}$','FontSize',18);
% title('Classic Scaling with Upstream $u_\tau$','FontSize',20);


% Finally with IBL Scaling
subplot(1,2,2);
for i = 1:Np60
    
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_UinftyIBL(i) - thisVel)./P60toP24_utau2(i);
    
    semilogx(P60toP24_yibl{i},thisDefect,...
        'o','MarkerSize',8,'Color',color12(i)); hold on
end
% ylim([0 5]);
semilogx(P60toP24_YDi,P60toP24_wakeIBL,'k--','LineWidth',3);
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$','FontSize',20);



%% New Velocity Profile Figures

%%%%%% P60 to P24

figure();
subplot(2,1,1);

%%%% SEMILOG
for i = 1:Np60
    semilogx(P60toP24_y{i}./P60toP24_dnu(i),...
        P60toP24_U{i}./P60toP24_utau2(i),'o','MarkerSize',8,...
        'Color',myColors(i)); hold on;
end
set(gca,'FontSize',16);
xlabel('$z^+$','FontSize',18);
ylabel('$U^+$','FontSize',18);

% for i = 1:Np60
%     plot(P60toP24_U{i}./P60toP24_utau2(i),P60toP24_y{i}./P60toP24_dnu(i),...
%         'o','MarkerSize',8,...
%         'Color',myColors(i)); hold on;
% end
% set(gca,'FontSize',16);
% ylabel('$z^+$','FontSize',18);
% xlabel('$U^+$','FontSize',18);
legend('$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','Interpreter','Latex',...
    'Location','NorthWest','NumColumns',2);

subplot(2,1,2)
for i = 1:Np60
    plot(P60toP24_U{i}./P60toP24_Uinf,P60toP24_ydelta{i}...
        ,'o','MarkerSize',8,...
        'Color',myColors(i)); hold on;
end
set(gca,'FontSize',16);
ylabel('$z/\delta$','FontSize',18);
xlabel('$U/U_\infty$','FontSize',18);


figure();
for i = 1:Np60
    plot(P60toP24_U{i}./P60toP24_UinftyIBL(i),P60toP24_yibl{i}...
        ,'o','MarkerSize',8,...
        'Color',myColors(i)); hold on;
end
set(gca,'FontSize',16);
ylabel('$z/\delta_i$','FontSize',18);
xlabel('$U/U_i$','FontSize',18);

figure();
for i = 1:Np60
    plot(P60toP24_U{i}./P60toP24_UinftyIBL(i),P60toP24_yibl{i}...
        ,'o','MarkerSize',8,...
        'Color',myColors(i)); hold on;
end
set(gca,'FontSize',16);
ylabel('$z/\delta_i$','FontSize',18);
xlabel('$U/U_i$','FontSize',18);
xlim([0 1])
ylim([0 1])


%% ZS Scaling

% Usual Values
for i = 2:Np60
    this_i = find(P60toP24_ydelta{i} >= 1.0, 1);
    this_z = P60toP24_y{i}(1:this_i);
   
    this_integrand = 1 - (P60toP24_U{i}./P60toP24_Uinf);
    this_integrand = this_integrand(1:this_i);
    
    delta_star(i) = trapz(this_z,this_integrand);
    
    u0(i) = P60toP24_Uinf*(delta_star(i)/P60toP24_deltaloc(i));
    
end

% Now with IBL Parameters
for i = 2:Np60
    this_zibl = P60toP24_yibl{i};
    
    this_i = find(this_zibl >= 1.0, 1);
    this_z = P60toP24_y{i}(1:this_i);
    
    this_integrand = 1 - (P60toP24_U{i}./P60toP24_Uinf);
    this_integrand = this_integrand(1:this_i);
    
    delta_i_star(i) = trapz(this_z,this_integrand);
    
    u0_i(i) = P60toP24_Uinf*(delta_i_star(i)/P60toP24_dibl(i));
    
end    

%% Plot

close all;

figure();
subplot(2,2,1);
for i = 2:Np60
    semilogx(P60toP24_ydelta{i},...
        P60toP24_velDefect{i}.*(P60toP24_utau1/P60toP24_utau2(i)),...
        'o','MarkerSize',8,'Color',color12(i)); hold on
end
ylim([0 15]);
set(gca,'FontSize',16);
xlabel('$z/\delta$','FontSize',18);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18);
title('Classic Scaling with Local $u_\tau$','FontSize',20);
legend('$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','Interpreter','Latex',...
    'Location','NorthEast','NumColumns',2);


% Finally with IBL Scaling
subplot(2,2,2);
for i = 2:Np60
    
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_UinftyIBL(i) - thisVel)./P60toP24_utau2(i);
    
    semilogx(P60toP24_yibl{i},thisDefect,...
        'o','MarkerSize',8,'Color',color12(i)); hold on
end
xlim([0 1]);
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$','FontSize',20);

subplot(2,2,3);
for i = 2:Np60
    
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_Uinf - thisVel)./u0(i);
    
    semilogx(P60toP24_ydelta{i},thisDefect,...
        'o','MarkerSize',8,'Color',color12(i)); hold on
end
set(gca,'FontSize',16);
xlabel('$z/\delta$','FontSize',18);
ylabel('$(U_\infty - U)/U_\infty\delta^*/\delta$','FontSize',18);
title('ZS Scaling','FontSize',20);

subplot(2,2,4);
for i = 2:Np60
    
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_UinftyIBL(i) - thisVel)./u0_i(i);
    
    semilogx(P60toP24_yibl{i},thisDefect,...
        'o','MarkerSize',8,'Color',color12(i)); hold on
end
xlim([0 1]);
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/U_i\delta_i^*/\delta$','FontSize',18);
title('ZS Scaling with IBL Parameters','FontSize',20);


%%%%%% Linear Plotting

figure();
subplot(2,2,1);
for i = 2:Np60
    plot(P60toP24_ydelta{i},...
        P60toP24_velDefect{i}.*(P60toP24_utau1/P60toP24_utau2(i)),...
        'o','MarkerSize',8,'Color',color12(i)); hold on
end
ylim([0 15]);
set(gca,'FontSize',16);
xlabel('$z/\delta$','FontSize',18);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18);
title('Classic Scaling with Local $u_\tau$','FontSize',20);
legend('$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','Interpreter','Latex',...
    'Location','NorthEast','NumColumns',2);


% Finally with IBL Scaling
subplot(2,2,2);
for i = 2:Np60
    
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_UinftyIBL(i) - thisVel)./P60toP24_utau2(i);
    
    plot(P60toP24_yibl{i},thisDefect,...
        'o','MarkerSize',8,'Color',color12(i)); hold on
end
xlim([0 1]);
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$','FontSize',20);

subplot(2,2,3);
for i = 2:Np60
    
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_Uinf - thisVel)./u0(i);
    
    plot(P60toP24_ydelta{i},thisDefect,...
        'o','MarkerSize',8,'Color',color12(i)); hold on
end
set(gca,'FontSize',16);
xlabel('$z/\delta$','FontSize',18);
ylabel('$(U_\infty - U)/U_\infty\delta^*/\delta$','FontSize',18);
title('ZS Scaling','FontSize',20);

subplot(2,2,4);
for i = 2:Np60
    
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_UinftyIBL(i) - thisVel)./u0_i(i);
    
    plot(P60toP24_yibl{i},thisDefect,...
        'o','MarkerSize',8,'Color',color12(i)); hold on
end
xlim([0 1]);
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/U_i\delta_i^*/\delta$','FontSize',18);
title('ZS Scaling with IBL Parameters','FontSize',20);

%% Functions

function wake = defectLaw(kappa,myPI,yd)
    wake =  (1/kappa)*( -1*log(yd) + myPI*(2 - 2*(sin((pi()/2)*yd).^2)));
end


%% END