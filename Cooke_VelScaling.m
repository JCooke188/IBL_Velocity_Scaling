%% Start

% Created by Justin Cooke to identify a structure of the mean velocity
% profile after a roughness transition

clc;
clear;
close all;

set(0,'defaultTextInterpreter','latex');

%% Load Dune Field Data

myDir = dir('./ExperimentalData/Cooke et al 2024/Sept13/x*');

Ndir = length(myDir);

ixyz = 0;
iu = 0;
iw = 0;
iuw = 0;


for i = 1:Ndir
    

    myName = myDir(i).name;
    myFolder = myDir(i).folder;
    
    myPath = strcat(myFolder,'/',myName);
    
    data = load(myPath);
    
    if contains(myPath,'README')
        ixyz = ixyz + 1;
        x(ixyz) = round(data(1,2));
        z_local{ixyz} = data(:,end);
    elseif contains(myPath,'comp(u,0)')
        iu = iu + 1;
        u{iu} = data(:,4:end);
    elseif contains(myPath,'comp(u,2)')
        iw = iw + 1;
        w{iw} = data(:,4:end);
    elseif contains(myPath,'comp(u_rey,1)')
        iuw = iuw + 1;
        uw{iuw} = data(:,4:end);
    end
        
end

clear my* Ndir i* data

myDir = dir('./ExperimentalData/Cooke et al 2024/SWSSData/SWSS*');

surf_N = length(myDir);

cx = 0;
cr = 0;

for i = 1:surf_N
   
    myName = myDir(i).name;
    myFolder = myDir(i).folder;
    
    loadMe = strcat(myFolder,'/',myName);
    
    if contains(loadMe,'README')
        swss_xyz = load(loadMe);
        swss_xyz = swss_xyz(:,3:end);
        cx = cx + 1;
        swss_y{cx} = swss_xyz;
    else
        swss_tau = load(loadMe);
        swss_tau = swss_tau(:,4:end);
        cr = cr + 1;
        swss_tauw{cr} = swss_tau;
    end
    
    
end

clear myName myFolder loadMe cr cx swss_xyz swss_tau

% Note that df_uxdata is set up as the cell # represents the x location,
% the size of each cell is Nt x Nz



%% Relevant data/parameters

delta = 300;
delta_ibl_withAF = [0.1,0.06356,0.1665,0.1665,0.1373,0.1372,0.2968,0.4362,...
    0.4362,0.5289].*delta;
nu = 0.32e-4; % Kin. viscosity of air
rho = 1.23; % Density of air 

z = linspace(1.5,400,100);
zdelta = z./delta;
xhat = x - 1850;

cookeCorr = 0.29.*((xhat(2:end)).^0.71);
cookeCorr(1) = 30;
cookeCorr(2) = 30;

Nx = length(xhat);
myColors = [ "black", "#daf8e3", "#97ebdb", "#00c2c7", "#0086ad", "#005582", ...
    "#ffc100", "#ff9a00", "#ff7400", "#bf0000" ];

%% Wall Stress

% Want to create a spanwise averaged u_tau for each xhat

% Spanwise station based averaging so now we have time and space average

for i = 1:Nx-1
   
    swss_taubar(i) = mean(swss_tauw{i},'All');
    utau_taubar(i) = sqrt(swss_taubar(i)/rho);
    dnu_taubar(i)  = nu/utau_taubar(i);
        
end

utau = [utau_taubar(1), utau_taubar];

%% Create Averages

N_u = length(u);
N_w = length(w);
N_uw = length(uw);

for i = 1:N_u

    tempu = u{i};
    U{i} = mean(tempu);
    up{i} = tempu - U{i};
    urms{i} = rms(up{i});
    
    tempw = w{i};
    W{i} = mean(tempw);
    wp{i} = tempw - W{i};
    wrms{i} = rms(up{i});
    
    tempuw = uw{i};
    UW{i} = mean(tempuw);
    
end
    
clear temp*

%% Find U_infty at z(\delta_i)

for i = 1:Nx-1
    
    ii = find(z >= delta_ibl_withAF(i),1);
    U_infty_i(i) = U{i}(ii);
    
end

%% Plot with U normalized

close all;

myPlot = zeros(10,1) + 1;

figure();
for i = 1:N_u-1
    plot(U{i}/U{i}(end),z/delta,'-','LineWidth',2,...
        'Color',myColors(i)); hold on
    plot(myPlot(i),delta_ibl_withAF(i)./delta,'ko',...
        'MarkerFaceColor',myColors(i),'MarkerSize',8); hold on
end
xlabel('$\langle U \rangle$/$U_\infty$','FontSize',18,'Interpreter','latex');
ylabel('$z/\delta$','FontSize',18,'Interpreter','latex');


figure();
for i = 1:N_u-1
    plot(U{i}/U_infty_i(i),z./delta_ibl_withAF(i),'-','LineWidth',2,...
        'Color',myColors(i)); hold on
end
yline(1,'LineWidth',2);
xlabel('$\langle U \rangle$/$U_\infty(z=\delta_i)$','FontSize',18,'Interpreter','latex');
ylabel('$z/\delta_i$','FontSize',18,'Interpreter','latex');
ylim([0 5]);

% figure();
% subplot(1,2,1);

%% Three Panel

myColors = ["#31b800", "#7b9e00", "#a57e00", "#c75000", "#d9000c", ...
    "#c2004a", "#a90069", "#88008e", "#2800c7"];

myColors = fliplr(myColors);

color12 = [ "#292f56", "#4a3664", "#6c3c6c", "#8d426d", ...
    "#bf5b5c", "#cc7050", "#d08a43", "#caa73d", "#9ddf62", ...
    "#70fa8e"];


close all;


figure();
subplot(1,3,1);
for i = 1:N_u-1
    if i == 1
        semilogx(z./(nu/utau(i)),U{i}./utau(i),'^','MarkerSize',8,...
        'Color','black'); hold on  
    else
        semilogx(z./(nu/utau(i)),U{i}./utau(i),'^','MarkerSize',8,...
        'Color',color12(i-1)); hold on        
    end
end
set(gca,'FontSize',16);
xlabel('$z^+$','FontSize',18);
ylabel('$\langle U\rangle^+$','FontSize',18);
title('Viscous Scaling','FontSize',20)
legend({'AF','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$',...
    '$\hat{x}_4$','$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$',...
    '$\hat{x}_8$','$\hat{x}_9$'},'Interpreter','Latex',...
    'Location','NorthWest','NumColumns',2);

subplot(1,3,2);
for i = 1:N_u-1
    if i == 1
        plot(U{i}./U{i}(end),zdelta,'^','MarkerSize',8,...
        'Color','black'); hold on
    else
        plot(U{i}./U{i}(end),zdelta,'^','MarkerSize',8,...
        'Color',color12(i-1)); hold on
    end
    
end
xlim([0 1]);
set(gca,'FontSize',16);
ylabel('$z/\delta_{99}$','FontSize',18);
xlabel('$\langle U\rangle/U_\infty$','FontSize',18);
title('Outer Scaling','FontSize',20)

subplot(1,3,3)
for i = 1:N_u-1
    if i == 1
        plot(U{i}./U_infty_i(i),z./delta_ibl_withAF(i),'^','MarkerSize',8,...
        'Color','black'); hold on
    else
        plot(U{i}./U_infty_i(i),z./delta_ibl_withAF(i),'^','MarkerSize',8,...
        'Color',color12(i-1)); hold on
    end
    
end
% xlim([0,1]);
set(gca,'FontSize',16);
ylabel('$z/\delta_i$','FontSize',18);
xlabel('$\langle U\rangle/U_i$','FontSize',18);
title('IBL Scaling','FontSize',20)

%% Smaller Velocity Profile

close all;


figure();
subplot(1,3,1);
% for i = 2:N_u-1   
%         semilogx(z./(nu/utau(i)),U{i}./utau(i),'^','MarkerSize',8,...
%         'Color',color12(i-1)); hold on        
% end
% set(gca,'FontSize',16);
% xlabel('$z^+$','FontSize',18);
% ylabel('$U^+$','FontSize',18);

for i = 2:N_u-1   
        plot(U{i}./utau(i),z./(nu/utau(i)),'^','MarkerSize',8,...
        'Color',color12(i-1)); hold on        
end
set(gca,'FontSize',16);
ylabel('$z^+$','FontSize',18);
xlabel('$U^+$','FontSize',18);
title('Viscous Scaling','FontSize',20)
legend({'$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$',...
    '$\hat{x}_4$','$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$',...
    '$\hat{x}_8$','$\hat{x}_9$'},'Interpreter','Latex',...
    'Location','NorthWest','NumColumns',3,'FontSize',10);

subplot(1,3,2);
for i = 2:N_u-1
    
        plot(U{i}./U{i}(end),zdelta,'^','MarkerSize',8,...
        'Color',color12(i-1)); hold on
    
end
xlim([0 1]);
set(gca,'FontSize',16);
ylabel('$z/\delta$','FontSize',18);
xlabel('$U/U_\infty$','FontSize',18);
title('Outer Scaling','FontSize',20)

subplot(1,3,3)
for i = 2:N_u-1
    
        plot(U{i}./U_infty_i(i),z./delta_ibl_withAF(i),'^','MarkerSize',8,...
        'Color',color12(i-1)); hold on
    
end
% xlim([0,1]);
set(gca,'FontSize',16);
ylabel('$z/\delta_i$','FontSize',18);
xlabel('$U/U_i$','FontSize',18);
title('IBL Scaling','FontSize',20)

%% Velocity Defect

close all;

delta_ibl = [0.1,0.06356,0.1665,0.1665,0.1373,0.1372,0.2968,0.4362,...
    0.4362,0.5289].*delta;



% myColors_AF = ["black", "#daf8e3", "#97ebdb", "#00c2c7", "#0086ad", "#005582", ...
%     "#ffc100", "#ff9a00", "#ff7400", "#bf0000" ];

figure();
subplot(1,2,1);
for i = 1:N_u-1
   
    semilogx(z./delta,(U{i}(end) - U{i})./utau(i),...
        'LineWidth',2,'Color',color12(i)); hold on
end
set(gca,'FontSize',16);
title('WMLES Data','FontSize',20);
xlabel('$z/\delta$','FontSize',18);
ylabel('$(U_\infty - \langle U\rangle) / u_{\tau,local}$','FontSize',18);

   
subplot(1,2,2);
for i = 1:N_u-1
   
    semilogx(z./delta_ibl(i),(U_infty_i(i) - U{i})./utau(i),...
        'LineWidth',2,'Color',color12(i)); hold on
end
set(gca,'FontSize',16);
title('WMLES Data','FontSize',20);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$U_\infty(z = \delta_i) - \langle U\rangle / u_{\tau,local}$',...
    'FontSize',18);

%% Three Panel

myColors = [ "black", "#daf8e3", "#97ebdb", "#00c2c7", "#0086ad", "#005582", ...
    "#ffc100", "#ff9a00", "#ff7400", "#bf0000" ];

close all;

figure();
subplot(1,3,1);
for i = 1:N_u-1
   
    semilogx(z./delta,(U{i}(end) - U{i})./utau(i),...
        '^','MarkerSize',8,'Color',myColors(i)); hold on
end
set(gca,'FontSize',16);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);
title('Classic Scaling with Local $u_\tau$');
legend({'AF','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$',...
    '$\hat{x}_4$','$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$',...
    '$\hat{x}_8$','$\hat{x}_9$'},'Interpreter','Latex',...
    'Location','NorthEast','NumColumns',2);

subplot(1,3,2);
for i = 1:N_u-1
   
    semilogx(z./delta,(U{i}(end) - U{i})./utau(1),...
        '^','MarkerSize',8,'Color',myColors(i)); hold on
end
set(gca,'FontSize',16);
% ylim([0 40]);
ylabel('$(U_\infty - U)/u_{\tau,0}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);
title('Classic Scaling with Upstream $u_\tau$');

subplot(1,3,3);
for i = 1:N_u-1
   
    semilogx(z./delta_ibl(i),(U_infty_i(i) - U{i})./utau(i),...
        '^','MarkerSize',8,'Color',myColors(i)); hold on
end
set(gca,'FontSize',16);
% ylim([0 40]);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$');

%% Plotting with Wake Function Regular

myPI = 0.55;
% myPI = 0.25;
kappa = 0.4;

thisYD = zdelta;

wake = defectLaw(kappa,myPI,thisYD);
    
close all;

figure();
for i = 1:N_u-1
   
    semilogx(z./delta,(U{i}(end) - U{i})./utau(i),...
        'k^','MarkerSize',8,'MarkerFaceColor',myColors(i)); hold on
end
semilogx(zdelta,wake,'k--','LineWidth',2);
set(gca,'FontSize',16);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);
title('Classic Scaling with Local $u_\tau$');
legend({'AF','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$',...
    '$\hat{x}_4$','$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$',...
    '$\hat{x}_8$','$\hat{x}_9$'},'Interpreter','Latex',...
    'Location','NorthEast','NumColumns',2);

%% Plotting with Wake Function IBL Scaling

myPI = 0.7;
% myPI = 0.25;
kappa = 0.4;

thisYD = z./delta_ibl(end);

ii = find(thisYD >= 1.2, 1);
thisYD = thisYD(1:ii);

wake = defectLaw(kappa,myPI,thisYD);
    
close all;

figure();
for i = 1:N_u-1
   
    semilogx(z./delta_ibl(i),(U_infty_i(i) - U{i})./utau(i),...
        'k^','MarkerSize',8,'MarkerFaceColor',myColors(i)); hold on
end
semilogx(thisYD,wake,'k--','LineWidth',3);
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$');
legend({'AF','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$',...
    '$\hat{x}_4$','$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$',...
    '$\hat{x}_8$','$\hat{x}_9$'},'Interpreter','Latex',...
    'Location','NorthEast','NumColumns',2);

%% Combined 

smoothPI = 0.55;
roughPI = 0.7;

thisYD = zdelta;
thisYDi = z./delta_ibl(end);

ii = find(thisYDi >= 1.2, 1);
thisYDi = thisYDi(1:ii);

wake_ibl = defectLaw(kappa,roughPI,thisYDi);
wake = defectLaw(kappa,roughPI,thisYD);
    
% close all;

figure();
subplot(1,2,1);
for i = 2:N_u-1
   
    semilogx(z./delta,(U{i}(end) - U{i})./utau(i),...
        '^','MarkerSize',8,'Color',color12(i-1)); hold on
end
semilogx(zdelta,wake,'k--','LineWidth',3);
set(gca,'FontSize',16);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);
title('Classic Scaling with Local $u_\tau$');
legend({'$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$',...
    '$\hat{x}_4$','$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$',...
    '$\hat{x}_8$','$\hat{x}_9$'},'Interpreter','Latex',...
    'Location','NorthEast','NumColumns',3,'FontSize',10);

subplot(1,2,2);
for i = 2:N_u-1
   
    semilogx(z./delta_ibl(i),(U_infty_i(i) - U{i})./utau(i),...
        '^','MarkerSize',8,'Color',color12(i-1)); hold on
end
semilogx(thisYD,wake,'k--','LineWidth',3);
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$');

%% New Velocity Profile Figures

figure();
subplot(2,1,1);
for i = 2:N_u-1   
        semilogx(z./(nu/utau(i)),U{i}./utau(i),'^','MarkerSize',8,...
        'Color',color12(i-1)); hold on        
end
set(gca,'FontSize',16);
xlabel('$z^+$','FontSize',18);
ylabel('$U^+$','FontSize',18);

% for i = 2:N_u-1   
%         plot(U{i}./utau(i),z./(nu/utau(i)),'^','MarkerSize',8,...
%         'Color',color12(i-1)); hold on        
% end
% set(gca,'FontSize',16);
% ylabel('$z^+$','FontSize',18);
% xlabel('$U^+$','FontSize',18);
legend({'$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$',...
    '$\hat{x}_4$','$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$',...
    '$\hat{x}_8$','$\hat{x}_9$'},'Interpreter','Latex',...
    'Location','NorthWest','NumColumns',3,'FontSize',10);

subplot(2,1,2);
for i = 2:N_u-1
    
        plot(U{i}./U{i}(end),zdelta,'^','MarkerSize',8,...
        'Color',color12(i-1)); hold on
    
end
xlim([0 1]);
set(gca,'FontSize',16);
ylabel('$z/\delta$','FontSize',18);
xlabel('$U/U_\infty$','FontSize',18);

figure();
for i = 2:N_u-1
    
        plot(U{i}./U_infty_i(i),z./delta_ibl_withAF(i),'^','MarkerSize',8,...
        'Color',color12(i-1)); hold on
    
end
% xlim([0,1]);
set(gca,'FontSize',16);
ylabel('$z/\delta_i$','FontSize',18);
xlabel('$U/U_i$','FontSize',18);

figure();
for i = 2:N_u-1
    
        plot(U{i}./U_infty_i(i),z./delta_ibl_withAF(i),'^','MarkerSize',8,...
        'Color',color12(i-1)); hold on
    
end
xlim([0,1]);
ylim([0,1]);
set(gca,'FontSize',16);
ylabel('$z/\delta_i$','FontSize',18);
xlabel('$U/U_i$','FontSize',18);


%% ZS Scaling

% Usual Values
for i = 2:N_u-1
    this_i = find(z./delta >= 1.0, 1);
    this_z = z(1:this_i);
    
    U_infty(i) = max(U{i});
    this_integrand = 1 - (U{i}./U_infty(i));
    this_integrand = this_integrand(1:this_i);
    
    delta_star(i) = trapz(this_z,this_integrand);
    
    u0(i) = U_infty(i)*(delta_star(i)/delta);
    
end

% Now with IBL Parameters
for i = 2:N_u-1
    thisZ = z;
    this_zibl = thisZ./delta_ibl(i);
    
    this_i = find(this_zibl >= 1.0, 1);
    this_z = thisZ(1:this_i);
    
    this_integrand = 1 - (U{i}./U_infty_i(i));
    this_integrand = this_integrand(1:this_i);
    
    delta_i_star(i) = trapz(this_z,this_integrand);
    
    u0_i(i) = U_infty_i(i)*(delta_i_star(i)/delta_ibl(i));
    
end    

%%

close all;

figure();
subplot(2,2,1);
for i = 2:N_u-1
   
    semilogx(z./delta,(U{i}(end) - U{i})./utau(i),...
        '^','MarkerSize',8,'Color',color12(i-1)); hold on
end
semilogx(zdelta,wake,'k--','LineWidth',3);
set(gca,'FontSize',16);
% xlim([0 1]);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);
title('Classic Scaling with Local $u_\tau$');
legend({'$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$',...
    '$\hat{x}_4$','$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$',...
    '$\hat{x}_8$','$\hat{x}_9$'},'Interpreter','Latex',...
    'Location','NorthEast','NumColumns',3,'FontSize',10);

subplot(2,2,2);
for i = 2:N_u-1
   
    semilogx(z./delta_ibl(i),(U_infty_i(i) - U{i})./utau(i),...
        '^','MarkerSize',8,'Color',color12(i-1)); hold on
end
semilogx(thisYD,wake,'k--','LineWidth',3);
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$');
% xlim([0 1]);

subplot(2,2,3);
for i = 2:N_u-1
    thisZS = (U_infty(i) - U{i})./u0(i);
    
    semilogx(zdelta,thisZS,...
        '^','MarkerSize',8,'Color',color12(i-1)); hold on
end
set(gca,'FontSize',16);
xlabel('$z/\delta$','FontSize',18);
ylabel('$(U_\infty - U)/U_\infty\delta^*/\delta$','FontSize',18);
title('ZS Scaling');
% xlim([0 1]);


subplot(2,2,4);
for i = 2:N_u-1
    thisZS = (U_infty_i(i) - U{i})./u0_i(i);
    thisZ = z;
    this_zibl = thisZ./delta_ibl(i);
    
    semilogx(this_zibl,thisZS,'^','MarkerSize',8,...
        'Color',color12(i-1)); hold on
end
set(gca,'FontSize',16);
% xlim([0 1]);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/U_i\delta_i^*/\delta_i$','FontSize',18);
title('ZS Scaling with IBL Parameters');


%%%%%% Linear Plotting

figure();
subplot(2,2,1);
for i = 2:N_u-1
   
    plot(z./delta,(U{i}(end) - U{i})./utau(i),...
        '^','MarkerSize',8,'Color',color12(i-1)); hold on
end
plot(zdelta,wake,'k--','LineWidth',3);
set(gca,'FontSize',16);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);
xlim([0 1]);
title('Classic Scaling with Local $u_\tau$');
legend({'$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$',...
    '$\hat{x}_4$','$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$',...
    '$\hat{x}_8$','$\hat{x}_9$'},'Interpreter','Latex',...
    'Location','NorthEast','NumColumns',3,'FontSize',10);

subplot(2,2,2);
for i = 2:N_u-1
   
    plot(z./delta_ibl(i),(U_infty_i(i) - U{i})./utau(i),...
        '^','MarkerSize',8,'Color',color12(i-1)); hold on
end
plot(thisYD,wake,'k--','LineWidth',3);
set(gca,'FontSize',16);
xlim([0 1]);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$');

subplot(2,2,3);
for i = 2:N_u-1
    thisZS = (U_infty(i) - U{i})./u0(i);
    
    plot(zdelta,thisZS,...
        '^','MarkerSize',8,'Color',color12(i-1)); hold on
end
set(gca,'FontSize',16);
xlim([0 1]);
xlabel('$z/\delta$','FontSize',18);
ylabel('$(U_\infty - U)/U_\infty\delta^*/\delta$','FontSize',18);
title('ZS Scaling');


subplot(2,2,4);
for i = 2:N_u-1
    thisZS = (U_infty_i(i) - U{i})./u0_i(i);
    thisZ = z;
    this_zibl = thisZ./delta_ibl(i);
    
    plot(this_zibl,thisZS,'^','MarkerSize',8,...
        'Color',color12(i-1)); hold on
end
set(gca,'FontSize',16);
xlim([0 1]);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/U_i\delta_i^*/\delta_i$','FontSize',18);
title('ZS Scaling with IBL Parameters');

%%
% figure();
% for i = 1:N_u-1
%    
%     plot((U_infty_i(i) - U{i})./utau(i),z./delta_ibl(i),...
%         'LineWidth',2,'Color',myColors_AF(i)); hold on
% end
% set(gca,'FontSize',16);
% title('WMLES Data','FontSize',20);
% xlabel('$z/\delta_i$','FontSize',18);
% ylabel('$U_\infty(z = \delta_i) - \langle U\rangle / u_{\tau,local}$',...
%     'FontSize',18);



%% Plot RSS

% close all;

Cooke_Colors = ["White","#babd00","#92f240","#00de69","#00b995","#00999c",...
    "#007c92","#00627f","#004a70","#0b1b84","White"]; % Yellow to Blue Gradient

% myColors = ["#daf8e3", "#97ebdb", "#00c2c7", "#0086ad", "#005582", ...
%     "#ffc100", "#ff9a00", "#ff7400", "#bf0000" ];

delta_ibl = [0.06356,0.1665,0.1665,0.1373,0.1372,0.2968,0.4362,...
    0.4362,0.5289].*delta;

uw_plot = cell2mat(UW');

figure();
tiledlayout(1,2);
p1 = nexttile;
%plot(uw_plot(1,:),z,'k--','LineWidth',2); hold on;
for i = 2:10
        plot(uw_plot(i,:),z,'Color',Cooke_Colors(i),'LineWidth',2); hold on;
end
set(gca,'FontName','SansSerif','FontSize',20);
ylabel('$z$','FontName','SansSerif','FontSize',36);
xlabel('$\langle u^\prime w^\prime \rangle$','FontName','SansSerif','FontSize',36);
ylim([0 200]);
legend({'$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$',...
    '$\hat{x}_4$','$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$',...
    '$\hat{x}_8$','$\hat{x}_9$'},'Interpreter','Latex',...
    'Location','NorthWest','NumColumns',3,'FontSize',30);


uw_plot2 = uw_plot(2:end,:);

p2 = nexttile;
% figure();
%plot(uw_plot(1,:)./(0.12^2),z./30,'k--','LineWidth',2); hold on;
for i = 2:10
    if i == 2
        plot(uw_plot2(i,:)./(0.12^2),z./30,'Color',Cooke_Colors(i),'LineWidth',2); hold on;
    else
        plot(uw_plot2(i,:)./(0.12^2),z./cookeCorr(i),'Color',Cooke_Colors(i),'LineWidth',2); hold on;
    end
end
set(gca,'FontName','SansSerif','FontSize',20);
ylabel('$z/\hat{\delta}$','FontName','SansSerif','FontSize',36);
xlabel('$\langle u^\prime w^\prime \rangle/u^2_{\tau,0}$','FontName','SansSerif','FontSize',36);
ylim([0 4]);

%%

figure();
for i = 2:10
    if i == 2
        plot(uw_plot(i,:)./(0.12^2),z/300,'Color',Cooke_Colors(i),'LineWidth',2); hold on;
    else
        plot(uw_plot(i,:)./(0.12^2),z/300,'Color',Cooke_Colors(i),'LineWidth',2); hold on;
    end
end
set(gca,'FontName','SansSerif','FontSize',20);
ylabel('$z/\delta$','FontName','SansSerif','FontSize',36);
xlabel('$\langle u^\prime w^\prime \rangle/u^2_{\tau,1}$','FontName','SansSerif','FontSize',36);
ylim([0 1]);


%% 

[m,n] = size(uw{1});

thisInd = 6;

dt = 0.0065*50; %0.0065*100;
tend = dt*m;
t = linspace(0,tend,m);
lett = t*U_infty(thisInd)/delta;


[T,Z] = meshgrid(lett,zdelta);

%% Kymograph 

close all;



figure();
contourf(T,Z,uw{thisInd}','LineColor','None'); hold on;
%yline(delta_ibl(thisInd-1)/delta,'--','LineWidth',3,'Color','#FA70CB');
colorbar;
colormap(flipud(jet));
caxis([-0.05 0]);
set(gca,'FontSize',20);
% set(gca,'XScale','log');
thisTitle = strcat('$\hat{x}_',num2str(thisInd-1),'$');
title(thisTitle,'FontSize',28);
xlabel('$T$','FontSize',36); % T = tU_\infty/\delta
ylabel('$z/\delta$','FontSize',36);
xlim([0 25]);
ylim([zdelta(1) 0.3]);


%% Functions

function wake = defectLaw(kappa,myPI,yd)
    wake =  (1/kappa)*( -1*log(yd) + myPI*(2 - 2*(sin((pi()/2)*yd).^2)));
end

%% End