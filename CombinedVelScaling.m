%% START

clc;
clear;
close all;

set(0,'defaultTextInterpreter','latex');

%% Load Gul Data

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
P24toS_data = P24toS_data(2:end,:); % Removes x = -1.5 data point
P24toS_xhat = P24toS_data(:,1);
P24toS_dibl = P24toS_data(:,2).*P24toS_delta0;
P24toS_utau2 = P24toS_data(:,3).*P24toS_utau1;

myDir = dir('./ExperimentalData/Gul 2022/P24toS/xhat*');

P24toS_N = length(P24toS_xhat);

for i = 1:P24toS_N
   
    thisName = myDir(i).name;
    thisFolder = myDir.folder;
    
    thisFile = strcat(thisFolder,'/',thisName);
    thisData = load(thisFile);
    
    P24toS_ydelta{i} = thisData(:,1);
    P24toS_velDefect{i} = thisData(:,2);
    
end

clear this*

P60toP24_data = load('./ExperimentalData/Gul 2022/P60toP24/P60toP24.txt');
P60toP24_data = P60toP24_data(2:end,:); % remove first data point at x=-1.5
P60toP24_xhat = P60toP24_data(:,1);
P60toP24_dibl = P60toP24_data(:,2).*P60toP24_delta0;
P60toP24_utau2 = P60toP24_data(:,3).*P60toP24_utau1;

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

%% Load Li Data

BL_Data = load('./ExperimentalData/Li et al 2021/Re07ks16/Re07ks16_BL');

Li_xhat = BL_Data(:,1);
Li_Uinfty = BL_Data(:,2);
Li_utau = BL_Data(:,3);
Li_nu = BL_Data(:,4);
Li_delta99 = BL_Data(:,end);

myDir = dir('./ExperimentalData/Li et al 2021/Re07ks16/Re07ks16_xhat*');

Li_N = length(myDir);

for i = 1:Li_N
   
    myName = myDir(i).name;
    myFolder = myDir(i).folder;
    
    data = load(strcat(myFolder,'/',myName));
    
    Li_zplus{i} = flipud(data(:,1));
    Li_zdel99{i} = flipud(data(:,2));
    Li_Uplus{i} = flipud(data(:,3));
    Li_UvelDef{i} = flipud(data(:,4));
    Li_uuplus{i} = flipud(data(:,end));
    
    Li_Nx(i) = length(Li_zplus{i});
        
end

Li_Nx = Li_Nx';

clear myFolder myDir myName

for i = 1:Li_N
   
    Li_uu{i} = Li_uuplus{i}.*(Li_utau(i)*Li_utau(i));
    Li_uuinfty{i} = Li_uu{i}./(Li_Uinfty(i)*Li_Uinfty(i));
    
end

Li_delta0 = 0.11;

Li_xhatdelta0 = Li_xhat./Li_delta0;
Li_logxhatdelta0 = log10(Li_xhatdelta0);

Li_deltai = (Li_delta0*0.094).*((Li_xhatdelta0).^0.77);

%% Li Velocity

for i = 1:Li_N
   
    thisZ = Li_zdel99{i};
    
    Li_U{i} = Li_Uplus{i}.*Li_utau(i);
    
    ii = find(thisZ > Li_deltai(i)/Li_delta99(i),1);
    Li_Udeltai(i) = Li_U{i}(ii);
    
    
end

%% Gul Velocity

for i = 1:P24toS_N
   
    thisP24toS_Defect = P24toS_velDefect{i};
    
    P24toS_U{i} = P24toS_Uinf - (thisP24toS_Defect.*P24toS_utau1);
    P24toS_ydelta0{i} = P24toS_ydelta{i}.*P24toS_delta0;
    
    
end

for i = 1:Np60
   
    thisDefect = P60toP24_velDefect{i};
    
    P60toP24_U{i} = P60toP24_Uinf - (thisDefect.*P60toP24_utau1);
    P60toP24_ydelta0{i} = P60toP24_ydelta{i}.*P60toP24_delta0;
    
    
end

clear this* ii


%% Need to calculate an approximated local BL height for Gul Velocity

nu = 1.52e-5; %Nu for air at T = 20C (just a guess)
P24toS_dnu = nu./P24toS_utau2;
P60toP24_dnu = nu./P60toP24_utau2;

% using delta/x = 0.16/(Re_x)^1/7

for i = 1:P24toS_N
    
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
% Raw Calculated Values
delta_ibl_withAF = [0.1,0.06356,0.1665,0.1665,0.1373,0.1372,0.2968,0.4362,...
    0.4362,0.5289].*delta;
nu = 0.32e-4; % Kin. viscosity of air
rho = 1.23; % Density of air 

z = linspace(1.5,400,100);
zdelta = z./delta;
xhat = x - 1850;

% Using the correlation delta_i = a*xhat^b which is fit to calc'd data
% Here a = 0.29 and b = 0.71 -- found prev. 
cookeCorr = 0.29.*((xhat(2:end)).^0.71);
% Because the IBL underpredicts the relevant length-scale at the first two
% xhat, we replace with ASL height.
cookeCorr(1) = 30; 
cookeCorr(2) = 30;

Nx = length(xhat);
myColors = [ "black", "#daf8e3", "#97ebdb", "#00c2c7", "#0086ad", "#005582", ...
    "#ffc100", "#ff9a00", "#ff7400", "#bf0000" ];

%% Cooke Data Wall Stress

% Want to create a spanwise averaged u_tau for each xhat

% Spanwise station based averaging so now we have time and space average

for i = 1:Nx-1
   
    swss_taubar(i) = mean(swss_tauw{i},'All');
    utau_taubar(i) = sqrt(swss_taubar(i)/rho);
    dnu_taubar(i)  = nu/utau_taubar(i);
        
end

utau = [utau_taubar(1), utau_taubar];

%% Create Averages for Cooke Velocity

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

%% Find U_infty at z(\delta_i) for Cooke data

%for i = 1:Nx-1
for i = 2:Nx-1
    
    %ii = find(z >= delta_ibl_withAF(i),1);
    ii = find(z >= cookeCorr(i-1),1);
    U_infty_i(i) = U{i}(ii);
    
end


%% Colors for plotting

% Orig starting color for y-b gradient - wasn't showing - "#fbff2b" 

% Gul Colors
P60toP24_Colors = ["#babd00","#3ce758","#00b398","#008195","#005978",...
    "#0b1b84"]; % Yellow to Blue Gradient 

%["white","#e8ffe0",'#bae1ac','#8cc27a','#5ea548','#288708']; % Green
%Gradient


% Li Colors
Li_7k_Colors = ["#babd00","#b0f736","#57ea52","#00d776","#00bd93",...
    "#00a49c","#008e9a","#007a91","#006783","#005575","#004370",...
    "#0b1b84"]; % Yellow to Blue Gradient 

        %["#b8dbff",'#a5c8f1','#92b5e2','#81a2d4','#708fc5',...
    %'#607db7','#516ba8','#435999','#36488b','#29377c','#1c266c','#0e165d']; % Blue Gradient

% Cooke Colors
Cooke_Colors = ["White","#babd00","#92f240","#00de69","#00b995","#00999c",...
    "#007c92","#00627f","#004a70","#0b1b84","White"]; % Yellow to Blue Gradient

        % ["#bbaeef",'#b2a1ea','#a994e4','#a187de','#997bd8','#926dd1',...
   % '#8a60ca','#8353c3','#7c44bc','#7635b4','#6f24ac','#690aa4']; % Purple Gradient




%% Plot Velocity Profiles in Sub-Plots

% 3x2 figure like:
% Cooke Viscous; Cooke Outer
% Gul Viscous; Gul Outer
% Li Viscous; Li Outer

close all;

figure();
tiledlayout(2,3)
p1 = nexttile;
for i = 2:N_u-1
    semilogx(z./(nu/utau(i)),U{i}./utau(i),'k^','MarkerSize',8,...
        'MarkerFaceColor',Cooke_Colors(i)); hold on
end
set(gca,'FontSize',18);
grid on;
xlabel('$z^+$');
ylabel('$U^+$');
title('Cooke24');

p3 = nexttile;
for i = 1:Np60
    semilogx(P60toP24_y{i}./P60toP24_dnu(i),...
        P60toP24_U{i}./P60toP24_utau2(i),'ksquare','MarkerSize',8,...
        'MarkerFaceColor',P60toP24_Colors(i)); hold on;
end
set(gca,'FontSize',18);
grid on;
xlabel('$z^+$');
title('Gul22');
%ylabel('$U^+$','FontSize',18);

p5 = nexttile;
for i = 1:Li_N
    semilogx(Li_zplus{i},Li_Uplus{i},'ko','MarkerSize',8,...
        'MarkerFaceColor',Li_7k_Colors(i)); hold on;
end
set(gca,'FontSize',18);
grid on;
xlabel('$z^+$');
title('Li21');
%ylabel('$U^+$','FontSize',18);

p2 = nexttile;
for i = 2:N_u-1
    plot(U{i}./max(U{i}),zdelta,'k^','MarkerSize',8,...
        'MarkerFaceColor',Cooke_Colors(i)); hold on
end
set(gca,'FontSize',18);
grid on;
xlim([0 1]);
ylim([0 1]);
xlabel('$z/\delta$');
ylabel('$U/U_\infty$','FontSize',18);



p4 = nexttile;
for i = 1:Np60
    thisP60U = P60toP24_U{i};
    plot(thisP60U./P60toP24_Uinf,P60toP24_ydelta{i},'ksquare','MarkerSize',8,...
        'MarkerFaceColor',P60toP24_Colors(i)); hold on;
end
set(gca,'FontSize',18);
grid on;
xlim([0 1]);
ylim([0 1]);
xlabel('$z/\delta$');
%ylabel('$U/U_\infty$','FontSize',18);



p6 = nexttile;
for i = 1:Li_N
    thisZ = Li_zdel99{i};
    plot(Li_U{i}./Li_Uinfty(i),thisZ,...
        'ko','MarkerSize',8,'MarkerFaceColor',Li_7k_Colors(i)); hold on;
end
set(gca,'FontSize',18);
xlim([0 1]);
ylim([0 1]);
grid on;
xlabel('$z/\delta$');
%ylabel('$U/U_\infty$','FontSize',18);

%% IBL Figures to Create Insets

thisLine = [0,1];
theseOnes = [1,1];

close all;

figure();
tiledlayout(1,3);
p1 = nexttile;
for i = 2:N_u-1
    plot(U{i}./U_infty_i(i),z./cookeCorr(i-1),'k^','MarkerSize',8,...
        'MarkerFaceColor',Cooke_Colors(i)); hold on
end
set(gca,'FontSize',24);
grid on;
xlim([0 1]);
ylim([0 1]);
ylabel('$z/\delta_i$','FontSize',36);
xlabel('$\langle U \rangle/U_i$','FontSize',36);

p2 = nexttile;
for i = 1:Np60
    thisP60U = P60toP24_U{i};
    plot(thisP60U./P60toP24_UinftyIBL(i),P60toP24_yibl{i},'ksquare','MarkerSize',8,...
        'MarkerFaceColor',P60toP24_Colors(i)); hold on;
end
set(gca,'FontSize',24);
grid on;
xlim([0 1]);
ylim([0 1]);
% ylabel('$z/\delta_i$','FontSize',18);
xlabel('$\langle U \rangle/U_i$','FontSize',36);

p3 = nexttile;
for i = 1:Li_N
    thisZ = Li_zdel99{i};
    plot(Li_U{i}./Li_Udeltai(i),thisZ.*Li_delta99(i)./Li_deltai(i),...
        'ko','MarkerSize',8,'MarkerFaceColor',Li_7k_Colors(i)); hold on;
end
set(gca,'FontSize',24);
grid on;
xlim([0 1]);
ylim([0 1]);
% ylabel('$z/\delta_i$','FontSize',18);
xlabel('$\langle U \rangle/U_i$','FontSize',36);


figure();
for i = 1:Np60
    thisP60U = P60toP24_U{i};
    plot(thisP60U./P60toP24_UinftyIBL(i),P60toP24_yibl{i},'ksquare','MarkerSize',10,...
        'MarkerFaceColor',P60toP24_Colors(i)); hold on;
end
plot(theseOnes,thisLine,'r-','LineWidth',2);
plot(thisLine,theseOnes,'r-','LineWidth',2);
set(gca,'FontSize',24);
grid on;
ylabel('$z/\delta_i$','FontSize',36);
xlabel('$\langle U \rangle/U_i$','FontSize',36);

figure();
for i = 1:Li_N
    thisZ = Li_zdel99{i};
    plot(Li_U{i}./Li_Udeltai(i),thisZ.*Li_delta99(i)./Li_deltai(i),...
        'ko','MarkerSize',10,'MarkerFaceColor',Li_7k_Colors(i)); hold on;
end
plot(theseOnes,thisLine,'r-','LineWidth',2);
plot(thisLine,theseOnes,'r-','LineWidth',2);
set(gca,'FontSize',24);
grid on;
xlim([0 1.5]);
ylabel('$z/\delta_i$','FontSize',36);
xlabel('$\langle U \rangle/U_i$','FontSize',36);

figure();
for i = 2:N_u-1
    plot(U{i}./U_infty_i(i),z./cookeCorr(i-1),'k^','MarkerSize',10,...
        'MarkerFaceColor',Cooke_Colors(i)); hold on
end
plot(theseOnes,thisLine,'r-','LineWidth',2);
plot(thisLine,theseOnes,'r-','LineWidth',2);
set(gca,'FontSize',24);
grid on;
ylabel('$z/\delta_i$','FontSize',36);
xlabel('$\langle U \rangle/U_i$','FontSize',36);



%% Velocity defect 

% Wake Parameter PI - Chosen based on downstream roughness
smoothPI = 0.55;
roughPI = 0.70; %%%% Could also do 0.70 

% Von Karman Constant
cookeKappa = 0.41;
gulKappa = 0.39;
liKappa = 0.4;

% Define y/delta input for wake function
cookeYD = zdelta;
cookeYDi = z./cookeCorr(end);

gulYD = P60toP24_ydelta{end};
gulYDi = P60toP24_yibl{end};

liYD = Li_zdel99{end};
liYDi = Li_zdel99{end}.*Li_delta99(end)./Li_deltai(end);

% Cutoff YDi to only have just beyond the IBL height
ii = find(cookeYDi >= 1.2, 1);
cookeYDi = cookeYDi(1:ii);

ii = find(gulYDi >= 1.2,1);
gulYDi = gulYDi(1:ii);

ii = find(liYDi >= 1.2,1);
liYDi = liYDi(1:ii);

% Calculate wake defect law functions
cookeWakeIBL = defectLaw(cookeKappa,roughPI,cookeYDi);
cookeWake = defectLaw(cookeKappa,roughPI,cookeYD);

gulWakeIBL = defectLaw(gulKappa,0.9,gulYDi);
gulWake = defectLaw(gulKappa,0.9,gulYD);

liWakeIBL = defectLaw(liKappa,smoothPI,liYDi);
liWake = defectLaw(liKappa,smoothPI,liYD);
    
%% Calc velocity defect for Li with IBL

for j = 1:Li_N
    thisVel = Li_U{j};
    Li_velDefect{j} = (Li_Udeltai(j) - thisVel)./Li_utau(j);
end


%% Individual Velocity Defect Plots

% 3x2 again
% Cooke: Classic; IBL
% Gul:   Classic; IBL
% Li:    Classic; IBL

close all;

figure();
tiledlayout(3,2)
p1 = nexttile;
for i = 2:N_u-1
   
    semilogx(z./delta,(U{i}(end) - U{i})./utau(i),...
        'k^','MarkerSize',8,'MarkerFaceColor',Cooke_Colors(i)); hold on
end
semilogx(zdelta,cookeWake,'k--','LineWidth',3);
set(gca,'FontSize',20);
grid on;
ylabel('$(U_\infty - \langle U \rangle)/u_{\tau,2}$','FontSize',24)
xline(0.2,'k','LineWidth',2);

p2 = nexttile;
for i = 2:N_u-1
   
    semilogx(z./cookeCorr(i-1),(U_infty_i(i) - U{i})./utau(i),...
        'k^','MarkerSize',8,'MarkerFaceColor',Cooke_Colors(i)); hold on
end
semilogx(cookeYDi,cookeWakeIBL,'k--','LineWidth',3);
set(gca,'FontSize',20);
grid on;
ylabel('$(U_i - \langle U \rangle)/u_{\tau,2}$','FontSize',24);
xlim([0 1]);
xline(0.2,'k','LineWidth',2);

p3 = nexttile;
for i = 1:Np60
    semilogx(P60toP24_ydelta{i},...
        P60toP24_velDefect{i}.*(P60toP24_utau1/P60toP24_utau2(i)),...
        'ksquare','MarkerSize',8,'MarkerFaceColor',P60toP24_Colors(i)); hold on
end
semilogx(gulYD,gulWake,'k--','LineWidth',3);
set(gca,'FontSize',20);
grid on;
ylabel('$(U_\infty - \langle U \rangle)/u_{\tau,2}$','FontSize',24);
xline(0.2,'k','LineWidth',2);

p4 = nexttile;
for i = 1:Np60
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_UinftyIBL(i) - thisVel)./P60toP24_utau2(i);
    
    semilogx(P60toP24_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,'MarkerFaceColor',P60toP24_Colors(i)); hold on
end
% ylim([0 5]);
semilogx(gulYDi,gulWakeIBL,'k--','LineWidth',3);
set(gca,'FontSize',20);
grid on;
ylabel('$(U_i - \langle U \rangle)/u_{\tau,2}$','FontSize',24);
xlim([10^-1 1]);
xline(0.2,'k','LineWidth',2);

p5 = nexttile;
for i = 1:Li_N
    p1 = semilogx(Li_zdel99{i},Li_UvelDef{i},'ko','MarkerSize',8,...
        'MarkerFaceColor',Li_7k_Colors(i)); hold on
end
semilogx(liYD,liWake,'k--','LineWidth',3);
set(gca,'FontSize',20);
grid on;
ylabel('$(U_\infty - \langle U \rangle)/u_{\tau,2}$','FontSize',24)
xlabel('$z/\delta$','FontSize',36);
xline(0.2,'k','LineWidth',2);

p6 = nexttile;
for i = 1:Li_N
    thisZ = Li_zdel99{i}.*Li_delta99(i)./Li_deltai(i);
    thisDefect = (Li_Udeltai(i) - Li_U{i})./Li_utau(i);
    semilogx(thisZ,thisDefect,'ko','MarkerSize',8,...
        'MarkerFaceColor',Li_7k_Colors(i)); hold on
end
semilogx(liYDi,liWakeIBL,'k--','LineWidth',3);
set(gca,'FontSize',20);
xlim([0 1]);
grid on;
xlabel('$z/\delta_i$','FontSize',36);
ylabel('$(U_i - \langle U \rangle)/u_{\tau,2}$','FontSize',24);
xline(0.2,'k','LineWidth',2);

%% 

% 2x3 Now 
% Classic: Cooke; Gul; Li
% IBL: Cooke; Gul; Li


close all;

figure();
tiledlayout(2,3)
p1 = nexttile;
for i = 2:N_u-1
   
    semilogx(z./delta,(U{i}(end) - U{i})./utau(i),...
        'k^','MarkerSize',8,'MarkerFaceColor',Cooke_Colors(i)); hold on
end
semilogx(zdelta,cookeWake,'k--','LineWidth',3);
set(gca,'FontSize',20);
grid on;
ylabel('$(U_\infty - \langle U \rangle)/u_{\tau,2}$','FontSize',24);
xlabel('$z/\delta$','FontSize',24);
xlim([10^(-1) 2.5]);
%xline(0.2,'k','LineWidth',2);

p3 = nexttile;
for i = 1:Np60
    semilogx(P60toP24_ydelta{i},...
        P60toP24_velDefect{i}.*(P60toP24_utau1/P60toP24_utau2(i)),...
        'ksquare','MarkerSize',8,'MarkerFaceColor',P60toP24_Colors(i)); hold on
end
semilogx(gulYD,gulWake,'k--','LineWidth',3);
set(gca,'FontSize',20);
grid on;
xlabel('$z/\delta$','FontSize',24);
xlim([.75*10^(-1) 2.5]);
%ylabel('$(U_\infty - \langle U \rangle)/u_{\tau,2}$','FontSize',24);
%xline(0.2,'k','LineWidth',2);

p5 = nexttile;
for i = 1:Li_N
    p1 = semilogx(Li_zdel99{i},Li_UvelDef{i},'ko','MarkerSize',8,...
        'MarkerFaceColor',Li_7k_Colors(i)); hold on
end
semilogx(liYD,liWake,'k--','LineWidth',3);
set(gca,'FontSize',20);
grid on;
%ylabel('$(U_\infty - \langle U \rangle)/u_{\tau,2}$','FontSize',24)
xlabel('$z/\delta$','FontSize',24);
xlim([5*10^(-4) 2.5]);
%xline(0.2,'k','LineWidth',2);

p2 = nexttile;
for i = 2:N_u-1
   
    semilogx(z./cookeCorr(i-1),(U_infty_i(i) - U{i})./utau(i),...
        'k^','MarkerSize',8,'MarkerFaceColor',Cooke_Colors(i)); hold on
end
semilogx(cookeYDi,cookeWakeIBL,'k--','LineWidth',3);
set(gca,'FontSize',20);
grid on;
ylabel('$(U_i - \langle U \rangle)/u_{\tau,2}$','FontSize',24);
xlabel('$z/\delta_i$','FontSize',24);
xlim([0 1]);
%xline(0.2,'k','LineWidth',2);

p4 = nexttile;
for i = 1:Np60
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_UinftyIBL(i) - thisVel)./P60toP24_utau2(i);
    
    semilogx(P60toP24_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,'MarkerFaceColor',P60toP24_Colors(i)); hold on
end
% ylim([0 5]);
semilogx(gulYDi,gulWakeIBL,'k--','LineWidth',3);
set(gca,'FontSize',20);
grid on;
xlabel('$z/\delta_i$','FontSize',24);
%ylabel('$(U_i - \langle U \rangle)/u_{\tau,2}$','FontSize',24);
xlim([10^-1 1]);
%xline(0.2,'k','LineWidth',2);

p6 = nexttile;
for i = 1:Li_N
    thisZ = Li_zdel99{i}.*Li_delta99(i)./Li_deltai(i);
    thisDefect = (Li_Udeltai(i) - Li_U{i})./Li_utau(i);
    semilogx(thisZ,thisDefect,'ko','MarkerSize',8,...
        'MarkerFaceColor',Li_7k_Colors(i)); hold on
end
semilogx(liYDi,liWakeIBL,'k--','LineWidth',3);
set(gca,'FontSize',20);
xlim([0 1]);
grid on;
xlabel('$z/\delta_i$','FontSize',24);
%ylabel('$(U_i - \langle U \rangle)/u_{\tau,2}$','FontSize',24);
%xline(0.2,'k','LineWidth',2);

%%

close all;

tiledlayout(1,3);
p2 = nexttile;
for i = 2:N_u-1
    semilogx(z./cookeCorr(i-1),(U_infty_i(i) - U{i})./utau(i),...
        'k^','MarkerSize',8,'MarkerFaceColor',Cooke_Colors(i)); hold on
end
semilogx(cookeYDi,cookeWakeIBL,'k--','LineWidth',3);
set(gca,'FontSize',20);
grid on;
ylabel('$(U_i - \langle U \rangle)/u_{\tau,2}$','FontSize',24);
xlabel('$z/\delta_i$','FontSize',24);
xlim([0 1]);
title('Cooke24');
%xline(0.2,'k','LineWidth',2);

p4 = nexttile;
for i = 1:Np60
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_UinftyIBL(i) - thisVel)./P60toP24_utau2(i);
    
    semilogx(P60toP24_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,'MarkerFaceColor',P60toP24_Colors(i)); hold on
end
% ylim([0 5]);
semilogx(gulYDi,gulWakeIBL,'k--','LineWidth',3);
set(gca,'FontSize',20);
grid on;
xlabel('$z/\delta_i$','FontSize',24);
title('Gul22');
%ylabel('$(U_i - \langle U \rangle)/u_{\tau,2}$','FontSize',24);
xlim([10^-1 1]);
%xline(0.2,'k','LineWidth',2);

p6 = nexttile;
for i = 1:Li_N
    thisZ = Li_zdel99{i}.*Li_delta99(i)./Li_deltai(i);
    thisDefect = (Li_Udeltai(i) - Li_U{i})./Li_utau(i);
    semilogx(thisZ,thisDefect,'ko','MarkerSize',8,...
        'MarkerFaceColor',Li_7k_Colors(i)); hold on
end
semilogx(liYDi,liWakeIBL,'k--','LineWidth',3);
set(gca,'FontSize',20);
xlim([0 1]);
grid on;
xlabel('$z/\delta_i$','FontSize',24);
title('Li21');
%ylabel('$(U_i - \langle U \rangle)/u_{\tau,2}$','FontSize',24);
%xline(0.2,'k','LineWidth',2);

%% Functions

function wake = defectLaw(kappa,myPI,yd)
    wake =  (1/kappa)*( -1*log(yd) + myPI*(2 - 2*(sin((pi()/2)*yd).^2)));
end


%%