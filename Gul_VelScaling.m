%% START

clc;
clear;
close all;

set(0,'defaultTextInterpreter','latex');

%% Load Data

genData = load('./ExperimentalData/Gul 2022/generalData.txt');

% Rows:
% P24 to S 
% P36 to S 
% P60 to S
% P24 to P60
% P36 to P60 
% P24 to P36 %%%%% Useless since no b0 or A0 values give
% P60 to P24
% P60 to P36
% P36 to P24 %%%%% Useless since no b0 or A0 values give

% P24 to Smooth
P24toS_Uinf = genData(1,1);
P24toS_utau1 = genData(1,2);
P24toS_delta0 = genData(1,3);

% P36 to Smooth
P36toS_Uinf = genData(2,1);
P36toS_utau1 = genData(2,2);
P36toS_delta0 = genData(2,3);

% P60 to Smooth
P60toS_Uinf = genData(3,1);
P60toS_utau1 = genData(3,2);
P60toS_delta0 = genData(3,3);

% P24 to P60
P24toP60_Uinf = genData(4,1);
P24toP60_utau1 = genData(4,2);
P24toP60_delta0 = genData(4,3);

% P36 to P60
P36toP60_Uinf = genData(5,1);
P36toP60_utau1 = genData(5,2);
P36toP60_delta0 = genData(5,3);

% P60 to P24
P60toP24_Uinf = genData(7,1);
P60toP24_utau1 = genData(7,2);
P60toP24_delta0 = genData(7,3);

% P60 to P36
P60toP36_Uinf = genData(8,1);
P60toP36_utau1 = genData(8,2);
P60toP36_delta0 = genData(8,3);

%%%%%%%%% 
 
% *to*.txt files 

% | xhat/delta_0 | delta_i/delta_0 | utau2/utau1 |
% First row has x -1.5 data so can just remove 

%%%%%%% P24 to Smooth

P24toS_data = load('./ExperimentalData/Gul 2022/P24toS/P24toS.txt');
P24toS_data = P24toS_data(2:end,:);
P24toS_xhat = P24toS_data(:,1);
P24toS_dibl = P24toS_data(:,2).*P24toS_delta0;
P24toS_utau2 = P24toS_data(:,3).*P24toS_utau1;

myDir = dir('./ExperimentalData/Gul 2022/P24toS/xhat*');

N = length(P24toS_xhat);

for i = 2:N+1
   
    thisName = myDir(i).name;
    thisFolder = myDir.folder;
    
    thisFile = strcat(thisFolder,'/',thisName);
    thisData = load(thisFile);
    
    P24toS_ydelta{i-1} = thisData(:,1);
    P24toS_velDefect{i-1} = thisData(:,2);
    
    thisDefect = P24toS_velDefect{i-1};
    
    P24toS_U{i-1} = P24toS_Uinf - (thisDefect.*P24toS_utau1);
    P24toS_ydelta0{i-1} = P24toS_ydelta{i-1}.*P24toS_delta0;
    
end

clear this*

%%%%%%% P36 to Smooth

P36toS_data = load('./ExperimentalData/Gul 2022/P36toS/P36toS.txt');
P36toS_data = P36toS_data(2:end,:);
P36toS_xhat = P36toS_data(:,1);
P36toS_dibl = P36toS_data(:,2).*P36toS_delta0;
P36toS_utau2 = P36toS_data(:,3).*P36toS_utau1;

myDir = dir('./ExperimentalData/Gul 2022/P36toS/xhat*');

P36toS_N = length(P36toS_xhat);

for i = 1:P36toS_N
   
    thisName = myDir(i).name;
    thisFolder = myDir.folder;
    
    thisFile = strcat(thisFolder,'/',thisName);
    thisData = load(thisFile);
    
    P36toS_ydelta{i} = thisData(:,1);
    P36toS_velDefect{i} = thisData(:,2);
    
    thisDefect = P36toS_velDefect{i};
    
    P36toS_U{i} = P36toS_Uinf - (thisDefect.*P36toS_utau1);
    P36toS_ydelta0{i} = P36toS_ydelta{i}.*P36toS_delta0;
    
end

clear this*

%%%%%%% P60 to Smooth

P60toS_data = load('./ExperimentalData/Gul 2022/P60toS/P60toS.txt');
P60toS_data = P60toS_data(2:end,:);
P60toS_xhat = P60toS_data(:,1);
P60toS_dibl = P60toS_data(:,2).*P60toS_delta0;
P60toS_utau2 = P60toS_data(:,3).*P60toS_utau1;

myDir = dir('./ExperimentalData/Gul 2022/P60toS/xhat*');

P60toS_N = length(P60toS_xhat);

for i = 1:P60toS_N
   
    thisName = myDir(i).name;
    thisFolder = myDir.folder;
    
    thisFile = strcat(thisFolder,'/',thisName);
    thisData = load(thisFile);
    
    P60toS_ydelta{i} = thisData(:,1);
    P60toS_velDefect{i} = thisData(:,2);
    
    thisDefect = P60toS_velDefect{i};
    
    P60toS_U{i} = P60toS_Uinf - (thisDefect.*P60toS_utau1);
    P60toS_ydelta0{i} = P60toS_ydelta{i}.*P60toS_delta0;
    
end

clear this*


%%%%%%% P24 to P60

P24toP60_data = load('./ExperimentalData/Gul 2022/P24toP60/P24toP60.txt');
P24toP60_data = P24toP60_data(2:end,:);
P24toP60_xhat = P24toP60_data(:,1);
P24toP60_dibl = P24toP60_data(:,2).*P24toP60_delta0;
P24toP60_utau2 = P24toP60_data(:,3).*P24toP60_utau1;

myDir = dir('./ExperimentalData/Gul 2022/P24toP60/xhat*');

P24toP60_N = length(P24toP60_xhat);

for i = 1:P24toP60_N
   
    thisName = myDir(i).name;
    thisFolder = myDir.folder;
    
    thisFile = strcat(thisFolder,'/',thisName);
    thisData = load(thisFile);
    
    P24toP60_ydelta{i} = thisData(:,1);
    P24toP60_velDefect{i} = thisData(:,2);
    
    thisDefect = P24toP60_velDefect{i};
    
    P24toP60_U{i} = P24toP60_Uinf - (thisDefect.*P24toP60_utau1);
    P24toP60_ydelta0{i} = P24toP60_ydelta{i}.*P24toP60_delta0;
    
end

clear this*


%%%%%%% P36 to P60

P36toP60_data = load('./ExperimentalData/Gul 2022/P36toP60/P36toP60.txt');
P36toP60_data = P36toP60_data(2:end,:);
P36toP60_xhat = P36toP60_data(:,1);
P36toP60_dibl = P36toP60_data(:,2).*P36toP60_delta0;
P36toP60_utau2 = P36toP60_data(:,3).*P36toP60_utau1;

myDir = dir('./ExperimentalData/Gul 2022/P36toP60/xhat*');

P36toP60_N = length(P36toP60_xhat);

for i = 1:P36toP60_N
   
    thisName = myDir(i).name;
    thisFolder = myDir.folder;
    
    thisFile = strcat(thisFolder,'/',thisName);
    thisData = load(thisFile);
    
    P36toP60_ydelta{i} = thisData(:,1);
    P36toP60_velDefect{i} = thisData(:,2);
    
    thisDefect = P36toP60_velDefect{i};
    
    P36toP60_U{i} = P36toP60_Uinf - (thisDefect.*P36toP60_utau1);
    P36toP60_ydelta0{i} = P36toP60_ydelta{i}.*P36toP60_delta0;
    
end

clear this*

%%%%%%% P60 to P24

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
    
    thisDefect = P60toP24_velDefect{i};
    
    P60toP24_U{i} = P60toP24_Uinf - (thisDefect.*P60toP24_utau1);
    P60toP24_ydelta0{i} = P60toP24_ydelta{i}.*P60toP24_delta0;
    
end

%%%%%%% P60 to P36

P60toP36_data = load('./ExperimentalData/Gul 2022/P60toP36/P60toP36.txt');
P60toP36_data = P60toP36_data(2:end,:); % remove first data point at x=-1.5
P60toP36_xhat = P60toP36_data(:,1);
P60toP36_dibl = P60toP36_data(:,2).*P60toP36_delta0;
P60toP36_utau2 = P60toP36_data(:,3).*P60toP36_utau1;

myDir = dir('./ExperimentalData/Gul 2022/P60toP36/xhat*');

P60toP36_N = length(P60toP36_xhat);

for i = 1:P60toP36_N
   
    thisName = myDir(i).name;
    thisFolder = myDir.folder;
    
    thisFile = strcat(thisFolder,'/',thisName);
    thisData = load(thisFile);
    
    P60toP36_ydelta{i} = thisData(:,1);
    P60toP36_velDefect{i} = thisData(:,2);
    
    thisDefect = P60toP36_velDefect{i};
    
    P60toP36_U{i} = P60toP36_Uinf - (thisDefect.*P60toP36_utau1);
    P60toP36_ydelta0{i} = P60toP36_ydelta{i}.*P60toP36_delta0;
    
end

clear this* myDir ii

% Gul Colors
Gul_Colors = ["#babd00","#3ce758","#00b398","#008195","#005978",...
    "#0b1b84"]; % Yellow to Blue Gradient 



%% Need to calculate an approximated local BL height

nu = 1.52e-5; %Nu for air at T = 20C (just a guess)

P24toS_dnu = nu./P24toS_utau2;
P36toS_dnu = nu./P36toS_utau2;
P60toS_dnu = nu./P60toS_utau2;
P24toP60_dnu = nu./P24toP60_utau2;
P36toP60_dnu = nu./P36toP60_utau2;
P60toP24_dnu = nu./P60toP24_utau2;
P60toP36_dnu = nu./P60toP36_utau2;

% using delta/x = 0.16/(Re_x)^1/7

%%%%%%%%% P24 to Smooth

for i = 1:N
    
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

%%%%%%%%% P36 to Smooth

for i = 1:N
    
    P36toS_x(i) = P36toS_xhat(i)*P36toS_delta0;
    P36toS_Rex(i) = P36toS_Uinf*P36toS_x(i)/nu;
    P36toS_deltaloc(i) = (0.16*P36toS_x(i))/(P36toS_Rex(i)^(1/7)) ...
        + P36toS_delta0;
    thisYdelta = P36toS_ydelta{i};
    P36toS_y{i} = thisYdelta.*P36toS_deltaloc(i);
    
    P36toS_yibl{i} = P36toS_y{i}./P36toS_dibl(i);
    
    ii = find(P36toS_yibl{i} > 1.0, 1);
    P36toS_UinftyIBL(i) = P36toS_U{i}(ii);
    
end

%%%%%%%%% P60 to Smooth

for i = 1:N
    
    P60toS_x(i) = P60toS_xhat(i)*P60toS_delta0;
    P60toS_Rex(i) = P60toS_Uinf*P60toS_x(i)/nu;
    P60toS_deltaloc(i) = (0.16*P60toS_x(i))/(P60toS_Rex(i)^(1/7)) ...
        + P60toS_delta0;
    thisYdelta = P60toS_ydelta{i};
    P60toS_y{i} = thisYdelta.*P60toS_deltaloc(i);
    
    P60toS_yibl{i} = P60toS_y{i}./P60toS_dibl(i);
    
    ii = find(P60toS_yibl{i} > 1.0, 1);
    P60toS_UinftyIBL(i) = P60toS_U{i}(ii);
    
end

%%%%%%%%% P24 to P60

for i = 1:N
    
    P24toP60_x(i) = P24toP60_xhat(i)*P24toP60_delta0;
    P24toP60_Rex(i) = P24toP60_Uinf*P24toP60_x(i)/nu;
    P24toP60_deltaloc(i) = (0.16*P24toP60_x(i))/(P24toP60_Rex(i)^(1/7)) ...
        + P24toP60_delta0;
    thisYdelta = P24toP60_ydelta{i};
    P24toP60_y{i} = thisYdelta.*P24toP60_deltaloc(i);
    
    P24toP60_yibl{i} = P24toP60_y{i}./P24toP60_dibl(i);
    
    ii = find(P24toP60_yibl{i} > 1.0, 1);
    P24toP60_UinftyIBL(i) = P24toP60_U{i}(ii);
    
end

%%%%%%%%% P36 to P60

for i = 1:N
    
    P36toP60_x(i) = P36toP60_xhat(i)*P36toP60_delta0;
    P36toP60_Rex(i) = P36toP60_Uinf*P36toP60_x(i)/nu;
    P36toP60_deltaloc(i) = (0.16*P36toP60_x(i))/(P36toP60_Rex(i)^(1/7)) ...
        + P36toP60_delta0;
    thisYdelta = P36toP60_ydelta{i};
    P36toP60_y{i} = thisYdelta.*P36toP60_deltaloc(i);
    
    P36toP60_yibl{i} = P36toP60_y{i}./P36toP60_dibl(i);
    
    ii = find(P36toP60_yibl{i} > 1.0, 1);
    P36toP60_UinftyIBL(i) = P36toP60_U{i}(ii);
    
end

%%%%%%%%% P60 to P24

for i = 1:N
    
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

%%%%%%%%% P60 to P36

for i = 1:N
    
    P60toP36_x(i) = P60toP36_xhat(i)*P60toP36_delta0;
    P60toP36_Rex(i) = P60toP36_Uinf*P60toP36_x(i)/nu;
    P60toP36_deltaloc(i) = (0.16*P60toP36_x(i))/(P60toP36_Rex(i)^(1/7)) ...
        + P60toP36_delta0;
    thisYdelta = P60toP36_ydelta{i};
    P60toP36_y{i} = thisYdelta.*P60toP36_deltaloc(i);
    
    P60toP36_yibl{i} = P60toP36_y{i}./P60toP36_dibl(i);
    
    ii = find(P60toP36_yibl{i} > 1.0, 1);
    P60toP36_UinftyIBL(i) = P60toP36_U{i}(ii);
    
end


%%  Mean Velocity Profiles

close all;

%%%%% Top Row Rough to Smooth/Less Rough Datasets Viscous Scaling

figure();
tiledlayout(2,5);
p1 = nexttile;
for i = 1:N
    semilogx(P24toS_y{i}./P24toS_dnu(i),...
        P24toS_U{i}./P24toS_utau2(i),...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z^+$','FontSize',28);
ylabel('$U^+$','FontSize',28);
title('$P24 \rightarrow S$','FontSize',24);

p2 = nexttile;
for i = 1:N
    semilogx(P36toS_y{i}./P36toS_dnu(i),...
        P36toS_U{i}./P36toS_utau2(i),...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z^+$','FontSize',28);
title('$P36 \rightarrow S$','FontSize',24);

p3 = nexttile;
for i = 1:N
    semilogx(P60toS_y{i}./P60toS_dnu(i),...
        P60toS_U{i}./P60toS_utau2(i),...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z^+$','FontSize',28);
title('$P60 \rightarrow S$','FontSize',24);

p4 = nexttile;
for i = 1:N
    semilogx(P24toP60_y{i}./P24toP60_dnu(i),...
        P24toP60_U{i}./P24toP60_utau2(i),...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z^+$','FontSize',28);
title('$P24 \rightarrow P60$','FontSize',24);

p5 = nexttile;
for i = 1:N
    semilogx(P36toP60_y{i}./P36toP60_dnu(i),...
        P36toP60_U{i}./P36toP60_utau2(i),...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z^+$','FontSize',28);
title('$P36 \rightarrow P60$','FontSize',24);

%%%% Bottom Row Rough to Smooth/Less Rough Datasets Outer Scaling


p6 = nexttile;
for i = 1:N
    plot(P24toS_U{i}./P24toS_Uinf,P24toS_ydelta{i},...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
ylabel('$\langle U \rangle/U_\infty$','FontSize',28);

p7 = nexttile;
for i = 1:N
    plot(P36toS_U{i}./P36toS_Uinf,P36toS_ydelta{i},...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);

p8 = nexttile;
for i = 1:N
    plot(P60toS_U{i}./P60toS_Uinf,P60toS_ydelta{i},...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);

p9 = nexttile;
for i = 1:N
    plot(P24toP60_U{i}./P24toP60_Uinf,P24toP60_ydelta{i},...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);

p10 = nexttile;
for i = 1:N
    plot(P36toP60_U{i}./P36toP60_Uinf,P36toP60_ydelta{i},...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);

%% Group Less Rough to More Rough

figure();

tiledlayout(2,2);
%%%% P60 to P24

p1 = nexttile;
for i = 1:N
    semilogx(P60toP24_y{i}./P60toP24_dnu(i),...
        P60toP24_U{i}./P60toP24_utau2(i),...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z^+$','FontSize',28);
ylabel('$U^+$','FontSize',28);
title('$P60 \rightarrow P24$','FontSize',24);

p2 = nexttile;
for i = 1:N
    semilogx(P60toP36_y{i}./P60toP36_dnu(i),...
        P60toP36_U{i}./P60toP36_utau2(i),...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z^+$','FontSize',28);
title('$P60 \rightarrow P36$','FontSize',24);

p4 = nexttile;
for i = 1:N
    plot(P60toP24_U{i}./P60toP24_Uinf,P60toP24_ydelta{i},...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
ylabel('$\langle U \rangle/U_\infty$','FontSize',28);

p4 = nexttile;
for i = 1:N
    plot(P60toP36_U{i}./P60toP36_Uinf,P60toP36_ydelta{i},...
        'ksquare','MarkerSize',8,'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);


%% IBL Scaling

thisLine = [0,1];
theseOnes = [1,1];

close all;

figure();
tiledlayout(2,2);
p1 = nexttile;
for i = 1:N
    plot(P24toS_U{i}./P24toS_UinftyIBL(i),P24toS_yibl{i},...
            'ksquare','MarkerSize',8,...
            'MarkerFaceColor',Gul_Colors(i)); hold on;
end
plot(theseOnes,thisLine,'r-','LineWidth',2);
plot(thisLine,theseOnes,'r-','LineWidth',2);
set(gca,'FontSize',20);
ylabel('$z/\delta_i$','FontSize',28);
xlabel('$U/U_i$','FontSize',28);

p2 = nexttile;
for i = 1:N
    plot(P24toS_U{i}./P24toS_UinftyIBL(i),P24toS_yibl{i},...
            'ksquare','MarkerSize',8,...
            'MarkerFaceColor',Gul_Colors(i)); hold on;
end
xlim([0 1]);
ylim([0 1]);
set(gca,'FontSize',20);
ylabel('$z/\delta_i$','FontSize',28);
xlabel('$U/U_i$','FontSize',28);

p3 = nexttile;
for i = 1:N
    plot(P60toP24_U{i}./P60toP24_UinftyIBL(i),P60toP24_yibl{i},...
            'ksquare','MarkerSize',8,...
            'MarkerFaceColor',Gul_Colors(i)); hold on;
end
plot(theseOnes,thisLine,'r-','LineWidth',2);
plot(thisLine,theseOnes,'r-','LineWidth',2);
set(gca,'FontSize',20);
ylabel('$z/\delta_i$','FontSize',28);
xlabel('$U/U_i$','FontSize',28);

p4 = nexttile;
for i = 1:N
    plot(P60toP24_U{i}./P60toP24_UinftyIBL(i),P60toP24_yibl{i},...
            'ksquare','MarkerSize',8,...
            'MarkerFaceColor',Gul_Colors(i)); hold on;
end
xlim([0 1]);
ylim([0 1]);
set(gca,'FontSize',20);
ylabel('$z/\delta_i$','FontSize',28);
xlabel('$U/U_i$','FontSize',28);

%% IBL For Supp Mat

close all;

figure();
tiledlayout(1,5);

p1 = nexttile;
for i = 1:N
    plot(P24toS_U{i}./P24toS_UinftyIBL(i),P24toS_yibl{i},...
            'ksquare','MarkerSize',8,...
            'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$\langle U \rangle/U_i$','FontSize',28);
title('$P24 \rightarrow S$','FontSize',24);
ylabel('$z/\delta_i$','FontSize',28);
xlim([0 1]);
ylim([0 1]);

p2 = nexttile;
for i = 1:N
    plot(P36toS_U{i}./P36toS_UinftyIBL(i),P36toS_yibl{i},...
            'ksquare','MarkerSize',8,...
            'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$\langle U \rangle/U_i$','FontSize',28);
title('$P36 \rightarrow S$','FontSize',24);
xlim([0 1]);
ylim([0 1]);

p3 = nexttile;
for i = 1:N
    plot(P60toS_U{i}./P60toS_UinftyIBL(i),P60toS_yibl{i},...
            'ksquare','MarkerSize',8,...
            'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$\langle U \rangle/U_i$','FontSize',28);
title('$P60 \rightarrow S$','FontSize',24);
xlim([0 1]);
ylim([0 1]);

p4 = nexttile;
for i = 1:N
    plot(P24toP60_U{i}./P24toP60_UinftyIBL(i),P24toP60_yibl{i},...
            'ksquare','MarkerSize',8,...
            'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$\langle U \rangle/U_i$','FontSize',28);
title('$P24 \rightarrow P60$','FontSize',24);
xlim([0 1]);
ylim([0 1]);

p5 = nexttile;
for i = 1:N
    plot(P36toP60_U{i}./P36toP60_UinftyIBL(i),P36toP60_yibl{i},...
            'ksquare','MarkerSize',8,...
            'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$\langle U \rangle/U_i$','FontSize',28);
title('$P36 \rightarrow P60$','FontSize',24);
xlim([0 1]);
ylim([0 1]);


%%%%% Less Rough to More Rough 

figure();
tiledlayout(1,2);

p1 = nexttile;
for i = 1:N
    plot(P60toP24_U{i}./P60toP24_UinftyIBL(i),P60toP24_yibl{i},...
            'ksquare','MarkerSize',8,...
            'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$\langle U \rangle/U_i$','FontSize',28);
title('$P60 \rightarrow P24$','FontSize',24);
ylabel('$z/\delta_i$','FontSize',28);
xlim([0 1]);
ylim([0 1]);

p2 = nexttile;
for i = 1:N
    plot(P60toP36_U{i}./P60toP36_UinftyIBL(i),P60toP36_yibl{i},...
            'ksquare','MarkerSize',8,...
            'MarkerFaceColor',Gul_Colors(i)); hold on;
end
set(gca,'FontSize',20);
xlabel('$\langle U \rangle/U_i$','FontSize',28);
title('$P60 \rightarrow P36$','FontSize',24);
xlim([0 1]);
ylim([0 1]);

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

P24toS_wakeIBL = defectLaw(kappa,0.3,P24toS_YDi);
P24toS_wake = defectLaw(kappa,0.3,P24toS_YD);

P60toP24_wakeIBL = defectLaw(kappa,roughPI,P60toP24_YDi);
P60toP24_wake = defectLaw(kappa,roughPI,P60toP24_YD);


%% Velocity Defect

close all;


figure();

%%%%%%% P24 to Smooth
tiledlayout(2,2);
p1 = nexttile;
for i = 1:N
    semilogx(P24toS_ydelta{i},...
        P24toS_velDefect{i}.*(P24toS_utau1/P24toS_utau2(i)),...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
semilogx(P24toS_YD,P24toS_wake,'k--','LineWidth',3);
ylim([0 15]);
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
ylabel('$(U_\infty - \langle U \rangle)/u_{\tau,2}$','FontSize',28);


p2 = nexttile;
for i = 1:N
    
    thisVel = P24toS_U{i};
    thisDefect = (P24toS_UinftyIBL(i) - thisVel)./P24toS_utau2(i);
    
    semilogx(P24toS_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
semilogx(P24toS_YDi,P24toS_wakeIBL,'k--','LineWidth',3);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);
ylabel('$(U_i - \langle U \rangle)/u_{\tau,2}$','FontSize',28);

%%%%%%% P60 to P24
p3 = nexttile;
for i = 1:N
    semilogx(P60toP24_ydelta{i},...
        P60toP24_velDefect{i}.*(P60toP24_utau1/P60toP24_utau2(i)),...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
semilogx(P60toP24_YD,P60toP24_wake,'k--','LineWidth',3);
ylim([0 15]);
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
ylabel('$(U_\infty - \langle U \rangle)/u_{\tau,2}$','FontSize',28);



% Finally with IBL Scaling
p4 = nexttile;
for i = 1:N
    
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_UinftyIBL(i) - thisVel)./P60toP24_utau2(i);
    
    semilogx(P60toP24_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
semilogx(P60toP24_YDi,P60toP24_wakeIBL,'k--','LineWidth',3);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);
ylabel('$(U_i - \langle U \rangle)/u_{\tau,2}$','FontSize',28);


%% Defect Sup Mat

%%%%%% Rough to Less Rough Top: Classical Scaling Bot: IBL Scaling

close all;

figure();

tiledlayout(2,5);
p1 = nexttile;
for i = 1:N
    semilogx(P24toS_ydelta{i},...
        P24toS_velDefect{i}.*(P24toS_utau1/P24toS_utau2(i)),...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
ylim([0 15]);
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P24 \rightarrow S$','FontSize',24);
ylabel('$(U_\infty - \langle U \rangle)/u_{\tau,2}$','FontSize',28);

p2 = nexttile;
for i = 1:N
    semilogx(P36toS_ydelta{i},...
        P36toS_velDefect{i}.*(P36toS_utau1/P36toS_utau2(i)),...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
ylim([0 15]);
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P36 \rightarrow S$','FontSize',24);

p3 = nexttile;
for i = 1:N
    semilogx(P60toS_ydelta{i},...
        P60toS_velDefect{i}.*(P60toS_utau1/P60toS_utau2(i)),...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
ylim([0 15]);
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P60 \rightarrow S$','FontSize',24);

p4 = nexttile;
for i = 1:N
    semilogx(P24toP60_ydelta{i},...
        P24toP60_velDefect{i}.*(P24toP60_utau1/P24toP60_utau2(i)),...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
ylim([0 15]);
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P24 \rightarrow P60$','FontSize',24);

p5 = nexttile;
for i = 1:N
    semilogx(P36toP60_ydelta{i},...
        P36toP60_velDefect{i}.*(P36toP60_utau1/P36toP60_utau2(i)),...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
ylim([0 15]);
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P36 \rightarrow P60$','FontSize',24);

%%%%% Bottom Row

p6 = nexttile;
for i = 1:N
    
    thisVel = P24toS_U{i};
    thisDefect = (P24toS_UinftyIBL(i) - thisVel)./P24toS_utau2(i);
    
    semilogx(P24toS_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);
ylabel('$(U_i - \langle U \rangle)/u_{\tau,2}$','FontSize',28);

p7 = nexttile;
for i = 1:N
    
    thisVel = P36toS_U{i};
    thisDefect = (P36toS_UinftyIBL(i) - thisVel)./P36toS_utau2(i);
    
    semilogx(P36toS_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);

p8 = nexttile;
for i = 1:N
    
    thisVel = P60toS_U{i};
    thisDefect = (P60toS_UinftyIBL(i) - thisVel)./P60toS_utau2(i);
    
    semilogx(P60toS_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);

p9 = nexttile;
for i = 1:N
    
    thisVel = P24toP60_U{i};
    thisDefect = (P24toP60_UinftyIBL(i) - thisVel)./P24toP60_utau2(i);
    
    semilogx(P24toP60_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);

p10 = nexttile;
for i = 1:N
    
    thisVel = P36toP60_U{i};
    thisDefect = (P36toP60_UinftyIBL(i) - thisVel)./P36toP60_utau2(i);
    
    semilogx(P36toP60_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);



%%%%%% Group Less to More Rough
figure();

tiledlayout(2,2);
p1 = nexttile;
for i = 1:N
    semilogx(P60toP24_ydelta{i},...
        P60toP24_velDefect{i}.*(P60toP24_utau1/P60toP24_utau2(i)),...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
ylim([0 15]);
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P60 \rightarrow P24$','FontSize',24);
ylabel('$(U_\infty - \langle U \rangle)/u_{\tau,2}$','FontSize',28);

p2 = nexttile;
for i = 1:N
    semilogx(P60toP36_ydelta{i},...
        P60toP36_velDefect{i}.*(P60toP36_utau1/P60toP36_utau2(i)),...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
ylim([0 15]);
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P60 \rightarrow P36$','FontSize',24);

% Finally with IBL Scaling
p3 = nexttile;
for i = 1:N
    
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_UinftyIBL(i) - thisVel)./P60toP24_utau2(i);
    
    semilogx(P60toP24_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);
ylabel('$(U_i - \langle U \rangle)/u_{\tau,2}$','FontSize',28);

p4 = nexttile;
for i = 1:N
    
    thisVel = P60toP36_U{i};
    thisDefect = (P60toP36_UinftyIBL(i) - thisVel)./P60toP36_utau2(i);
    
    semilogx(P60toP36_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);

%% ZS Scaling

% Usual Values
for j = 1:7 % number of different datasets
    
    for i = 1:N
        
        if j == 1
            this_i = find(P24toS_ydelta{i} >= 1.0, 1);
            this_z = P24toS_y{i}(1:this_i);

            this_integrand = 1 - (P24toS_U{i}./P24toS_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_star{j}(i) = trapz(this_z,this_integrand);

            u0{j}(i) = P24toS_Uinf*(delta_star{j}(i)/P24toS_deltaloc(i));
        elseif j == 2
            this_i = find(P36toS_ydelta{i} >= 1.0, 1);
            this_z = P36toS_y{i}(1:this_i);

            this_integrand = 1 - (P36toS_U{i}./P36toS_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_star{j}(i) = trapz(this_z,this_integrand);

            u0{j}(i) = P36toS_Uinf*(delta_star{j}(i)/P36toS_deltaloc(i));
        elseif j == 3
            this_i = find(P60toS_ydelta{i} >= 1.0, 1);
            this_z = P60toS_y{i}(1:this_i);

            this_integrand = 1 - (P60toS_U{i}./P60toS_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_star{j}(i) = trapz(this_z,this_integrand);

            u0{j}(i) = P60toS_Uinf*(delta_star{j}(i)/P60toS_deltaloc(i));
        elseif j == 4
            this_i = find(P24toP60_ydelta{i} >= 1.0, 1);
            this_z = P24toP60_y{i}(1:this_i);

            this_integrand = 1 - (P24toP60_U{i}./P24toP60_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_star{j}(i) = trapz(this_z,this_integrand);

            u0{j}(i) = P24toP60_Uinf*(delta_star{j}(i)/P24toP60_deltaloc(i));
        elseif j == 5
            this_i = find(P36toP60_ydelta{i} >= 1.0, 1);
            this_z = P36toP60_y{i}(1:this_i);

            this_integrand = 1 - (P36toP60_U{i}./P36toP60_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_star{j}(i) = trapz(this_z,this_integrand);

            u0{j}(i) = P36toP60_Uinf*(delta_star{j}(i)/P36toP60_deltaloc(i));
        elseif j == 6
            this_i = find(P60toP24_ydelta{i} >= 1.0, 1);
            this_z = P60toP24_y{i}(1:this_i);

            this_integrand = 1 - (P60toP24_U{i}./P60toP24_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_star{j}(i) = trapz(this_z,this_integrand);

            u0{j}(i) = P60toP24_Uinf*(delta_star{j}(i)/P60toP24_deltaloc(i));
        else
            this_i = find(P60toP36_ydelta{i} >= 1.0, 1);
            this_z = P60toP36_y{i}(1:this_i);

            this_integrand = 1 - (P60toP36_U{i}./P60toP36_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_star{j}(i) = trapz(this_z,this_integrand);

            u0{j}(i) = P60toP36_Uinf*(delta_star{j}(i)/P60toP36_deltaloc(i));
        end

    end
    
end

%%

% Now with IBL Parameters
for j = 1:7
    for i = 2:N
        
        if j == 1 
            this_zibl = P24toS_yibl{i};

            this_i = find(this_zibl >= 1.0, 1);
            this_z = P24toS_y{i}(1:this_i);

            this_integrand = 1 - (P24toS_U{i}./P24toS_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_i_star{j}(i) = trapz(this_z,this_integrand);

            u0_i{j}(i) = P24toS_Uinf*(delta_i_star{j}(i)/P24toS_dibl(i));
        elseif j == 2
            this_zibl = P36toS_yibl{i};

            this_i = find(this_zibl >= 1.0, 1);
            this_z = P36toS_y{i}(1:this_i);

            this_integrand = 1 - (P36toS_U{i}./P36toS_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_i_star{j}(i) = trapz(this_z,this_integrand);

            u0_i{j}(i) = P36toS_Uinf*(delta_i_star{j}(i)/P36toS_dibl(i));
        elseif j == 3
            this_zibl = P60toS_yibl{i};

            this_i = find(this_zibl >= 1.0, 1);
            this_z = P60toS_y{i}(1:this_i);

            this_integrand = 1 - (P60toS_U{i}./P60toS_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_i_star{j}(i) = trapz(this_z,this_integrand);

            u0_i{j}(i) = P60toS_Uinf*(delta_i_star{j}(i)/P60toS_dibl(i));
        elseif j == 4 
            this_zibl = P24toP60_yibl{i};

            this_i = find(this_zibl >= 1.0, 1);
            this_z = P24toP60_y{i}(1:this_i);

            this_integrand = 1 - (P24toP60_U{i}./P24toP60_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_i_star{j}(i) = trapz(this_z,this_integrand);

            u0_i{j}(i) = P24toP60_Uinf*(delta_i_star{j}(i)/P24toP60_dibl(i));
        elseif j == 5
            this_zibl = P36toP60_yibl{i};

            this_i = find(this_zibl >= 1.0, 1);
            this_z = P36toP60_y{i}(1:this_i);

            this_integrand = 1 - (P36toP60_U{i}./P36toP60_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_i_star{j}(i) = trapz(this_z,this_integrand);

            u0_i{j}(i) = P36toP60_Uinf*(delta_i_star{j}(i)/P36toP60_dibl(i));
        elseif j == 6
            this_zibl = P60toP24_yibl{i};

            this_i = find(this_zibl >= 1.0, 1);
            this_z = P60toP24_y{i}(1:this_i);

            this_integrand = 1 - (P60toP24_U{i}./P60toP24_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_i_star{j}(i) = trapz(this_z,this_integrand);

            u0_i{j}(i) = P60toP24_Uinf*(delta_i_star{j}(i)/P60toP24_dibl(i));
        elseif j == 7 && i ~= 1
            this_zibl = P60toP36_yibl{i};

            this_i = find(this_zibl >= 1.0, 1);
            this_z = P60toP36_y{i}(1:this_i);

            this_integrand = 1 - (P60toP36_U{i}./P60toP36_Uinf);
            this_integrand = this_integrand(1:this_i);

            delta_i_star{j}(i) = trapz(this_z,this_integrand);

            u0_i{j}(i) = P60toP36_Uinf*(delta_i_star{j}(i)/P60toP36_dibl(i));
        end
    end
end    

%% Plot

close all;

figure();
%%%%%%%% ZS
tiledlayout(2,5);
nexttile;
thisIND = 1;
for i = 2:N
    
    thisVel = P24toS_U{i};
    thisDefect = (P24toS_Uinf - thisVel)./u0{thisIND}(i);
    
    semilogx(P24toS_ydelta{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P24 \rightarrow S$','FontSize',24);
ylabel('$(U_\infty - \langle U \rangle)/U_\infty\delta^*/\delta$','FontSize',28);

nexttile;
thisIND = 2;
for i = 2:N
    
    thisVel = P36toS_U{i};
    thisDefect = (P36toS_Uinf - thisVel)./u0{thisIND}(i);
    
    semilogx(P36toS_ydelta{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P36 \rightarrow S$','FontSize',24);

nexttile;
thisIND = 3;
for i = 2:N
    
    thisVel = P60toS_U{i};
    thisDefect = (P60toS_Uinf - thisVel)./u0{thisIND}(i);
    
    semilogx(P60toS_ydelta{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P60 \rightarrow S$','FontSize',24);

nexttile;
thisIND = 4;
for i = 2:N
    
    thisVel = P24toP60_U{i};
    thisDefect = (P24toP60_Uinf - thisVel)./u0{thisIND}(i);
    
    semilogx(P24toP60_ydelta{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P24 \rightarrow P60$','FontSize',24);

nexttile;
thisIND = 5;
for i = 2:N
    
    thisVel = P36toP60_U{i};
    thisDefect = (P36toP60_Uinf - thisVel)./u0{thisIND}(i);
    
    semilogx(P36toP60_ydelta{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P36 \rightarrow P60$','FontSize',24);

% Finally with IBL Scaling
nexttile;
thisIND = 1;
for i = 2:N
    
    thisVel = P24toS_U{i};
    thisDefect = (P24toS_UinftyIBL(i) - thisVel)./u0_i{thisIND}(i);
    
    semilogx(P24toS_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);
ylabel('$(U_i - \langle U \rangle)/U_\infty\delta_i^*/\delta_i$','FontSize',28);

nexttile;
thisIND = 2;
for i = 2:N
    
    thisVel = P36toS_U{i};
    thisDefect = (P36toS_UinftyIBL(i) - thisVel)./u0_i{thisIND}(i);
    
    semilogx(P36toS_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);

nexttile;
thisIND = 3;
for i = 2:N
    
    thisVel = P60toS_U{i};
    thisDefect = (P60toS_UinftyIBL(i) - thisVel)./u0_i{thisIND}(i);
    
    semilogx(P60toS_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);

nexttile;
thisIND = 4;
for i = 2:N
    
    thisVel = P24toP60_U{i};
    thisDefect = (P24toP60_UinftyIBL(i) - thisVel)./u0_i{thisIND}(i);
    
    semilogx(P24toP60_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);

nexttile;
thisIND = 5;
for i = 2:N
    
    thisVel = P36toP60_U{i};
    thisDefect = (P36toP60_UinftyIBL(i) - thisVel)./u0_i{thisIND}(i);
    
    semilogx(P36toP60_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);



%%%% Less to More Rough

figure();
%%%%%%%% ZS
tiledlayout(2,2);
nexttile;
thisIND = 6;
for i = 2:N
    
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_Uinf - thisVel)./u0{thisIND}(i);
    
    semilogx(P60toP24_ydelta{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P60 \rightarrow P24$','FontSize',24);
ylabel('$(U_\infty - \langle U \rangle)/U_i\delta^*/\delta$','FontSize',28);

nexttile;
thisIND = 7;
for i = 2:N
    
    thisVel = P60toP36_U{i};
    thisDefect = (P60toP36_Uinf - thisVel)./u0{thisIND}(i);
    
    semilogx(P60toP36_ydelta{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
set(gca,'FontSize',20);
xlabel('$z/\delta$','FontSize',28);
title('$P60 \rightarrow P36$','FontSize',24);

nexttile;
thisIND = 6;
for i = 2:N
    
    thisVel = P60toP24_U{i};
    thisDefect = (P60toP24_UinftyIBL(i) - thisVel)./u0_i{thisIND}(i);
    
    semilogx(P60toP24_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);
ylabel('$(U_i - \langle U \rangle)/U_i\delta_i^*/\delta_i$','FontSize',28);

nexttile;
thisIND = 7;
for i = 2:N
    
    thisVel = P60toP36_U{i};
    thisDefect = (P60toP36_UinftyIBL(i) - thisVel)./u0_i{thisIND}(i);
    
    semilogx(P60toP36_yibl{i},thisDefect,...
        'ksquare','MarkerSize',8,...
        'MarkerFaceColor',Gul_Colors(i)); hold on
end
xlim([10^-1 1]);
set(gca,'FontSize',20);
xlabel('$z/\delta_i$','FontSize',28);


%% Functions

% function wake = defectLaw(kappa,myPI,yd)
%     wake =  (1/kappa)*( -1*log(yd) + myPI*(2 - 2*(sin((pi()/2)*yd).^2)));
% end


%% END