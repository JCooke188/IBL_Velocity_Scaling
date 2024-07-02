%% START

clc;
clear;
close all;

set(0,'defaultTextInterpreter','latex');

%% Load Gul Data

genData = load('./Gul 2022/generalData.txt');

P24toS_Uinf = genData(1,1);
P24toS_utau1 = genData(1,2);
P24toS_delta0 = genData(1,3);

P24toS_data = load('./Gul 2022/P24toS/P24toS.txt');
P24toS_xhat = P24toS_data(:,1);
P24toS_dibl = P24toS_data(:,2).*P24toS_delta0;
P24toS_utau2 = P24toS_data(:,3).*P24toS_utau1;

myDir = dir('./Gul 2022/P24toS/xhat*');

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

%% Load Li Data

BL_Data = load('./Li et al 2021/Re07ks16/Re07ks16_BL');

Li_xhat = BL_Data(:,1);
Li_Uinfty = BL_Data(:,2);
Li_utau = BL_Data(:,3);
P24toS_nu = BL_Data(:,4);
Li_delta99 = BL_Data(:,end);

myDir = dir('./Li et al 2021/Re07ks16/Re07ks16_xhat*');

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

clear this* ii


P24toS_nu = 1.52e-5; %Nu for air at T = 20C (just a guess)
P24toS_dnu = P24toS_nu./P24toS_utau2;

% using delta/x = 0.16/(Re_x)^1/7

for i = 2:P24toS_N
    
    P24toS_x(i) = P24toS_xhat(i)*P24toS_delta0;
    P24toS_Rex(i) = P24toS_Uinf*P24toS_x(i)/P24toS_nu;
    % Attempt to find a local BL height -- based on correlation above,
    % which is added to delta_0 
    P24toS_deltaloc(i) = (0.16*P24toS_x(i))/(P24toS_Rex(i)^(1/7)) ...
        + P24toS_delta0;
    thisYdelta = P24toS_ydelta{i};
    P24toS_y{i} = thisYdelta.*P24toS_deltaloc(i);
    
    P24toS_yibl{i} = P24toS_y{i}./P24toS_dibl(i);
    P24toS_yplus{i} = P24toS_y{i}./P24toS_dnu(i);
    
    ii = find(P24toS_yibl{i} > 1.0, 1);
    P24toS_UinftyIBL(i) = P24toS_U{i}(ii);
    
end

%% Load Dune Field Data

myDir = dir('/home/jpcooke/Desktop/Research/Dunes/Quadrant Analysis/Data/Sept13/x*');

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

myDir = dir('/home/jpcooke/Desktop/Research/Dunes/Amplitude Modulation/DuneField/WSS/SWSS*');

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

%% Plot everything

% Viscous Scaling
figure();
subplot(1,2,1);
for i = 2:P24toS_N
    semilogx(P24toS_yplus{i},...
        P24toS_U{i}./P24toS_utau2(i),'-','LineWidth',2,...
        'Color','#293256'); hold on;
end
for i = 1:Li_N
    semilogx(Li_zplus{i},Li_Uplus{i},'-.','LineWidth',2,...
        'Color','#0094ab'); hold on;
end
for i = 1:N_u-1
    semilogx(z./(nu/utau(i)),U{i}./utau(i),'--','LineWidth',2,...
        'Color','#70fa97'); hold on
end
set(gca,'FontSize',16);
xlabel('$z^+$','FontSize',18);
ylabel('$U^+$','FontSize',18);

subplot(1,2,2);
for i = 2:P24toS_N
    thisP24U = P24toS_U{i};
    semilogx(P24toS_yibl{i},thisP24U./P24toS_UinftyIBL(i),...
        '-','LineWidth',2,'Color','#293256'); hold on;
end
for i = 1:Li_N
    thisZ = Li_zdel99{i};
    semilogx(thisZ.*Li_delta99(i)./Li_deltai(i),Li_U{i}./Li_Udeltai(i),...
        '-.','LineWidth',2,'Color','#0094ab'); hold on;
end
for i = 1:N_u-1
    semilogx(z./delta_ibl_withAF(i),U{i}./U_infty_i(i),'--','LineWidth',2,...
        'Color','#70fa97'); hold on
end
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$U/U_\infty(z=\delta_i)$','FontSize',18);




%% Velocity Defect Plot

close all;

figure();
subplot(1,2,1);
for i = 2:P24toS_N
    semilogx(P24toS_ydelta{i},P24toS_velDefect{i},'-','LineWidth',2,...
        'Color','#293256'); hold on;
end
for i = 1:Li_N
    semilogx(Li_zdel99{i},Li_UvelDef{i},'-.','LineWidth',2,...
        'Color','#0094ab'); hold on;
end
for i = 1:N_u-1
    semilogx(z./delta,(U{i}(end) - U{i})./utau(i),'--','LineWidth',2,...
        'Color','#70fa97'); hold on
end
set(gca,'FontSize',16);
xlabel('$z/\delta$','FontSize',18);
ylabel('$(U_\infty - \langle U \rangle)/u_\tau)$','FontSize',18);

subplot(1,2,2);
for i = 2:P24toS_N
    thisP24U = P24toS_U{i};
    thisP24toS_Defect = (P24toS_UinftyIBL(i) - thisP24U)./P24toS_utau2(i);
    
    p1 = semilogx(P24toS_yibl{i},thisP24toS_Defect,...
        '-','LineWidth',2,'Color','#293256'); hold on;
end
for i = 1:Li_N
    thisLi_vel = Li_U{i};
    thisLi_Defect = (Li_Udeltai(i) - thisLi_vel)./Li_utau(i);

    thisZ = Li_zdel99{i};
    p2 = semilogx(thisZ.*Li_delta99(i)./Li_deltai(i),thisLi_Defect,...
        '-.','LineWidth',2,'Color','#0094ab'); hold on;
end
for i = 1:N_u-1
    p3 = semilogx(z./delta_ibl_withAF(i),(U_infty_i(i) - U{i})./utau(i),'--','LineWidth',2,...
        'Color','#70fa97'); hold on
end
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_\infty(z=\delta_i) - \langle U \rangle)/u_\tau)$','FontSize',18);
legend([p1 p2 p3],{'Gul 2022 Expt. $R \rightarrow S$',...
    'Li 2021 Expt $R \rightarrow S$','Cooke WMLES 2024 $S \rightarrow R$'},...
    'Interpreter','Latex');

%%