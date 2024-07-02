%% START

clc;
clear;
close all;

set(0,'defaultTextInterpreter','latex');

%% Load Data

BL_Data = load('./Experimental Data/Li et al 2021/Re07ks16/Re07ks16_BL');

xhat = BL_Data(:,1);
Uinfty = BL_Data(:,2);
utau = BL_Data(:,3);
nu = BL_Data(:,4);
delta99 = BL_Data(:,end);

myDir = dir('./Experimental Data/Li et al 2021/Re07ks16/Re07ks16_xhat*');

myN = length(myDir);

for i = 1:myN
   
    myName = myDir(i).name;
    myFolder = myDir(i).folder;
    
    data = load(strcat(myFolder,'/',myName));
    
    zplus{i} = flipud(data(:,1));
    zdel99{i} = flipud(data(:,2));
    Uplus{i} = flipud(data(:,3));
    UvelDef{i} = flipud(data(:,4));
    uuplus{i} = flipud(data(:,end));
    
    N(i) = length(zplus{i});
        
end

N = N';

clear myFolder myDir myName

%% Calculate IBL Height Based on Power Law Fit

for i = 1:myN
   
    uu{i} = uuplus{i}.*(utau(i)*utau(i));
    uuinfty{i} = uu{i}./(Uinfty(i)*Uinfty(i));
    
end

delta0 = 0.11;

xhatdelta0 = xhat./delta0;
logxhatdelta0 = log10(xhatdelta0);

deltai = (delta0*0.094).*((xhatdelta0).^0.77);

%% Velocity

for i = 1:myN
   
    thisZ = zdel99{i};
    
    U{i} = Uplus{i}.*utau(i);
    
    ii = find(thisZ > deltai(i)/delta99(i),1);
    Udeltai(i) = U{i}(ii);
    
    
end


%% Mean Velocity Profiles

close all;

figure();
subplot(1,2,1);
for i = 1:myN
    semilogx(zplus{i},Uplus{i},'^','MarkerSize',8); hold on
end
set(gca,'FontSize',16);
xlabel('$z^+$','FontSize',18);
ylabel('$\langle U\rangle^+$','FontSize',18);
title('Li $\textit{et al.}$ $R \rightarrow S$','FontSize',20)

subplot(1,2,2);
for i = 1:myN
    
    thisZ = zdel99{i};
    
    semilogx(thisZ.*delta99(i)./deltai(i),U{i}./Udeltai(i),'^',...
        'MarkerSize',8); hold on
end
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$\langle U\rangle/U_\infty(z = \delta_i)$','FontSize',18);
title('Li $\textit{et al.}$ $R \rightarrow S$','FontSize',20)

%%

for i = 1:myN

    thisvel = U{i};
    velDef{i} = (Udeltai(i) - thisvel)./utau(i);

end

%%

close all;

figure();
subplot(1,2,1);
for i = 1:myN 
    semilogx(zdel99{i},UvelDef{i},'^',...
        'MarkerSize',8); hold on
end
set(gca,'FontSize',16);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);

subplot(1,2,2);
for i = 1:myN
    thisZ = zdel99{i};   
    semilogx(thisZ.*delta99(i)./deltai(i),velDef{i},'^',...
        'MarkerSize',8); hold on
end
set(gca,'FontSize',16);
ylabel('$(U_\infty(z=\delta_i) - U)/u_{\tau,2}$','FontSize',18)
xlabel('$z/\delta_i$','FontSize',18);


%% End