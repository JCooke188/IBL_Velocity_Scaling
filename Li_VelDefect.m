%% START

clc;
clear;
close all;

set(0,'defaultTextInterpreter','latex');

%% Load BL Based Data

rootdir = './Experimental Data/Li et al 2021/';
filelist = dir(fullfile(rootdir,'Re*/*BL'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list

N_bl = length(filelist);

for i = 1:N_bl
   
    thisname = filelist(i).name;
    thisdir = filelist(i).folder;
    thisfile = strcat(thisdir,'/',thisname);
    
    BL_Data = load(thisfile);
    
    xhat{i} = BL_Data(:,1);
    Uinfty{i} = BL_Data(:,2);
    utau{i} = BL_Data(:,3);
    nu{i} = BL_Data(:,4);
    delta99{i} = BL_Data(:,end);
    
end

utau_upstream = [1.0114 1.0437 1.0572 1.0191];

clear this* BL_Data filelist

%%

filelist = dir(fullfile(rootdir,'Re*/*xhat*'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list

thisN = length(filelist);

i7 = 0;
i10 = 0;
i14 = 0;
i21 = 0;

for i = 1:thisN
   
    thisname = filelist(i).name;
    thisfolder = filelist(i).folder;
    thisfile = strcat(thisfolder,'/',thisname);
    
    thisData = load(thisfile);
    
    if contains(thisname,'Re07')
        
        i7 = i7 + 1;
        
        Re07_zplus{i7} = flipud(thisData(:,1));
        Re07_zdel99{i7} = flipud(thisData(:,2));
        Re07_Uplus{i7} = flipud(thisData(:,3));
        Re07_UvelDef{i7} = flipud(thisData(:,4));
        Re07_uuplus{i7} = flipud(thisData(:,end));

        Re07_N(i7) = length(Re07_zplus{i7});
        
    elseif contains(thisname,'Re10')
        
        i10 = i10 + 1;
        
        Re10_zplus{i10} = flipud(thisData(:,1));
        Re10_zdel99{i10} = flipud(thisData(:,2));
        Re10_Uplus{i10} = flipud(thisData(:,3));
        Re10_UvelDef{i10} = flipud(thisData(:,4));
        Re10_uuplus{i10} = flipud(thisData(:,end));

        Re10_N(i10) = length(Re10_zplus{i10});
    
    elseif contains(thisname,'Re14')
        
        i14 = i14 + 1;
        
        Re14_zplus{i14} = flipud(thisData(:,1));
        Re14_zdel99{i14} = flipud(thisData(:,2));
        Re14_Uplus{i14} = flipud(thisData(:,3));
        Re14_UvelDef{i14} = flipud(thisData(:,4));
        Re14_uuplus{i14} = flipud(thisData(:,end));

        Re14_N(i14) = length(Re14_zplus{i14});
    
    elseif contains(thisname,'Re21')
        
        i21 = i21 + 1;
        
        Re21_zplus{i21} = flipud(thisData(:,1));
        Re21_zdel99{i21} = flipud(thisData(:,2));
        Re21_Uplus{i21} = flipud(thisData(:,3));
        Re21_UvelDef{i21} = flipud(thisData(:,4));
        Re21_uuplus{i21} = flipud(thisData(:,end));

        Re21_N(i21) = length(Re21_zplus{i21});
    
    end
end

clear this* filelist i7 i10 i14 i21


%% Calculate IBL Height Based on Power Law Fit

delta0 = [0.11 0.15 0.22 0.32];

for i = 1:N_bl
   
    xhatdelta0 = xhat{i}./delta0(i);
    
    delta_ibl{i} = ((delta0(i)*0.094).*((xhatdelta0).^0.77));
    
end


%% Velocity

for i = 1:N_bl
    
    if i == 1
        thisN = length(Re07_N);
        for j = 1:thisN
            thisZ = Re07_zdel99{j};
            
            Re07_U{j} = Re07_Uplus{j}.*utau{i}(j);
            
            ii = find(thisZ > delta_ibl{i}(j)/delta99{i}(j),1);
            
            Re07_Udelta_ibl(j) = Re07_U{j}(ii);
        end
    
    elseif i == 2 
        thisN = length(Re10_N);
        for j = 1:thisN
            thisZ = Re10_zdel99{j};
            
            Re10_U{j} = Re10_Uplus{j}.*utau{i}(j);
            
            ii = find(thisZ > delta_ibl{i}(j)/delta99{i}(j),1);
            
            Re10_Udelta_ibl(j) = Re10_U{j}(ii);
        end
        
    elseif i == 3 
        thisN = length(Re14_N);
        for j = 1:thisN
            thisZ = Re14_zdel99{j};
            
            Re14_U{j} = Re14_Uplus{j}.*utau{i}(j);
            
            ii = find(thisZ > delta_ibl{i}(j)/delta99{i}(j),1);
            
            Re14_Udelta_ibl(j) = Re14_U{j}(ii);
        end
        
    else 
        thisN = length(Re21_N);
        for j = 1:thisN
            thisZ = Re21_zdel99{j};
            
            Re21_U{j} = Re21_Uplus{j}.*utau{i}(j);
            
            ii = find(thisZ > delta_ibl{i}(j)/delta99{i}(j),1);
            
            Re21_Udelta_ibl(j) = Re21_U{j}(ii);
        end
        
    end
       
    
end


%% Mean Velocity Profiles

close all;

theseLengths = [length(Re07_N) length(Re10_N) length(Re14_N) length(Re21_N)];
theseColors = ["#effa70","#f39163","#a54b74","#322956"];

figure();
subplot(1,2,1);
for i = 1:theseLengths(1)
    p1 = semilogx(Re07_zplus{i},Re07_Uplus{i},'^','MarkerSize',8,...
    'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    p2 = semilogx(Re10_zplus{i},Re10_Uplus{i},'^','MarkerSize',8,...
    'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    p3 = semilogx(Re14_zplus{i},Re14_Uplus{i},'^','MarkerSize',8,...
    'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    p4 = semilogx(Re21_zplus{i},Re21_Uplus{i},'^','MarkerSize',8,...
    'Color',theseColors(4)); hold on
end
set(gca,'FontSize',16);
xlabel('$z^+$','FontSize',18);
ylabel('$\langle U\rangle^+$','FontSize',18);
title('Viscous Scaling','FontSize',20)
legend([p1 p2 p3 p4],...
    {'$Re_\tau = 7000$','$Re_\tau = 10000$','$Re_\tau = 14000$','$Re_\tau = 21000$'},...
    'Interpreter','Latex','Location','Northwest');

subplot(1,2,2);
for i = 1:theseLengths(1)
    thisZ = Re07_zdel99{i};
    semilogx(thisZ.*delta99{1}(i)./delta_ibl{1}(i),...
        Re07_U{i}./Re07_Udelta_ibl(i),'^','MarkerSize',8,...
        'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    thisZ = Re10_zdel99{i};
    semilogx(thisZ.*delta99{2}(i)./delta_ibl{2}(i),...
        Re10_U{i}./Re10_Udelta_ibl(i),'^','MarkerSize',8,...
        'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    thisZ = Re14_zdel99{i};
    semilogx(thisZ.*delta99{3}(i)./delta_ibl{3}(i),...
        Re14_U{i}./Re14_Udelta_ibl(i),'^','MarkerSize',8,...
        'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    thisZ = Re21_zdel99{i};
    semilogx(thisZ.*delta99{4}(i)./delta_ibl{4}(i),...
        Re21_U{i}./Re21_Udelta_ibl(i),'^','MarkerSize',8,...
        'Color',theseColors(4)); hold on
end
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$\langle U\rangle/U_\infty(z = \delta_i)$','FontSize',18);
title('IBL Scaling','FontSize',20)

%% Non-Log Scaling

figure();
subplot(1,2,1);
for i = 1:theseLengths(1)
    p1 = plot(Re07_U{i}./Uinfty{1}(i),Re07_zdel99{i},'^','MarkerSize',8,...
    'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    p2 = plot(Re10_U{i}./Uinfty{2}(i),Re10_zdel99{i},'^','MarkerSize',8,...
    'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    p3 = plot(Re14_U{i}./Uinfty{3}(i),Re14_zdel99{i},'^','MarkerSize',8,...
    'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    p4 = plot(Re21_U{i}./Uinfty{4}(i),Re21_zdel99{i},'^','MarkerSize',8,...
    'Color',theseColors(4)); hold on
end
xlim([0 1]);
set(gca,'FontSize',16);
ylabel('$z/\delta_{99}$','FontSize',18);
xlabel('$\langle U\rangle$','FontSize',18);
title('Viscous Scaling','FontSize',20)
legend([p1 p2 p3 p4],...
    {'$Re_\tau = 7000$','$Re_\tau = 10000$','$Re_\tau = 14000$','$Re_\tau = 21000$'},...
    'Interpreter','Latex','Location','Northwest');

subplot(1,2,2);
for i = 1:theseLengths(1)
    thisZ = Re07_zdel99{i};
    plot(Re07_U{i}./Re07_Udelta_ibl(i),thisZ.*delta99{1}(i)./delta_ibl{1}(i),...
        '^','MarkerSize',8,'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    thisZ = Re10_zdel99{i};
    plot(Re10_U{i}./Re10_Udelta_ibl(i),thisZ.*delta99{2}(i)./delta_ibl{2}(i),...
        '^','MarkerSize',8,'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    thisZ = Re14_zdel99{i};
    plot(Re14_U{i}./Re14_Udelta_ibl(i),thisZ.*delta99{3}(i)./delta_ibl{3}(i),...
        '^','MarkerSize',8,'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    thisZ = Re21_zdel99{i};
    plot(Re21_U{i}./Re21_Udelta_ibl(i),thisZ.*delta99{4}(i)./delta_ibl{4}(i),...
        '^','MarkerSize',8,'Color',theseColors(4)); hold on
end
xlim([0,1]);
set(gca,'FontSize',16);
ylabel('$z/\delta_i$','FontSize',18);
xlabel('$\langle U\rangle/U_\infty(z = \delta_i)$','FontSize',18);
title('IBL Scaling','FontSize',20)

%% Three Window Plot

close all;

theseLengths = [length(Re07_N) length(Re10_N) length(Re14_N) length(Re21_N)];
% theseColors = ["#effa70","#f39163","#a54b74","#322956"];
theseColors = ["#31b800","#bc6200","#ba0055","#2800c7"];


figure();
subplot(1,3,1);
for i = 1:theseLengths(1)
    p1 = semilogx(Re07_zplus{i},Re07_Uplus{i},'^','MarkerSize',8,...
    'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    p2 = semilogx(Re10_zplus{i},Re10_Uplus{i},'square','MarkerSize',8,...
    'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    p3 = semilogx(Re14_zplus{i},Re14_Uplus{i},'+','MarkerSize',8,...
    'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    p4 = semilogx(Re21_zplus{i},Re21_Uplus{i},'o','MarkerSize',8,...
    'Color',theseColors(4)); hold on
end
set(gca,'FontSize',16);
xlabel('$z^+$','FontSize',18);
ylabel('$\langle U\rangle^+$','FontSize',18);
title('Viscous Scaling','FontSize',20)
legend([p1 p2 p3 p4],...
    {'$Re_\tau = 7000$','$Re_\tau = 10000$','$Re_\tau = 14000$','$Re_\tau = 21000$'},...
    'Interpreter','Latex','Location','Northwest');

subplot(1,3,2);
for i = 1:theseLengths(1)
    plot(Re07_U{i}./Uinfty{1}(i),Re07_zdel99{i},'^','MarkerSize',8,...
    'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    plot(Re10_U{i}./Uinfty{2}(i),Re10_zdel99{i},'square','MarkerSize',8,...
    'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    plot(Re14_U{i}./Uinfty{3}(i),Re14_zdel99{i},'+','MarkerSize',8,...
    'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    plot(Re21_U{i}./Uinfty{4}(i),Re21_zdel99{i},'o','MarkerSize',8,...
    'Color',theseColors(4)); hold on
end
xlim([0 1]);
set(gca,'FontSize',16);
ylabel('$z/\delta$','FontSize',18);
xlabel('$\langle U\rangle/U_\infty$','FontSize',18);
title('Outer Scaling','FontSize',20)

% Non-Log
subplot(1,3,3);
for i = 1:theseLengths(1)
    thisZ = Re07_zdel99{i};
    plot(Re07_U{i}./Re07_Udelta_ibl(i),thisZ.*delta99{1}(i)./delta_ibl{1}(i),...
        '^','MarkerSize',8,'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    thisZ = Re10_zdel99{i};
    plot(Re10_U{i}./Re10_Udelta_ibl(i),thisZ.*delta99{2}(i)./delta_ibl{2}(i),...
        'square','MarkerSize',8,'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    thisZ = Re14_zdel99{i};
    plot(Re14_U{i}./Re14_Udelta_ibl(i),thisZ.*delta99{3}(i)./delta_ibl{3}(i),...
        '+','MarkerSize',8,'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    thisZ = Re21_zdel99{i};
    plot(Re21_U{i}./Re21_Udelta_ibl(i),thisZ.*delta99{4}(i)./delta_ibl{4}(i),...
        'o','MarkerSize',8,'Color',theseColors(4)); hold on
end
% xlim([0,1]);
set(gca,'FontSize',16);
ylabel('$z/\delta_i$','FontSize',18);
xlabel('$\langle U\rangle/U_\infty(z = \delta_i)$','FontSize',18);
title('IBL Scaling','FontSize',20)

% Semilog
% subplot(1,3,3);
% for i = 1:theseLengths(1)
%     thisZ = Re07_zdel99{i};
%     semilogx(thisZ.*delta99{1}(i)./delta_ibl{1}(i),...
%         Re07_U{i}./Re07_Udelta_ibl(i),'^','MarkerSize',8,...
%         'Color',theseColors(1)); hold on
% end
% for i = 1:theseLengths(2)
%     thisZ = Re10_zdel99{i};
%     semilogx(thisZ.*delta99{2}(i)./delta_ibl{2}(i),...
%         Re10_U{i}./Re10_Udelta_ibl(i),'^','MarkerSize',8,...
%         'Color',theseColors(2)); hold on
% end
% for i = 1:theseLengths(3)
%     thisZ = Re14_zdel99{i};
%     semilogx(thisZ.*delta99{3}(i)./delta_ibl{3}(i),...
%         Re14_U{i}./Re14_Udelta_ibl(i),'^','MarkerSize',8,...
%         'Color',theseColors(3)); hold on
% end
% for i = 1:theseLengths(4)
%     thisZ = Re21_zdel99{i};
%     semilogx(thisZ.*delta99{4}(i)./delta_ibl{4}(i),...
%         Re21_U{i}./Re21_Udelta_ibl(i),'^','MarkerSize',8,...
%         'Color',theseColors(4)); hold on
% end
% set(gca,'FontSize',16);
% xlabel('$z/\delta_i$','FontSize',18);
% ylabel('$\langle U\rangle/U_i$','FontSize',18);
% title('IBL Scaling','FontSize',20)

%% Calculate new defect

for i = 1:N_bl
    
    thisN = theseLengths(i);
    thisUtau = utau{i};
    
    if i == 1
        for j = 1:thisN
            thisVel = Re07_U{j};
            Re07_velDefect{j} = (Re07_Udelta_ibl(j) - thisVel)./thisUtau(j);
        end
    elseif i == 2
        for j = 1:thisN
            thisVel = Re10_U{j};
            Re10_velDefect{j} = (Re10_Udelta_ibl(j) - thisVel)./thisUtau(j);
        end
    elseif i == 3
        for j = 1:thisN
            thisVel = Re14_U{j};
            Re14_velDefect{j} = (Re14_Udelta_ibl(j) - thisVel)./thisUtau(j);
        end
    else
        for j = 1:thisN
            thisVel = Re21_U{j};
            Re21_velDefect{j} = (Re21_Udelta_ibl(j) - thisVel)./thisUtau(j);
        end
    end
end

% Normalizing old defect with upstream utau
for i = 1:N_bl
    
    thisN = theseLengths(i);
    thisUtau = utau{i};
    thisUtauUP = utau_upstream(i);
    
    if i == 1
        for j = 1:thisN
            thisVel = Re07_UvelDef{j};
            Re07_vdUP{j} = (thisVel.*thisUtau(j))./thisUtauUP(i);
        end
    elseif i == 2
        for j = 1:thisN
            thisVel = Re10_UvelDef{j};
            Re10_vdUP{j} = (thisVel.*utau{i}(j))./utau_upstream(i);
        end
    elseif i == 3
        for j = 1:thisN
            thisVel = Re14_UvelDef{j};
            Re14_vdUP{j} = (thisVel.*utau{i}(j))./utau_upstream(i);
        end
    else
        for j = 1:thisN
            thisVel = Re21_UvelDef{j};
            Re21_vdUP{j} = (thisVel.*utau{i}(j))./utau_upstream(i);
        end
    end
end

%% Plot Velocity Defects

close all;

figure();
subplot(1,2,1);
for i = 1:theseLengths(1)
    p1 = semilogx(Re07_zdel99{i},Re07_UvelDef{i},'^','MarkerSize',8,...
    'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    p2 = semilogx(Re10_zdel99{i},Re10_UvelDef{i},'^','MarkerSize',8,...
    'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    p3 = semilogx(Re14_zdel99{i},Re14_UvelDef{i},'^','MarkerSize',8,...
    'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    p4 = semilogx(Re21_zdel99{i},Re21_UvelDef{i},'^','MarkerSize',8,...
    'Color',theseColors(4)); hold on
end
set(gca,'FontSize',16);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);
title('Classic Scaling with Local $u_\tau$');
legend([p1 p2 p3 p4],...
    {'$Re_\tau = 7k$','$Re_\tau = 10k$','$Re_\tau = 14k$','$Re_\tau = 21k$'},...
    'Interpreter','Latex','Location','Northeast');

subplot(1,2,2);
for i = 1:theseLengths(1)
    thisZ = Re07_zdel99{i};
    semilogx(thisZ.*delta99{1}(i)./delta_ibl{1}(i),...
        Re07_velDefect{i},'^','MarkerSize',8,...
        'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    thisZ = Re10_zdel99{i};
    semilogx(thisZ.*delta99{2}(i)./delta_ibl{2}(i),...
        Re10_velDefect{i},'^','MarkerSize',8,...
        'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    thisZ = Re14_zdel99{i};
    semilogx(thisZ.*delta99{3}(i)./delta_ibl{3}(i),...
        Re14_velDefect{i},'^','MarkerSize',8,...
        'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    thisZ = Re21_zdel99{i};
    semilogx(thisZ.*delta99{4}(i)./delta_ibl{4}(i),...
        Re21_velDefect{i},'^','MarkerSize',8,...
        'Color',theseColors(4)); hold on
end
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_\infty(z=\delta_i) - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$');

%% Three Panel Defect

close all;

figure();
subplot(1,3,1);
for i = 1:theseLengths(1)
    p1 = semilogx(Re07_zdel99{i},Re07_UvelDef{i},'^','MarkerSize',8,...
    'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    p2 = semilogx(Re10_zdel99{i},Re10_UvelDef{i},'square','MarkerSize',8,...
    'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    p3 = semilogx(Re14_zdel99{i},Re14_UvelDef{i},'+','MarkerSize',8,...
    'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    p4 = semilogx(Re21_zdel99{i},Re21_UvelDef{i},'o','MarkerSize',8,...
    'Color',theseColors(4)); hold on
end
set(gca,'FontSize',16);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);
title('Classic Scaling with Local $u_\tau$');
legend([p1 p2 p3 p4],...
    {'$Re_\tau = 7k$','$Re_\tau = 10k$','$Re_\tau = 14k$','$Re_\tau = 21k$'},...
    'Interpreter','Latex','Location','Northeast');

subplot(1,3,2);
for i = 1:theseLengths(1)
    p1 = semilogx(Re07_zdel99{i},Re07_vdUP{i},'^','MarkerSize',8,...
    'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    p2 = semilogx(Re10_zdel99{i},Re10_vdUP{i},'square','MarkerSize',8,...
    'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    p3 = semilogx(Re14_zdel99{i},Re14_vdUP{i},'+','MarkerSize',8,...
    'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    p4 = semilogx(Re21_zdel99{i},Re21_vdUP{i},'o','MarkerSize',8,...
    'Color',theseColors(4)); hold on
end
set(gca,'FontSize',16);
ylim([0 40]);
ylabel('$(U_\infty - U)/u_{\tau,0}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);
title('Classic Scaling with Upstream $u_\tau$');

subplot(1,3,3);
for i = 1:theseLengths(1)
    thisZ = Re07_zdel99{i};
    semilogx(thisZ.*delta99{1}(i)./delta_ibl{1}(i),...
        Re07_velDefect{i},'^','MarkerSize',8,...
        'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    thisZ = Re10_zdel99{i};
    semilogx(thisZ.*delta99{2}(i)./delta_ibl{2}(i),...
        Re10_velDefect{i},'square','MarkerSize',8,...
        'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    thisZ = Re14_zdel99{i};
    semilogx(thisZ.*delta99{3}(i)./delta_ibl{3}(i),...
        Re14_velDefect{i},'+','MarkerSize',8,...
        'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    thisZ = Re21_zdel99{i};
    semilogx(thisZ.*delta99{4}(i)./delta_ibl{4}(i),...
        Re21_velDefect{i},'o','MarkerSize',8,...
        'Color',theseColors(4)); hold on
end
set(gca,'FontSize',16);
xlim([0 1]);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$');

%% Non-Log Scaling

figure();
subplot(1,2,1);
for i = 1:theseLengths(1)
    p1 = plot(Re07_zdel99{i},Re07_UvelDef{i},'^','MarkerSize',8,...
    'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    p2 = plot(Re10_zdel99{i},Re10_UvelDef{i},'^','MarkerSize',8,...
    'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    p3 = plot(Re14_zdel99{i},Re14_UvelDef{i},'^','MarkerSize',8,...
    'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    p4 = plot(Re21_zdel99{i},Re21_UvelDef{i},'^','MarkerSize',8,...
    'Color',theseColors(4)); hold on
end
set(gca,'FontSize',16);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);
title('Classic Scaling with Local $u_\tau$');
legend([p1 p2 p3 p4],...
    {'$Re_\tau = 7k$','$Re_\tau = 10k$','$Re_\tau = 14k$','$Re_\tau = 21k$'},...
    'Interpreter','Latex','Location','Northeast');

subplot(1,2,2);
for i = 1:theseLengths(1)
    thisZ = Re07_zdel99{i};
    plot(thisZ.*delta99{1}(i)./delta_ibl{1}(i),...
        Re07_velDefect{i},'^','MarkerSize',8,...
        'Color',theseColors(1)); hold on
end
for i = 1:theseLengths(2)
    thisZ = Re10_zdel99{i};
    plot(thisZ.*delta99{2}(i)./delta_ibl{2}(i),...
        Re10_velDefect{i},'^','MarkerSize',8,...
        'Color',theseColors(2)); hold on
end
for i = 1:theseLengths(3)
    thisZ = Re14_zdel99{i};
    plot(thisZ.*delta99{3}(i)./delta_ibl{3}(i),...
        Re14_velDefect{i},'^','MarkerSize',8,...
        'Color',theseColors(3)); hold on
end
for i = 1:theseLengths(4)
    thisZ = Re21_zdel99{i};
    plot(thisZ.*delta99{4}(i)./delta_ibl{4}(i),...
        Re21_velDefect{i},'^','MarkerSize',8,...
        'Color',theseColors(4)); hold on
end
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_\infty(z=\delta_i) - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$');

%% Imposing One defect with Wake Function

color12 = [ "#292f56", "#4a3664", "#6c3c6c", "#8d426d", "#a94b67", ...
    "#bf5b5c", "#cc7050", "#d08a43", "#caa73d", "#b9c346", "#9ddf62", ...
    "#70fa8e"];

smoothPI = 0.55;
weakPI = 0.25;
kappa = 0.4;

% thisYD = Re07_zdel99{end};
% 
% myIND = theseLengths(1);
% thisYDi = Re07_zdel99{myIND}.*delta99{1}(myIND)./delta_ibl{1}(myIND);

thisYD = Re14_zdel99{end};

myIND = theseLengths(3);
thisYDi = Re14_zdel99{myIND}.*delta99{1}(myIND)./delta_ibl{1}(myIND);

ii = find(thisYD >= 1.2,1);
thisYDi = thisYDi(1:ii);

wake_ibl = defectLaw(kappa,smoothPI,thisYDi);
wake = defectLaw(kappa,smoothPI,thisYD);
    

close all;

figure();
subplot(1,2,1);
for i = 1:myIND
%     p1 = semilogx(Re07_zdel99{i},Re07_UvelDef{i},'square','MarkerSize',8,...
%         'Color',color12(i)); hold on
    p1 = semilogx(Re14_zdel99{i},Re14_UvelDef{i},'square','MarkerSize',8,...
        'Color',color12(i)); hold on
end
semilogx(thisYD,wake,'k--','LineWidth',3);
set(gca,'FontSize',16);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);
title('Classic Scaling with Local $u_\tau$');
legend('$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$','$\hat{x}_8$',...
    '$\hat{x}_9$','$\hat{x}_{10}$','$\hat{x}_{11}$','$\hat{x}_{12}$',...
    'Interpreter','Latex',...
    'Location','NorthEast','NumColumns',3,'FontSize',10);

subplot(1,2,2);
for i = 1:myIND
%     thisZ = Re07_zdel99{i};
%     semilogx(thisZ.*delta99{1}(i)./delta_ibl{1}(i),...
%         Re07_velDefect{i},'square','MarkerSize',8,...
%         'Color',color12(i)); hold on
    thisZ = Re14_zdel99{i};
    semilogx(thisZ.*delta99{1}(i)./delta_ibl{1}(i),...
        Re14_velDefect{i},'square','MarkerSize',8,...
        'Color',color12(i)); hold on
end
semilogx(thisYDi,wake_ibl,'k--','LineWidth',3);
set(gca,'FontSize',16);
xlim([0 1]);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$');

%% Re 7k only mean velocity

figure();
subplot(1,3,1);
% for i = 1:theseLengths(1)
%     p1 = semilogx(Re07_zplus{i},Re07_Uplus{i},'square','MarkerSize',8,...
%     'Color',color12(i)); hold on
% end
% set(gca,'FontSize',16);
% xlabel('$z^+$','FontSize',18);
% ylabel('$U^+$','FontSize',18);

for i = 1:theseLengths(1)
    p1 = semilogy(Re07_Uplus{i},Re07_zplus{i},'square','MarkerSize',8,...
    'Color',color12(i)); hold on
end
set(gca,'FontSize',16);
ylabel('$z^+$','FontSize',18);
xlabel('$U^+$','FontSize',18);
legend('$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$','$\hat{x}_8$',...
    '$\hat{x}_9$','$\hat{x}_{10}$','$\hat{x}_{11}$','$\hat{x}_{12}$',...
    'Interpreter','Latex',...
    'Location','NorthWest','NumColumns',3);
% title('Viscous Scaling','FontSize',20)

subplot(1,3,2);
for i = 1:theseLengths(1)
    plot(Re07_U{i}./Uinfty{1}(i),Re07_zdel99{i},'square','MarkerSize',8,...
    'Color',color12(i)); hold on
end
xlim([0 1]);
set(gca,'FontSize',16);
ylabel('$z/\delta$','FontSize',18);
xlabel('$U/U_\infty$','FontSize',18);
% title('Outer Scaling','FontSize',20)

subplot(1,3,3);
for i = 1:theseLengths(1)
    thisZ = Re07_zdel99{i};
    plot(Re07_U{i}./Re07_Udelta_ibl(i),thisZ.*delta99{1}(i)./delta_ibl{1}(i),...
        'square','MarkerSize',8,'Color',color12(i)); hold on
    
end
set(gca,'FontSize',16);
% ylabel('$z/\delta_i$','FontSize',18);
% xlabel('$U/U_i$','FontSize',18);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$U/U_i$','FontSize',18);
xlim([0 1]);
ylim([0 1]);
% title('IBL Scaling','FontSize',20)

%% New Velocity Profile Figures

thisChoice = 4;

figure();
subplot(2,1,1);
for i = 1:theseLengths(thisChoice)
    p1 = semilogx(Re21_zplus{i},Re21_Uplus{i},'square','MarkerSize',8,...
    'Color',color12(i)); hold on
end
set(gca,'FontSize',16);
xlabel('$z^+$','FontSize',18);
ylabel('$U^+$','FontSize',18);

% for i = 1:theseLengths(thisChoice)
%     p1 = semilogy(Re21_Uplus{i},Re21_zplus{i},'square','MarkerSize',8,...
%     'Color',color12(i)); hold on
% end
% set(gca,'FontSize',16);
% ylabel('$z^+$','FontSize',18);
% xlabel('$U^+$','FontSize',18);
legend('$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$','$\hat{x}_8$',...
    '$\hat{x}_9$','$\hat{x}_{10}$',...'$\hat{x}_{11}$','$\hat{x}_{12}$',...
    'Interpreter','Latex',...
    'Location','NorthWest','NumColumns',3);
% title('Viscous Scaling','FontSize',20)

subplot(2,1,2);
for i = 1:theseLengths(thisChoice)
    plot(Re21_U{i}./Uinfty{thisChoice}(i),Re21_zdel99{i},'square','MarkerSize',8,...
    'Color',color12(i)); hold on
end
xlim([0 1]);
set(gca,'FontSize',16);
ylabel('$z/\delta$','FontSize',18);
xlabel('$U/U_\infty$','FontSize',18);
% title('Outer Scaling','FontSize',20)

figure();
for i = 1:theseLengths(thisChoice)
    thisZ = Re21_zdel99{i};
    plot(Re21_U{i}./Re21_Udelta_ibl(i),thisZ.*delta99{thisChoice}(i)./delta_ibl{thisChoice}(i),...
        'square','MarkerSize',8,'Color',color12(i)); hold on
    
end
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$U/U_i$','FontSize',18);
% xlim([0 1]);
% ylim([0 1]);
% title('IBL Scaling','FontSize',20)

figure();
for i = 1:theseLengths(thisChoice)
    thisZ = Re21_zdel99{i};
    plot(Re21_U{i}./Re21_Udelta_ibl(i),thisZ.*delta99{thisChoice}(i)./delta_ibl{thisChoice}(i),...
        'square','MarkerSize',8,'Color',color12(i)); hold on
    
end
set(gca,'FontSize',16);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$U/U_i$','FontSize',18);
xlim([0 1]);
ylim([0 1]);

%% ZS Scaling

% Usual Values
for i = 1:theseLengths(4)
    this_i = find(Re21_zdel99{i} >= 1.0, 1);
    this_z = Re21_zdel99{i}.*delta99{4}(i);
    this_z = this_z(1:this_i);
    
    this_integrand = 1 - (Re21_U{i}./Uinfty{4}(i));
    this_integrand = this_integrand(1:this_i);
    
    delta_star(i) = trapz(this_z,this_integrand);
    
    u0(i) = Uinfty{4}(i)*(delta_star(i)/delta99{4}(i));
    
end

% Now with IBL Parameters
for i = 1:theseLengths(4)
    thisZ = Re21_zdel99{i};
    this_zibl = thisZ.*delta99{4}(i)./delta_ibl{4}(i);
    
    this_i = find(this_zibl >= 1.0, 1);
    this_z = Re21_zdel99{i}.*delta99{4}(i);
    this_z = this_z(1:this_i);
    
    this_integrand = 1 - (Re21_U{i}./Re21_Udelta_ibl(i));
    this_integrand = this_integrand(1:this_i);
    
    delta_i_star(i) = trapz(this_z,this_integrand);
    
    u0_i(i) = Re21_Udelta_ibl(i)*(delta_i_star(i)/delta_ibl{4}(i));
    
end    


%% Plot ZS Scaling

myIND = theseLengths(4);

close all;

figure();
subplot(2,2,1);
for i = 1:myIND
    p1 = semilogx(Re21_zdel99{i},Re21_UvelDef{i},'square','MarkerSize',8,...
        'Color',color12(i)); hold on
%     p1 = semilogx(Re14_zdel99{i},Re14_UvelDef{i},'square','MarkerSize',8,...
%         'Color',color12(i)); hold on
end
% semilogx(thisYD,wake,'k--','LineWidth',3);
set(gca,'FontSize',16);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);
title('Classic Scaling with Local $u_\tau$');
legend('$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$','$\hat{x}_8$',...
    '$\hat{x}_9$','$\hat{x}_{10}$','$\hat{x}_{11}$','$\hat{x}_{12}$',...
    'Interpreter','Latex',...
    'Location','NorthEast','NumColumns',3,'FontSize',10);

subplot(2,2,2);
for i = 1:myIND
    thisZ = Re21_zdel99{i};
    semilogx(thisZ.*delta99{4}(i)./delta_ibl{4}(i),...
        Re21_velDefect{i},'square','MarkerSize',8,...
        'Color',color12(i)); hold on
%     thisZ = Re14_zdel99{i};
%     semilogx(thisZ.*delta99{1}(i)./delta_ibl{1}(i),...
%         Re14_velDefect{i},'square','MarkerSize',8,...
%         'Color',color12(i)); hold on
end
% semilogx(thisYDi,wake_ibl,'k--','LineWidth',3);
set(gca,'FontSize',16);
xlim([0 1]);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$');

subplot(2,2,3);
for i = 1:myIND
    thisZS = (Uinfty{4}(i) - Re21_U{i})./u0(i);
    
    semilogx(Re21_zdel99{i},thisZS,'square','MarkerSize',8,...
        'Color',color12(i)); hold on
end
set(gca,'FontSize',16);
xlabel('$z/\delta$','FontSize',18);
ylabel('$(U_\infty - U)/U_\infty\delta^*/\delta$','FontSize',18);
title('ZS Scaling');


subplot(2,2,4);
for i = 1:myIND
    thisZS = (Re21_Udelta_ibl(i) - Re21_U{i})./u0_i(i);
    thisZ = Re21_zdel99{i};
    this_zibl = thisZ.*delta99{4}(i)./delta_ibl{4}(i);
    
    semilogx(this_zibl,thisZS,'square','MarkerSize',8,...
        'Color',color12(i)); hold on
end
set(gca,'FontSize',16);
xlim([0 1]);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/U_i\delta_i^*/\delta_i$','FontSize',18);
title('ZS Scaling');

% Non Semilog



figure();
subplot(2,2,1);
for i = 1:myIND
    p1 = plot(Re21_zdel99{i},Re21_UvelDef{i},'square','MarkerSize',8,...
        'Color',color12(i)); hold on
%     p1 = semilogx(Re14_zdel99{i},Re14_UvelDef{i},'square','MarkerSize',8,...
%         'Color',color12(i)); hold on
end
xlim([0 1]);
% semilogx(thisYD,wake,'k--','LineWidth',3);
set(gca,'FontSize',16);
ylabel('$(U_\infty - U)/u_{\tau,2}$','FontSize',18)
xlabel('$z/\delta$','FontSize',18);
title('Classic Scaling with Local $u_\tau$');
legend('$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$','$\hat{x}_8$',...
    '$\hat{x}_9$','$\hat{x}_{10}$','$\hat{x}_{11}$','$\hat{x}_{12}$',...
    'Interpreter','Latex',...
    'Location','NorthEast','NumColumns',3,'FontSize',10);

subplot(2,2,2);
for i = 1:myIND
    thisZ = Re21_zdel99{i};
    plot(thisZ.*delta99{4}(i)./delta_ibl{4}(i),...
        Re21_velDefect{i},'square','MarkerSize',8,...
        'Color',color12(i)); hold on
%     thisZ = Re14_zdel99{i};
%     semilogx(thisZ.*delta99{1}(i)./delta_ibl{1}(i),...
%         Re14_velDefect{i},'square','MarkerSize',8,...
%         'Color',color12(i)); hold on
end
% semilogx(thisYDi,wake_ibl,'k--','LineWidth',3);
set(gca,'FontSize',16);
xlim([0 1]);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/u_{\tau,2}$','FontSize',18);
title('IBL Scaling with Local $u_\tau$');

subplot(2,2,3);
for i = 1:myIND
    thisZS = (Uinfty{4}(i) - Re21_U{i})./u0(i);
    
    plot(Re21_zdel99{i},thisZS,'square','MarkerSize',8,...
        'Color',color12(i)); hold on
end
xlim([0 1]);
set(gca,'FontSize',16);
xlabel('$z/\delta$','FontSize',18);
ylabel('$(U_\infty - U)/U_\infty\delta^*/\delta$','FontSize',18);
title('ZS Scaling');


subplot(2,2,4);
for i = 1:myIND
    thisZS = (Re21_Udelta_ibl(i) - Re21_U{i})./u0_i(i);
    thisZ = Re21_zdel99{i};
    plot(thisZ.*delta99{4}(i)./delta_ibl{4}(i),thisZS,...
        'square','MarkerSize',8,...
        'Color',color12(i)); hold on
end
set(gca,'FontSize',16);
xlim([0 1]);
xlabel('$z/\delta_i$','FontSize',18);
ylabel('$(U_i - U)/U_i\delta_i^*/\delta_i$','FontSize',18);
title('ZS Scaling');


%% Functions

function wake = defectLaw(kappa,myPI,yd)
    wake =  (1/kappa)*( -1*log(yd) + myPI*(2 - 2*(sin((pi()/2)*yd).^2)));
end

%% End