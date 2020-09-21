%%% 1D test best on rotating angles
close all;
clear
clc

angles = (0:0.5:90)/180*pi;

%%% Butchered this section to get correct RI values,
    %%% Rethink in other
ePerp = 1.53^2; % Short axis
ePara = 1.47^2; (1.53-0.0024)^2; % Long axis   
%Swapped - I think...
eDiff = ePerp - ePara;

reflectedComponent = zeros(length(angles),2);

RIs = zeros(length(angles),4);

for i = 1:length(angles)
    % Starts in parallel fibers, goes to perendicular
    n1 = (ePara + eDiff*cos(angles(i)).^2.*cos(pi/2).^2).^0.5;
    n2 = (ePara + eDiff*cos(pi/2-angles(i)).^2.*cos(0).^2).^0.5;
    
    RIs(i,1:2) = [n1, n2];
    
    if n1 > n2
        critAngle = asin(n2/n1);
    else
        critAngle = pi/2;
    end
    
    transAngle = asin(n1/n2*sin(angles(i)));
    
    cosT = sqrt(1-(n1/n2*sin(angles(i)))^2);
    
    if angles(i) <= critAngle
        Rperp = ((n1*cos(angles(i)) - n2*cosT)/(n1*cos(angles(i)) + n2*cosT))^2;

        % Possible invert first part
        Rpara = ((n2*cos(angles(i)) - n1*cosT)/(n2*cos(angles(i)) + n1*cosT))^2;
        
        reflectedComponent(i,1) = (Rperp + Rpara)/2;
    else
        reflectedComponent(i,1) = 1;
    end
    
    % Opposite, start in perendicular, goes to parallel
    n1 = (ePara + eDiff*cos(pi/2-angles(i)).^2.*cos(0).^2).^0.5;
    n2 = (ePara + eDiff*cos(angles(i)).^2.*cos(pi/2).^2).^0.5;
    
    RIs(i,3:4) = [n1, n2];
    
    if n1 > n2
        critAngle = asin(n2/n1);
    else
        critAngle = pi/2;
    end
    
    transAngle = asin(n1/n2*sin(angles(i)));
    
    cosT = sqrt(1-(n1/n2*sin(angles(i)))^2);
    
    if angles(i) <= critAngle
        Rperp = ((n1*cos(angles(i)) - n2*cosT)/(n1*cos(angles(i)) + n2*cosT))^2;

        % Possible invert first part
        Rpara = ((n2*cos(angles(i)) - n1*cosT)/(n2*cos(angles(i)) + n1*cosT))^2;
        
        reflectedComponent(i,2) = (Rperp + Rpara)/2;
    else
        reflectedComponent(i,2) = 1;
    end
end

figure; subplot(1,2,1); hold on;
plot(angles/pi*180, reflectedComponent(:,1), 'r')
plot(angles/pi*180, 1-reflectedComponent(:,1), 'b')
xlim([0 90])
ylabel('Power'); xlabel('Short-axis angle to vertical')
title('Ray starts in long axis fibers (low RI)')


subplot(1,2,2); hold on;
plot(angles/pi*180, reflectedComponent(:,2), 'r')
plot(angles/pi*180, 1-reflectedComponent(:,2), 'b')
xlim([0 90])
ylabel('Power'); xlabel('Short-axis angle to vertical')
title('Ray starts in short axis fibers (high RI)')
legend('Reflected', 'Transmitted')

figure; subplot(1,2,1); hold on;
plot(angles/pi*180, RIs(:,1), 'g')
plot(angles/pi*180, RIs(:,2), 'c')

subplot(1,2,2); hold on;
plot(angles/pi*180, RIs(:,3), 'g')
plot(angles/pi*180, RIs(:,4), 'c')
legend('N1', 'N2')
plot(angles/pi*180, ePara+sin(angles*2-pi/2)*eDiff/2 + eDiff/2, 'r')

%% Test out RI profile at border

lamellaeThickness = 1;

layerThickness = 0.005;

layerTurn = 180/(lamellaeThickness/layerThickness);

distances = 0:layerThickness:2*lamellaeThickness;

angles = (0:layerTurn:360)/180*pi;

figure; hold on

nValues = (ePara+sin(angles*2-pi/2)*eDiff/2 + eDiff/2).^0.5;

plot(distances, nValues, 'g')

nBottomInds = find(nValues < (ePara + eDiff*0.1)^0.5);

nTopInds = find(nValues > (ePerp - eDiff*0.1)^0.5);

plot(distances(nBottomInds), nValues(nBottomInds), 'bx')

plot(distances(nTopInds), nValues(nTopInds), 'rx')

%% Sanity check
% Flat - theta 0
the = 0;
[ePerp + eDiff*cos(the).^2.*cos(0).^2, ...
ePerp + eDiff*cos(the).^2.*sin(0).^2, ...
ePerp + eDiff*sin(the).^2].^0.5

[ePerp + eDiff*cos(the).^2.*cos(pi/2).^2, ...
ePerp + eDiff*cos(the).^2.*sin(pi/2).^2, ...
ePerp + eDiff*sin(the).^2].^0.5

[ePerp + eDiff*cos(the).^2.*cos(pi).^2, ...
ePerp + eDiff*cos(the).^2.*sin(pi).^2, ...
ePerp + eDiff*sin(the).^2].^0.5

% Extra sanity checking
% % 45 degree incline
% the = pi/4;
% [ePerp + eDiff*cos(the).^2.*cos(0).^2, ...
% ePerp + eDiff*cos(the).^2.*sin(0).^2, ...
% ePerp + eDiff*sin(the).^2].^0.5
% 
% [ePerp + eDiff*cos(the).^2.*cos(pi/2).^2, ...
% ePerp + eDiff*cos(the).^2.*sin(pi/2).^2, ...
% ePerp + eDiff*sin(the).^2].^0.5
% 
% [ePerp + eDiff*cos(the).^2.*cos(pi).^2, ...
% ePerp + eDiff*cos(the).^2.*sin(pi).^2, ...
% ePerp + eDiff*sin(the).^2].^0.5
% 
% % Vertical - theta pi/2
% the = pi/2;
% [ePerp + eDiff*cos(the).^2.*cos(0).^2, ...
% ePerp + eDiff*cos(the).^2.*sin(0).^2, ...
% ePerp + eDiff*sin(the).^2].^0.5
% 
% [ePerp + eDiff*cos(the).^2.*cos(pi/2).^2, ...
% ePerp + eDiff*cos(the).^2.*sin(pi/2).^2, ...
% ePerp + eDiff*sin(the).^2].^0.5
% 
% [ePerp + eDiff*cos(the).^2.*cos(pi).^2, ...
% ePerp + eDiff*cos(the).^2.*sin(pi).^2, ...
% ePerp + eDiff*sin(the).^2].^0.5
% 
% % Backflip - theta pi
% the = pi;
% [ePerp + eDiff*cos(the).^2.*cos(0).^2, ...
% ePerp + eDiff*cos(the).^2.*sin(0).^2, ...
% ePerp + eDiff*sin(the).^2].^0.5
% 
% [ePerp + eDiff*cos(the).^2.*cos(pi/2).^2, ...
% ePerp + eDiff*cos(the).^2.*sin(pi/2).^2, ...
% ePerp + eDiff*sin(the).^2].^0.5
% 
% [ePerp + eDiff*cos(the).^2.*cos(pi).^2, ...
% ePerp + eDiff*cos(the).^2.*sin(pi).^2, ...
% ePerp + eDiff*sin(the).^2].^0.5