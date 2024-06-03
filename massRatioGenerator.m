%EOSV9 - Updated: 7/13/21
clear
clc
massList = zeros(95,1);
warning('off')

    %% Choose material 

    rho0F = 8430; %% kg/m^3 
    rho0M = 4176; 

w = waitbar(0, 'Calculating...');

for i=1:199 %%radius CMB, out of 200 
    %% initial conditions 
    prop = i*4;
    rE = 1; %earth radii
    r = (6356.752*10^3)*rE;
    V = (4/3)*pi*r^3; %volume (m^3)
    m = V*rho0M; %mass kg assuming constant density
    rList = linspace(0.01*6356.752*10^3, r, 800); %List of radii, rList(1) is bottom edge of first layer startin at 0.01 earth radius 
    rhoList(1:prop) = rho0F; %assume constant density, kg/m^3 
    rhoList(prop:800) = rho0M;

    PLayer1 = gravPresCalc(rList, rhoList); %pressure for first iteration 
    rhoList1 = rhoList; %constant density first iteration 
    PLayer2 = zeros(1,length(PLayer1));
    rhoList2 = zeros(1,length(rhoList1)); 
  
    %% iterate until pressure and density match 
    while(any(abs(PLayer1 - PLayer2) > 0.01) || any(abs(rhoList1 - rhoList2) > 10))
        if(rhoList2(2) > 0) %does not effect first cycle
            PLayer1 = PLayer2; %update old pressure and density to match new 
            rhoList1 = rhoList2;
        end
        prop = i*4;
            rhoList2((prop+1):800) = calcDensity(rhoList1((prop+1):800), PLayer1((prop+1:800)), 4176, 265.5, 4.16,1); %calculate density of mantle
            PLayer2 = gravPresCalc(rList, rhoList2);
            rhoList2(1:prop) = calcDensity(rhoList1(1:prop), PLayer2(1:prop), 7874, 162.5, 5.5,1);
            [PLayer2, m] = gravPresCalc(rList, rhoList2);        
%         rhoList2 = calcDensity(rhoList1, PLayer1); %new density under updated pressure
%        (abs(mean(rhoList1(1:length(rhoList1)),'all') - mean(rhoList2(1:length(rhoList2)),'all')))
%        
%         [PLayer2,m] = gravPresCalc(rList,rhoList2);%new pressure under updated density 
%        (abs(mean(PLayer1(2:length(rhoList1)),'all') - mean(PLayer2(2:length(rhoList2)),'all')))
        %a = "-------------------------------"
    end 
    massList(i,1) = sum(m);%(mean(rhoList2(1,1:length(rhoList2)),'all'))*(4/3)*pi*r^3;
    coreMass(i) = sum(m(1:prop));
    mantleMass(i) = sum(m(prop:800));
    massRatio(i) = mantleMass(i)/coreMass(i)
    CMF = coreMass(i)/(mantleMass(i)+coreMass(i))
    
    PList(i,1) = PLayer2(1);
    waitbar(i/195)
end
close(w)
% %% Create Plot 
 eMass = 5.9724e24;
 earthMassList = massList/eMass;
 
% earthRadiusList = 0.02:0.02:2;
% plot(earthMassList, earthRadiusList);



%% Calculate Pressure
function [PLayer, mlayer] = gravPresCalc(rList, rhoList)
        G = 6.67408*10^-11; %gravitational constant 
    %gravity and mass calculations 
    mlayer = zeros(1,length(rList));
     %on the surface the layer mass is zero

    for i = 1:length(rList)-1 %calculate layer masses 
        mlayer(i) = rhoList(i)*(rList(i+1)^3 - rList(i)^3)*(4/3)*pi; %layer 1 goes from 0.01 to 0.03 earth radii 
    end
    mlayer(end) = 0;
    
    minside = zeros(1,length(rList));
    grav = zeros(1,length(rList)-1);
    dp = zeros(1,length(rList)-1);
    for i = 1:length(rList)-1 %calculate mass inside 
        minside(i) = sum(mlayer(1:i));
        grav(i) = G*minside(i)/((rList(i+1))^2);
        dp(i) = rhoList(i)*grav(i)*(rList(end)/length(rList));
    end 
    %pressure calculation 

    for i = length(dp):-1:1
        PLayer(i) = sum(dp(i:length(dp)))*10^-9;
    end
    PLayer(length(rList)) = 0; %no pressure on the surface
end

%% Solver 
function newRhoList = calcDensity(rhoList, PLayer, rho0, k0, k0p, EOS)
xList = zeros(1,length(rhoList)); 

%solve each layer 
guess = 0.9;
options = optimset('Display','off'); %, 'Tolx', 1e-8);
 for i=1:length(rhoList)
    if(EOS == 1)
    Peqn = @(x) 3*k0 * (x^(-2)) * (1-x)*exp((1.5*k0p-1.5)*(1-x)) - PLayer(i); %Vinet
    %else  
    %else other EOS
    end
    xList(i) = fsolve(Peqn, guess,options);
    guess = xList(i);

 end 
    newRhoList = rho0./(xList.^3); %rho = rho0*n
end