%EOSV9 - Updated: 7/13/21
clear
clc
warning('off')
undifGPEList = importdata("undifGPEList_109_1xR.mat");
undifGravList = importdata("undifGravList_109_1xR.mat");

    %% Choose material 

    rho0F = 7874*0.875; %% kg/m^3 
    rho0M = 4103; 
    

    GPIncrement = zeros(1,799);
    massList = zeros(1,199);
    deltaRList = zeros(1,199);
    coreMass = zeros(1,199);
    mantleMass = zeros(1,199);
    GPEofLayer = zeros (1,199);

w = waitbar(0, 'Calculating...');
 
for i=109
    %% initial conditions 
    prop = 4*i;
    rE = 1; %earth radii
    r = (6371*10^3)*rE;
    V = (4/3)*pi*r^3; %volume (m^3)
    m = V*rho0M; %mass kg assuming constant density
    rList = linspace(0.01*r, r, 800); %List of radii, rList(1) is bottom edge of first layer startin at 0.01 earth radius *
    rhoList(1:prop) = rho0F; %assume constant density, kg/m^3 
    rhoList(prop:800) = rho0M;


    PLayer1 = gravPresCalc(rList, rhoList, prop); %pressure for first iteration
    rhoList1 = rhoList; %constant density first iteration 
    PLayer2 = zeros(1,length(PLayer1));
    rhoList2 = zeros(1,length(rhoList1)); 
    
    prop = i*4;
    rhoList2((prop+1):800) = calcDensity(rhoList1((prop+1):800), PLayer1(prop+1:800), 4103, 265.5, 4.16,1); %calculate density of mantle 4176
    rhoList2(1:prop) = calcDensity(rhoList1(1:prop), PLayer1(1:prop), 78740.875, 162.5, 5.5,1); %8430
    [PLayer2, finTotalMassLayer, gravList, minside, rlistSquared] = gravPresCalc(rList, rhoList2, prop); 
    x=0
    %% iterate until pressure and density match 
    while(any(abs(PLayer1 - PLayer2) > 0.01) || any(abs(rhoList1 - rhoList2) > 10))
    r = (6371*10^3)*rE;
        if(rhoList2(2) > 0) %does not effect first cycle
            PLayer1 = PLayer2; %update old pressure and density to match new 
            rhoList1 = rhoList2;
        end
        prop = i*4;
        rhoList2((prop+1):800) = calcDensity(rhoList1((prop+1):800), PLayer1(prop+1:800), 4103, 265.5, 4.16,1); %calculate density of mantle, 0.98 accounts for thermal expansion
        rhoList2(1:prop) = calcDensity(rhoList1(1:prop), PLayer1(1:prop), 78740.875, 162.5, 5.5,1); %0.94 accounts for thermal expansion, melting, and light elements
        [PLayer2, finTotalMassLayer, gravList, minside, rlistSquared] = gravPresCalc(rList, rhoList2, prop);   
        massList(i) = sum(finTotalMassLayer);
%         rhoList2 = calcDensity(rhoList1, PLayer1); %new density under updated pressure
%        (abs(mean(rhoList1(1:length(rhoList1)),'all') - mean(rhoList2(1:length(rhoList2)),'all')))
%        
%         [PLayer2,m] = gravPresCalc(rList,rhoList2);%new pressure under updated density 
%        (abs(mean(PLayer1(2:length(rhoList1)),'all') - mean(PLayer2(2:length(rhoList2)),'all')))
        %a = "-------------------------------"


    end 
    massList(i) = sum(finTotalMassLayer);
        %checks if masses of undifferentiated and differentiated planets match
   
    %checks if masses of undifferentiated and differentiated planets match
    deltaR = r - 6356.752*10^3;
    FeMassTotal = sum(finTotalMassLayer(1:prop));
    MgMassTotal = sum(finTotalMassLayer((prop+1):end));
    innerMass = sum(finTotalMassLayer(1:400));
    outerMass = sum(finTotalMassLayer(401:800));
    massTotal = innerMass + outerMass;
    
    
    % gravity calculations- only done once masses match and pressure /density converge
    PList(i,1) = PLayer2(1);
    gravList(800) = 0;
    j = 798;
    while j > 0 
        k = 799;
        while k > 1 
            GPIncrement(k) = finTotalMassLayer(j)*gravList(k)*(rList(k)-rList(k-1));
            k = k - 1;
        end
        GPEofLayer(j) = sum(GPIncrement);
        j = j - 1;
    end
    tGravPotential = sum(GPEofLayer);
    FeGravPotential = sum(GPEofLayer(1:prop));
    mgGravPotential = sum(GPEofLayer((prop+1):end));

    difGPEList(i) = -tGravPotential;
    coreMass(i) = sum(finTotalMassLayer(1:prop));
    mantleMass(i) = sum(finTotalMassLayer((prop+1):end));
    CorePressure_50(i) = PLayer2(prop/2);
    CorePressure_25(i) = PLayer2(prop/4);
    CorePressure_75(i) = PLayer2(prop*(3/4));
    centerPressureList(i) = PLayer2(prop);
    
    
    waitbar(x/5)
     
end
close(w)
% %% Create Plot 
 PEofDif = undifGPEList - difGPEList;
 Me = massList(i)/(5.973*10^24)
 eMass = 5.9724e24;
 earthMassList = massList/eMass;
 earthRadiusList = 1:800;
 plot(earthRadiusList,rhoList2);
 plot(earthRadiusList,PLayer2);

%% Calculate Pressure
function [PLayer, mlayer, grav, minside, rlistSquared] = gravPresCalc(rList, rhoList, prop)
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
    for j = 1:(length(rList))-1 %calculate mass inside 
        minside(j) = sum(mlayer(1:j));
        rlistSquared(j) = ((rList(j+1))^2);
        grav(j) = G*minside(j)/((rList(j+1))^2);% gravity is highest at r = 500, so layers there have the most potential energy. Really not sure what to do about this :(
        dp(j) = rhoList(j)*grav(j)*(rList(end)/length(rList));
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
guess = 1;
options = optimset('fsolve'); %, 'Tolx', 1e-8); 
 for i=1:length(rhoList)
    Peqn = @(x) 3*k0 * (x^(-2)) * (1-x)*exp((1.5*k0p-1.5)*(1-x)) - PLayer(i); %Vinet
    xList(i) = fsolve(Peqn, guess,options);
    guess = xList(i);

 end 
    newRhoList = rho0./(xList.^3); %rho = rho0*n
end