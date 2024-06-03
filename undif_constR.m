   %EOSV9 - Updated: 7/13/21
clear
clc
warning('off')
massRatio = importdata("massRatio_109_1xR.mat");
difMassList = importdata("difMassList_109_1xR.mat");
massRatio(94) = 4;
    %% Choose material 

    rho0F = 8430; %% kg/m^3 
    rho0M = 4176; 

w = waitbar(0, 'Calculating...');
q = 109
for i=109
    %% initial conditions 
    prop = i; %out of 200
    rE = 1; %earth radii
    r = ((6371*10^3)*rE);
    V = (4/3)*pi*r^3; %volume (m^3)
    m = V*rho0M; %mass kg assuming constant density
    rList = linspace(0.01*r, r, 800); %List of radii, rList(1) is bottom edge of first layer startin at 0.01 earth radius 
    rhoList = zeros(1,800);

    rhoList(2:2:end) = rho0F; %Iron = even
    rhoList(1:2:end) = rho0M; %MgSiO = odd
    
    PLayer1 = gravPresCalc(rList, rhoList, i, massRatio); %pressure for first iteration 
    rhoList1 = rhoList; %constant density first iteration 
    PLayer2 = zeros(1,length(PLayer1));
    rhoList2 = zeros(1,length(rhoList1)); 
    finTotalMassLayer = zeros(1,length(rList));
    k = 1;

    rhoList2(2:2:end) = calcDensity(rhoList1(2:2:end), PLayer1(2:2:end), 7874*0.895, 162.5, 5.5, 1, k); %8430
    rhoList2(1:2:end) = calcDensity(rhoList1(1:2:end), PLayer1(1:2:end), 3910*1.01, 265.5, 4.16, 1, k); %4176
    [PLayer2, finTotalMassLayer, gravList, thicknessRatio, thicknessPercentage] = gravPresCalc(rList, rhoList, prop, massRatio);
    %% iterate until pressure and density match 
    while(any(abs(PLayer1 - PLayer2) > 0.01) || any(abs(rhoList1 - rhoList2) > 10))
        if(rhoList2(2) > 0) %does not effect first cycle
            PLayer1 = PLayer2; %update old pressure and density to match new 
            rhoList1 = rhoList2;
        end

        for p = 2:length(rList) - 2
            rList(p) = rList(p-1)+(thicknessPercentage(p)*(2*(rE*6356.752*10^3)/800));
        end
        rhoListDif = abs(mean(rhoList1(1:length(rhoList1)),'all') - mean(rhoList2(1:length(rhoList2)),'all'));
        PLayerDif = abs(mean(PLayer1(1:length(PLayer1)),'all') - mean(PLayer2(1:length(PLayer2)),'all'));

        rhoList2(2:2:end) = calcDensity(rhoList2(2:2:end), PLayer1(2:2:end), 7874*0.895, 162.5, 5.5, 1, k); 
        rhoList2(1:2:end) = calcDensity(rhoList2(1:2:end), PLayer1(1:2:end), 3910*1.01, 265.5, 4.16, 1, k); %calculate density of mantle
        [PLayer2, finTotalMassLayer, gravList, thicknessRatio, thicknessPercentage] = gravPresCalc(rList, rhoList2, i, massRatio);%calculate pressure @ each layer
       
%         rhoList2 = calcDensity(rhoList1, PLayer1); %new density under updated pressure
%        (abs(mean(rhoList1(1:length(rhoList1)),'all') - mean(rhoList2(1:length(rhoList2)),'all')))
%        
%         [PLayer2,m] = gravPresCalc(rList,rhoList2);%new pressure under updated density 
%        (abs(mean(PLayer1(2:length(rhoList1)),'all') - mean(PLayer2(2:length(rhoList2)),'all')))
        %a = "-------------------------------"
    
    end


        
    massList(i) = sum(finTotalMassLayer);%(mean(rhoList2(1,1:length(rhoList2)),'all'))*(4/3)*pi*r^3;
    mgMassTotal = sum(finTotalMassLayer(1:2:end)); 
    FeMassTotal = sum(finTotalMassLayer(2:2:end));
    innerMass = sum(finTotalMassLayer(1:400));
    outerMass = sum(finTotalMassLayer(401:800));
    massTotal = innerMass + outerMass;
    PList(i,1) = PLayer2(1);
    gravList(800) = 0;
    %%old version gravPotentialofLayer = finTotalMassLayer.*gravList.*rList; 

    
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
    mgGPE = sum(GPEofLayer(1:2:end)); 
    FeGPE = sum(GPEofLayer(2:2:end));
    undifGPEList(i) = -tGravPotential;
    q = q+1;
    waitbar(q/800)

end
    if (massList(i) < difMassList(i))
        while ((massList(i) < difMassList(i)))
            r = r*1.005
            rList = linspace(0.01*r, r, 800);
            prop = i;
            rhoList2(2:2:end) = calcDensity(rhoList1(2:2:end), PLayer1(2:2:end), 7874*0.895, 162.5, 5.5, 1, k); %8430
            rhoList2(1:2:end) = calcDensity(rhoList1(1:2:end), PLayer1(1:2:end), 3910*1.01, 265.5, 4.16, 1, k); %4176
            [PLayer2, finTotalMassLayer, gravList, thicknessRatio, thicknessPercentage] = gravPresCalc(rList, rhoList2, prop, massRatio);  
            massList(i) = sum(finTotalMassLayer)
        end
    end
    if (massList(i) > difMassList(i))
        while ((massList(i) > difMassList(i)))
            r = r*0.995
            rList = linspace(0.01*r, r, 800);
            prop = i;
            rhoList2(2:2:end) = calcDensity(rhoList1(2:2:end), PLayer1(2:2:end), 7874*0.895, 162.5, 5.5, 1, k); %8430
            rhoList2(1:2:end) = calcDensity(rhoList1(1:2:end), PLayer1(1:2:end), 3910*1.01, 265.5, 4.16, 1, k); %4176
            [PLayer2, finTotalMassLayer, gravList, thicknessRatio, thicknessPercentage] = gravPresCalc(rList, rhoList2, prop, massRatio);  
            massList(i) = sum(finTotalMassLayer)
        end
    end
% %% Create Plot 
 eMass = 5.9724e24;
 earthMassList = massList/eMass;
 earthRadiusList = 0.01:0.01:1.99;
 feRhoList = rhoList2(2:2:end);
 mgRhoList = rhoList2(1:2:end);

 plot(rList(1:800),rhoList2)
 

%% Calculate Pressure
function [PLayer, totalMassLayer, grav, thicknessRatio, thicknessPercentage] = gravPresCalc(rList, rhoList, coreRadius, massRatio)
        G = 6.67408*10^-11; %gravitational constant 
    %gravity and mass calculations 
    thicknessRatio = zeros(1,length(rhoList)/2);
    for k=1:2:length(rhoList)
        thicknessRatio(k) = (massRatio(coreRadius))/(rhoList(k)/rhoList(k+1)); %%should be thickness of Mg over thicknessFe does this iteration converge?
    end
    feMassLayer = zeros(1,length(rList));
    mgMassLayer = zeros(1,length(rList));
    totalMassLayer = zeros(1,length(rList));
    for j = 1:2:(length(rList)-2)
        %calculate layer masses 
        thicknessPercentage(j+1) = ((1/thicknessRatio(j))/((1/thicknessRatio(j))+1)) %%this is a volume percentage for Fe
        thicknessPercentage(j) = (thicknessRatio(j)/(thicknessRatio(j)+1)) %%this is a volume percentage for Mg
        feMassLayer(j+1) = rhoList(j+1)*thicknessPercentage(j+1)*((rList(j+3)^3 - rList(j+1)^3)*(4/3)*pi); %layer 1 goes from 0.01 to 0.03 earth radii 
        mgMassLayer(j) = rhoList(j)*thicknessPercentage(j)*((rList(j+2)^3 - rList(j)^3)*(4/3)*pi);
    end

    for q = 1:800
        totalMassLayer(q) = feMassLayer(q) + mgMassLayer(q);
    end
    
    minside = zeros(1,length(rList));
    grav = zeros(1,length(rList)-1);
    dp = zeros(1,length(rList)-1);
    for j = 1:(length(rList))-3 %calculate mass inside 
        minside(j) = sum(totalMassLayer(1:j));
        grav(j) = G*minside(j)/((rList(j+1))^2);
        weightedRhoAvg(j) = rhoList(j)*thicknessPercentage(j) + rhoList(j+1) * thicknessPercentage(j+1);
        dp(j) = weightedRhoAvg(j)*grav(j)*(rList(end)/length(rList)); % redo this
    end 
    %pressure calculation - this is all off, do tomorrow

    for p = length(dp):-1:1
        PLayer(p) = sum(dp(p:length(dp)))*10^-9;
    end
    PLayer(length(rList)) = 0; %no pressure on the surface
end

%% Solver 
function newRho = calcDensity(rhoList, PLayer, rho0, k0, k0p, EOS, k)
xList = zeros(1,length(rhoList)); 
%solve each layer 
options = optimset('fsolve'); %, 'Tolx', 1e-8);
guess = 1

 for i=1:length(rhoList)
    Peqn = @(x) 3*k0 * (x^(-2)) * (1-x)*exp((1.5*k0p-1.5)*(1-x)) - PLayer(i); %Vinet
    xList(i) = fsolve(Peqn, guess,options);
    guess = xList(i);
 end 
 
newRho = rho0./(xList.^3) %rho = rho0*n
end