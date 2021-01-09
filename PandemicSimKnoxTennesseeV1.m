% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Here is a model for the spread of a pandemic. The spread of a disease
%%%%% in a cell will go as dI/dt = Co*exp(-r(t-to)^2). It might not,
%%%%% though. Pick a cell size between 30 and 60. Let MATLAB decide how to
%%%%% size the cells. Cells should overlap. If someone in your cell is a 2,
%%%%% you should be much more likely to become a 2 than if no one in your
%%%%% cell is a 2. How likely you are to become a 2 will have an Arrhenius
%%%%% dependence. Each cell will have a temperature. The spread condition
%%%%% will be an individual event. A value of '0' will imply that the
%%%%% individual does not have covid. A value of '1' will indicate that
%%%%% the individual has covid and is contagious. A value of '2' will
%%%%% indicate that the individual is removed from the stream. Base this
%%%%% spread on the Arrhenius rate law :-)

% function [dIdt7Sim] = PandemicSimFunc
%%%%% Partition the U.S. into a number of cells. America is a matrix of
%%%%% cells with somewhere between 25 and 100 individuals.
% clear variables; close all; clc;

GoodFit = 0; R2D2Best = 1e9; 

Threshold = 0.0675; % Threshold = 0.0675;
CloseVal = 9; FarVal = 13; % CloseVal = 10; FarVal = 16;
CellOdds = 5; 
while GoodFit == 0
% RNGSeed = randi([-5, 5]); rng(RNGSeed);
% nIndividuals counts how many individuals the simulation thinks live in
% the county. The cells are square. This block builds the county. 
    nIndividuals = 0; k1Max = 76; k2Max = 76;
    America = cell(k1Max, k2Max); % cell(2715, 2715);
    for i = 1:numel(America)
        nSize = randi([5, 10]);
        America{i} = zeros(nSize, nSize, 2, 'single');
        %     Duration{i} = America{i};
        nIndividuals = nIndividuals + nSize*nSize;
    end
fprintf('Number of People in Simulation: %i \n', nIndividuals)
% pause

%%%%% Seed the infection. Right now, dictate where infection will start %%%
    CellOne = America{2};
    CellOne(2,2,1) = 1;
    CellOne(2,2,2) = 1;
    America{2} = CellOne;

    CellTwo = America{3};
    CellTwo(3,2,1) = 1;
    America{3} = CellTwo;

    CellThree = America{5};
    CellThree(2,3,1) = 1;
    America{5} = CellThree;

% Duration is the memory of how long the individual has been infected.
% Duration = America;

%%%%% What happens if a cell has an infection? %%%%%
% Is the cell infected? If 1, then yes. If '0', then no. If '2', then cell
% cannot be infected because it was infected. rngMax is the max value of
% the random number generator. 'Contagious' is a measure of how contagious
% the virus is. The greater the number, the more likely an infection.
% Instead of using a random number generator to determine infection, I am
% using the volume diffusion probability from the sintering model. Rho is a
% density, and will be the fraction of those infected in the cell. C is a
% constant that will be set to 1.00001. Do is a constant with value of
% 6.8668. This was determined from previous work on sintering. RT is a
% measure of the kinetics of the community. The nonzero function is used to
% determine the number of people who had the virus in the cell.

rngMax = 1e6; Contagious = 1.5*2.314497e-16; 
gbContagious = 1.5*2.314497e-16;
ConTime = 14; Do = 6.8668; C = 1e-8; R = 8.314; Tmin = 1380; Tmax = 1570;
nDays = 253; time = (linspace(1, 250, 250))'; T = zeros(nDays, 1);

QVol = 1.7555*475e3*ones(nDays, 1); % 1.75*475e3
Qgb = 1.25*532e3*ones(nDays, 1); % 1.2*532e3

DayVec = zeros(nDays, 1, 'single');
dIdtVec = zeros(nDays, 1, 'single');
InterdIdtVec = zeros(nDays, 1, 'single');
IntradIdtVec = zeros(nDays, 1, 'single');

InterEx = 0; IntraEx = 0; InterTrans = 0; IntraTrans = 0;
InterExpAtt = 0; IntraExpAtt = 0; 
InterTransProbVec = zeros(nDays, 1); IntraTransProbVec = zeros(nDays, 1);
InterExProbVec = zeros(nDays, 1); IntraExProbVec = zeros(nDays, 1); 
for day = 1:nDays
    DailyDiff = 100; ECount = 0;
    if day > 1
        QVol(day) = QVol(day - 1); Qgb(day) = Qgb(day - 1);
    end
    while DailyDiff > 10
    
%     Threshold = 0.0675; % Threshold = 0.0675;
%     CloseVal = 9; FarVal = 13; % CloseVal = 10; FarVal = 16;
%     CellOdds = 5; 
        
% From day to day, it is reasonable to expect activation energy to change.
% The values that control the spread are the exposure values, the contagion
% values, and the threshold for intercell spread. In populations that tend
% to mix more, the threshold will be lower. In more segregated communities,
% threshold will be higher. 
    
%     if mod(day, 3) == 0
%         CellOdds = 4;
%     else
%         CellOdds = 4;
%     end
%     if day > 170
%         CellOdds = 3; QVol = 1.7*475e3; Qgb = 1.325*532e3;
%         
%     end
    fprintf('Day %i \n', day)
%     if day > 15 && day <= 209
%         QVol = 1.2*441e3;
%         Qgb  = 1.2*544e3;
%     end
%     if day > 210
%         QVol = 0.8*441.5e3;
%         Qgb = 0.8*544e3;
%     end        
        
    InterdIdt = zeros(1, 'single'); IntradIdt = zeros(1, 'single');
%     k1Max = 250; k2Max = 250;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Have to infect neighboring cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k1 = 1:k1Max % 15
        for k2 = 1:k2Max % 15
% Pick a cell
            InterStack = America{k1,k2};
            InterMatInf = InterStack(:,:,1);
            InterMatRec = InterStack(:,:,2);
            
% Check to see how many are infected in the cell
            GateOne = InterMatInf == 1;
            
% Check to see how many are contagious in the cell
            GateTwo = InterMatRec >= 2;
            
% Look at how many are actually contagious
            SumInf = 0;
            for iii = 1:numel(GateOne)
                if GateOne(iii) == 1 && GateTwo(iii) == 1
                    SumInf = SumInf + 1;
                end
            end
            
% Find fraction of cell infected and contagious
            FracInf = SumInf/numel(InterMatInf);
            
% If there are infected individuals in a cell, then the
% infection can spread to a neighboring cell. This is to see if
% there is exposure.
% Do you interact outside your cell?
% Cell temp is a measure of how a community moves. It is the kinetics of
% the community.
            

            if SumInf > 0
                % Look at neighboring cells
            for ii = -1:1
                for jj = -1:1
                    %                         if randi([1,10]) == randi([1, 10])
                    
                    if day < 651
                        T(day) = 2850; % randi([2200, 2800]);
%                         T(day) = Tmin + 1/32*randi([64, 512])*(Tmax-Tmin);
%                     elseif day >= 150 && day <= 202
%                         T(day) = randi([2300, 2900]); % randi([2600, 3200]);
% %                         T(day) = Tmin + 1/32*randi([64, 1024])*(Tmax-Tmin);
%                     elseif day > 202
%                         T(day) = randi([6800, 8800]);
% %                         T(day) = Tmin + 1/32*randi([64, 2048])*(Tmax-Tmin);
                    end
                    
% The greater the contagious fraction, the more likely the spread.
                    InterExpAtt = InterExpAtt + 1;
                    if FracInf > Threshold &&...
                            1/(1 - FracInf)*exp(-Qgb(day)/(R*T(day))) > gbContagious
                        InterEx = InterEx + 1;
% if SumInf > 0 && 1/(1+C-FracInf)*exp(-Qgb/(R*T)) > gbContagious
                        
% dx, dy: to make sure I am not changing variables within a loop. Also need
% to stay in bounds
                        dx = ii; dy = jj;
                        if (k1 + dx) > k1Max || (k1 + dx) == 0
                            dx = 0;
                        end
                        if (k2 + dy) > k2Max || (k2 + dy) == 0
                            dy = 0;
                        end
                        
% Each neighbor should have a chance to get infected
%                         if (dx ~= 0 && dy ~= 0)
                        if (abs(dx) + abs(dy)) == 1
                            Enya = CloseVal;
                        elseif (abs(dx) + abs(dy)) == 2
                            Enya = FarVal;   
%                         elseif (abs(dx) + abs(dy)) == 0
%                             Enya = 785;
                            
                            % This is to see if the exposure leads to infection
                            if randi([1, Enya]) == randi([1, Enya])
                                InterTrans = InterTrans + 1;
%                                 InterCellSpread = InterCellSpread + 1;
                                % 40 40
                                
                                % Pick the cell to be infected, and modify the appropriate record matrix
                                NeighborStack = America{k1+dx, k2+dy};
                                NeighborMat = NeighborStack(:,:,1);
                                NeighborRec = NeighborStack(:,:,2);
                                
                                % Get the size to see who gets infected, according to the RNG
                                [NeiRow, NeiCol] = size(NeighborMat);
                                NeiX = randi([1,NeiRow]);
                                NeiY = randi([1,NeiCol]);
                                
                                % if the randomly selected individual in the neighboring matrix is healthy,
                                % they will now get sick.
                                if NeighborMat(NeiX, NeiY) == 0
                                    NeighborMat(NeiX, NeiY) = 1;
                                    NeighborRec(NeiX, NeiY) = 1;
                                    %fprintf('Intercell Infection \n')
                                    % Daily infection rate increases
                                    InterdIdt = InterdIdt + 1;
%                                     InterTrans = InterTrans + 1;
                                    
                                    % Matrices tracking infected are modified and cell is modified.
                                    NeighborStack(:,:,1) = NeighborMat;
                                    NeighborStack(:,:,2) = NeighborRec;
                                    America{k1+dx, k2+dy} = NeighborStack;
%                                     fprintf('Intracell Infection: (%i, %i) cell (%i, %i) \n',...
%                                         NeiX, NeiY, k1, k2)
                                end
                            end
                        end
                        clear NeighborMat; clear NeighborStack;
                        %                         end
                    end
                end
            end
            end
        end
    end
    
    %%%%% This loop is the intercell virus spread %%%%%
    clear i
    for k = 1:numel(America)
        %     fprintf('Loop %i \n', i)
        TheMatrix = America{k}; CelldIdt = 0;
        DisMat = TheMatrix(:,:,1);
        
        % This checks to see if someone in the cell is infected. If so, check to
        % see if they were exposed.
        SumOne = sum(sum(TheMatrix(:,:,1) == 1));
        
        %%% Enter this loop if there is a contagious individual in the cell %%%
        if SumOne > 0
%             fprintf('Virus detected in Cell %i \n', k)
            %             InfInCell = SumOne; Cont = sum(sum(TheMatrix(:,:,2) > 1));
            %             ContInCell = (TheMatrix(:,:,1) == 1) == (TheMatrix(:,:,2) > 1);
            %             ContFrac = sum(sum(ContInCell))/numel(ContInCell);
            %             fprintf('Fraction of cell %i infected and contagious: %5.4f \n', k, ContFrac)
            [Col, Row] = size(TheMatrix(:,:,1));
            
            for i = 1:Col
                for j = 1:Row
                    %%% If an individual is infected, then there is a chance someone else in
                    %%% the cell will be infected %%% You can only be infected if you have not
                    %%% been infected
                    if TheMatrix(i,j,1) == 1
                        
                        % Modify how long each individual has been inected
                        TheMatrix(i,j,2) = TheMatrix(i,j,2) + 1;
                        
                        %%% No longer contagious after 14 days %%%
                        if TheMatrix(i,j,2) > ConTime
                            TheMatrix(i,j,1) = 2;
                            
                            % fprintf('Cell %i Modified \n', k)
                        end
                        %%% Condition for NN and NNN infection %%%
                        % Redefine the cell to account for those no longer infected
                        America{k} = TheMatrix;
                    end
                end
            end
% Define the matrix to be modified
%                         DisMat = TheMatrix(:,:,1);
% See how many are infected
                    
                    CellStart = sum(sum(TheMatrix(:,:,1) == 1));
% See how many can be infected, maybe???
%                     Remaining = numel(TheMatrix(:,:,1)) - nnz(TheMatrix(:,:,1));
                    Rho = SumOne/numel(TheMatrix(:,:,1));
%                         fprintf('Fraction of cell %i infected: %5.4f \n', k, Rho)
%                         if Remaining > 0
                    if CellStart >= round(numel(TheMatrix(:,:,1))/2)
                        nEnd = numel(TheMatrix(:,:,1));
                    else 
                        nEnd = 2*CellStart;
                    end
                    
                    for n = 1:nEnd % numel(TheMatrix(:,:,1))
% Will the infected individual be exposed to the ones with CoVid
%                         if day < 150
%                             T = Tmin + 1/32*randi([-128, 128])*(Tmax-Tmin);
%                         elseif day >= 150 && day <= 202
%                             T = Tmin + 1/32*randi([-256, 256])*(Tmax-Tmin);
%                         elseif day > 202 
%                             T = Tmin + 1/32*randi([-1024, 1024])*(Tmax-Tmin);
%                         end
%                         T = Tmin + (Rho)*(Tmax - Tmin) + randi([-512, 128]);
                        IntraExpAtt = IntraExpAtt + 1;
                        if 1/(Rho)*exp(-QVol(day)/(R*T(day))) > Contagious
                            IntraEx = IntraEx + 1;
%%% if loop in here to simulate the random nature of catching a disease %%%
%                                 if randi([-10, 10]) < 20
%                             fprintf('Exposure \n')
% Individual immunity to covid??? %
                            if randi([1,CellOdds]) == randi([1,CellOdds])
                                IntraTrans = IntraTrans + 1;
                                xInf = randi([1,Row]);
                                yInf = randi([1,Col]);
                                if DisMat(xInf, yInf) == 0
                                    %                                 figure; spy(DisMat)
%                                     fprintf('Individual (%i, %i) now infected \n',...
%                                         xInf, yInf)

%                                     if (randi([1, 10]) + randi([1, 10])) < 21 % randi([2, 20])
%                                     IntraTrans = IntraTrans + 1;
                                    DisMat(xInf, yInf) = 1;
                                    IntradIdt = IntradIdt + 1;
                                    CelldIdt = CelldIdt + 1;
%                                     end
%                                         fprintf('Number Infected on Day %i in Cell %i: %i \n',...
%                                             day, k, CelldIdt)
                                end
                            end
                            %                                 end
                        end
                    end
                    if CelldIdt > 2*CellStart
                        fprintf('Too many infections \n')
                    end
%                 end
                %                         end
                %                     end
%             end
        end
        dIdt = InterdIdt + IntradIdt;
        if dIdt > 0
            TheMatrix(:,:,1) = DisMat; America{k} = TheMatrix;
%             fprintf('Cell %i modified \n', k)
        end
    end
    clear TheMatrix; clear dIdt; clear DisMat;

    dIdt = InterdIdt + IntradIdt;
    DayVec(day) = day;
    InterdIdtVec(day) = InterdIdt;
    IntradIdtVec(day) = IntradIdt;
    dIdtVec(day) = dIdt;
%     CelldIdtVec(day) = CelldIdt;
% fprintf('dIdt: %i \n', dIdt)


% For each day, check daily infection rate of simulation against daily
% infection rate of data. 
% Threshold = 0.0675; % Threshold = 0.0675;
% CloseVal = 9; FarVal = 13; % CloseVal = 10; FarVal = 16;
% CellOdds = 5; 

% QVol = 1.7555*475e3; % 1.75*475e3
% Qgb = 1.25*532e3; % 1.2*532e3

if day > 7
% Get a 7-day average
    dIdt7 = mean(dIdtVec(day-7:day));
    DailyDiff = abs(dIdt7 - didt7Data(day));
% elseif day <= 7
%     DailyDiff = abs(dIdtVec(day) - didt7Data(day);
    if DailyDiff > 10
%         RNGSeed = randi([0, 24]); rng(RNGSeed);
%         rng default
        ECount = ECount + 1;
        if ECount < 16
% Suppose simulation is too high.
            if dIdt7 > didt7Data(day)
% There should be a critical number where the spread will be controlled
% by grain-boundary activation energy. Before that, volume diffusion...
                if sum(dIdtVec) < 215
                    QVol(day) = 1.0125*QVol(day);
%                     Threshold = 1.01*Threshold;
                else
                    QVol(day) = 1.0125*QVol(day); Qgb(day) = 1.0125*Qgb(day);
%                     Threshold = 1.01*Threshold;

                end
                if ECount > 3
                        fprintf('Threshold too low: Increasing Threshold \n')
                        Threshold = 1.01*Threshold;
                        if ECount > 7
                            CellOdds = CellOdds + 1;
                            fprintf('Cell Contagiousness Decreased \n')
                        end
                end
                
% Suppose simulation is too low :-)
            elseif dIdt7 < didt7Data(day)
                if sum(dIdtVec) < 215
                    QVol(day) = 0.9875*QVol(day);
%                     Threshold = 0.99*Threshold;
                else
                    QVol(day) = 0.9875*QVol(day); Qgb(day) = 0.9875*Qgb(day);
%                     Threshold = 0.99*Threshold;
                end
                if ECount > 3
                        fprintf('Threshold too high: Decreasing Threshold \n')
                        Threshold = 0.99*Threshold;
                        if ECount > 7
                            CellOdds = CellOdds - 1;
                            fprintf('Cell Contagiousness Increased \n')
                        end
                end
                
            end
            if CellOdds < 1
                CellOdds = 1;
            end
        else 
            QVol(day) = QVol(day - 1); Qgb(day) = Qgb(day - 1);
            fprintf('Energy change too extreme \n')
            DailyDiff = 5;
        end
    end
else
    DailyDiff = 5;
end



    end
    
% figure; plot(DayVec, dIdtVec, 'bo')
% xlabel('Time (Days)'); ylabel('Daily Infection Rate')
% set(gca, 'fontsize', 14, 'box', 'off')
% 
% TotInf = cumtrapz(DayVec, dIdtVec);
% figure; plot(DayVec, TotInf, 'gs')
% xlabel('Time (Days)'); ylabel('Total Infected')
% set(gca, 'fontsize', 14, 'box', 'off')

% for i = 1:10000
%     SpyMat = America{i};
%     if sum(sum(SpyMat)) > 1
%         figure; spy(SpyMat(:,:,1))
%     end
% end

% figure; plot(InterdIdtVec, 'bo')
% hold on; plot(IntradIdtVec, 'rp'); hold off;

dIdt7Sim = zeros(numel(dIdtVec) - 7, 1);
for n = 8:numel(dIdtVec)
    dIdt7Sim(n-7) = mean(dIdtVec(n-7:n));
end

SimNum = numel(dIdt7Sim);
if SimNum >= numel(didt7Data)
    SimComp = dIdt7Sim(1:numel(didt7Data));
    DataComp = didt7Data;
    R2 = (dIdt7Sim(1:numel(didt7Data)) - didt7Data);
    R2D2 = R2'*R2;
elseif SimNum < numel(didt7Data)
    SimComp = dIdt7Sim;
    DataComp = didt7Data(1:SimNum);
    R2 = (dIdt7Sim - DataComp);
    R2D2 = R2'*R2;
end

if R2D2 <= R2D2Best
    R2D2Best = R2D2; 
    
    InterTransProb = InterTrans/InterEx;
    IntraTransProb = IntraTrans/IntraEx;
    InterExProb = InterEx/InterExpAtt;
    IntraExProb = IntraEx/IntraExpAtt;
    
    if day == nDays
        figure; plot(dIdt7Sim, 'bo')
        hold on; plot(TimeYo(12:end)-TimeYo(12), didt7Data, 'k.')
        TheLeg = legend('Simulation', 'Data');
        set(TheLeg, 'Location', 'Best', 'box', 'off')
        xlabel('Time (Days)'); ylabel('Daily Infection Rate')
        set(gca, 'fontsize', 14, 'box', 'off')
        
        figure; plot(cumtrapz(dIdt7Sim), 'bo')
        hold on; plot(TimeYo(12:end)-TimeYo(12), cumtrapz(didt7Data), 'k.')
        TheLeg = legend('Simulation', 'Data');
        set(TheLeg, 'Location', 'Best', 'box', 'off')
        xlabel('Time (Days)'); ylabel('Total Infected')
        set(gca, 'fontsize', 14, 'box', 'off')
    end
    dIdt7SimBest = dIdt7Sim;
    fprintf('Between cell Transmission Probability: %5.4f \n', InterTransProb)
    fprintf('Inside cell Transmission Probability: %5.4f \n', IntraTransProb)
    fprintf('Between cell Exposure Probability: %5.4f \n', InterExProb)
    fprintf('Inside cell Exposure Probability: %5.4f \n', IntraExProb)
    fprintf('R-Squared: %5.4f \n', R2D2Best)
    if day < nDays
        GoodFit = 0;
    elseif day == nDays
        GoodFit = input('Is Fit Good \n 1: Yes \n 0: No\n');
    end
% elseif R2D2 > R2D2Best
%     if mean(dIdt7Sim - didt7Data(1:SimNum)) < 0
%         fprintf('Threshold too high: Decreasing Threshold \n')
%         Threshold = 0.999*Threshold;
%     elseif mean(dIdt7Sim - didt7Data(1:SimNum)) > 0
%         fprintf('Threshold too low: Increasing Threshold \n')
%         Threshold = 1.001*Threshold;
%     end
    
    if day == nDays
        figure; plot(dIdt7Sim, 'bo')
        hold on; plot(TimeYo(12:end)-TimeYo(12), didt7Data, 'k.')
        TheLeg = legend('Simulation', 'Data');
        set(TheLeg, 'Location', 'Best', 'box', 'off')
        xlabel('Time (Days)'); ylabel('Daily Infection Rate')
        set(gca, 'fontsize', 14, 'box', 'off')
        
        figure; plot(cumtrapz(dIdt7Sim), 'bo')
        hold on; plot(TimeYo(12:end)-TimeYo(12), cumtrapz(didt7Data), 'k.')
        TheLeg = legend('Simulation', 'Data');
        set(TheLeg, 'Location', 'Best', 'box', 'off')
        xlabel('Time (Days)'); ylabel('Total Infected')
        set(gca, 'fontsize', 14, 'box', 'off')
    end
    
end

    InterTransProbVec(day) = InterTrans/InterEx;
    IntraTransProbVec(day) = IntraTrans/IntraEx;
    InterExProbVec(day) = InterEx/InterExpAtt;
    IntraExProbVec(day) = IntraEx/IntraExpAtt;

for i = 1:numel(T)-7
    T7(i) = mean(T(i:i+7));
end
% figure; plot(T7, 'k.')
end
end

figure; plot(dIdt7SimBest, 'bo')
hold on; plot(TimeYo(12:end)-TimeYo(12), didt7Data, 'g^')
TheLeg = legend('Simulation', 'Data');
set(TheLeg, 'Location', 'Best', 'box', 'off')
xlabel('Time (Days)'); ylabel('Daily Infection Rate')
set(gca, 'fontsize', 14, 'box', 'off')

figure; plot(cumtrapz(dIdt7SimBest), 'bo')
hold on; plot(TimeYo(12:end)-TimeYo(12), cumtrapz(didt7Data), 'g^')
TheLeg = legend('Simulation', 'Data');
set(TheLeg, 'Location', 'Best', 'box', 'off')
xlabel('Time (Days)'); ylabel('Total Infected')
set(gca, 'fontsize', 14, 'box', 'off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extrapolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Set time of extrapolation as about half the length of outbreak %%%%%%
nExDays = nDays + floor(nDays/2);

%%%%% Get activation energies; or exposure and contagiousness %%%%%
QVol7 = zeros(numel(QVol) - 7); Qgb7 = QVol7; 
for i = 1:numel(QVol)-7
QVol7(i) = mean(QVol(i:i+7)); Qgb7(i) = mean(Qgb(i:i+7));
end

QVolEx = mean(QVol7); QgbEx = mean(Qgb7);

%%%%% Linear fit of the outbreak %%%%%
% use findpeaks to get the high points :-) Those peaks can be fit with 
% ln(dI/dt) = ln(co) - r(t - to)^2

InterFish = InterTransProbVec.*InterExProbVec;
IntraFish = IntraTransProbVec.*IntraExProbVec;

InterProb = InterFish(end); IntraProb = IntraFish(end);
ProbMean = (InterProb + IntraProb)/2;

%%% For each point, keep the resulting probability mean around where it
%%% started



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Garbage Can %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if day == 1
% %     figure; plot(day, dIdt, 'bo'); hold on;
% else
% %     plot(day, dIdt, 'bo');
% end
% America{k} = TheMatrix;
% load gong.mat; soundsc(y)

%                     for ii = -1:1
%                         for jj = -1:1
%                             %%% INSTEAD of changing the index, set the
%                             %%% index equal to a variable. Then change the
%                             %%% variable.
%                             dx = ii; dy = jj;
%                             if (i + dx) > Col || (i + dx) == 0
%                                 dx = 0;
%                             end
%                             if (j + dy) > Row || (j + dy) == 0
%                                 dy = 0;
%                             end
%                             if DisMat(i+dx, j+dy) == 0 &&...
%                                 (ii ~= 0 && jj ~= 0)
% %                                 disp('If 3 entered')
%                                 Rho = nnz(DisMat(:,:,1))/...
%                                     numel(DisMat(:,:,1));
%                                 fprintf('Fraction of cell %i infected: %5.4f \n', k, Rho)
%
% %                                     America{k} = TheMatrix;

%     if i == 1
%         figure; plot(day, Rho, 'rp'); hold on;
%     else
%         plot(day, Rho, 'rp');
%     end
%     if size(DurationMat) ~= size(TheMatrix)
%         fprintf('Something is Wrong Day: %i Loop %i \n', day, i)
%         pause
%     end

%     clear TheMatrix; clear DurationMat;
%     %%%%% This is the loop for intracell virus spread %%%%%
%     [ColCell, RowCell] = size(America);
%     for i = 1:ColCell
%             for j = 1:RowCell
%                 if sum(sum(America{i,j})) > 1 %% Check the nearest neighbors %%
%                     T = randi([1380, 1570]);  %% Set the kinetics
%
%
%
%                     ('If 2 entered')
%                     if TheMatrix(i,j,2) > 14
%                         TheMatrix(i,j,1) = 2;
%                     end
%                     %%% Condition for NN and NNN infection %%%
%                     for ii = -1:1
%                         for jj = -1:1
%                             if (i + ii) > Col || (i + ii) == 0
%                                 ii = 0;
%                             end
%                             if (j + jj) > Row || (j + jj) == 0
%                                 jj = 0;
%                             end
%                             if TheMatrix(i+ii,j+jj) == 0 && ((ii + jj) > 0)
%                                 disp('If 3 entered')
%                                 if 1/(sqrt(abs(ii+jj)))*1e-6*...
%                                         randi([1, 1e6]) > Contagious
%                                     TheMatrix(i+ii,j+jj,1) = 1;
%                                     dIdt = dIdt + 1;
%                                     disp('Individual Infected')
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
