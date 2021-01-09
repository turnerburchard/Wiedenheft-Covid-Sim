%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sorting data is difficult. So sort the data from the state and county
% names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [timeYo, NInf, dIdt7, FitTime, FitData, XtrapTime, XtrapFit, Duration, LogisticCap, ...
    CapTime, yBest, FishFit, ParaDataX, ParaDataY] = CountyFitsFunk
% clear variables; clear functions; close all; clc

%%%%% This is just the names of the counties %%%%%
filenameYo = 'AllCounties17Oct.xlsx';
[~, StateCounty] = xlsread(filenameYo, 'B2:C784213');

%%%%% Get the dates %%%%%
FishStyx1 = xlsread(filenameYo, 'A2:A784213');

%%%%% Get the cases %%%%%
FishStyx2 = xlsread(filenameYo, 'E2:E784213');

%%%%% Select state and county of interest %%%%%
StateName = input('Enter state name of interest in single quotes: ');
County = input('\nEnter county of interest in single quotes: ');

%%%%% This loop determines the size of the data set for a %%%%%
%%%%% particular county %%%%%
TheCount = 0;
for i = 1:length(StateCounty)
    if strcmp(StateName, StateCounty{i,2}) == 1 &&...
            strcmp(County, StateCounty{i,1}) == 1
        TheCount = TheCount + 1;
    end % End of the if statement
end % End of for loop for counting number of entries for vector

%%%%% This loop selects data associated with the county of interest %%%%%
k = 0;
FishData = zeros(TheCount, 2); StartDate = min(FishStyx1);
for i = 1:length(StateCounty)
    if strcmp(StateName, StateCounty{i,2}) == 1 &&...
            strcmp(County, StateCounty{i,1}) == 1
        k = k + 1;
        FishData(k,1) = FishStyx1(i,1) - StartDate;
        FishData(k,2) = FishStyx2(i,1);
    end % End loop populating county data
end % End for loop for selecting county data
clear StateCounty; clear TheCount; clear FishStyx1; clear FishStyx2;
timeYo = FishData(:,1); NInf = FishData(:,2);

%%%%% Figures of interest (1) NInf and dI/dt (2) time and dI/dt (3) time
%%%%% and NInf %%%%%

%%%%% Get dI/dt %%%%%
dIdt = zeros(numel(NInf), 1);
for j = 2:numel(NInf)
    dIdt(j) = (NInf(j) - NInf(j-1))/(timeYo(j) - timeYo(j-1));
    if isfinite(dIdt) == 0
        fprintf('dIdt NaN %i \n', j)
    end
end
dIdt(1) = dIdt(2);

%%%%% Time averate the data to some number of days given by nDays %%%%%
nDays = 7;
dIdt7 = zeros(numel(dIdt), 1);
for i = nDays+1:numel(dIdt)
    dIdt7(i) = mean(dIdt(i-nDays:i));
end

figure; plot(timeYo, NInf, 'k.');
xlabel('Time Yo!'); ylabel('Number Infected N')
figure; plot(timeYo, dIdt7, 'k.');
xlabel('Time Yo!'); ylabel('Daily Rate dN/dt')
figure; plot(NInf, dIdt7, 'k.');
xlabel('Number Infected N'); ylabel('Daily Rate dN/dt')
% fprintf('Run to line 73')

%%%%% Stuff to get meeting together %%%%%
% Fishy1 = NInf(6:20); Fishy2 = dIdt7(6:20);
% NinjaFit = polyfit(Fishy1, Fishy2, 2); display(NinjaFit);
% CarCap = -NinjaFit(2)/NinjaFit(1); display(CarCap);
% NinjaVal = polyval(NinjaFit, NInf(6:46));
%
% figure; plot(NInf(6:46), dIdt7(6:46), 'k.')
% hold on; plot(NInf(6:46), NinjaVal, 'bo');
% plot(NInf(6:46), dIdt7(6:46)-NinjaVal, 'rp');
% hold off;
% xlabel('Number Infected'); ylabel('Daily Infection Rate')
% legend('Data', 'Fit based on 14 days')
%
% TheLeg = legend('Data', 'Fit based on 14 days', 'Difference',...
%     'Orientation', 'horizontal');
% set(TheLeg, 'box', 'off'); axis([0 225 -2 10]); box off

%%%% Select data for fit of dN/dt and N %%%%%

% LoopLime1 = floor(0.25*numel(NInf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% This is a weird section of code. Take it out. See what it does

for i = 1:10
    if dIdt7(i) == dIdt7(i+1)
%         NInf = NInf(i+1:end);
%         timeYo = timeYo(i+1:end);
        dIdt7 = dIdt7(i+1:end);
    end
end
SumCounty = cumtrapz(dIdt7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
KeepGoing = 1; dIdtCheck = dIdt7; TimeFitCells = cell(24); kFit = 0;
dIdtTimeCheck = dIdt7;  tEnd = 6; close all;
% CapTime = zeros(24, 1);
while KeepGoing == 1
    %%%%% Find the peaks with the N - dN/dt plots %%%%%
    [FishPeaks, Locs] = findpeaks(dIdtCheck);
    
    figure; findpeaks(dIdtCheck);
    SUp = round(max(dIdt7)/8);
    text(Locs, FishPeaks+SUp, num2str((1:numel(FishPeaks))'))
    
    %%%%% Get data associated with first peak %%%%%
    
    % fprintf('From figure, enter number associated with first peak of interest. \n')
    nPeaks = 1; % input('Enter number of peak of interest: \n');
    
    PeakVec = zeros(nPeaks, 1); % yKeep = zeros(nPeaks, 3);
    for i = 1:nPeaks
        fprintf('Enter number associated with peak %i: ', i)
        PeakVec(i,1) = input('\nPeak Number: ');
    end
    
    %%%%% This is the vector of peak locations of interest %%%%%
    FishVec = Locs(PeakVec);
    
    %%%%% Some number of points to the left of the peak will be below some
    %%%%% specifiable threshold or above the value of the peak. Once that is
    %%%%% achieved, a sufficient number of data points have been collected
    
    DragonVec = cell(nPeaks, 1); InfVec = DragonVec; TimeVec = InfVec;
    for j = 1:numel(FishVec)
        k = 0;
        
        while dIdtCheck(FishVec(j) - k) > 1 && ...
                dIdtCheck(FishVec(j) - k) <= dIdtCheck(FishVec(j)) + 2 && ...
                k < FishVec(j)-1
            k = k + 1;
            %         HalfPab(k,1) = dIdt7(FishVec(j) - k);
        end % End while loop
        
        kk = 0;
        
        while dIdtCheck(FishVec(j) + kk + 1) < dIdtCheck(FishVec(j) + kk) + 0.3 || kk < 4
            kk = kk + 1;
            %         HalfPab(k,1) = dIdt7(FishVec(j) - k);
        end % End while loop
        
        DragonVec{j} = dIdtCheck(FishVec(j)-k:FishVec(j)+kk);
        InfVec{j} = SumCounty(FishVec(j)-k:FishVec(j)+kk);
        TimeVec{j} = timeYo(FishVec(j)-k:FishVec(j)+kk);
    end % End for loop through peaks
    
    figure; plot(InfVec{j}, DragonVec{j}, 'k.');
    
    %%%%% From the selected data points, get a curve fit %%%%%
    %     FitTimeYo = InfVec{1}; FitDataYo = DragonVec{1};
    FishFit = polyfit(InfVec{j}, DragonVec{j}, 2);
    
    jYo = 0; FishyVal = polyval(FishFit, SumCounty(FishVec(j)-k:FishVec(j)+1));
    
%     if FishVec(j) < 4
%         %         CapParabola = -FishFit(2)/FishFit(1) - SumCounty(FishVec(j)-k);
%         CapParabola = SumCounty(FishVec(j)+1) - SumCounty(FishVec(j)-k);
%     else
%         CapParabola = SumCounty(FishVec(j)+1) - SumCounty(FishVec(j)-k);
%     end
    while FishyVal(end) > 1
        FishyVal = polyval(FishFit, SumCounty(FishVec(j)-k:FishVec(j)+jYo));
        jYo = jYo + 1;
    end
    
    %%%%% %%%%% %%%%%
    RootsParabola = roots(FishFit);
    % LogisticCap = (-FishFit(2)/FishFit(1) - 2*SumCounty(FishVec(j)-k));
    LogisticCap = abs(RootsParabola(1) - RootsParabola(2));
    
    %%%%% Determine the difference between prediction and reality %%%%%
    Reality = dIdtCheck(FishVec(j)-k:FishVec(j)+jYo-1);
    RealDiff = Reality - FishyVal;
    figure; plot(SumCounty(FishVec(j)-k:FishVec(j)+jYo-1), RealDiff, 'k.')
    
    %%%%% Use that difference to make a prediction %%%%%
    fprintf('Max of difference: %5.4f \n', max(RealDiff))
    
    %%
    
    %     while sum(OutYo) == 0
    if max(RealDiff) > 0.0
        FitTime = TimeVec{j}; FitData = DragonVec{j};
        
        %         [OutBest, yBest] = TimeCurveFit(FitTime, FitData);
        %%%%% Define the function for the fitting
        
        %         Funk = @(y) y(1)*exp(y(2)*(FitTime(1:tEnd) - y(3)).^2) - FitData(1:tEnd);
        Funk = @(y) y(1)*exp(y(2)*(FitTime - y(3)).^2) - FitData;
        
        %%%%% Define the lower and upper bounds
        LoB = [0, -1, FitTime(1)];
        UpB = [round(1.8*max(FitData)), 0, FitTime(end)];
        R2Best = 0.005;       %%%%% Initial R-Squared value for fitting
        clear OutBest;
        %%%%% Run fitting loop %%%%%
        OutYo = 0;
        while sum(OutYo) < 1
%             for kj = 1:3
                
                %%%%% Randomized initial condition %%%%%
                y0 = [randi([LoB(1), UpB(1)]), (randi([100*LoB(2), 100*UpB(2)]))/100,...
                    randi([LoB(3), UpB(3)])];
                
%                 y = lsqnonlin(Funk, y0, LoB, UpB);
                
                Range = 0.5;
                Lo = 1 - Range; Hi = 1 + Range;
                v1 = linspace(Lo*y0(1), Hi*y0(1), 5);
                v2 = linspace(Lo*y0(2), Hi*y0(2), 5);
                v3 = linspace(Lo*y0(3), Hi*y0(3), 5);
                
                %                 OutYo = 0;
                %                 while sum(OutYo) == 0
                for i = 1:numel(v1)
                    %     fprintf('Step 1 \n')
                    for ii = 1:numel(v2)
                        %         fprintf('Step 2 \n')
                        for iii = 1:numel(v3)
                            %             fprintf('Step 3 \n')
                            
                            y0 = [v1(i), v2(ii), v3(iii)]';
                            
                            y = lsqnonlin(Funk, y0, LoB, UpB);
                            
                            OutYo = y(1)*exp(y(2)*(FitTime - y(3)).^2);
                            
                            R2 = corrcoef(FitData, OutYo); R2D2 = R2(2,1)^2;
                            % fprintf('Made it here \n')
                            if R2D2 > R2Best && R2D2 < 1
                                R2Best = R2D2;
                                yBest(1) = y(1); yBest(2) = y(2); yBest(3) = y(3);
                                
                                OutBest = OutYo;
                                
                                fprintf('R2D2: %5.4f \n', R2D2)
                                %                             figure;
                                %                             plot(TimeDog{ij}, DataDog{ij}, 'k.', 'MarkerSize', 14); hold on;
                                %                             plot(FitTime, OutYo, 'bo'); hold off;
                                %                             pause
                            end % End if loop for tracking best fits
                            
                        end % End shift parameter loop
                    end % End exponential multiplication factor loop
                end % End pre-exponential factor loop
%             end
            clear y0;
        end % End loop that goes through the random exploration of parameter space
        kFit = kFit + 1;
    end
    jYoFish = 0; FishyTimeVal = OutBest;
    while numel(FishyTimeVal) < numel(FishyVal)
        FishyTimeVal = yBest(1)*exp(yBest(2)*...
            (timeYo(FishVec(j)-k:FishVec(j)+jYoFish) - yBest(3)).^2);
        
        %             polyval(FishFit, SumCounty(FishVec(j)-k:FishVec(j)+jYoFish));
        jYoFish = jYoFish + 1;
    end
    %
    XtrapTime = timeYo(FishVec(j)-k:FishVec(j)+jYoFish - 1);
    XtrapFit = FishyTimeVal;
%     if numel(FishyTimeVal) > numel(FitTime)
%         XtrapTime = timeYo(FishVec(j)-k:FishVec(j)+jYoFish - 1);
%         XtrapFit = FishyTimeVal;
%     else
%         XtrapTime = FitTime; XtrapFit = OutBest;
%     end
    
    
    %%%%%% Plot the data %%%%%%%%
    figure; plot(InfVec{j}, DragonVec{j}, 'k.'); hold on;
    %%%%%% Plot the polynomial fit of the data %%%%%
    plot(SumCounty(FishVec(j)-k:FishVec(j)+jYo-1), FishyVal, 'bo'); 
    %%%%%% Plot the Gaussian fit in Infected - Infection Rate Space %%%%%%%
    nInfGauss = cumtrapz(XtrapTime, XtrapFit);
    plot(nInfGauss+SumCounty(FishVec(j)-k), XtrapFit, 'rp')
    GaussPoly = polyfit(nInfGauss+SumCounty(FishVec(j)-k), XtrapFit, 2);
    GaussPolyVal = polyval(GaussPoly,...
        nInfGauss+SumCounty(FishVec(j)-k));
    %%%%%% Plot the parabolic fit of the Gaussian Outbreak %%%%%%%%%
    plot(nInfGauss+SumCounty(FishVec(j)-k), GaussPolyVal, 'c^'); 
    hold off;
    xlabel('People Infected'); ylabel('Daily Infection Rate');
    set(gca, 'fontsize', 14);
    TheLeg = legend('Data', 'Poly Fit', 'Gauss Fit', 'Gauss Poly Fit');
    set(TheLeg, 'Orientation', 'Vertical', 'box', 'off', 'Location', 'best')
    R2Polys = corrcoef(GaussPolyVal, FishyVal);
    R2D2Polys = R2Polys(2,1);
    if R2D2Polys < 0.93
        RootsGaussPara = roots(GaussPoly);
        LogisticCap = abs(RootsGaussPara(2) - RootsGaussPara(1));
        FishFit = GaussPoly;
        fprintf('Switched Poly Fit \n');
%         axis([RootsGaussPara(2) RootsGaussPara(1) 0 1.07*max(XtrapFit)])
    else
%         axis([RootsParabola(2) RootsParabola(1) 0 1.07*max(XtrapFit)])
    end    
        
    CapTime = nInfGauss(end); 
    Duration = abs(timeYo(FishVec(j)-k) - timeYo(FishVec(j)+jYoFish-1));
    fprintf('Capacity from Parabola: %i \n', round(LogisticCap))
    fprintf('Capacity from Gaussian: %i \n', round(CapTime))
    fprintf('Duration from Gaussian: %i \n', round(Duration))
    TimeFitCells{kFit} = OutBest;
    
    %%%%% Plot time fit results %%%%%
    figure; plot(FitTime-FitTime(1), OutBest, 'k.', 'MarkerSize', 14);
    hold on;
    plot(FitTime-FitTime(1), FitData, 'bo');
    hold off;
    R2C = corrcoef(FitData, OutBest); R2C2 = R2C(2,1)^2;
    fprintf('R-Squared for extrapolation and Data Set: %5.4f \n', R2C2)
    pause
    %%%%% End Plot time fit results section %%%%%
    
    for i = (FishVec(j)-k):numel(dIdt7)
        if i-(FishVec(j)-k)+1 <= numel(FishyVal)
            dIdtTimeCheck(i) = dIdtTimeCheck(i) - FishyTimeVal(i+1-(FishVec(j)-k));
        else
            dIdtTimeCheck(i) = dIdtTimeCheck(i);
        end
    end
    
    figure; plot(dIdtTimeCheck, 'k.')
    %
    %         %         [dNdtFit, NFit] = TimeCurveFunc(TimeVec{j}, InfVec{j});
    % end
    %     end
    
    %%
    %     pause
    
    % dIdtCheck = zeros(numel(dIdt7), 1);
    for i = (FishVec(j)-k):numel(dIdt7)
        if i-(FishVec(j)-k)+1 <= numel(FishyVal)
            dIdtCheck(i) = dIdtCheck(i) - FishyVal(i+1-(FishVec(j)-k));
        else
            dIdtCheck(i) = dIdtCheck(i);
        end
    end
    
    figure; plot(dIdtCheck, 'k.')
    
    KeepGoing = input('\nWant to keep going? 1: Yes 2: No. ');
end

ParaDataX = InfVec{j};
ParaDataY = DragonVec{j};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is the stuff I like to keep :-)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     if max(RealDiff) > 2.0
%         [~, MaxDex] = max(RealDiff);
%         FishyFit = polyfit(SumCounty(MaxDex-5:MaxDex), dIdtCheck(MaxDex-5:MaxDex), 2);
%         Cap = -FishyFit(2)/FishyFit(1);
%         fprintf('Number of infected: %5.1f \n', Cap);
%
%         %%%%% Also perform the fit in time with a Gaussian %%%%%
%         %%%%% Define the function for the fitting
%         Funk = @(y) y(1)*exp(y(2)*(TimeYo - y(3)).^2) - FitData;
%
%         %%%%% Define the lower and upper bounds
%         LoB = [0, -1, 0];
%         UpB = [800, 0, 120];
%
%         %%%%% Run fitting loop %%%%%
%         for k = 1:3
%
%             %%%%% Randomized initial condition %%%%%
%             y0 = [randi([LoB(1), UpB(1)]), (randi([100*LoB(2), 100*UpB(2)]))/100,...
%                 randi([LoB(3), UpB(3)])];
%
%             y = lsqnonlin(Funk, y0, LoB, UpB);
%
%             Range = 0.25;
%             Lo = 1 - Range; Hi = 1 + Range;
%             v1 = linspace(Lo*y(1), Hi*y(1), 3);
%             v2 = linspace(Lo*y(2), Hi*y(2), 3);
%             v3 = linspace(Lo*y(3), Hi*y(3), 3);
%
%             for i = 1:numel(v1)
%                 %     fprintf('Step 1 \n')
%                 for ii = 1:numel(v2)
%                     %         fprintf('Step 2 \n')
%                     for iii = 1:numel(v3)
%                         %             fprintf('Step 3 \n')
%
%                         y0 = [v1(i), v2(ii), v3(iii)]';
%
%                         y = lsqnonlin(Funk, y0, LoB, UpB);
%
%                         OutYo = y(1)*exp(y(2)*(FitTime - y(3)).^2);
%
%                         R2 = corrcoef(FitData, OutYo); R2D2 = R2(2,1)^2;
%                         % fprintf('Made it here \n')
%                         if R2D2 > R2Best && R2D2 < 1
%                             R2Best = R2D2;
%                             yBest(1) = y(1); yBest(2) = y(2); yBest(3) = y(3);
%
%                             OutBest = OutYo;
%
%                             fprintf('R2D2: %5.4f \n', R2D2)
%                             figure
%                             plot(TimeDog{ij}, DataDog{ij}, 'k.', 'MarkerSize', 14); hold on;
%                             plot(FitTime, OutYo, 'bo'); hold off;
%                             pause
%                         end % End if loop for tracking best fits
%
%                     end % End shift parameter loop
%                 end % End exponential multiplication factor loop
%             end % End pre-exponential factor loop
%         end % End loop that goes through the random exploration of parameter space
%         yKeep(ij,:) = yBest;
%         OutTwo = yBest(1)*exp(yBest(2)*(TimeDog{ij} - yBest(3)).^2);
%         figure;
%         plot(TimeDog{ij}, DataDog{ij}, 'k.', 'MarkerSize', 14); hold on;
%         plot(TimeDog{ij}, OutTwo, 'bo')
%         R2C = corrcoef(DataDog{ij}, OutTwo); R2C2 = R2C(2,1)^2;
%         fprintf('R-Squared for extrapolation and Data Set: %5.4f \n', R2C2)
%     end % End loop that goes through the number of outbreaks of a community




