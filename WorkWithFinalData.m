clear; clc
%% Get the data from excel file
for GetDataFromExcel=1:1
% Specify the file path to your Excel file
filePath = 'C:\Users\netas\Documents\Lammel lab\Obesity paper\Nature 12 2023\Revision March 2024\All resubmission materials\Statistics info\combinedStatisticsTable7.xlsx';
% Read the data from the Excel file into a table
dataTable = readtable(filePath, 'Sheet', 'Final Data');
% Clean the table:  Loop through each variable in the table and remove '''
% and ' from string variables, convert all to lowercase
for i = 1:width(dataTable)
    if iscellstr(dataTable{:, i}) % Check if the column contains strings
        % Remove all instances of both '''' and single quote '
        dataTable{:, i} = erase(dataTable{:, i}, '''');
        dataTable{:, i} = erase(dataTable{:, i}, '"');
                % Convert to lowercase
        dataTable{:, i} = lower(dataTable{:, i});
    end
end
% get the values you want to analyse
IsResponding = dataTable(:,18);
IsTagged = dataTable(:,30);
DietType = dataTable(:,3);
ClassType = dataTable(:,25);
UnitNumber= dataTable(:,1);
Condition = dataTable(:,8);
MouseName = dataTable(:,4);
TimeSpent= dataTable(:,34);
MakeName=table2cell(dataTable(:,4:7));
MakeName=string(MakeName);
UnitName = join(MakeName, ' ', 2); % Concatenates along columns with a space separator
NameAndNumber.UnitName=UnitName;
NameAndNumber.UnitNumber=UnitNumber;
%% Get information about cell classification


%% hendle the variables
% Convert condition to numeric index for easier handling
ConditionIndex = nan(height(Condition),1);
for i = 1:height(Condition)
    switch char(Condition{i,1})
        case 'chow'
            ConditionIndex(i) = 1;
        case 'empty'
            ConditionIndex(i) = 2;
        case 'jelly'
            ConditionIndex(i) = 3;
        case 'laser'
            ConditionIndex(i) = 4;
        case 'piezo_chow'
            ConditionIndex(i) = 5;
        case 'piezo_jelly'
            ConditionIndex(i) = 6;
        case 'rearing'
            ConditionIndex(i) = 7;
        case 'running'
            ConditionIndex(i) = 8;
        case 'stopping'
            ConditionIndex(i) = 9;
        case 'trotting'
            ConditionIndex(i) = 10;
        case 'turnleft'
            ConditionIndex(i) = 11;
        case 'turnright'
            ConditionIndex(i) = 12;
        case 'walking'
            ConditionIndex(i) = 13;
    end
end
% Convert variables to logical
IsReg = strcmp(DietType{:,:}, 'reg');
IsTagged = logical(table2array(IsTagged));
IsResponding = logical(table2array(IsResponding));
% Convert columns to arrays if necessary
ClassType = double(table2array(ClassType));
ClassType(isnan(ClassType))=0;
% Get unique values for each criterion
uniqueDiets = unique(IsReg);
uniqueTypes = unique(ClassType);
uniqueTypes(uniqueTypes==0)= [];
uniqueConditionIndex = unique(ConditionIndex);
uniqueDecisions = unique(IsResponding);
uniqueConditionNames=unique(Condition);
uniqueTagged=unique(IsTagged);
% make nan array to store the data
RespondingClass=nan(size(IsTagged));
for k=1:length(RespondingClass)
    TempResponse=IsResponding(k);
    TempClass=ClassType(k);
    if TempResponse % statistically significant
        if TempClass==1
            RespondingClass(k)=1;
        elseif TempClass==2
        RespondingClass(k)=-1;
        else
       RespondingClass(k)=0;
        end
    else % Not statistically significant
               RespondingClass(k)=0;
    end
end
   


end
%% fix issues
for FindIssues=1:1
FindMissingCondition=cell(height(Condition),3);
FindMissingCondition(:,1)=num2cell(ConditionIndex);
FindMissingCondition(:,2)=Condition.Var8;
FindMissingCondition(:,3)=num2cell(UnitNumber.Var1);

% Define the required set of conditions
requiredConditions = {'chow', 'empty', 'jelly', 'laser', 'piezo_chow', ...
    'piezo_jelly', 'rearing', 'running', 'stopping', 'trotting', ...
    'turnleft', 'turnright', 'walking'};
% Convert cell array to table for easier manipulation
FindMissingTable = cell2table(FindMissingCondition, 'VariableNames', {'ID', 'Condition', 'Unit'});
% Get unique unit numbers
uniqueUnits = unique(FindMissingTable.Unit);
% Initialize results to store missing and duplicate conditions for each unit
missingResults = [];
duplicateResults = [];
% Loop through each unit and find missing and duplicate conditions
for i = 1:length(uniqueUnits)
    currentUnit = uniqueUnits(i);
 
    % Extract conditions for the current unit
    unitConditions = FindMissingTable.Condition(FindMissingTable.Unit == currentUnit);
    
    % Find missing conditions
    missingConditions = setdiff(requiredConditions, unitConditions);
    
    % Store missing conditions for this unit
    for j = 1:length(missingConditions)
        missingResults = [missingResults; {currentUnit, missingConditions{j}}];
    end
    
    % Find duplicate conditions by counting occurrences
    [uniqueConditions, ~, conditionIndices] = unique(unitConditions);
    conditionCounts = histc(conditionIndices, 1:numel(uniqueConditions));
    duplicateConditions = uniqueConditions(conditionCounts > 1);
    
    % Store duplicate conditions for this unit
    for j = 1:length(duplicateConditions)
        duplicateResults = [duplicateResults; {currentUnit, duplicateConditions{j}}];
    end
end
% % Convert results to tables
% missingConditionsTable = cell2table(missingResults, 'VariableNames', {'Unit', 'MissingCondition'});
% duplicateConditionsTable = cell2table(duplicateResults, 'VariableNames', {'Unit', 'DuplicateCondition'});
% % Display the tables
% disp('Missing Conditions:');
% disp(missingConditionsTable);
% disp('Duplicate Conditions:');
% disp(duplicateConditionsTable);
end
%% Task 1: Collect counts based on conditions and perform chi-square tests
for ChiTestDiet=1:1
results = [];
for cond = 1:length(uniqueConditionIndex)
    for typ = 1:length(uniqueTypes)
        for Tag=1:length(uniqueTagged)
%             disp([uniqueConditionNames(1,cond),uniqueConditionIndex(cond)]);
        % Select rows based on the current condition and type
        reg_Responding = sum((IsReg == 1) & (ClassType == uniqueTypes(typ)) & (ConditionIndex == uniqueConditionIndex(cond))& (IsTagged == uniqueTagged(Tag)) & (IsResponding == 1));
        hfd_Responding = sum((IsReg == 0) & (ClassType == uniqueTypes(typ)) & (ConditionIndex == uniqueConditionIndex(cond))& (IsTagged == uniqueTagged(Tag))  & (IsResponding == 1));
        reg_non_Responding = sum((IsReg == 1) & (ConditionIndex == uniqueConditionIndex(cond))& (IsTagged == uniqueTagged(Tag))  & (IsResponding == 0));
        hfd_non_Responding = sum((IsReg == 0) & (ConditionIndex == uniqueConditionIndex(cond))& (IsTagged == uniqueTagged(Tag))  & (IsResponding == 0));
regPrecent=100*(reg_Responding/(reg_Responding+reg_non_Responding));
hfdPrecent=100*(hfd_Responding/(hfd_Responding+hfd_non_Responding));
regTotals=(reg_Responding+reg_non_Responding);
hfdTotals=(hfd_Responding+hfd_non_Responding);

        % Construct observed table
        observed = [reg_Responding, hfd_Responding; reg_non_Responding, hfd_non_Responding];
        
        % Calculate expected frequencies
        total = sum(observed(:));
        row_totals = sum(observed, 2);
        col_totals = sum(observed, 1);
        expected = (row_totals * col_totals) / total;

        % Perform chi-square test manually
        if all(expected(:) > 0)
            chi2_stat = sum((observed(:) - expected(:)).^2 ./ expected(:));
            p = 1 - chi2cdf(chi2_stat, 1);  % Degrees of freedom = 1 for 2x2 table
        else
            chi2_stat = NaN;
            p = NaN;
        end
        
        % Determine significance level
        if p < 0.001
            decision = '***';
        elseif p < 0.01
            decision = '**';
        elseif p < 0.05
            decision = '*';
        else
            decision = 'ns';
        end
        
        % Store results for each condition and type
        results = [results; {uniqueConditionIndex(cond), uniqueTypes(typ), uniqueTagged(Tag),reg_non_Responding,reg_Responding,regPrecent, regTotals, hfd_non_Responding, hfd_Responding,hfdPrecent, hfdTotals,chi2_stat, p, decision}];
        end
    end
end
% Convert results to a table
finalResultsTable = cell2table(results, 'VariableNames', {'Condition', 'Type','Tagged', 'REG-', 'REG+','REG%','REG total','HFD-','HFD+','HFD%','HFD total', 'Chi_Square_Statistic', 'p_value', 'Decision'});
% Replace numbers in the 'condition' column with corresponding names
for i = 1:length(uniqueConditionIndex)
    conditionNum = uniqueConditionIndex(i);
    conditionName = char(uniqueConditionNames{i,1}); % Get the corresponding name
    % Find rows with the current condition number and replace with the name
    finalResultsTable.ConditionName(finalResultsTable.Condition == conditionNum) = {conditionName};
end
% Save to an Excel file (optional)
writetable(finalResultsTable, 'C:\Users\netas\Documents\Lammel lab\Obesity paper\Nature 12 2023\Revision March 2024\All resubmission materials\Statistics info\finalResultsTable5.xlsx', 'Sheet', 'ProportionsChiTest');




end
%% Task 2: Collect counts based on conditions and perform chi-square tests
for ChiTestPiezo=1:1
results = [];
% find the piezo jelly condition
% Find the locations of 'piezo_chow' and 'piezo_jelly'
uniqueConditionNames=table2array(uniqueConditionNames);
location_piezo_chow = find(strcmp(uniqueConditionNames, 'piezo_chow'));
location_piezo_jelly = find(strcmp(uniqueConditionNames, 'piezo_jelly'));
uniquePiezo=[location_piezo_chow,location_piezo_jelly];
for Diets = 1:length(uniqueDiets)
    for typ = 1:length(uniqueTypes)
        for Tag=1:length(uniqueTagged)
        % Select rows based on the current condition and type
        jelly_Responding = sum((ConditionIndex == uniqueConditionIndex(location_piezo_jelly)) & (IsResponding == 1)&(IsReg == uniqueDiets(Diets))& (IsTagged == uniqueTagged(Tag)) & (ClassType == uniqueTypes(typ)));
        chow_Responding = sum((ConditionIndex == uniqueConditionIndex(location_piezo_chow)) & (IsResponding == 1)&(IsReg == uniqueDiets(Diets)) & (IsTagged == uniqueTagged(Tag))  & (ClassType == uniqueTypes(typ)));
        jelly_non_Responding = sum((ConditionIndex == uniqueConditionIndex(location_piezo_jelly)) & (IsResponding == 0)&(IsReg == uniqueDiets(Diets))& (IsTagged == uniqueTagged(Tag))  & (ClassType == uniqueTypes(typ)));
        chow_non_Responding = sum((ConditionIndex == uniqueConditionIndex(location_piezo_chow)) & (IsResponding == 0)&(IsReg == uniqueDiets(Diets))& (IsTagged == uniqueTagged(Tag))  & (ClassType == uniqueTypes(typ)));
    jellyPrecent=100*(jelly_Responding/(jelly_Responding+jelly_non_Responding));
    chowPrecent=100*(chow_Responding/(chow_Responding+chow_non_Responding));
        % Construct observed table
        observed = [jelly_Responding, chow_Responding; jelly_non_Responding, chow_non_Responding];
        % Calculate expected frequencies
        total = sum(observed(:));
        row_totals = sum(observed, 2);
        col_totals = sum(observed, 1);
        expected = (row_totals * col_totals) / total;
        % Perform chi-square test manually
        if all(expected(:) > 0)
            chi2_stat = sum((observed(:) - expected(:)).^2 ./ expected(:));
            p = 1 - chi2cdf(chi2_stat, 1);  % Degrees of freedom = 1 for 2x2 table
        else
            chi2_stat = NaN;
            p = NaN;
        end       
        % Determine significance level
        if p < 0.001
            decision = '***';
        elseif p < 0.01
            decision = '**';
        elseif p < 0.05
            decision = '*';
        else
            decision = 'ns';
        end        
         % Store results for each condition and type
        results = [results; {uniqueDiets(Diets), uniqueTypes(typ), uniqueTagged(Tag),jelly_non_Responding,jelly_Responding,jellyPrecent,chow_non_Responding,chow_Responding,chowPrecent, chi2_stat, p, decision}];
        end
    end
end
% Convert results to a table
finalResultsTable = cell2table(results, 'VariableNames', {'IsReg', 'Type','Tagged', 'Jelly-', 'Jelly+','Jelly%','Chow-','Chow+','Chow%', 'Chi_Square_Statistic', 'p_value', 'Decision'});
% Save to an Excel file (optional)
writetable(finalResultsTable, 'C:\Users\netas\Documents\Lammel lab\Obesity paper\Nature 12 2023\Revision March 2024\All resubmission materials\Statistics info\finalResultsTable.xlsx', 'Sheet', 'PiezoChiTestbyDiet');
end
%% Task 3: Collect counts based on conditions and perform chi-square tests
for ChiTestDLC=1:1 % revise
results = [];
% find the piezo jelly condition
% Find the locations of 'piezo_chow' and 'piezo_jelly'
try uniqueConditionNames=table2array(uniqueConditionNames);catch;end
location_piezo_chow = find(strcmp(uniqueConditionNames, 'piezo_chow'));
location_piezo_jelly = find(strcmp(uniqueConditionNames, 'piezo_jelly'));
uniquePiezo=[location_piezo_chow,location_piezo_jelly];
for Diets = 1:length(uniqueDiets)
    for typ = 1:length(uniqueTypes)
        for Tag=1:length(uniqueTagged)
        % Select rows based on the current condition and type
        jelly_Responding = sum((ConditionIndex == uniqueConditionIndex(location_piezo_jelly)) & (IsResponding == 1)&(IsReg == uniqueDiets(Diets))& (IsTagged == uniqueTagged(Tag)) & (ClassType == uniqueTypes(typ)));
        chow_Responding = sum((ConditionIndex == uniqueConditionIndex(location_piezo_chow)) & (IsResponding == 1)&(IsReg == uniqueDiets(Diets)) & (IsTagged == uniqueTagged(Tag))  & (ClassType == uniqueTypes(typ)));
        jelly_non_Responding = sum((ConditionIndex == uniqueConditionIndex(location_piezo_jelly)) & (IsResponding == 0)&(IsReg == uniqueDiets(Diets))& (IsTagged == uniqueTagged(Tag))  & (ClassType == uniqueTypes(typ)));
        chow_non_Responding = sum((ConditionIndex == uniqueConditionIndex(location_piezo_chow)) & (IsResponding == 0)&(IsReg == uniqueDiets(Diets))& (IsTagged == uniqueTagged(Tag))  & (ClassType == uniqueTypes(typ)));
    jellyPrecent=100*(jelly_Responding/(jelly_Responding+jelly_non_Responding));
    chowPrecent=100*(chow_Responding/(chow_Responding+chow_non_Responding));
        % Construct observed table
        observed = [jelly_Responding, chow_Responding; jelly_non_Responding, chow_non_Responding];
        % Calculate expected frequencies
        total = sum(observed(:));
        row_totals = sum(observed, 2);
        col_totals = sum(observed, 1);
        expected = (row_totals * col_totals) / total;
        % Perform chi-square test manually
        if all(expected(:) > 0)
            chi2_stat = sum((observed(:) - expected(:)).^2 ./ expected(:));
            p = 1 - chi2cdf(chi2_stat, 1);  % Degrees of freedom = 1 for 2x2 table
        else
            chi2_stat = NaN;
            p = NaN;
        end       
        % Determine significance level
        if p < 0.001
            decision = '***';
        elseif p < 0.01
            decision = '**';
        elseif p < 0.05
            decision = '*';
        else
            decision = 'ns';
        end        
         % Store results for each condition and type
        results = [results; {uniqueDiets(Diets), uniqueTypes(typ), uniqueTagged(Tag),jelly_non_Responding,jelly_Responding,jellyPrecent,chow_non_Responding,chow_Responding,chowPrecent, chi2_stat, p, decision}];
        end
    end
end
% Convert results to a table
finalResultsTable = cell2table(results, 'VariableNames', {'IsReg', 'Type','Tagged', 'Jelly-', 'Jelly+','Jelly%','Chow-','Chow+','Chow%', 'Chi_Square_Statistic', 'p_value', 'Decision'});
% Save to an Excel file (optional)
writetable(finalResultsTable, 'C:\Users\netas\Documents\Lammel lab\Obesity paper\Nature 12 2023\Revision March 2024\All resubmission materials\Statistics info\finalResultsTable.xlsx', 'Sheet', 'PiezoChiTestbyDiet');
end
%% Task 4: make a heatmap
for MakeHeatmap = 1:1
    % Prepare data table with relevant columns and additional variables
    Data2Plot = cell(height(dataTable), 4);
    Data2Plot(:, 1) = num2cell(dataTable.Var1); % Unit number
    Data2Plot(:, 2) = dataTable.Var8;           % Condition names
    Data2Plot(:, 3) = num2cell(ConditionIndex);
    Data2Plot(:, 4) = num2cell(dataTable.Var15);%DATA TO PLOT in heatmap!!
    Data2Plot(:, 5) = num2cell(dataTable.Var34);

    Data2PlotTable = cell2table(Data2Plot, 'VariableNames', ...
        {'Unit', 'Condition', 'ConditionIndex', 'DataValue','TimeSpent'});

    % Add IsReg, IsTagged, ClassType, and UnitName to the table
    Data2PlotTable.IsReg = IsReg;       % Logical variable (true/false) indicating REG or non-REG
    Data2PlotTable.IsTagged = IsTagged; % Logical variable (true/false) indicating Tagged or non-Tagged
    Data2PlotTable.ClassType = ClassType; % Integer values (0, 1, 2) indicating class type
    Data2PlotTable.UnitName = UnitName;   % UnitName for each unit

    % Define the new condition order
    conditions = {'jelly', 'chow', 'empty', 'rearing', 'stopping', ...
        'walking', 'trotting', 'running', 'turnleft', 'turnright', ...
        'laser', 'piezo_chow', 'piezo_jelly'};

    % Define combinations of IsReg and IsTagged
    combinations = {
        true,  true,  'REG+ Tagged+';
        false, true,  'REG- Tagged+';
        true,  false, 'REG+ Tagged-';
        false, false, 'REG- Tagged-'
    };

    % Initialize output cell array to store data for each heatmap
    heatmapOutputs = cell(size(combinations, 1), 5); % Adding a 5th column for UnitName

    % Loop through each condition combination and generate a heatmap
    for i = 1:size(combinations, 1)
        % Filter data based on IsReg and IsTagged
        isReg = combinations{i, 1};
        isTagged = combinations{i, 2};
        subset = Data2PlotTable(Data2PlotTable.IsReg == isReg & Data2PlotTable.IsTagged == isTagged, :);
        
        % Get "jelly" condition rows and sort by ClassType in 2, 0, 1 order
        jellySubset = subset(strcmp(subset.Condition, 'jelly'), :);
        [~, sortIdx] = sort(jellySubset.ClassType, 'ascend'); % Sort by ClassType

        % Extract sorted unit list and UnitName, keeping only unique entries
        sortedUnits = unique(jellySubset.Unit(sortIdx), 'stable');      % Unique sorted units for the "jelly" condition
        sortedUnitNames = unique(jellySubset.UnitName(sortIdx), 'stable'); % Unique sorted unit names
% 
        % Filter subset to include only rows with these unique sorted units
        subset = subset(ismember(subset.Unit, sortedUnits), :);

        % Initialize heatmap data matrix
        heatmapData = NaN(length(conditions), length(sortedUnits));

        % Populate heatmap data matrix
        for j = 1:height(subset)
            unitIdx = find(sortedUnits == subset.Unit(j));
            conditionIdx = find(strcmp(conditions, subset.Condition{j}));
            heatmapData(conditionIdx, unitIdx) = subset.DataValue(j);
        end

        % Store the heatmap data, units, conditions, ClassType information, and UnitName
        heatmapOutputs{i, 1} = heatmapData;
        heatmapOutputs{i, 2} = sortedUnits;
        heatmapOutputs{i, 3} = conditions;
        heatmapOutputs{i, 4} = jellySubset.ClassType(sortIdx); % Store sorted ClassType info for "jelly" condition units
% Find indices in Data2PlotTable for the sorted unique units
[~, unitIndices] = ismember(sortedUnits, Data2PlotTable.Unit);
% Retrieve UnitName in the sorted order and store in column 6 of heatmapOutputs
heatmapOutputs{i, 5} = Data2PlotTable.UnitName(unitIndices);
heatmapOutputs{i, 6} = Data2PlotTable.TimeSpent(unitIndices);


        % Create heatmap using unit indices as labels on the x-axis
        figure;
        unitIndices = 1:length(sortedUnits); % Using indices instead of names for plotting
        h = heatmap(unitIndices, conditions, heatmapData, 'MissingDataColor', [0.8 0.8 0.8]);
        h.Title = sprintf('Heatmap for %s', combinations{i, 3});
        h.XLabel = 'Unit Index (sorted by ClassType for "jelly")';
        h.YLabel = 'Condition';

        % Apply custom colormap and color range
        colormap(h, 'jet');
h.ColorLimits = [-1, 1]; % Set the color range from -0.5 to 5
        colorbar; % Add colorbar for reference
    end

    % Display or access heatmap data as needed
    % Each row of `heatmapOutputs` now contains:
    % - heatmapOutputs{i, 1}: Data matrix for the heatmap
    % - heatmapOutputs{i, 2}: Units on the x-axis, sorted by ClassType for the "jelly" condition
    % - heatmapOutputs{i, 3}: Conditions on the y-axis in the specified order
    % - heatmapOutputs{i, 4}: ClassType information for sorted units based on the "jelly" condition
    % - heatmapOutputs{i, 5}: UnitName array for each unit (stored, not displayed)
end



