function [data] = readExcel(fileName)

% Reading top row of excel data
excelData = readcell(fileName);

% Number of input cases
data.numCases = (size(excelData, 2) - 1) / 2;

% Determining selected case
for ii = 3:2:size(excelData, 2)
    bool = excelData{1, ii};
    if bool == 1
        data.selectedCase = (ii - 1) / 2;
        break
    end
end

% Getting value names from excelData
data.valNames = excelData(:, 1);

% Getting values from excelData
data.vals = excelData(:, data.selectedCase * 2);

% Getting units from excelData
data.units = excelData(:, (data.selectedCase * 2) + 1);


end
