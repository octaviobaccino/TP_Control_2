function [VarName1, VarName2, VarName3, VarName4, VarName5] = importfile_MOTOR(workbookFile, sheetName, startRow, endRow)
%IMPORTFILE Import data from a spreadsheet
%  [VARNAME1, VARNAME2, VARNAME3, VARNAME4, VARNAME5] = IMPORTFILE(FILE)
%  reads data from the first worksheet in the Microsoft Excel
%  spreadsheet file named FILE.  Returns the data as column vectors.
%
%  [VARNAME1, VARNAME2, VARNAME3, VARNAME4, VARNAME5] = IMPORTFILE(FILE,
%  SHEET) reads from the specified worksheet.
%
%  [VARNAME1, VARNAME2, VARNAME3, VARNAME4, VARNAME5] = IMPORTFILE(FILE,
%  SHEET, STARTROW, ENDROW) reads from the specified worksheet for the
%  specified row interval(s). Specify STARTROW and ENDROW as a pair of
%  scalars or vectors of matching size for dis-contiguous row intervals.
%
%  Example:
%  [VarName1, VarName2, VarName3, VarName4, VarName5] = importfile("D:\Facultad Octavio\Control 2\Curvas_Medidas_Motor_2024.xls", "Hoja1", 1, 32370);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 04-Apr-2024 16:38:23

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 3
    startRow = 1;
    endRow = 32370;
end

%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + startRow(1) + ":E" + endRow(1);

% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Import the data
tbl = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:length(startRow)
    opts.DataRange = "A" + startRow(idx) + ":E" + endRow(idx);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    tbl = [tbl; tb]; %#ok<AGROW>
end

%% Convert to output type
VarName1 = tbl.VarName1;
VarName2 = tbl.VarName2;
VarName3 = tbl.VarName3;
VarName4 = tbl.VarName4;
VarName5 = tbl.VarName5;
end