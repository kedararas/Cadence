function T = cadence_read_table(file, text_cols)
%CADENCE_READ_TABLE  readtable for the batch CSVs with the text columns kept as text.
%
%   T = cadence_read_table(file)              % every non-numeric-looking column as detected
%   T = cadence_read_table(file, text_cols)   % these columns are read as strings, always
%
%   readtable guesses each column's type from its contents.  A column of time
%   stamps ("2025-11-10 09:41:07") is turned into datetime, and one that is
%   mostly empty into double, which silently drops the few values it holds when
%   the table is written back.  Naming the text columns avoids both.

    if nargin < 2, text_cols = {}; end
    opts = detectImportOptions(char(file), 'Delimiter', ',', 'TextType', 'string');
    cols = intersect(cellstr(text_cols), opts.VariableNames, 'stable');
    if ~isempty(cols), opts = setvartype(opts, cols, 'string'); end
    T = readtable(char(file), opts);
    for k = 1:numel(cols)
        v = T.(cols{k});  v(ismissing(v)) = "";  T.(cols{k}) = v;
    end
end
