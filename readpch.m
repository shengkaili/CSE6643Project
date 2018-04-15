function [ K, F ] = readpch(filename)

fid = fopen(filename);

% dimension of problem?
n = 2;

tline = fgets(fid);

while ~isnumeric(tline)
   
    if strcmp(tline(1:5), 'DMIG*')
        % get row we are working on
        x = str2double(strrep(tline(33:40),' ', ''));
        y = str2double(strrep(tline(49:56),' ', ''));
        col = (x*n-1)+(y-1);
        % skip
    elseif strcmp(tline(1:4), 'DMIG')
        % initialize appropriate matrix on header
        flag = tline(9:12);
        if strcmp(flag, 'KAAX')
            m = str2double(strrep(tline(65:72),' ', ''));
            K = zeros(m,m);
        elseif strcmp(flag, 'PAX ')
            F = zeros(m,1);
        elseif strcmp(flag, 'VAX ')
            V = zeros(m,1);
        else
            keyboard
        end
    elseif strcmp(tline(1), '*')
        % get column and value of we are working on
        x = str2double(strrep(tline(17:24),' ', ''));
        y = str2double(strrep(tline(33:40),' ', ''));
        row = (x*n-1)+(y-1);
        value = str2double(strrep(strrep(tline(41:56),' ', ''),'D','E'));
        % assign to the appropriate matrix
        if strcmp(flag, 'KAAX')
            K(col,row) = value;
        elseif strcmp(flag, 'PAX ')
            F(row,1) = value;
        elseif strcmp(flag, 'VAX ')
            V(row,1) = value;
        else
            keyboard
        end
    else
        keyboard
    end
    
    % get next line
    tline = fgets(fid);
    
end

% remove entries that are not being used
rmv = ~V;
K(rmv,:) = [];
K(:,rmv) = [];
F(rmv,:) = [];

% reflect K matrix 
K = K + tril(K,-1)';