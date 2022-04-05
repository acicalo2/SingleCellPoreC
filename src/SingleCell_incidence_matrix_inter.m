%% Load dataset and convert to matrix
load('ES_Barcode11_Pairs_reads.mat') % Fibroblast experiment V2 pairs table
binSize = 1e6;
ES_Barcode11_Pairs.posA = ceil(ES_Barcode11_Pairs.posA/binSize);
ES_Barcode11_Pairs.posB = ceil(ES_Barcode11_Pairs.posB/binSize);
ES_Barcode11_Pairs = sortrows(ES_Barcode11_Pairs, 'read_id', 'ascend');
porecMatrix = table2array(ES_Barcode11_Pairs);

%% Add chromosome offsets
chrLength = [196, 182, 160, 157, 152, 150, 145, 131, 125, 131, 122,...
    121, 121, 126, 105, 99, 96, 91, 62, 170, 92];
for chr = 2:21
    binAdd = sum(chrLength(1:chr-1));
    binIdxL = porecMatrix(:, 2)==chr;
    porecMatrix(binIdxL, 3) = porecMatrix(binIdxL, 3)+binAdd; 
    binIdxR = porecMatrix(:, 4)==chr;
    porecMatrix(binIdxR, 5) = porecMatrix(binIdxR, 5)+binAdd;
end

%% Threshold pairwise contacts removing noise
%{
[~, uniqReadContacts, ~] = unique(porecMatrix(:, [1 3 5]), 'rows');
porecMatrix = porecMatrix(uniqReadContacts, :);
pairedContacts = porecMatrix(:, [3 5]);
[~, ~, contactIdx] = unique(pairedContacts, 'rows', 'stable');  
contactHist = accumarray(contactIdx, 1);
contactCount = contactHist(contactIdx);
eps = prctile(contactCount, 85); % 85th percentile 
porecMatrix = porecMatrix(contactCount>=eps, :);
%}

%% Extract two-way inter-chromosome contacts
[~, ~, readIDReIdx] = unique(porecMatrix(:, 1)); % Reindex the ReadID
porecMatrix(:, 1) = readIDReIdx;
allChrContacts = cell(readIDReIdx(end), 1);
for i = 1:readIDReIdx(end)
    chrContact = unique(porecMatrix(porecMatrix(:, 1)==i, [2 4]));
    allChrContacts{i} = chrContact(:)';
end
twoWayInter = [1 11];
interChrLength = find((cellfun(@(v) all(ismember(v, twoWayInter)) & length(v)==length(twoWayInter), allChrContacts))>0);
interChrIntersect = intersect(porecMatrix(:, 1), interChrLength);
porecMatrixInterChr = porecMatrix(ismember(porecMatrix(:, 1), interChrIntersect), :);

%% Create an incidence matrix
[~, ~, readIDReIdx] = unique(porecMatrixInterChr(:, 1)); % Reindex the ReadID
porecMatrixInterChr(:, 1) = readIDReIdx;
incidenceMatrixA = sparse(porecMatrixInterChr(:, 3), porecMatrixInterChr(:, 1), 1,...
    sum(chrLength), porecMatrixInterChr(end, 1)); 
incidenceMatrixB = sparse(porecMatrixInterChr(:, 5), porecMatrixInterChr(:, 1), 1,...
    sum(chrLength), porecMatrixInterChr(end, 1));
incidenceMatrix = (incidenceMatrixA+incidenceMatrixB)>0;
incidenceMatrix = unique(incidenceMatrix', 'rows')'; % Remove duplicate edge

%% Extract all higer-order hyperedges
highOrderContacts = cell(size(incidenceMatrix, 2), 1);
for i = 1:size(incidenceMatrix, 2)
    highOrderContacts{i} = find(incidenceMatrix(:, i)>0);    
end
highOrderContacts = flipud(highOrderContacts(cellfun('length', highOrderContacts)>=3)); % Order >=3

%% Create CSV file for PAOHvis
CSVIncidenceMatrix = [];
label = [];
for i = 1:length(highOrderContacts)
    contact = highOrderContacts{i};
    for j = 1:length(contact)
        CSVIncidenceMatrix = [CSVIncidenceMatrix; [i contact(j)]]; %#ok<AGROW>
        if contact(j) >= 1 && contact(j) <= 196
            label{end+1} = 'Chr 1'; %#ok<*SAGROW>
        elseif contact(j) >= 1530 && contact(j) <= 1651
            label{end+1} = 'Chr 11';
        else 
            label{end+1} = 'Other';
        end
    end 
end 
        
CSVIncidenceMatrixTable = table(CSVIncidenceMatrix(:, 1), num2str(CSVIncidenceMatrix(:, 2), '%04d'), ...
    cell(size(CSVIncidenceMatrix, 1), 1), cell(size(CSVIncidenceMatrix, 1), 1), label');
writetable(CSVIncidenceMatrixTable, 'Inter_Chr_1_11.csv')

