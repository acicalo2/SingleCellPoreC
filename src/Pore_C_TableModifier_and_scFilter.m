%% Cut the table down that comes from the snakemake output
FT =  FT(:,[1 3 4 5 13 14 15 16 17 18 19 20 28 29 30 31 32 33 34]);

%IR = FT(:,[1 5 6 7 8 11 12 13 14]);
clearvars -except all_cells

%% Rename Variables
%IR = renamevars(IR,["RunID", "align1_chrom", "align1_mapping_quality", "align2_chrom", "align2_mapping_quality"],...
%    ["read_id", "chrA", "mapA", "chrB", "mapB"]);
FT = renamevars(FT,["align1_mapping_quality", "align2_mapping_quality"...
   "align1_align_score", "align2_align_score", "align1_align_base_qscore" ,"align2_align_base_qscore","align1_chrom","align2_chrom"],...
   ["mapA", "mapB", "alignA", "alignB", "qscoreA", "qscoreB","chrA","chrB"]);
%% PosA
%start_end_frag_A = IR(:,[3 4]);
start_end_frag_A = FT(:,[7 8]);

posA = mean(table2array(start_end_frag_A),2);
%% PosB
start_end_frag_B = FT(:,[14 15]);
posB = mean(table2array(start_end_frag_B),2);

%% Remove the fragment positions
FT = removevars(FT,[7 8 14 15]);

%%% Append PosA and PosB
FT.posA = posA;
FT.posB = posB;

FT = movevars(FT,'posA','Before','chrA');
FT = movevars(FT,'posB','Before','chrB');


clearvars -except FT

%% Map Chromosome Number to ID
%chromosomes = readtable('C:\Pore-C\data_cleaning\Human Chromosomes');
addpath C:\Pore-C\GiantScript
addpath C:\Pore-C\data_cleaning
chromosomes = readtable('Human Chromosomes');
chromosomes = chromosomes(:,[1 3]);
%ES_Barcode11 = renamevars(ES_Barcode11,["align1_chrom", "align2_chrom"],...
%   ["chrA", "chrB"]);
c2n = containers.Map(chromosomes.Var3,chromosomes.Var1); 
FT = FT(ismember(FT.chrA,keys(c2n)),:);
FT = FT(ismember(FT.chrB,keys(c2n)),:);
FT.chrA = cell2mat(values(c2n,cellstr(FT.chrA)));
FT.chrB = cell2mat(values(c2n,cellstr(FT.chrB)));

clearvars -except Vari068
%% Single Cell Filter
FT = FT(FT.mapA >=30 & FT.mapB >=30,:);
FT = FT(FT.qscoreA >=9 & FT.qscoreB >=9,:);
FT = FT(FT.align1_frag_repeats<=4 & FT.align2_frag_repeats<=4,:);

%% Remove Vars
FT=removevars(FT,{'mapA','mapB','qscoreA','qscoreB','alignA','alignB'});                         % as to not upset anybodies hard coded indices
FT=removevars(FT,{'read_name','read_length'});
FT=removevars(FT,{'align1_align_idx','RunID','align1_strand','align1_frag_repeats','align2_frag_repeats','align2_strand'});

Vari068_pairs = FT(:,[2 5 4 8 7]);