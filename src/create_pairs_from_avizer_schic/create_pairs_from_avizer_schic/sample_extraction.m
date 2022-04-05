%   Inputs
%   Decoder Table
%   'adj' file for single cell data      
%
%   Outputs
%   Pairs File From the Cells Analyzed in 
%   Nagano, Takashi, et al. "Single-cell Hi-C reveals cell-to-cell variability in 
%   chromosome structure." Nature 502.7469 (2013): 59-64. 


disp('sample_extraction');
clearvars -except decoder_table
clc
close all

if ~exist('decoder_table','var')
    load('decoder_table.mat');
end

% load an 'adj' file for single cell data and covert to pairs format
adj=readtable('adj','FileType','text','Delimiter','\t');
adj.fend1=adj.fend1+1;                                                     % one based indexing in matlab
adj.fend2=adj.fend2+1;                                                     % one based indexing in matlab
adj.chrA=decoder_table.chr(adj.fend1);
adj.posA=decoder_table.coord(adj.fend1);
adj.chrB=decoder_table.chr(adj.fend2);
adj.posB=decoder_table.coord(adj.fend2);
adj=removevars(adj,{'fend1','fend2','count'});
head(adj)