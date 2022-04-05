disp('fancy stuff');

% PoreC_Hypergraphs
% Inpect hypergraphs and generate data for Paohviz
clear
clc
close all


% USER SETTINGS
%typ = "scPore-C";      % option are "sc" (single-cell) or "pop" (population)
chr = [1 11];         
res = 1e6;       % resolution (bp)
%hoc = 2;         % number of higher order contacts
%runID = 11;
%use_filter=true; % filter contacts
%qc_filter=true;

load('n2s.mat');
chrom_sizes=readtable('GRCm39_Chrom_Sizes.txt','ReadVariableNames',false);
chrom_sizes.Properties.VariableNames={'Chrom','Size'};

if(res==1e6)
    resname="1Mb";
elseif(res==1e5)
    resname="100kb";
else
    resname=string(res);
end
load('unique_pairs_barcode11.mat');
load('Pairs_HiC.mat');

close all

%{
idx1=ismember(pairs.chrA,chr);
idx2=ismember(pairs.chrB,chr);
pairs=pairs(idx1&idx2,:);
pairs.posA=floor(pairs.posA./res)+1 ;
pairs.posB=floor(pairs.posB./res)+1;

chrom_sizes=chrom_sizes(ismember(chrom_sizes.Chrom,chr),:);
chrom_sizes.Size=floor(chrom_sizes.Size./res)+1;
ne=height(chrom_sizes);
chrom_sizes.Offset=zeros(ne,1);
cs=cumsum(chrom_sizes.Size);
chrom_sizes.Offset(2:end)=cs(1:end-1);

c2o=containers.Map(chrom_sizes.Chrom,chrom_sizes.Offset);
pairs.posA=pairs.posA+cell2mat(values(c2o,num2cell(pairs.chrA)));
pairs.posB=pairs.posB+cell2mat(values(c2o,num2cell(pairs.chrB)));

sp=sparse(pairs.posA,pairs.posB,ones(height(pairs),1),cs(end),cs(end));
f=full(sp);
%}
[f1,cs1]=fancy_stuff(uniquepairs,chr,chrom_sizes,res);
%[f2,cs2]=fancy_stuff(Pairs_HiC,chr,chrom_sizes,res);
%f=f+triu(f,1)';
f1=double(f1>0);
f1 = f1 + f1';
f1=double(f1>0);
f1 = triu(f1);
%f2=double(f2>0);
%f2=triu(f2)';
fcmb=f1+ f1';
diag_idx=logical(diag(ones(size(fcmb,1),1)));
fcmb(diag_idx)=fcmb(diag_idx)./2;
imagesc(fcmb);
axis image;

if(isequal(cs1,cs1))
    cs=cs1;
    for i=1:length(cs)
        xline(cs(i),'Color','red');
        yline(cs(i),'Color','red');
    end
else
    warn("Something went wrong");
    exit
end


ax=gca;
tickvalues=calculate_tick_location(cs);
ax.XTick=tickvalues;
ax.YTick=tickvalues;
%ax.XTick=[mean([0 196]) mean([196 378])];
%ax.YTick=[mean([0 196]) mean([196 378])];
ticklabels=values(n2s,num2cell(chr));
%ax.XTickLabel={'Chr1','Chr2'};
%ax.YTickLabel={'Chr1','Chr2'};
ax.XTickLabels=ticklabels;
ax.YTickLabels=ticklabels;

function [f,cs]=fancy_stuff(pairs,chr,chrom_sizes,res)
    idx1=ismember(pairs.chrA,chr);
    idx2=ismember(pairs.chrB,chr);
    pairs=pairs(idx1&idx2,:);
    pairs.posA=floor(pairs.posA./res)+1 ;
    pairs.posB=floor(pairs.posB./res)+1;

    chrom_sizes=chrom_sizes(ismember(chrom_sizes.Chrom,chr),:);
    chrom_sizes.Size=floor(chrom_sizes.Size./res)+1;
    ne=height(chrom_sizes);
    chrom_sizes.Offset=zeros(ne,1);
    cs=cumsum(chrom_sizes.Size);
    chrom_sizes.Offset(2:end)=cs(1:end-1);

    c2o=containers.Map(chrom_sizes.Chrom,chrom_sizes.Offset);
    pairs.posA=pairs.posA+cell2mat(values(c2o,num2cell(pairs.chrA)));
    pairs.posB=pairs.posB+cell2mat(values(c2o,num2cell(pairs.chrB)));

    sp=sparse(pairs.posA,pairs.posB,ones(height(pairs),1),cs(end),cs(end));
    f=full(sp);
end

function tickvalues=calculate_tick_location(cs)
   cs=[0;cs];
   ncs=length(cs);
   tickvalues=zeros(ncs-1,1);
   for i=1:length(tickvalues)
       tickvalues(i)=mean([cs(i) cs(i+1)]);
   end
end

