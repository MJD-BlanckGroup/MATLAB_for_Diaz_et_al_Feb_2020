%% discovery-based program correlating gene expression to KM survival distinction
% inputs: RNASeq is the exome file from cBioPortal, CBIO contains the read
% groups depicted by KM analysis, polif_apop_genes and immune_genes
% contains lists of well-researched effector and immune genes for sample
% comparison

% the shell of the program -- for fitting many chem features in single file
function Watson(RNASeq,CBIO,prolif_apop_genes, immune_genes)
for i=1:4
    gothroughallgenesincolumn(RNASeq,i,CBIO,prolif_apop_genes,immune_genes)
end
end

% need to automatically index through thousands of genes
function gothroughallgenesincolumn(RNASeq,i,CBIO, pa, immune)

for r= 2:20532
    gene= RNASeq(r,1);
    sigtesting(gene,i,RNASeq,CBIO,pa,immune);
end
end

% need to compute whether group gene expression is significant
function sigtesting(gene,column, RNASeq, CBIO, pa, immune)

% "i" refers to the excel column, which refers to the CDR3 chemical feature
e_row=RNASeq(:,1);e_col=RNASeq(1,:);
indapop=strfind(e_row,gene);
 
for n =1:numel(indapop)
    if indapop{n,1}==1
        genelocation=n;
    end
    if exist('genelocation') == 1
        break
    end
end

% setting base parameters
hcolval=[]; lcolval=[];
z= CBIO(:,column); 
Upper_percentile=[];Lower_percentile=[];
 
if column==1
    for o=1:numel(z)
        if o<9
            hstr= z(o,1);

            fileind = strfind(e_col,hstr);
             for x=1:numel(fileind)
                if fileind{1,x}==1
                hcolval(end+1)=x;
                break
                end
                end
        end
                      
        if o>9
            lstr = z(o,1);
            fileind2=strfind(e_col,lstr);
            for x=1:numel(fileind2)
                if fileind2{1,x}==1
                    lcolval(end+1)= x;
                    break
                end
            end
        end
    end
    
        for u=1:numel(lcolval)
            Upper_percentile(end+1)= RNASeq(genelocation, lcolval(u));
        end
        for l=1:numel(hcolval)
            Lower_percentile(end+1)= RNASeq(genelocation, hcolval(l));
        end
    if isempty(Upper_percentile) || isempty(Lower_percentile)
        disp('no rows returned')
    else
    [~,p] = ttest2(Upper_percentile,Lower_percentile);
    end
end

if column==2
 
    for o=1:numel(z)
        if o<29
            hstr= z(o,1);
            fileind= strfind(e_col,hstr);
            for x =1:numel(fileind)
                if fileind{1,x}==1
                  hcolval(end+1)= x;
                  break  
                end
            end
        end
        if o>29
            lstr = z(o,1);
            fileind2=strfind(e_col,lstr);
            for x=1:numel(fileind2)
                if fileind2{1,x}==1
                    lcolval(end+1)= x;
                    break
                end
            end
        end
    end
    for u=1:numel(hcolval)
        Upper_percentile(end+1) = RNASeq(genelocation, hcolval(u));
    end
    for l=1:numel(lcolval)
        Lower_percentile(end+1) = RNASeq(genelocation, lcolval(l));
       
    end
    if isempty(Upper_percentile) || isempty(Lower_percentile)
        disp('no rows returned')
    else
    [~,p] = ttest2(Upper_percentile,Lower_percentile);
    end
end
 
if column==3
 
    for o=1:numel(z)
        if o<16
            hstr= z(o,1);
            fileind= strfind(e_col,hstr);
            for x =1:numel(fileind)
                if fileind{1,x}==1
                  hcolval(end+1)= x;
                  break
                end
            end
        end
        if o>16
            lstr = z(o,1);
            fileind2=strfind(e_col,lstr);
            for x=1:numel(fileind2)
                if fileind2{1,x}==1
                    lcolval(end+1)= x;
                    break
                end
            end
        end
    end
        for u=1:numel(hcolval)
        Upper_percentile(end+1) = RNASeq(genelocation, hcolval(u));
        end
        for l=1:numel(lcolval)
        Lower_percentile(end+1) = RNASeq(genelocation, lcolval(l));
        end
        if isempty(Upper_percentile) || isempty(Lower_percentile)
        disp('no rows returned')
        else
            [~,p] = ttest2(Upper_percentile,Lower_percentile);
        end
end
 
if column==4

    for o=1:numel(z)
        if o<17
            hstr= z(o,1);
            fileind= strfind(e_col,hstr);
            for x =1:numel(fileind)
                if fileind{1,x}==1
                  hcolval(end+1)= x;
                  break
                end
            end
        end
        if o>17
            lstr = z(o,1);
            fileind2=strfind(e_col,lstr);
            for x=1:numel(fileind2)
                if fileind2{1,x}==1
                    lcolval(end+1)= x;
                    break
                end
            end
        end
    end
for u=1:numel(hcolval)
    Upper_percentile(end+1) = RNASeq(genelocation, hcolval(u));
end

for l=1:numel(lcolval)
    Lower_percentile(end+1) = RNASeq(genelocation, lcolval(l));
end
if isempty(Upper_percentile) || isempty(Lower_percentile)
        disp('no rows returned')
else
    [~,p] = ttest2(Upper_percentile,Lower_percentile);
end
end

if p<0.05

% organizing our results by CDR3 chemical feature  
if column==1
    Attribute = 'fraction_charged_TRA_TRB';
elseif column==2 
    Attribute = 'isolectric_point_IGH_IGK';
elseif column==3
    Attribute = 'ncpr_IGH_IGL';
elseif column==4
    Attribute = 'isoelectric_point_IGH_IGL';
end

% indexing location of the prolif, apop, and immune genes
proliffiltered=pa(:,1);
apoptofiltered=pa(:,2);
cp=contains(proliffiltered,gene);ca=contains(apoptofiltered,gene);ci=contains(immune,gene);

% determining whether our gene falls into one of the above categories
if (sum(cp) == 1 && (p<0.05))
     PROLIFERATION='Y'; APOPTOSIS_EFF='N';IMMUNE='N';
elseif (sum(ca) == 1 && (p<0.05))
     APOPTOSIS_EFF='Y'; PROLIFERATION='N'; IMMUNE='N';
elseif (sum(ci) == 1 && (p<0.05))
    IMMUNE='Y'; APOPTOSIS_EFF='N'; PROLIFERATION='N';
elseif p<0.05
     PROLIFERATION='N'; APOPTOSIS_EFF='N'; IMMUNE='N';
end

%calculating top&bottom half averages
Top_50_Avg= mean(Upper_percentile, 'all'); %AVERAGE GENE EXPRESSION FOR TOP SURVIVING GROUP
Bot_50_Avg = mean(Lower_percentile, 'all'); %AVERAGE GENE EXPRESSION FOR BOTTOM SURVIVING GROUP

% collects esults into a text file which can then be imported into excel
results = [gene, Attribute, p, Top_50_Avg, Bot_50_Avg, PROLIFERATION, APOPTOSIS_EFF, IMMUNE];
fileID = fopen('COADresults.txt','a+');
fprintf(fileID,'%s %s %f %f %f %s %s %s\n', results);
fclose(fileID); %closes file after result has been reported

end

end