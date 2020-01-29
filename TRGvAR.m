%% discovery-based program correlating gene expression to KM survival distinction
% inputs: RNASeq is the exome file from cBioPortal, CBIO contains the read
% groups depicted by KM analysis, polif_apop_genes and immune_genes
% contains lists of well-researched effector and immune genes for sample
% comparison
   
function TRGvAR(RNASeq,CBIO,prolif_apop_genes, immune_genes)
for i=1
    gothroughallgenesincolumn(RNASeq,i,CBIO,prolif_apop_genes,immune_genes)
end
end

function gothroughallgenesincolumn(RNASeqUPD,i,updCBIO, pa, immune)

for r= 2:20532
    gene= RNASeqUPD(r,1);
    sigtesting(gene,i,RNASeqUPD,updCBIO, pa, immune);
end
end

function sigtesting(gene,column, RNASeq, updCBIO, pa, immune)
 
e_row=RNASeq(:,1);e_col=RNASeq(1,:);
indapop=strfind(e_row,gene);
 
for i =1:numel(indapop)
    if indapop{i,1}==1
        genelocation=i;
    end
    if exist('genelocation') == 1
        break
    end
end

TRGcolval=[]; ARcolval=[];
z= updCBIO(:,column); 
TRG=[];AR=[];
 
if column==1 %always true//
    for o=1:numel(z)
        if o<20
            hstr= z(o,1);
%checking for proliferation genes
            fileind = strfind(e_col,hstr);
             for i=1:numel(fileind)
                if fileind{1,i}==1
                TRGcolval(end+1)=i;
                break
                end
                end
        end
                      
        if o>20
            lstr = z(o,1);
            fileind2=strfind(e_col,lstr);
            for i=1:numel(fileind2)
                if fileind2{1,i}==1
                    ARcolval(end+1)= i;
                    break
                end
            end
        end
    end
    
        for u=1:numel(ARcolval)
            TRG(end+1)= RNASeq(genelocation, ARcolval(u));
        end
        for l=1:numel(TRGcolval)
            AR(end+1)= RNASeq(genelocation, TRGcolval(l));
        end
    if isempty(TRG) || isempty(AR)
        disp('no rows returned')
    else
    [~,p] = ttest2(TRG,AR);
    end
end


if p<0.05

Attribute = 'TRGvAR';


%checking for proliferation genes
proliffiltered=pa(:,1);
apoptofiltered=pa(:,2);
cp=contains(proliffiltered,gene);ca=contains(apoptofiltered,gene);ci=contains(immune,gene);

if (sum(cp) == 1 && (p<0.05))
     PROLIFERATION='Y'; APOPTOSIS_EFF='N';IMMUNE='N';
%checking for apoptosis effector genes
elseif (sum(ca) == 1 && (p<0.05))
     APOPTOSIS_EFF='Y'; PROLIFERATION='N'; IMMUNE='N';
elseif (sum(ci) == 1 && (p<0.05))
    IMMUNE='Y'; APOPTOSIS_EFF='N'; PROLIFERATION='N';
elseif p<0.05
     PROLIFERATION='N'; APOPTOSIS_EFF='N'; IMMUNE='N';
end

%calculating TRG&AR averages
TRG_Avg = mean(TRG, 'all');
AR_Avg = mean(AR, 'all');

results = [gene, Attribute, p, TRG_Avg, AR_Avg, PROLIFERATION, APOPTOSIS_EFF, IMMUNE];
fileID = fopen('TRGvAR.txt','a+');
fprintf(fileID,'%s %s %f %f %f %s %s %s\n', results);
fclose(fileID);

end

end


