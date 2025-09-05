function [CS, UCS] = CACK(Problem, Con, Dec)
[mc,nc]  = size(Con);       
FR  = zeros(1, nc);
for i = 1:nc
    FR(i) = sum(Con(:,i) > 1e-6) / mc;      
end
%-------------------------------------------------------------------------%
R = zeros(1,Problem.D);
for i = 1:Problem.D
    CC = zeros(1,nc);
    for j = 1:nc
        [m,~] = find(Con(:,j) > 1e-6);
        if ~isempty(m)
           D1    = Dec(m,:);
           C1    = Con(m,:);
           CC(j) = corr(D1(:,i),C1(:,j),'Type','Pearson');
        else
           CC(j) = 0;
        end
    end
    R(i) = sum(FR.*CC);
end
weights = abs(R); 
CS = knee(weights);
UCS = setdiff(1:Problem.D, CS);  
%-------------------------------------------------------------------------%
P1            = Dec;
P1(:,UCS)     = 0;           
CVP1           = Problem.CalCon(P1);   
CV1           = sum(max(0,CVP1),2);
CI            = zeros(1,length(UCS));
for i  = 1:length(UCS)
    U         = UCS;
    P2        = Dec;
    U(:,i)    = [];
    P2(:,U)   = 0;
    CVP2      = Problem.CalCon(P2);   
    CV2       = sum(max(0,CVP2),2);
    CI(i)     = mean(CV2)-mean(CV1);
end
%-------------------------------------------------------------------------%
norCI = (CI - min(CI)) / (max(CI) - min(CI));
IS = knee(norCI);
CS=[CS,IS];
UCS = setdiff(1:Problem.D, CS);
end