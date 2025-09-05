function result = is_update(Objvalues,gen,Population,N,change_threshold,M)
    result = 0;
    [FrontNo,~]=NDSort(Population.objs,size(Population.objs,1));
    NC=size(find(FrontNo==1),2);
    max_change = abs(Objvalues(gen)-Objvalues(gen-1));
    if NC == N
        change_threshold = change_threshold * abs(((Objvalues(gen) / N))/(M))*10^(M-2);
        if max_change <= change_threshold
            result = 1;
        end
    end
end