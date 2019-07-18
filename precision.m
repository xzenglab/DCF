function rate = precision(T,P,topR)
    % clear;
    % splits = [1 0 0 1; 0 0 0 1; 0 0 1 1; 0 1 1 0; 0 0 0 0];
    % ScoreMatrix = [8 3 1 6; 3 4 3 7; 1 2 5 1; 5 7 2 3; 2 1 1 4];
    [rank,index] = sort(P,'descend');
    for r = 1:topR  
        numSingl = 0 ;
         for i = 1:size(T,2);
            ind = find(T(:,i));
            if isempty(ind)
                numSingl = numSingl + 1;
                rate(i,r)= 0;
            else
                num = length(ind);
                count = 0;
                    for j = 1:num
                        for k = 1:r
                            if(index(k,i) == ind(j))
                                count = count + 1;
                            end
                        end
                    end
                    rate(i,r) = count / r;
            end
         end
    end

    rate = sum(rate) / (size(rate,1) - numSingl);
