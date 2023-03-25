
function [result] = concordanceindex(Pred, Time, Status)
n=length(Pred); count= 0; score=0;
for i = 1:(n-1)
    for j = (i+1):n 
        if ((Status(i) == 1) && (Time(j) > Time(i))) ||...
                ((Status(j) == 1) && (Time(i) > Time(j))) 
            count = count+1;
            if Pred(i) == Pred(j)
                score = score + 0.5;
            end
            if (Time(i) < Time(j)) && (Pred(i) > Pred(j))
                score = score + 1;
            end
            if Time(i) > Time(j) && (Pred(i) < Pred(j))
                score = score + 1;
            end
        end
    end
end
result=score/count;
end
