tmp = g0.';
for i = 1:size(g1,1)
    flag = true;
    for j = 1:size(g1,1)
        if abs(tmp(i,j)) > 1e-6
            flag = false;
        end
    end
    if flag
        i
    end
end