for i = 1:size(g1,1)
    flag = True;
    for j = 1:size(g1,1)
        if abs(g1(i,j)) > 1e-6
            flag = False;
        end
    end
    if flag
        print(i, j);
    end
end