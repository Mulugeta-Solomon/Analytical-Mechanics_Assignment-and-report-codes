function index = nearest_index ( vec, value )
% searching a specified value in ascending ordered column vector
% index number of which element in 'vec' is nearest to 'value'
% 昇順に並んでいる列ベクトルの検索
% 昇順に並んでいる列ベクトルvecの要素でvalueに最も近い要素の添字
    sz = size(vec); n = sz(1);
    i = 1; j = n;
    
    if vec(i) >= value
        index = i;
        return;
    end
    if vec(j) <= value;
        index = j;
        return;
    end
    
    while (i < j-1)
        %fprintf("%d %d\n", i, j);
        k = floor((i+j)/2);
        if vec(k) <= value
            i = k;
        end
        if vec(k) >= value
            j = k;
        end
    end
    
    if value - vec(i) < vec(j) - value
        index = i;
    else
        index = j;
    end
end
