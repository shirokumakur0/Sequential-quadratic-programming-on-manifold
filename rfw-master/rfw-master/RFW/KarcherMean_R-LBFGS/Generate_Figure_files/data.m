function data(Fig, position)

r = [37 9000 20 47 8500];
init_time = zeros(1,5);

for times = 1 : length(r)
    rand('state', r(times));
    randn('state', r(times));
    if(Fig == 1)
        k = 3;
        n = 3;
        f = 1;
    elseif(Fig == 2)
        if(strcmp(position, 'top'))
            k = 100;
            n = 3;
            f = 2;
        elseif(strcmp(position, 'bottom'))
            k = 100;
            n = 3;
            f = 5;
        end
    elseif(Fig == 3)
        if(strcmp(position, 'top'))
            k = 30;
            n = 100;
            f = 1;
        elseif(strcmp(position, 'bottom'))
            k = 30;
            n = 100;
            f = 5;
        end
    end

    
    [Kmean, As, cond_Kmean, cond_As] = constructed_data(k,n,f);
    
    tic;
    Ar = As{1};
    Ha = pinv(As{1});
    for i = 2:k
        Ar = Ar + As{i};
        Ha = Ha + pinv(As{i});
    end
    Ar = Ar/k;
    Ha = pinv(Ha/k);
    Arh = Ar^(0.5);
    Arinvh = Ar^(-0.5);
    temp = Arinvh * Ha * Arinvh;
    AH = Arh * temp^(0.5) * Arh;
    init_time(times) = toc;
    save(['Kmean',int2str(times)],'Kmean','As','cond_Kmean','cond_As','AH');
end

initial_time = mean(init_time);
save(['initial_time'],'initial_time')
end