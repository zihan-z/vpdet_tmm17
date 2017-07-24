function tr = pick (numb, k, par)
    ind = 0;
    %n = 20;
    n = max(size (numb));
    ind = 1;
    tr (1, 1:3) = [1 2 3];
    next = 1;
    prev = tr (1, :);
    while (next)
        [trip, next] = plus1 (prev, n);
        if (trip (1) ~= trip(2) && trip (2) ~= trip(3) && trip (1) ~= trip(3))
            ind = ind + 1;
            tr (ind, 1:3) = trip (:);
        end
        prev = trip;
    end
end
    

function [trip, next] = plus1 (tr, n)
    p = 1;
    i = 3;
    while (p && i > 0)
        if (tr (i) < n)
            tr(i) = tr(i) + 1;
            p = 0;
        else
            tr(i) = 1;
            i = i - 1;
        end
    end
    next = (i ~=0);
    trip = tr;
end
       