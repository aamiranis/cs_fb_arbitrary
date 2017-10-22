function beta = image_downsampling_fn(s_im,opt)
E = [];
O = [];
switch opt
    case 'rectangle'
        for j = 1:2:s_im
            E = [E ((j-1)*s_im + (1:2:s_im))];
            O = [O ((j-1)*s_im + (2:2:s_im))];
        end
        
        for j = 2:2:s_im
            E = [E ((j-1)*s_im + (2:2:s_im))];
            O = [O ((j-1)*s_im + (1:2:s_im))];
        end
    case 'diamond'
        for j = 1:2:s_im
            index = max(((j-1)*s_im + (1:1:s_im)));
            if index <= s_im*s_im
                E = [E ((j-1)*s_im + (1:1:s_im))];
            end
            index = max((j*s_im + (1:1:s_im)));
            if index <= s_im*s_im
                O = [O (j*s_im + (1:1:s_im))];
            end
        end
    case 'vertical'
        for j = 1:s_im
            E = [E ((j-1)*s_im + (1:2:s_im))];
            O = [O ((j-1)*s_im + (2:2:s_im))];
        end
    case 'horizontal'
        for j = 1:2:s_im
            E = [E ((j-1)*s_im + (1:s_im))];
            O = [O (j*s_im + (1:s_im))];
        end
end
beta = zeros(s_im*s_im,1);
beta(E) = 1;
beta(O) = -1;