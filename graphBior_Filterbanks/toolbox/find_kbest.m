%% Function for non-linear approximation
function f = find_kbest(f_w, nnz, d_level,Colorednodes)
if nnz < 1
    lin_vec = [];
    for level = 1:d_level
        lin_vec = [lin_vec; f_w{level}(Colorednodes{level,2},2); f_w{level}(Colorednodes{level,3},3); f_w{level}(Colorednodes{level,4},4)];
    end
    nCoeffs = floor(length(lin_vec)*nnz);
    lin_vec = sort(abs(lin_vec(:)),'descend');
    thresh = lin_vec(nCoeffs+1);
    for level = 1:d_level
        temp = f_w{level}(:,2:4);
        temp(abs(temp) <thresh) = 0;
        %     temp(isolated{level},:) = f_w{level}(isolated{level},2:4);
        f_w{level}(:,2:4) = temp;
    end
end
f = f_w;
end