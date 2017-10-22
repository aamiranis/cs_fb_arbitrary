%% Function for non-linear approximation
function [f thresh] = find_kbest_biorth(f_w, nnz, d_level,Colorednodes)
thresh = 0;
if nnz == 0
    for level = 1:(d_level-1)
        f_w{level} = 0*f_w{level};
    end
    f_w{d_level}(:,2:4) = 0;
else
    lin_vec = [];
    for level = 1:d_level
        lin_vec = [lin_vec; f_w{level}(Colorednodes{level,2},2); f_w{level}(Colorednodes{level,3},3); f_w{level}(Colorednodes{level,4},4)];
    end
    nCoeffs = floor(nnz*length(lin_vec));
    lin_vec = sort(abs(lin_vec(:)),'descend');
    thresh = lin_vec(nCoeffs+1);
    for level = 1:d_level
        for F = 2:4
            temp = f_w{level}(:,F);
            temp(abs(temp) <thresh) = 0;
            f_w{level}(:,F) = temp;
        end
        %     temp(isolated{level},:) = f_w{level}(isolated{level},2:4);
    end
end
f = f_w;
end