function [GC GC_inv] = compute_GC_Bior(Ln_bpt,dfactor,hi_d,lo_d,normtype)

%% compute square kernels and their Chebychev approximations
filterlen_hi = length(roots(hi_d));
filterlen_lo = length(roots(lo_d));
arange=  [0 2];
h0 = @(x)(polyval(lo_d,x));
h1 = @(x)(polyval(hi_d,x));
c_d{1}=sgwt_cheby_coeff(h0,filterlen_lo,filterlen_lo+1,arange);
c_d{2}=sgwt_cheby_coeff(h1,filterlen_hi,filterlen_hi+1,arange);

% p_lo= conv(lo_d,lo_d);
% p_hi = conv(hi_d,hi_d);
% p0 = @(x)(polyval(p_lo,x));
% p1 = @(x)(polyval(p_hi,x));
% c_p{1} = sgwt_cheby_coeff(p0,2*filterlen_lo,2*filterlen_lo+1,arange);
% c_p{2} = sgwt_cheby_coeff(p1,2*filterlen_hi,2*filterlen_hi+1,arange);

%% Estimate the diagonal of corresponding filters.
max_level = size(Ln_bpt,1); % # decomposition levels
theta = size(Ln_bpt,2); % # bpt graphs
GC = cell(max_level,theta);
GC_inv = cell(max_level,theta);
switch normtype
    case 'asym'
        for level = 1:max_level
            N = length(Ln_bpt{level,1});
            s_im = sqrt(N);
            for i = 1:theta
                fin = ones(N,1);
                fout = sgwt_cheby_op(fin,Ln_bpt{level,i},c_d{1},arange);
                GC_inv{level,i}(:,1) =fout;
                GC{level,i}(:,1)= 1./(fout+10^-16);
                if i ==1
                    fin= image_downsampling_fn(s_im,'rectangle');
                else
                    fin = image_downsampling_fn(s_im,'diamond');
                end
                fout = sgwt_cheby_op(fin,Ln_bpt{level,i},c_d{2},arange);
                GC_inv{level,i}(:,2) =fout./fin;
                GC{level,i}(:,2)= fin./(fout+10^-16);
                
            end
        end
    case 'sym'
        
        for level = 1:max_level
            N = length(Ln_bpt{level,1});
            s_im = sqrt(N);
            for i = 1:theta
                fin = dfactor{level}(:,i).*ones(N,1);
                fout = sgwt_cheby_op(fin,Ln_bpt{level,i},c_d{1},arange);
                GC_inv{level,i}(:,1) =fout./(fin +10^-16);
                GC{level,i}(:,1)= fin./(fout+10^-16);
                if i ==1
                    fin= image_downsampling_fn(s_im,'rectangle');
                else
                    fin = image_downsampling_fn(s_im,'diamond');
                end
                fin = dfactor{level}(:,i).*fin;
                fout = sgwt_cheby_op(fin,Ln_bpt{level,i},c_d{2},arange);
                GC_inv{level,i}(:,2) =fout./(fin + 10^-16);
                GC{level,i}(:,2)= fin./(fout + 10^-16);
                
            end
        end
    otherwise
        disp('unknown option')
        return;
end
