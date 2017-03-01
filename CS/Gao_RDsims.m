function y =...
    Gao_RDsims(Wvec, Rvec, Kvec, dvec, kvec, rllseq)
% Simulations of Random Demodulator(RD) using RLL sequences.


% Use Wvec, Rvec, Kvec, dvec, and kvec for the desired parameter
% W, R, K, d, and k to run.  Note: 'dvec' and 'kvec' must be of equal
% length.
%
% 'rllseq' specifies what type of random chipping sequence to use:
%   'general'
%   'repcode'


verbose = 0; % Produce extra logs to the MATLAB command line

% number of iterations for each 4-tuple (K,R,W,d)

% RD Parameters
W_vec = Wvec;  % DFT (signal) size
R_vec = Rvec;  % Sampling rate
K_vec = Kvec;  % Sparsity level

% RLL sequence parameters: [d,k]-code
d_vec = dvec;   % min parameter of consective zeros
k_vec = kvec;   % max parameter

% threshold = 10^(-4);  % Desired precision for reconstruction

% Multi-dim matrix for storing success probability values, or MSE for the
% noisy case, for (K,R,W,d)


for var_W=1:length(W_vec)
    W = W_vec(var_W);
    
    % Create Unitary DFT matrix (size W x W)
    n=0:(W-1);  m=0:(W-1);
    F=exp(2*pi*1i*n'*m/W)/sqrt(W);  % Normalized DFT matrix
    if (W < R_vec(end))
        R_ind = sum((R_vec==W).*(1:length(R_vec)));
        R_use = R_vec(1:R_ind);
    else
        R_use = R_vec;
    end
    for var_R=1:length(R_use)
        R = R_use(var_R);
        r = (W/R);
        if verbose
            disp(['Running now with (W,R)=(' num2str(W) ',' num2str(R) ')']);
        end
        % Create summation-operation matrix (size RxW)
        H = create_H_matrix(R,W);
        K_max = min(length(K_vec),floor(W));
        for var_K=1:K_max
            K = K_vec(var_K);
            if verbose
                disp(['K = ' num2str(K)]);
            end
            for var_d=1:length(d_vec)
                d = d_vec(var_d);
                k = k_vec(var_d);
                
                % ---------------------------------------------
                    
                    % -----Create Sensing Matrix-------------------
                    % Create random sequence matrix D (size W x W)
                    t = gen_rll_waveform(d,k,W,rllseq); % rll switching sequence
                    D = diag(t); % diag matrix with switching sequence t
                    
                    % Create measurement matrix
                    Phi = H*D*F;
                    Phi = Phi/sqrt(r); % scale to orthonormal rows (for yall1)
                    

            end
        end
    end
end
y = Phi;
end % function

