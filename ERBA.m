function dsites = ERBA(f,dsites,rbf,rho,mode,tol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The implementation of the Efficient Raduced Basis Algorithm (ERBA)
% Calls on: DistanceMatrix.m by G. Fasshauer
%
% Input
%   f : function handle (not needed in the power-based case)
%   dsites: Mxs matrix representing a set of M data sites in R^s
%              (i.e., each row contains one s-dimensional point)
%   rbf: function handle, the chosen RBF for the interpolation
%   rho: positive integer, the number of nodes to be removed at each step
%   mode: 0 (power-based case) or 1 (residual-based case)
%   tol: positive real number, the tolerance used in the algorithm
%
% Output
%   dsites: Nxs matrix representing the set of selected data sites (N <= M)
%
% To use this function, please cite:
% F. Marchetti, E. Perracchione, "Efficient Reduced Basis Algorithm (ERBA) 
% for kernel-based approximation", submitted.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('ERBA \n\n')

n_dsites = size(dsites,1);
n_folds = floor(n_dsites/rho);

rng(42);
folds = cvpartition(n_dsites,'KFold',n_folds);
w = 0;

if mode == 0
    
    while w < tol
    
        invIM = pinv(rbf(DistanceMatrix(dsites,dsites)));

        pow_array = zeros(n_folds,1);

        for i=1:n_folds

            test_ind = test(folds,i);
            powfun = real(sqrt(sum((invIM(:,test_ind)*pinv(invIM(test_ind,test_ind))).*rbf(DistanceMatrix(dsites,dsites(test_ind,:))),1)));
            pow_array(i) = sqrt(mean(powfun.^2));

        end    

        [w,i_star] = min(pow_array);
        
        if w < tol
            
            dsites(test(folds,i_star),:) = [];
            n_dsites = size(dsites,1);
            n_folds = floor(n_dsites/rho);
            
            if floor(n_dsites/rho) <= 1
                
                fprintf('Number of remaining data not sufficient to continue \n\n')
                return
                
            end

            rng(42);
            folds = cvpartition(n_dsites,'KFold',n_folds);
            
        end
        
        fprintf('Number of remaining data: %d \n', n_dsites)
        fprintf('Value of the weight: %e \n\n', w)        
    
    end
end

if mode == 1
    
    while w < tol
        
        rhs = f(dsites(:,1),dsites(:,2));
        invIM = pinv(rbf(DistanceMatrix(dsites,dsites)));
        coeffs = invIM*rhs;

        res_array = zeros(n_folds,1);

        for i=1:n_folds

            test_ind = test(folds,i);
            residual = invIM(test_ind,test_ind)\coeffs(test_ind);
            res_array(i) = sqrt(mean(residual.^2));

        end    

        [w,i_star] = min(res_array);
        
        if w < tol
            
            dsites(test(folds,i_star),:) = [];
            n_dsites = size(dsites,1);
            n_folds = floor(n_dsites/rho);
            
            if floor(n_dsites/rho) <= 1
                
                fprintf('Number of remaining data not sufficient to continue \n\n')
                return
                
            end

            rng(42);
            folds = cvpartition(n_dsites,'KFold',n_folds);
            
        end
        
        fprintf('Number of remaining data: %d \n', n_dsites)
        fprintf('Value of the weight: %e \n\n', w)     
    
    end
end

end
