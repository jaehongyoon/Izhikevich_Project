function [str, v, d, delta, tau] = type_(A)
% type of fixed point type for given Jacobian is calculated
tau = trace(A);
delta = det(A);
del_tau = tau^2 - 4*delta;

    if(delta < 0)
        % saddle points
        str = 'saddle point';
        % calculate eigenvector/eigenvalue
        [v,d]= eig(A);
    elseif((del_tau > 0) && (tau > 0) && (delta > 0))
        % unstable node
        str = 'unstable node';
        [v,d]= eig(A);
    elseif((del_tau > 0) && (tau < 0) && (delta > 0))
        % stable node
        str = 'stable node';
        [v,d]= eig(A);
    elseif((delta == 0))
        % non-isolated fixed node
        str = 'non-isolated fixed node';
        v = []; d = [];
    elseif((tau == 0))
        % center
        str = 'center';
        v = []; d = [];
    elseif((del_tau < 0) && (tau > 0) && (delta > 0))
        % unstable spirals
        str = 'unstable spirals';
        v = []; d = [];
    elseif((del_tau < 0) && (tau < 0) && (delta > 0))
        % stable spirals
        str = 'stable spirals';
        v = []; d = [];
    else
        str = 'ERROR';

    end
end