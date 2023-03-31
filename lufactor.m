function [L,U,P] = lufactor(A)
    U = A;
    [m,n] = size(A);
    L = eye(m);
    P = eye(m);
    for k = 1:1:m-1
        max = abs(U(k,k));
        maxindex = k;
        for i = k:1:m-1
            if (abs(U(i,k)) > max)
                max = abs(U(i,k));
                maxindex = i;
            end
        end    
        rowholder = U(k,k:m);
        U(k,k:m) = U(maxindex,k:m);
        U(maxindex,k:m) = rowholder;

        lrowholder = L(k,1:k-1);
        L(k,1:k-1) = L(maxindex,1:k-1);
        L(maxindex,1:k-1) = lrowholder;

        prowholder = P(k,1:m);
        P(k,1:m) = P(maxindex, 1:m);
        P(maxindex, 1:m) = prowholder;
        for j = k+1:1:m
            L(j,k) = U(j,k)/U(k,k);
            U(j,k:m) = U(j,k:m)-L(j,k)*U(k,k:m);
        end
    end
end
