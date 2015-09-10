function [A,B] = same(C,M,N, Is, Js)
%     Is
%     Js
    lvl = length(Is);
    p = localfactor(C);
    P = p(lvl);
    q = fliplr(cumprod(fliplr(p)));
    P2 = 2*P-1;
    Ct= prod(p(lvl:end));
    Mt= 2*Ct;
    Nt= N*Ct/C;
    
%     I = mod(bsxfun(@plus, 0:P-1, (0:P2-1)'), P2);
%     J = [P2*ones(P2,1), mod(bsxfun(@plus, P2-1:-1:P,  (0:P2-1)'*(P-1)), P2)];
    
    I = mod(bsxfun(@plus, 0:P-1, (0:P2-1)'*(P-1)), P2);
    J = [P2*ones(P2,1), mod(bsxfun(@plus, P2-1:-1:P,  (0:P2-1)'*(P-1)), P2)];
    
    for k=2:2:P2
        Jt = fliplr(I(k,2:P));
        I(k,2:P) = fliplr(J(k,2:P));
        J(k,2:P) = Jt;
    end
    
    A = zeros(2*(Mt-1), Ct);
    B = A;
    if P == Ct;
        % at the deepest level
        I0 = sum(Is.*q(1:lvl));
        J0 = sum(Js.*q(1:lvl));
        
        A = I0+[I; flipud(I)];
        B = J0+[J; flipud(J)];
    else
        % not the deepest level
        for k=1:P2-1
            ks = (k-1)*C/P+1:k*C/P;
            for l=1:P
                ls = (l-1)*C/P+1:l*C/P;
                [A(ks,ls),B(ks,ls)] = diff(C,M,N, [Is I(k,l)], [Js J(k,l)]);
            end
        end
        
        ks = (P2-1)*Ct/P+(1:2*(Mt/P-1));
        for l=1:P
            ls = (l-1)*Ct/P+1:l*Ct/P;
            [At,Bt] = same(Ct/P,Mt/P,Nt/P, 0,0);
            dIJ = abs(I(end,l)-J(end,l));
            mIJ = min(I(end,l),J(end,l));
            At(At >= Ct/P) = At(At >= Ct/P)+(dIJ-1)*Ct/P;
            Bt(Bt >= Ct/P) = Bt(Bt >= Ct/P)+(dIJ-1)*Ct/P;
            A(ks,ls) = mIJ*Ct/P+At;
            B(ks,ls) = mIJ*Ct/P+Bt;
        end
        
        for k=1:P2-1
            ks = 2*(Mt-1)+1-((k-1)*C/P+1:k*C/P);
            for l=1:P
                ls = (l-1)*C/P+1:l*C/P;
                [A(ks,ls),B(ks,ls)] = diff(C,M,N, [Is I(k,l)], [Js J(k,l)]);
            end
        end
    end
end

function [A,B] = diff(C,M,N, Is, Js)
    lvl = length(Is);
    p = localfactor(C);
    P = p(lvl);
    q = fliplr(cumprod(fliplr(p)));
    Ct= prod(p(lvl:end));
    
    I = kron(0:P-1, ones(P,1));
    J = mod(bsxfun(@plus, 0:P-1, (0:P-1)'), P);
    
    A = zeros(Ct, Ct);
    B = A;
    if P == Ct;
        I0 = sum(Is.*q(1:lvl));
        J0 = sum(Js.*q(1:lvl));
        % at the deepest level
        A = I0+I;
        B = J0+J;
    else
        for k=1:P
            ks = (k-1)*Ct/P+1:k*Ct/P;
            for l=1:P
                ls = (l-1)*Ct/P+1:l*Ct/P;
                [A(ks,ls),B(ks,ls)] = diff(C,M,N, [Is I(k,l)], [Js J(k,l)]);
            end
        end
    end
end

function p = localfactor(C)
    p = factor(C);
end