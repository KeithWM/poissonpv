% clc
clear all
close all

% for C = 2.^(4:7);
for C = 1:12;
    disp(C);
    M = C*2;
    N = M;

    [A,B] = same_new(C,M,N, 0,0, 1);
    A = uint8(A);
    B = uint8(B);

    save(sprintf('orders/order_fractal_new_C%d', C), 'A', 'B');

    [A,B] = same(C,M,N, 0,0);
    A = uint8(A);
    B = uint8(B);

    save(sprintf('orders/order_fractal_C%d', C), 'A', 'B');
    
    P = C;
    P2= 2*P-1;
    A = mod(bsxfun(@plus, 0:P-1, (0:P2-1)'), P2);
    B = [P2*ones(P2,1), mod(bsxfun(@plus, P2-1:-1:P,  (0:P2-1)'), P2)];
    A = uint8([A; flipud(A)]);
    B = uint8([B; flipud(B)]);

    save(sprintf('orders/order_naive_C%d', C), 'A', 'B');
    
    A = mod(bsxfun(@plus, 0:P-1, (0:P2-1)'*(P-1)), P2);
    B = [P2*ones(P2,1), mod(bsxfun(@plus, P2-1:-1:P,  (0:P2-1)'*(P-1)), P2)];
    
    for k=2:2:P2
        Bt = fliplr(A(k,2:P));
        A(k,2:P) = fliplr(B(k,2:P));
        B(k,2:P) = Bt;
    end
    A = uint8([A; flipud(A)]);
    B = uint8([B; flipud(B)]);    

    save(sprintf('orders/order_smart_C%d', C), 'A', 'B');
end