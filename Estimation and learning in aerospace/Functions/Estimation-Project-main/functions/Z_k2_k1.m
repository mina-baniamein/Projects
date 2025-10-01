function Z_out = Z_k2_k1(u, y, k2, k1)

if k1 > k2
    disp('error: k2 must be higher thah k1')
    return 
end

z = zeros(1,2*(k2 - k1 +1));
index = 1;


for i= k1:k2 
    if i == 0
        z(1,2*index-1) = 0;
        z(1,2*index) = 0;
    else
        z(1,2*index-1) = y(i);
        z(1,2*index) = u(i);
    end
    index = index + 1;
end

Z_out = (flip(z))';
end
