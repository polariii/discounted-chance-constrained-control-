function y=H_matrix(N,A,B)
y=[];

for i=N:-1:1
    y=[y,(A^(i-1))*B];
end

