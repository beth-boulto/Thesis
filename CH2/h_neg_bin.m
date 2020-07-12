function[h] = h_neg_bin(z)

%negative binomial function from Basic Isham Model calculating $h(z), h'(z)
%and h''(z)
%import p and r values from main code
global r p
k=r;
h=zeros(3,1);
%calculate h(z)
h(1) = (1-p)^k/((1-p*z)^k); 
%alculate h'
h(2) = k*(p)*(1-p)^(k)/((1-p*z)^(k+1));
%calculate h''
h(3) = k*(k+1)*p^2*(1-p)^(k)/((1-p*z)^(k+2));

end

