function [c,s] = givens(x1,x2)
%GIVENS gives elements in Givens Rotation [c,s;-s,c]. Follows Algorithm
% 5.1.3 in text book (Golub & Van Loan, 2013).
%    x1,x2: vecor to apply Givens rotation
%    c,s: elements in Givens rotation
    if abs(x2)>abs(x1)
        t = -x1/x2;
        s = 1/sqrt(1+t^2);
        c = s*t;
    else
        t = -x2/x1;
        c = 1/sqrt(1+t^2);
        s = c*t;
    end
end