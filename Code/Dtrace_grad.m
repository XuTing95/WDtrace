function [f, grad] = Dtrace_grad(S1, S2, Delta)

    b = S1 - S2;
    A = S1 * Delta * S2;
    B = S2 * Delta * S1;
    grad = 0.5 * (A+B) - b;
    f = 0.5 * sum(diag(Delta * A)) - sum(diag(Delta * b));
    
end