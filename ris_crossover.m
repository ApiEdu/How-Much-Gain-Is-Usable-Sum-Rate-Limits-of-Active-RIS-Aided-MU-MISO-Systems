function [child1, child2] = ris_crossover(parent1, parent2)
    % RIS_CROSSOVER Crossover for RIS optimization problem
    
    % Crossover for precoding vectors
    alpha = rand;
    child1.w = alpha * parent1.w + (1 - alpha) * parent2.w;
    child2.w = alpha * parent2.w + (1 - alpha) * parent1.w;
    
    % Crossover for RIS amplification matrix
    alpha = rand;
    p1 = diag(parent1.P);
    p2 = diag(parent2.P);
    child1.P = diag(alpha * p1 + (1 - alpha) * p2);
    child2.P = diag(alpha * p2 + (1 - alpha) * p1);
    
    % Crossover for RIS phase shifts
    alpha = rand;
    theta1 = angle(diag(parent1.Theta));
    theta2 = angle(diag(parent2.Theta));
    
    % Handle phase wrap-around
    diff = mod(theta2 - theta1 + pi, 2*pi) - pi;
    new_theta1 = theta1 + alpha * diff;
    new_theta2 = theta2 - (1 - alpha) * diff;
    
    child1.Theta = diag(exp(1j * new_theta1));
    child2.Theta = diag(exp(1j * new_theta2));
end