function selected_indices = tournament_selection(fitness_values, tournament_size)
 % TOURNAMENT_SELECTION Selects individuals using tournament selection
    %
    % Inputs:
    %   fitness_values   : Vector of fitness values for each individual
    %   tournament_size  : Number of individuals in each tournament
    %
    % Output:
    %   selected_indices : Indices of selected individuals
    
    population_size = length(fitness_values);
    selected_indices = zeros(population_size, 1);
    
    % Ensure tournament_size doesn't exceed population_size
    tournament_size = min(tournament_size, population_size);
    
    for i = 1:population_size
        % Randomly select unique individuals for the tournament
        tournament_indices = randperm(population_size, tournament_size);
        %assignin('base','tournament_indices_result',selected_indices)
         
        
        % Find the best individual in the tournament
        [~, best_idx] = max(fitness_values(tournament_indices));
        selected_indices(i) = tournament_indices(best_idx);
    end
end