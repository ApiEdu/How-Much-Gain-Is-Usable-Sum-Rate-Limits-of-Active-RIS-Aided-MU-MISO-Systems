function offspring = generate_offspring(selected_population, population_size, crossover_func)
    % GENERATE_OFFSPRING Creates offspring through crossover
    %
    % Inputs:
    %   selected_population : Cell array of selected individuals
    %   population_size     : Desired size of offspring population
    %   crossover_func      : Function handle for crossover operation
    %
    % Output:
    %   offspring          : Cell array of offspring individuals
    
    num_selected = length(selected_population);
    offspring = cell(population_size, 1);
    
    % Shuffle the selected population to ensure random pairing
    shuffled_indices = randperm(num_selected);
    shuffled_population = selected_population(shuffled_indices);
    
    % Create offspring through crossover
    for i = 1:2:population_size-1
        % Get parents (cycling through if needed)
        parent1_idx = mod(i-1, num_selected) + 1;
        parent2_idx = mod(i, num_selected) + 1;
        
        parent1 = shuffled_population{parent1_idx};
        parent2 = shuffled_population{parent2_idx};
        
        % Perform crossover (assuming it returns two children)
        [child1, child2] = crossover_func(parent1, parent2);
        
        offspring{i} = child1;
        offspring{i+1} = child2;
    end
    
    % Handle odd population size by cloning the last individual
    if mod(population_size, 2) == 1
        offspring{population_size} = shuffled_population{mod(population_size-1, num_selected) + 1};
    end
end