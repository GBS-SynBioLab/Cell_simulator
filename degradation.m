    % degradation(species), where:
    % -- species denotes the diluted molecular species (M) diluted
    function deg_rate = degradation(species)
        global  dm
        deg_rate= dm * species;
    end