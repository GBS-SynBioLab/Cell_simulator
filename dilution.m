    % dilution(species), where:
    % -- species denotes the diluted molecular species (M/sec) diluted
    function dil_rate = dilution(species)
        global lam;
        dil_rate= lam * species;
    end