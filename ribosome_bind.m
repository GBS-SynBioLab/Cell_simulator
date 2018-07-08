    % ribosome_bind(species), where:
    % -- species denotes the  molecular species (M) bound
    function ribocomplex = ribosome_bind(species)
        global kb r;
        ribocomplex= kb * r * species;
    end