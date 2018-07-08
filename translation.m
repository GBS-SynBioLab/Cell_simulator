    % translation(ribocomplex,codons) where:
    % -- ribocomplex denotes the ribosome/mRNA complex species (M)
    % -- codons denotes size of mRNA in bases triplets
    function p = translation(ribocomplex,codons)
        global gamma;
        p = ribocomplex/codons*gamma;
    end