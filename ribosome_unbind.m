    % ribosome_unbind(complex), where:
    % -- complex species denotes the  dissociating ribocomplex (M) 

    function r_or_mRNA = ribosome_unbind(complex)
        global ku;
        r_or_mRNA= ku * complex;
    end